// yield_curve.hpp
//
// Fixed-income layer:
//   * YieldCurve: zero curve with linear interpolation on zero rates
//   * Bond analytics: cash flows, pricing off the curve or a yield, YTM,
//     modified duration, convexity
//   * Ho-Lee short-rate lattice calibrated exactly to the curve, used to
//     price callable bonds (issuer call option)
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace fin {

// ---------------------------------------------------------------- Zero curve

class YieldCurve {
public:
    // tenors in years (strictly increasing, > 0); zero rates continuously compounded.
    YieldCurve(std::vector<double> tenors, std::vector<double> zero_rates)
        : tenors_(std::move(tenors)), zeros_(std::move(zero_rates)) {
        if (tenors_.size() < 2 || tenors_.size() != zeros_.size())
            throw std::invalid_argument("YieldCurve: matching tenor/rate arrays with >= 2 points required");
        if (tenors_.front() <= 0.0)
            throw std::invalid_argument("YieldCurve: tenors must be positive");
        for (std::size_t i = 1; i < tenors_.size(); ++i)
            if (tenors_[i] <= tenors_[i - 1])
                throw std::invalid_argument("YieldCurve: tenors must be strictly increasing");
    }

    // Linear interpolation on zero rates; flat extrapolation at both ends.
    double zero_rate(double t) const {
        if (t <= tenors_.front()) return zeros_.front();
        if (t >= tenors_.back())  return zeros_.back();
        const auto it = std::upper_bound(tenors_.begin(), tenors_.end(), t);
        const std::size_t i = static_cast<std::size_t>(it - tenors_.begin()) - 1;
        const double w = (t - tenors_[i]) / (tenors_[i + 1] - tenors_[i]);
        return (1.0 - w) * zeros_[i] + w * zeros_[i + 1];
    }

    double discount(double t) const {
        if (t <= 0.0) return 1.0;
        return std::exp(-zero_rate(t) * t);
    }

    // Continuously compounded forward rate between t1 and t2.
    double forward_rate(double t1, double t2) const {
        if (t1 < 0.0 || t2 <= t1)
            throw std::invalid_argument("forward_rate: need 0 <= t1 < t2");
        return -(std::log(discount(t2)) - std::log(discount(t1))) / (t2 - t1);
    }

private:
    std::vector<double> tenors_;
    std::vector<double> zeros_;
};

// ---------------------------------------------------------------- Bonds

struct BondSpec {
    double face = 100.0;
    double coupon_rate = 0.0;  // annual coupon as a fraction of face; paid at t = 1..maturity
    int maturity_years = 1;
    double call_price = 0.0;   // issuer call price (ex-coupon); <= 0 means non-callable
    int first_call_year = 1;   // calls permitted on coupon dates at years >= this
};

struct CashFlow {
    double t;
    double amount;
};

inline std::vector<CashFlow> bond_cashflows(const BondSpec& b) {
    if (b.maturity_years < 1)
        throw std::invalid_argument("bond maturity must be at least 1 year");
    std::vector<CashFlow> cfs;
    for (int k = 1; k <= b.maturity_years; ++k) {
        double amount = b.coupon_rate * b.face;
        if (k == b.maturity_years) amount += b.face;
        cfs.push_back({static_cast<double>(k), amount});
    }
    return cfs;
}

inline double price_from_curve(const std::vector<CashFlow>& cfs, const YieldCurve& curve) {
    double pv = 0.0;
    for (const auto& cf : cfs) pv += cf.amount * curve.discount(cf.t);
    return pv;
}

// Annual-compounding yield.
inline double price_from_yield(const std::vector<CashFlow>& cfs, double y) {
    double pv = 0.0;
    for (const auto& cf : cfs) pv += cf.amount / std::pow(1.0 + y, cf.t);
    return pv;
}

inline double yield_to_maturity(const std::vector<CashFlow>& cfs, double price) {
    double lo = -0.99, hi = 5.0;
    if (price_from_yield(cfs, lo) < price || price_from_yield(cfs, hi) > price)
        throw std::domain_error("yield_to_maturity: price outside solvable range");
    for (int it = 0; it < 300 && hi - lo > 1e-14; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (price_from_yield(cfs, mid) > price) lo = mid; else hi = mid;  // price falls as yield rises
    }
    return 0.5 * (lo + hi);
}

inline double modified_duration(const std::vector<CashFlow>& cfs, double y) {
    const double p = price_from_yield(cfs, y);
    double weighted = 0.0;
    for (const auto& cf : cfs) weighted += cf.t * cf.amount / std::pow(1.0 + y, cf.t);
    return weighted / (p * (1.0 + y));
}

inline double convexity(const std::vector<CashFlow>& cfs, double y) {
    const double p = price_from_yield(cfs, y);
    double c = 0.0;
    for (const auto& cf : cfs) c += cf.t * (cf.t + 1.0) * cf.amount / std::pow(1.0 + y, cf.t + 2.0);
    return c / p;
}

// ---------------------------------------------------------- Ho-Lee lattice

struct HoLeeLattice {
    int years = 0;
    int steps_per_year = 0;
    double dt = 0.0;
    // rate[i][j]: short rate at step i, node j (j up-moves), applied over [t_i, t_{i+1}].
    std::vector<std::vector<double>> rate;
};

// Ho-Lee (normal) model on a recombining lattice:
//     r(i, j) = alpha_i + sigma * sqrt(dt) * (2j - i)
// alpha_i is solved forward in time from Arrow-Debreu prices so that the lattice
// reprices every discount factor on the grid exactly. sigma is an absolute (normal)
// rate volatility, e.g. 0.01 = 100bp per sqrt(year).
inline HoLeeLattice build_ho_lee(const YieldCurve& curve, double sigma, int years, int steps_per_year) {
    if (sigma <= 0.0 || years < 1 || steps_per_year < 1)
        throw std::invalid_argument("build_ho_lee: sigma, years and steps_per_year must be positive");

    HoLeeLattice lat;
    lat.years = years;
    lat.steps_per_year = steps_per_year;
    lat.dt = 1.0 / steps_per_year;
    const int n_steps = years * steps_per_year;
    const double sd = sigma * std::sqrt(lat.dt);

    std::vector<double> q{1.0};  // Arrow-Debreu prices at the current time level
    for (int i = 0; i < n_steps; ++i) {
        const double target = curve.discount((i + 1) * lat.dt);
        auto pv = [&](double alpha) {
            double s = 0.0;
            for (int j = 0; j <= i; ++j)
                s += q[j] * std::exp(-(alpha + sd * (2.0 * j - i)) * lat.dt);
            return s;  // strictly decreasing in alpha
        };

        double lo = -1.0, hi = 1.0;
        if (pv(lo) < target || pv(hi) > target)
            throw std::runtime_error("build_ho_lee: calibration bracket failed");
        for (int it = 0; it < 300 && hi - lo > 1e-15; ++it) {
            const double mid = 0.5 * (lo + hi);
            if (pv(mid) > target) lo = mid; else hi = mid;
        }
        const double alpha = 0.5 * (lo + hi);

        std::vector<double> r(i + 1);
        for (int j = 0; j <= i; ++j) r[j] = alpha + sd * (2.0 * j - i);
        lat.rate.push_back(r);

        // Propagate Arrow-Debreu prices one step forward (risk-neutral p = 1/2).
        std::vector<double> next(i + 2, 0.0);
        for (int j = 0; j <= i + 1; ++j) {
            const double from_up   = (j >= 1) ? q[j - 1] * std::exp(-r[j - 1] * lat.dt) : 0.0;
            const double from_down = (j <= i) ? q[j]     * std::exp(-r[j]     * lat.dt) : 0.0;
            next[j] = 0.5 * (from_up + from_down);
        }
        q = std::move(next);
    }
    return lat;
}

// Price a fixed-rate bond on the lattice. If call_price > 0, the issuer may call at
// coupon dates from first_call_year onward, paying call_price plus that year's coupon.
// The lattice horizon must equal the bond maturity.
inline double lattice_bond_price(const HoLeeLattice& lat, const BondSpec& b) {
    if (b.maturity_years != lat.years)
        throw std::invalid_argument("lattice_bond_price: lattice horizon must equal bond maturity");

    const int spy = lat.steps_per_year;
    const int n = lat.years * spy;
    const bool callable = b.call_price > 0.0;
    const double coupon = b.coupon_rate * b.face;

    std::vector<double> v(n + 1, b.face + coupon);  // maturity node: final principal + coupon
    for (int i = n - 1; i >= 0; --i) {
        const bool coupon_date = (i > 0 && i % spy == 0);
        const double c = coupon_date ? coupon : 0.0;
        const bool call_date = callable && coupon_date && (i / spy) >= b.first_call_year;
        for (int j = 0; j <= i; ++j) {
            double cont = std::exp(-lat.rate[i][j] * lat.dt) * 0.5 * (v[j] + v[j + 1]);
            if (call_date) cont = std::min(cont, b.call_price);
            v[j] = cont + c;
        }
    }
    return v[0];
}

}  // namespace fin
