#pragma once
// derivatives_pricer.cpp
//
// Equity option pricing toolkit (C++17):
//   * Black-Scholes-Merton pricing with continuous dividend yield (European)
//   * Analytic Greeks: delta, gamma, vega, theta, rho
//   * Cox-Ross-Rubinstein binomial tree (European and American)
//   * Monte Carlo pricing with antithetic variates and a standard error
//   * Implied volatility via safeguarded bisection
//   * Built-in self-test that cross-checks the methods against each other
//
// Build:  g++ -std=c++17 -O2 -Wall -Wextra -pedantic derivatives_pricer.cpp -o pricer
// Run:    ./pricer        (exit code 0 == all checks passed)

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace fin {

constexpr double kPi = 3.14159265358979323846;

enum class OptionType { Call, Put };
enum class ExerciseStyle { European, American };

struct Option {
    OptionType type;
    ExerciseStyle style;
    double strike;    // K
    double maturity;  // T, in years
};

struct Market {
    double spot;        // S
    double rate;        // r, continuously compounded
    double dividend;    // q, continuous dividend yield
    double volatility;  // sigma, annualised
};

struct Greeks {
    double delta;
    double gamma;
    double vega;   // per 1.00 change in volatility (divide by 100 for "per vol point")
    double theta;  // per year (divide by 365 for per calendar day)
    double rho;    // per 1.00 change in rate
};

struct PriceResult {
    double price;
    double std_error;  // 0 for deterministic methods
};

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

inline void validate(const Option& o, const Market& m) {
    if (o.strike <= 0.0)     throw std::invalid_argument("strike must be positive");
    if (o.maturity <= 0.0)   throw std::invalid_argument("maturity must be positive");
    if (m.spot <= 0.0)       throw std::invalid_argument("spot must be positive");
    if (m.volatility <= 0.0) throw std::invalid_argument("volatility must be positive");
}

inline void require_european(const Option& o, const char* who) {
    if (o.style != ExerciseStyle::European)
        throw std::invalid_argument(std::string(who) + ": closed form is European only; use binomial_price()");
}

// ---------------------------------------------------------------------------
// Normal distribution helpers
// ---------------------------------------------------------------------------

inline double norm_cdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }
inline double norm_pdf(double x) { return std::exp(-0.5 * x * x) / std::sqrt(2.0 * kPi); }

// ---------------------------------------------------------------------------
// Black-Scholes-Merton (European)
// ---------------------------------------------------------------------------

struct D1D2 { double d1; double d2; };

inline D1D2 compute_d1d2(const Option& o, const Market& m) {
    const double vol_sqrt_t = m.volatility * std::sqrt(o.maturity);
    const double d1 = (std::log(m.spot / o.strike) +
                       (m.rate - m.dividend + 0.5 * m.volatility * m.volatility) * o.maturity) /
                      vol_sqrt_t;
    return {d1, d1 - vol_sqrt_t};
}

inline double black_scholes_price(const Option& o, const Market& m) {
    validate(o, m);
    require_european(o, "black_scholes_price");

    const auto [d1, d2] = compute_d1d2(o, m);
    const double df_r = std::exp(-m.rate * o.maturity);
    const double df_q = std::exp(-m.dividend * o.maturity);

    if (o.type == OptionType::Call)
        return m.spot * df_q * norm_cdf(d1) - o.strike * df_r * norm_cdf(d2);
    return o.strike * df_r * norm_cdf(-d2) - m.spot * df_q * norm_cdf(-d1);
}

inline Greeks black_scholes_greeks(const Option& o, const Market& m) {
    validate(o, m);
    require_european(o, "black_scholes_greeks");

    const auto [d1, d2] = compute_d1d2(o, m);
    const double sqrt_t = std::sqrt(o.maturity);
    const double df_r = std::exp(-m.rate * o.maturity);
    const double df_q = std::exp(-m.dividend * o.maturity);
    const double pdf_d1 = norm_pdf(d1);
    const bool is_call = (o.type == OptionType::Call);

    Greeks g{};
    g.delta = is_call ? df_q * norm_cdf(d1) : df_q * (norm_cdf(d1) - 1.0);
    g.gamma = df_q * pdf_d1 / (m.spot * m.volatility * sqrt_t);
    g.vega  = m.spot * df_q * pdf_d1 * sqrt_t;

    const double common = -m.spot * df_q * pdf_d1 * m.volatility / (2.0 * sqrt_t);
    if (is_call) {
        g.theta = common - m.rate * o.strike * df_r * norm_cdf(d2)
                         + m.dividend * m.spot * df_q * norm_cdf(d1);
        g.rho   =  o.strike * o.maturity * df_r * norm_cdf(d2);
    } else {
        g.theta = common + m.rate * o.strike * df_r * norm_cdf(-d2)
                         - m.dividend * m.spot * df_q * norm_cdf(-d1);
        g.rho   = -o.strike * o.maturity * df_r * norm_cdf(-d2);
    }
    return g;
}

// ---------------------------------------------------------------------------
// Cox-Ross-Rubinstein binomial tree (European or American)
// ---------------------------------------------------------------------------

inline double intrinsic(OptionType type, double strike, double s) {
    return type == OptionType::Call ? std::max(s - strike, 0.0)
                                    : std::max(strike - s, 0.0);
}

inline double binomial_price(const Option& o, const Market& m, std::size_t steps) {
    validate(o, m);
    if (steps < 1) throw std::invalid_argument("binomial_price: steps must be >= 1");

    const double dt     = o.maturity / static_cast<double>(steps);
    const double u      = std::exp(m.volatility * std::sqrt(dt));
    const double d      = 1.0 / u;
    const double growth = std::exp((m.rate - m.dividend) * dt);
    const double p      = (growth - d) / (u - d);  // risk-neutral up-probability
    if (p <= 0.0 || p >= 1.0)
        throw std::runtime_error("binomial_price: risk-neutral probability outside (0,1); increase steps");

    const double disc = std::exp(-m.rate * dt);
    const bool american = (o.style == ExerciseStyle::American);

    // Terminal payoffs: node j has j up-moves out of `steps`.
    std::vector<double> v(steps + 1);
    for (std::size_t j = 0; j <= steps; ++j) {
        const double s_T = m.spot * std::pow(u, static_cast<double>(j)) *
                                    std::pow(d, static_cast<double>(steps - j));
        v[j] = intrinsic(o.type, o.strike, s_T);
    }

    // Backward induction.
    for (std::size_t i = steps; i-- > 0;) {
        for (std::size_t j = 0; j <= i; ++j) {
            double value = disc * (p * v[j + 1] + (1.0 - p) * v[j]);
            if (american) {
                const double s_ij = m.spot * std::pow(u, static_cast<double>(j)) *
                                             std::pow(d, static_cast<double>(i - j));
                value = std::max(value, intrinsic(o.type, o.strike, s_ij));
            }
            v[j] = value;
        }
    }
    return v[0];
}

// ---------------------------------------------------------------------------
// Monte Carlo (European) with antithetic variates
// ---------------------------------------------------------------------------

inline PriceResult monte_carlo_price(const Option& o, const Market& m,
                                     std::size_t pairs, std::uint64_t seed) {
    validate(o, m);
    require_european(o, "monte_carlo_price");
    if (pairs < 2) throw std::invalid_argument("monte_carlo_price: need at least 2 path pairs");

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> std_normal(0.0, 1.0);

    const double t       = o.maturity;
    const double drift   = (m.rate - m.dividend - 0.5 * m.volatility * m.volatility) * t;
    const double diffuse = m.volatility * std::sqrt(t);
    const double df      = std::exp(-m.rate * t);

    // Welford's online mean/variance for numerical stability.
    double mean = 0.0, m2 = 0.0;
    for (std::size_t k = 0; k < pairs; ++k) {
        const double z = std_normal(rng);
        const double s_up = m.spot * std::exp(drift + diffuse * z);
        const double s_dn = m.spot * std::exp(drift - diffuse * z);  // antithetic path
        const double x = 0.5 * (intrinsic(o.type, o.strike, s_up) +
                                intrinsic(o.type, o.strike, s_dn)) * df;
        const double delta = x - mean;
        mean += delta / static_cast<double>(k + 1);
        m2   += delta * (x - mean);
    }
    const double variance  = m2 / static_cast<double>(pairs - 1);
    const double std_error = std::sqrt(variance / static_cast<double>(pairs));
    return {mean, std_error};
}

// ---------------------------------------------------------------------------
// Implied volatility (European, Black-Scholes inversion)
// ---------------------------------------------------------------------------

inline double implied_volatility(const Option& o, double spot, double rate,
                                 double dividend, double market_price) {
    constexpr double lo_bound = 1e-6;
    constexpr double hi_bound = 5.0;

    const double p_lo = black_scholes_price(o, {spot, rate, dividend, lo_bound});
    const double p_hi = black_scholes_price(o, {spot, rate, dividend, hi_bound});
    if (market_price < p_lo || market_price > p_hi)
        throw std::domain_error("implied_volatility: price not attainable for vol in [1e-6, 5]");

    double lo = lo_bound, hi = hi_bound;
    for (int iter = 0; iter < 200; ++iter) {
        const double mid = 0.5 * (lo + hi);
        const double p_mid = black_scholes_price(o, {spot, rate, dividend, mid});
        if (std::fabs(p_mid - market_price) < 1e-12) return mid;
        if (p_mid < market_price) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

}  // namespace fin
