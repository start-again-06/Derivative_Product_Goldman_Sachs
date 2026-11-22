// vol_surface.hpp
//
// Implied volatility surface: a rectangular grid of implied vols over
// (strike, maturity), with bilinear interpolation and flat extrapolation.
// Prices are Black-Scholes prices evaluated at the surface vol for each option.
#pragma once

#include "derivatives_pricer.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fin {

class VolSurface {
public:
    // vols[t][k]: implied volatility at maturities[t] and strikes[k].
    // Both grids must be strictly increasing.
    VolSurface(std::vector<double> strikes, std::vector<double> maturities,
               std::vector<std::vector<double>> vols)
        : strikes_(std::move(strikes)),
          maturities_(std::move(maturities)),
          vols_(std::move(vols)) {
        require_grid(strikes_, "strikes");
        require_grid(maturities_, "maturities");
        if (vols_.size() != maturities_.size())
            throw std::invalid_argument("VolSurface: one vols row per maturity required");
        for (const auto& row : vols_) {
            if (row.size() != strikes_.size())
                throw std::invalid_argument("VolSurface: one vols entry per strike required");
            for (double v : row)
                if (v <= 0.0) throw std::invalid_argument("VolSurface: vols must be positive");
        }
    }

    // Bilinear in (strike, maturity); held flat outside the grid.
    double vol(double strike, double maturity) const {
        if (strike <= 0.0 || maturity <= 0.0)
            throw std::invalid_argument("VolSurface::vol: strike and maturity must be positive");
        const auto [ki, kw] = locate(strikes_, strike);
        const auto [ti, tw] = locate(maturities_, maturity);
        const double v_lo = (1.0 - kw) * vols_[ti][ki]     + kw * vols_[ti][ki + 1];
        const double v_hi = (1.0 - kw) * vols_[ti + 1][ki] + kw * vols_[ti + 1][ki + 1];
        return (1.0 - tw) * v_lo + tw * v_hi;
    }

    // European Black-Scholes price using the surface vol at the option's strike and maturity.
    double price(const Option& o, double spot, double rate, double dividend) const {
        return black_scholes_price(o, Market{spot, rate, dividend, vol(o.strike, o.maturity)});
    }

private:
    static void require_grid(const std::vector<double>& g, const char* name) {
        if (g.size() < 2)
            throw std::invalid_argument(std::string("VolSurface: ") + name + " needs >= 2 points");
        for (std::size_t i = 1; i < g.size(); ++i)
            if (g[i] <= g[i - 1])
                throw std::invalid_argument(std::string("VolSurface: ") + name + " must be strictly increasing");
    }

    // Returns the lower cell index i (0..n-2) and weight w in [0,1] with x = (1-w) g[i] + w g[i+1].
    static std::pair<std::size_t, double> locate(const std::vector<double>& g, double x) {
        if (x <= g.front()) return {0, 0.0};
        if (x >= g.back())  return {g.size() - 2, 1.0};
        const auto it = std::upper_bound(g.begin(), g.end(), x);
        const std::size_t i = static_cast<std::size_t>(it - g.begin()) - 1;
        return {i, (x - g[i]) / (g[i + 1] - g[i])};
    }

    std::vector<double> strikes_;
    std::vector<double> maturities_;
    std::vector<std::vector<double>> vols_;
};

}  // namespace fin
