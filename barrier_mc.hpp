// barrier_mc.hpp
//
// Monte Carlo pricing of single-barrier European options (knock-in / knock-out).
//   * Discrete monitoring at `steps` points, with the Broadie-Glasserman-Kou
//     continuity correction (barrier shifted by 0.5826 * sigma * sqrt(dt)).
//   * Knock-in and knock-out are computed on the same paths, so
//     knock_in + knock_out == vanilla holds path by path.
#pragma once

#include "derivatives_pricer.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>

namespace fin {

enum class BarrierDirection { Up, Down };

struct BarrierResult {
    PriceResult knock_in;   // pays only if the barrier is touched at a monitoring date
    PriceResult knock_out;  // pays only if the barrier is never touched
    PriceResult vanilla;    // knock_in.price + knock_out.price == vanilla.price
};

// Broadie, Glasserman & Kou (1997) correction for discretely monitored barriers.
inline constexpr double kBroadieGlassermanKouBeta = 0.5826;

namespace detail {

// Welford online mean and variance.
struct Welford {
    double mean = 0.0;
    double m2 = 0.0;

    void add(double x, std::size_t n_before) {
        const double d = x - mean;
        mean += d / static_cast<double>(n_before + 1);
        m2 += d * (x - mean);
    }

    PriceResult result(std::size_t n) const {
        const double var = m2 / static_cast<double>(n - 1);
        return {mean, std::sqrt(var / static_cast<double>(n))};
    }
};

}  // namespace detail

inline BarrierResult barrier_monte_carlo(const Option& o, const Market& m, double barrier,
                                         BarrierDirection dir, std::size_t paths,
                                         std::size_t steps, std::uint64_t seed) {
    validate(o, m);
    require_european(o, "barrier_monte_carlo");
    if (barrier <= 0.0)
        throw std::invalid_argument("barrier must be positive");
    if (dir == BarrierDirection::Up && barrier <= m.spot)
        throw std::invalid_argument("up barrier must be above spot");
    if (dir == BarrierDirection::Down && barrier >= m.spot)
        throw std::invalid_argument("down barrier must be below spot");
    if (paths < 2 || steps < 1)
        throw std::invalid_argument("barrier_monte_carlo: need >= 2 paths and >= 1 step");

    const bool up = (dir == BarrierDirection::Up);
    const double T = o.maturity;
    const double dt = T / static_cast<double>(steps);
    const double sdt = m.volatility * std::sqrt(dt);
    const double drift = (m.rate - m.dividend - 0.5 * m.volatility * m.volatility) * dt;
    const double bgk = kBroadieGlassermanKouBeta * sdt;
    const double log_barrier = std::log(barrier) + (up ? bgk : -bgk);
    const double df = std::exp(-m.rate * T);
    const double log_s0 = std::log(m.spot);

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> std_normal(0.0, 1.0);
    detail::Welford ki_acc, ko_acc, v_acc;

    for (std::size_t k = 0; k < paths; ++k) {
        double x = log_s0;
        bool hit = false;
        for (std::size_t s = 0; s < steps; ++s) {
            x += drift + sdt * std_normal(rng);
            // Keep drawing after a hit so every barrier level consumes the same random
            // numbers. This makes comparisons across barrier levels exact path by path.
            if (up ? (x >= log_barrier) : (x <= log_barrier)) hit = true;
        }
        const double v = intrinsic(o.type, o.strike, std::exp(x)) * df;
        ki_acc.add(hit ? v : 0.0, k);
        ko_acc.add(hit ? 0.0 : v, k);
        v_acc.add(v, k);
    }

    return {ki_acc.result(paths), ko_acc.result(paths), v_acc.result(paths)};
}

}  // namespace fin
