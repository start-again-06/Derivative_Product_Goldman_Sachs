// extensions_test.cpp - demo and self-tests for the three extension modules:
//   1. Volatility surface   (vol_surface.hpp)
//   2. Yield curve / bonds  (yield_curve.hpp)  incl. Ho-Lee callable bond lattice
//   3. Barrier options MC   (barrier_mc.hpp)
//
// Build: g++ -std=c++17 -O2 -Wall -Wextra -pedantic extensions_test.cpp -o ext_test
// Run:   ./ext_test   (exit code 0 == all checks passed)

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "derivatives_pricer.hpp"
#include "vol_surface.hpp"
#include "yield_curve.hpp"
#include "barrier_mc.hpp"

namespace {

int g_failures = 0;

void check(bool ok, const std::string& label) {
    std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << label << '\n';
    if (!ok) ++g_failures;
}

}  // namespace

int main() {
    using namespace fin;

    try {
        std::cout << std::fixed << std::setprecision(6);

        // ------------------------------------------------------------------
        // 1. Volatility surface
        // ------------------------------------------------------------------
        std::cout << "1. Volatility surface\n";
        const std::vector<double> strikes = {80, 90, 100, 110, 120};
        const std::vector<double> mats    = {0.5, 1.0, 2.0};
        const std::vector<std::vector<double>> vols = {
            {0.26, 0.22, 0.20, 0.21, 0.24},
            {0.25, 0.21, 0.19, 0.20, 0.23},
            {0.24, 0.20, 0.18, 0.19, 0.22}};
        const VolSurface surf(strikes, mats, vols);

        const Market mkt{100.0, 0.05, 0.02, 0.20};
        const Option atm {OptionType::Call, ExerciseStyle::European, 100.0, 1.0};
        const Option wing{OptionType::Call, ExerciseStyle::European, 110.0, 1.0};

        check(std::fabs(surf.vol(100.0, 1.0) - 0.19) < 1e-12, "surface reproduces a grid node");
        check(std::fabs(surf.vol(95.0, 1.0) - 0.20) < 1e-12, "bilinear interpolation in strike");
        check(std::fabs(surf.vol(100.0, 0.75) - 0.195) < 1e-12, "bilinear interpolation in maturity");
        check(std::fabs(surf.vol(50.0, 5.0) - 0.24) < 1e-12, "flat extrapolation outside the grid");

        const double surf_price = surf.price(wing, mkt.spot, mkt.rate, mkt.dividend);
        const double surf_iv = implied_volatility(wing, mkt.spot, mkt.rate, mkt.dividend, surf_price);
        std::cout << "  surface price K=110,T=1: " << surf_price
                  << "   implied vol: " << surf_iv << '\n';
        check(std::fabs(surf_iv - surf.vol(110.0, 1.0)) < 1e-6, "surface price inverts to surface vol");

        const VolSurface flat({80.0, 120.0}, {0.5, 2.0}, {{0.2, 0.2}, {0.2, 0.2}});
        check(std::fabs(flat.price(atm, mkt.spot, mkt.rate, mkt.dividend) -
                        black_scholes_price(atm, mkt)) < 1e-12,
              "flat surface equals constant-vol Black-Scholes");

        // ------------------------------------------------------------------
        // 2. Yield curve and bonds
        // ------------------------------------------------------------------
        std::cout << "\n2. Yield curve and bonds\n";
        const YieldCurve curve({0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 30.0},
                               {0.035, 0.037, 0.039, 0.040, 0.042, 0.043, 0.044, 0.045});
        const BondSpec bond{100.0, 0.045, 10, 0.0, 1};
        const auto cfs = bond_cashflows(bond);

        const double p_curve = price_from_curve(cfs, curve);
        const double ytm = yield_to_maturity(cfs, p_curve);
        const double mod_dur = modified_duration(cfs, ytm);
        const double cvx = convexity(cfs, ytm);

        std::cout << "  10y 4.5% annual bond  price " << p_curve
                  << "  YTM " << ytm * 100.0 << "%\n";
        std::cout << "  modified duration " << mod_dur << "   convexity " << cvx << '\n';
        std::cout << "  1y-2y forward rate " << curve.forward_rate(1.0, 2.0) * 100.0 << "%\n";

        check(std::fabs(price_from_yield(cfs, ytm) - p_curve) < 1e-8, "YTM reprices the bond");

        const double h = 1e-6;
        const double mod_dur_fd =
            -(price_from_yield(cfs, ytm + h) - price_from_yield(cfs, ytm - h)) / (2.0 * h * p_curve);
        check(std::fabs(mod_dur_fd - mod_dur) < 1e-6, "modified duration matches numerical derivative");

        // ------------------------------------------------------------------
        // 3. Callable bond on a Ho-Lee lattice calibrated to the curve
        // ------------------------------------------------------------------
        std::cout << "\n3. Callable bond (Ho-Lee lattice calibrated to curve)\n";
        const auto lat_1pct = build_ho_lee(curve, 0.01, 10, 12);
        const auto lat_2pct = build_ho_lee(curve, 0.02, 10, 12);

        BondSpec callable = bond;
        callable.call_price = 100.0;
        callable.first_call_year = 3;

        const double straight_tree = lattice_bond_price(lat_1pct, bond);
        const double callable_tree = lattice_bond_price(lat_1pct, callable);
        const double callable_hi   = lattice_bond_price(lat_2pct, callable);

        std::cout << "  straight bond (lattice)      " << straight_tree << '\n';
        std::cout << "  callable, sigma 1%           " << callable_tree
                  << "   (issuer call value " << straight_tree - callable_tree << ")\n";
        std::cout << "  callable, sigma 2%           " << callable_hi << '\n';

        check(std::fabs(straight_tree - p_curve) < 1e-8,
              "lattice reprices straight bond to the curve");
        check(callable_tree <= straight_tree + 1e-10,
              "callable <= straight (issuer call is worth >= 0)");
        check(callable_hi <= callable_tree + 1e-10,
              "higher vol raises the issuer call value");

        // ------------------------------------------------------------------
        // 4. Barrier options (Monte Carlo)
        // ------------------------------------------------------------------
        std::cout << "\n4. Barrier options (Monte Carlo, 252 monitoring steps, 100k paths)\n";
        const std::size_t N = 100000, M = 252;
        const std::uint64_t seed = 7;
        const auto down90  = barrier_monte_carlo(atm, mkt, 90.0,  BarrierDirection::Down, N, M, seed);
        const auto up120   = barrier_monte_carlo(atm, mkt, 120.0, BarrierDirection::Up,   N, M, seed);
        const auto far1    = barrier_monte_carlo(atm, mkt, 1.0,   BarrierDirection::Down, N, M, seed);
        const auto down80  = barrier_monte_carlo(atm, mkt, 80.0,  BarrierDirection::Down, N, M, seed);
        const double bs_call = black_scholes_price(atm, mkt);

        std::cout << "  vanilla call (MC)       " << down90.vanilla.price
                  << " +/- " << down90.vanilla.std_error
                  << "   (Black-Scholes " << bs_call << ")\n";
        std::cout << "  down-and-out, H=90      " << down90.knock_out.price
                  << " +/- " << down90.knock_out.std_error << '\n';
        std::cout << "  down-and-in,  H=90      " << down90.knock_in.price
                  << " +/- " << down90.knock_in.std_error << '\n';
        std::cout << "  up-and-out,   H=120     " << up120.knock_out.price
                  << " +/- " << up120.knock_out.std_error << '\n';
        std::cout << "  up-and-in,    H=120     " << up120.knock_in.price
                  << " +/- " << up120.knock_in.std_error << '\n';

        check(std::fabs(down90.knock_in.price + down90.knock_out.price - down90.vanilla.price) < 1e-10,
              "down: knock-in + knock-out == vanilla (path by path)");
        check(std::fabs(up120.knock_in.price + up120.knock_out.price - up120.vanilla.price) < 1e-10,
              "up: knock-in + knock-out == vanilla (path by path)");
        check(std::fabs(far1.knock_out.price - bs_call) < 4.0 * far1.knock_out.std_error,
              "barrier far away converges to the vanilla price");
        check(down80.knock_out.price >= down90.knock_out.price - 1e-12,
              "lower barrier gives more knock-out value");
        check(down90.knock_out.price <= down90.vanilla.price + 1e-12,
              "knock-out <= vanilla");

        std::cout << '\n';
        if (g_failures == 0) {
            std::cout << "All extension checks passed.\n";
            return 0;
        }
        std::cout << g_failures << " check(s) FAILED.\n";
        return 1;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 2;
    }
}
