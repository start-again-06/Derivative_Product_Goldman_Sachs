#include "derivatives_pricer.hpp"
#include <iomanip>
#include <iostream>
#include <string>

// ===========================================================================
// Demo and self-test
// ===========================================================================

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
        const Market mkt{100.0, 0.05, 0.02, 0.20};
        const Option call{OptionType::Call, ExerciseStyle::European, 105.0, 1.0};
        const Option put {OptionType::Put,  ExerciseStyle::European, 105.0, 1.0};
        const Option am_put_100{OptionType::Put, ExerciseStyle::American, 100.0, 1.0};
        const Option eu_put_100{OptionType::Put, ExerciseStyle::European, 100.0, 1.0};

        std::cout << std::fixed << std::setprecision(6);

        // 1. Analytic prices and Greeks ------------------------------------
        const double c = black_scholes_price(call, mkt);
        const double p = black_scholes_price(put, mkt);
        const Greeks gc = black_scholes_greeks(call, mkt);

        std::cout << "Market: S=" << mkt.spot << " r=" << mkt.rate << " q=" << mkt.dividend
                  << " sigma=" << mkt.volatility << "   Option: K=105, T=1y\n\n";
        std::cout << "Black-Scholes call price : " << c << '\n';
        std::cout << "Black-Scholes put price  : " << p << '\n';
        std::cout << "Call delta / gamma       : " << gc.delta << " / " << gc.gamma << '\n';
        std::cout << "Call vega (per vol pt)   : " << gc.vega / 100.0 << '\n';
        std::cout << "Call theta (per day)     : " << gc.theta / 365.0 << '\n';
        std::cout << "Call rho (per rate pt)   : " << gc.rho / 100.0 << "\n\n";

        // 2. Self-tests ----------------------------------------------------
        std::cout << "Self-tests:\n";

        // Put-call parity: C - P = S e^{-qT} - K e^{-rT}
        const double parity_lhs = c - p;
        const double parity_rhs = mkt.spot * std::exp(-mkt.dividend * call.maturity) -
                                  call.strike * std::exp(-mkt.rate * call.maturity);
        check(std::fabs(parity_lhs - parity_rhs) < 1e-10, "put-call parity holds");

        // Binomial tree converges to Black-Scholes for European options.
        const double tree = binomial_price(call, mkt, 2000);
        check(std::fabs(tree - c) < 0.01, "binomial(2000 steps) matches Black-Scholes");

        // Monte Carlo lies within 4 standard errors of the analytic price.
        const PriceResult mc = monte_carlo_price(call, mkt, 250000, 42);
        check(std::fabs(mc.price - c) < 4.0 * mc.std_error, "Monte Carlo within 4 std errors");

        // Analytic delta matches a central finite difference.
        const double h = 1e-4;
        const double delta_fd =
            (black_scholes_price(call, {mkt.spot + h, mkt.rate, mkt.dividend, mkt.volatility}) -
             black_scholes_price(call, {mkt.spot - h, mkt.rate, mkt.dividend, mkt.volatility})) /
            (2.0 * h);
        check(std::fabs(delta_fd - gc.delta) < 1e-6, "analytic delta matches finite difference");

        // Implied volatility round-trips the input volatility.
        const double iv = implied_volatility(call, mkt.spot, mkt.rate, mkt.dividend, c);
        check(std::fabs(iv - mkt.volatility) < 1e-6, "implied vol recovers 20% exactly");

        // American put is worth at least its European counterpart (early exercise premium).
        const double am = binomial_price(am_put_100, mkt, 2000);
        const double eu = black_scholes_price(eu_put_100, mkt);
        check(am >= eu - 1e-9, "American put >= European put");

        std::cout << "\nAmerican vs European put (K=100, T=1y):\n";
        std::cout << "  American (binomial) : " << am << '\n';
        std::cout << "  European (BS)       : " << eu << '\n';
        std::cout << "  Early-exercise premium: " << (am - eu) << "\n\n";

        std::cout << "Monte Carlo call (250k antithetic pairs): " << mc.price
                  << " +/- " << mc.std_error << "\n";
        std::cout << "Implied vol from BS call price: " << iv << "\n\n";

        if (g_failures == 0) {
            std::cout << "All checks passed.\n";
            return 0;
        }
        std::cout << g_failures << " check(s) FAILED.\n";
        return 1;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 2;
    }
}
