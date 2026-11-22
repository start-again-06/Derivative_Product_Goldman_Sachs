# Equity Derivatives Pricer (C++17)

A self-contained options pricing library and test harness for European and American equity options. It implements analytic pricing with Greeks, a binomial lattice, Monte Carlo simulation with variance reduction, and implied volatility inversion, and checks the methods against each other in a built-in self-test.

## Contents

| File | Purpose |
|---|---|
| `derivatives_pricer.hpp` | Core library: Black-Scholes, Greeks, binomial tree, Monte Carlo, implied volatility |
| `derivatives_pricer.cpp` | Demo and self-test for the core library |
| `vol_surface.hpp` | Implied volatility surface with bilinear interpolation |
| `yield_curve.hpp` | Zero curve, bond analytics (YTM, duration, convexity), Ho-Lee lattice for callable bonds |
| `barrier_mc.hpp` | Knock-in / knock-out barrier options by Monte Carlo |
| `extensions_test.cpp` | Demo and self-tests for the three extension modules |
| `README.md` | This document |

All code is header-based C++17 with no dependencies beyond the standard library.

## Build and run

```
g++ -std=c++17 -O2 -Wall -Wextra -pedantic derivatives_pricer.cpp -o pricer && ./pricer
g++ -std=c++17 -O2 -Wall -Wextra -pedantic extensions_test.cpp -o ext_test && ./ext_test
echo $?    # 0 = all self-tests passed, 1 = a check failed, 2 = invalid input
```

Requirements: any C++17 compiler (GCC 7+, Clang 5+, MSVC 2017+). No external libraries.

## Model and conventions

- **Underlying**: geometric Brownian motion with constant volatility, continuous dividend yield `q`, and constant continuously compounded risk-free rate `r`.
- **Time**: maturity `T` is in years. Greeks are returned in raw units; the demo rescales them for display.
- **Sign conventions**: theta is per year (negative for long options as time decays); rho and vega are per unit change in rate and volatility.

## Public API (namespace `fin`)

### Types

- `OptionType { Call, Put }`
- `ExerciseStyle { European, American }`
- `Option { type, style, strike, maturity }`
- `Market { spot, rate, dividend, volatility }`
- `Greeks { delta, gamma, vega, theta, rho }`
- `PriceResult { price, std_error }`

### Functions

| Function | Returns | Notes |
|---|---|---|
| `black_scholes_price(opt, mkt)` | `double` | European only; throws for American |
| `black_scholes_greeks(opt, mkt)` | `Greeks` | European only |
| `binomial_price(opt, mkt, steps)` | `double` | European or American; throws if the risk-neutral probability leaves (0,1) |
| `monte_carlo_price(opt, mkt, pairs, seed)` | `PriceResult` | European; antithetic pairs; Welford mean and variance |
| `implied_volatility(opt, spot, rate, div, price)` | `double` | Bisection on [1e-6, 5]; throws if the price is outside the attainable range |

All inputs are validated. Invalid strike, maturity, spot, or volatility throws `std::invalid_argument`.

## Method details

**Black-Scholes-Merton.** With `d1 = [ln(S/K) + (r - q + σ²/2)T] / (σ√T)` and `d2 = d1 - σ√T`:

- Call = `S e^{-qT} N(d1) - K e^{-rT} N(d2)`
- Put = `K e^{-rT} N(-d2) - S e^{-qT} N(-d1)`

**Greeks.** Delta is `e^{-qT} N(d1)` for calls and `e^{-qT}(N(d1) - 1)` for puts. Gamma and vega share `φ(d1)`. Theta and rho follow the standard closed forms.

**Binomial tree (CRR).** `u = e^{σ√Δt}`, `d = 1/u`, `p = (e^{(r-q)Δt} - d)/(u - d)`. The tree is built from terminal payoffs by backward induction. For American options, each node takes `max(continuation, intrinsic)`.

**Monte Carlo.** Terminal price `S exp((r - q - σ²/2)T ± σ√T Z)`, with each draw paired with its antithetic mirror to reduce variance. The estimator's standard error is computed online with Welford's algorithm.

**Implied volatility.** Black-Scholes price is strictly increasing in volatility for European options, so bisection is guaranteed to converge.

## Self-test suite

The demo in `main()` runs six checks and prints PASS or FAIL for each:

1. Put-call parity: `C - P = S e^{-qT} - K e^{-rT}` to 1e-10.
2. Binomial with 2,000 steps agrees with Black-Scholes to 0.01.
3. Monte Carlo lies within 4 standard errors of the analytic call price.
4. Analytic delta agrees with a central finite difference to 1e-6.
5. Implied volatility recovers the input 20% to 1e-6.
6. American put price is at least the European put price.

## Reference output

Inputs: S = 100, r = 5%, q = 2%, σ = 20%, T = 1 year.

| Quantity | Value |
|---|---|
| BS call, K = 105 | 6.986920 |
| BS put, K = 105 | 8.846142 |
| Call delta / gamma | 0.492465 / 0.019552 |
| Call vega (per vol point) | 0.391036 |
| Call theta (per day) | -0.013804 |
| Call rho (per rate point) | 0.422596 |
| American put, K = 100 (binomial, 2,000 steps) | 6.660226 |
| European put, K = 100 (BS) | 6.330081 |
| Early-exercise premium | 0.330145 |
| MC call, K = 105 (250k pairs) | 6.995484 ± 0.014336 |
| Implied vol from BS call | 0.200000 |

## Limitations

- Constant volatility, rate, and dividend yield only. No volatility smile or term structure.
- Greeks are analytic for European options only. American Greeks are not provided.
- Monte Carlo covers European payoffs only.
- Implied volatility uses bisection, which is robust but slower than Newton's method.
- The binomial tree is O(N²) in time and uses `std::pow` per node, so very large step counts are slow.

## Possible extensions

- Barrier, Asian, and lookback payoffs in Monte Carlo.
- Local or stochastic volatility models (Heston, SABR).
- A volatility surface with interpolation across strikes and maturities.
- Newton-Raphson implied volatility with a bisection fallback.
- Python bindings via pybind11.
