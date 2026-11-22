# Equity Derivatives Pricer (C++17)

A self-contained derivatives pricing library and test harness. The core covers European and American equity options: analytic pricing with Greeks, a binomial lattice, Monte Carlo with variance reduction, and implied volatility inversion. Three extension modules add a volatility surface, a fixed-income layer (zero curve, bond analytics, and callable bonds on a Ho-Lee lattice), and barrier options by Monte Carlo. Every method is checked against the others in built-in self-tests, and a plotting pipeline renders the results.

## Contents

| File | Purpose |
|---|---|
| `derivatives_pricer.hpp` | Core library: Black-Scholes, Greeks, binomial tree, Monte Carlo, implied volatility |
| `derivatives_pricer.cpp` | Demo and self-test for the core library |
| `vol_surface.hpp` | Implied volatility surface with bilinear interpolation and flat extrapolation |
| `yield_curve.hpp` | Zero curve, bond analytics (YTM, duration, convexity), Ho-Lee lattice for callable bonds |
| `barrier_mc.hpp` | Knock-in and knock-out barrier options by Monte Carlo with a continuity correction |
| `extensions_test.cpp` | Demo and self-tests for the three extension modules |
| `make_data.cpp` | Writes CSV results from the core and extension modules to a `data/` folder |
| `plot_results.py` | Renders the CSV results as PNG plots (requires `numpy` and `matplotlib`) |
| `README.md` | This document |

All C++ code is header-based C++17 with no dependencies beyond the standard library. Only the plotting script needs Python packages.

## Build and run

```
g++ -std=c++17 -O2 -Wall -Wextra -pedantic derivatives_pricer.cpp -o pricer && ./pricer
g++ -std=c++17 -O2 -Wall -Wextra -pedantic extensions_test.cpp -o ext_test && ./ext_test
echo $?    # 0 = all self-tests passed, 1 = a check failed, 2 = invalid input
```

Requirements: any C++17 compiler (GCC 7+, Clang 5+, MSVC 2017+).

## Generate data and plots

```
g++ -std=c++17 -O2 make_data.cpp -o make_data && ./make_data data   # about 11 s
python3 plot_results.py data plots                                  # needs numpy, matplotlib
```

`make_data` writes these CSVs into `data/`: `core_bs.csv`, `vol_smile.csv`, `vol_surface.csv`, `curve.csv`, `forward.csv`, `callable_vs_vol.csv`, `barrier_down.csv`, and `barrier_up.csv`. `plot_results.py` turns them into the six figures in `plots/`, listed under Results.

## Model and conventions

- **Underlying**: geometric Brownian motion with constant volatility, continuous dividend yield `q`, and constant continuously compounded risk-free rate `r`. The volatility surface is the one exception: it supplies a volatility that varies with strike and maturity.
- **Time**: maturity `T` is in years. Greeks are returned in raw units; the demo rescales them for display.
- **Sign conventions**: theta is per year (negative for long options as time decays); rho and vega are per unit change in rate and volatility.
- **Rates**: the fixed-income module uses continuously compounded zero rates and annual-compounding bond yields. Rate volatility is absolute (normal), so `0.01` means 100 bp per square-root year.

## Public API

All names live in namespace `fin`.

### Core types

- `OptionType { Call, Put }`
- `ExerciseStyle { European, American }`
- `Option { type, style, strike, maturity }`
- `Market { spot, rate, dividend, volatility }`
- `Greeks { delta, gamma, vega, theta, rho }`
- `PriceResult { price, std_error }`

### Core functions

| Function | Returns | Notes |
|---|---|---|
| `black_scholes_price(opt, mkt)` | `double` | European only; throws for American |
| `black_scholes_greeks(opt, mkt)` | `Greeks` | European only |
| `binomial_price(opt, mkt, steps)` | `double` | European or American; throws if the risk-neutral probability leaves (0,1) |
| `monte_carlo_price(opt, mkt, pairs, seed)` | `PriceResult` | European; antithetic pairs; Welford mean and variance |
| `implied_volatility(opt, spot, rate, div, price)` | `double` | Bisection on [1e-6, 5]; throws if the price is outside the attainable range |

All inputs are validated. An invalid strike, maturity, spot, or volatility throws `std::invalid_argument`.

### Volatility surface (`vol_surface.hpp`)

| Member | Notes |
|---|---|
| `VolSurface(strikes, maturities, vols)` | `vols[t][k]`; both grids strictly increasing, at least 2 points each |
| `vol(strike, maturity)` | Bilinear interpolation; held flat outside the grid |
| `price(opt, spot, rate, div)` | Black-Scholes price at the surface vol for the option's strike and maturity |

### Fixed income (`yield_curve.hpp`)

| Type or function | Notes |
|---|---|
| `YieldCurve(tenors, zero_rates)` | Linear interpolation on continuously compounded zero rates; flat extrapolation |
| `curve.discount(t)`, `curve.zero_rate(t)`, `curve.forward_rate(t1, t2)` | Discount factor, zero rate, and continuously compounded forward |
| `BondSpec { face, coupon_rate, maturity_years, call_price, first_call_year }` | Annual coupons. `call_price <= 0` means non-callable |
| `bond_cashflows(bond)` | Returns `std::vector<CashFlow>` |
| `price_from_curve(cfs, curve)` | Present value off the curve |
| `price_from_yield(cfs, y)`, `yield_to_maturity(cfs, price)` | Annual-compounding yield and its inverse |
| `modified_duration(cfs, y)`, `convexity(cfs, y)` | Risk measures at a given yield |
| `build_ho_lee(curve, sigma, years, steps_per_year)` | Ho-Lee lattice calibrated exactly to the curve's discount factors |
| `lattice_bond_price(lattice, bond)` | Prices straight or callable bonds on the lattice; the lattice horizon must equal the bond maturity |

### Barrier options (`barrier_mc.hpp`)

| Function | Notes |
|---|---|
| `barrier_monte_carlo(opt, mkt, barrier, dir, paths, steps, seed)` | European options only; `dir` is `BarrierDirection::Up` or `Down`. Returns `BarrierResult { knock_in, knock_out, vanilla }`, each a `PriceResult` |

Up barriers must lie above spot and down barriers below spot. Knock-in and knock-out are computed on the same paths, so their sum equals the vanilla price exactly.

## Method details

**Black-Scholes-Merton.** With `d1 = [ln(S/K) + (r - q + σ²/2)T] / (σ√T)` and `d2 = d1 - σ√T`:

- Call = `S e^{-qT} N(d1) - K e^{-rT} N(d2)`
- Put = `K e^{-rT} N(-d2) - S e^{-qT} N(-d1)`

**Greeks.** Delta is `e^{-qT} N(d1)` for calls and `e^{-qT}(N(d1) - 1)` for puts. Gamma and vega share `φ(d1)`. Theta and rho follow the standard closed forms.

**Binomial tree (CRR).** `u = e^{σ√Δt}`, `d = 1/u`, `p = (e^{(r-q)Δt} - d)/(u - d)`. The tree is built from terminal payoffs by backward induction. For American options, each node takes `max(continuation, intrinsic)`.

**Monte Carlo.** Terminal price `S exp((r - q - σ²/2)T ± σ√T Z)`, with each draw paired with its antithetic mirror to reduce variance. The estimator's standard error is computed online with Welford's algorithm.

**Implied volatility.** The Black-Scholes price is strictly increasing in volatility for European options, so bisection is guaranteed to converge.

**Volatility surface.** Bilinear interpolation in (strike, maturity) on the implied-vol grid. The interpolation is linear in strike and does not enforce no-arbitrage conditions.

**Ho-Lee lattice.** Short rate `r(i, j) = α_i + σ√Δt (2j − i)` on a recombining tree with risk-neutral up-probability 1/2. Each `α_i` is solved by bisection so that the lattice reproduces the curve's discount factor at the next time step. Arrow-Debreu prices propagate forward, and bonds are valued by backward induction. A callable bond's value is `min(continuation, call price)` at each permitted coupon date, plus the coupon.

**Barrier Monte Carlo.** Log-price paths with `steps` monitoring dates. The barrier is shifted by the Broadie-Glasserman-Kou correction, `β = 0.5826`, times `σ√Δt`, which approximates continuous monitoring. Every path draws all of its random numbers even after a hit, so results for different barrier levels are comparable path by path.

## Self-test suites

**`derivatives_pricer.cpp`** runs six checks:

1. Put-call parity: `C − P = S e^{−qT} − K e^{−rT}` to 1e-10.
2. Binomial with 2,000 steps agrees with Black-Scholes to 0.01.
3. Monte Carlo lies within 4 standard errors of the analytic call price.
4. Analytic delta agrees with a central finite difference to 1e-6.
5. Implied volatility recovers the input 20% to 1e-6.
6. American put price is at least the European put price.

**`extensions_test.cpp`** runs 17 checks:

- Volatility surface: grid nodes are reproduced, bilinear interpolation in strike and in maturity is exact on the test points, flat extrapolation holds, the surface price inverts to the surface vol, and a flat surface matches constant-vol Black-Scholes.
- Fixed income: the YTM reprices the bond, and modified duration matches a numerical derivative. The Ho-Lee lattice reprices the straight bond to the curve, the callable bond is never worth more than the straight bond, and a higher rate volatility raises the issuer's call value.
- Barrier options: knock-in plus knock-out equals vanilla on both the down and up sides, a far-away barrier converges to the vanilla price, a lower down barrier gives more knock-out value, and knock-out never exceeds vanilla.

## Reference output

### Core (S = 100, r = 5%, q = 2%, σ = 20%, T = 1 year)

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

### Extensions

| Quantity | Value |
|---|---|
| Surface price, call K = 110, T = 1 | 5.188582 (implied vol 0.200000) |
| 10-year 4.5% bond, curve price | 100.295587 |
| Yield to maturity | 4.462713% |
| Modified duration / convexity | 7.918386 / 77.902566 |
| 1y–2y forward rate | 4.100000% |
| Callable bond, σ = 1% (lattice) | 96.182788 (issuer call value 4.112799) |
| Callable bond, σ = 2% (lattice) | 91.723067 |
| Vanilla call, K = 100 (MC, 100k paths, 252 steps) | 9.208327 ± 0.043693 (BS 9.227006) |
| Down-and-out call, H = 90 | 7.993708 ± 0.043414 |
| Down-and-in call, H = 90 | 1.214618 ± 0.014780 |
| Up-and-out call, H = 120 | 1.445397 ± 0.011358 |
| Up-and-in call, H = 120 | 7.762930 ± 0.044771 |

The extension values use the same inputs as `extensions_test.cpp`: the market above and the curve with pillars at 0.5, 1, 2, 3, 5, 7, 10, and 30 years.

## Results

Running `make_data` and `plot_results.py` produces these figures in `plots/`:

1. `01_core_black_scholes.png`: call and put price, and call delta, against spot.
2. `02_vol_smile.png`: the implied volatility smile at maturities of 0.5, 1, and 2 years.
3. `03_vol_surface_heatmap.png`: the full volatility surface over strike and maturity.
4. `04_yield_curve.png`: zero rates, discount factors, and forward rates.
5. `05_callable_bond_vs_vol.png`: straight and callable bond prices against rate volatility, with the issuer call value.
6. `06_barrier_prices.png`: knock-out and knock-in prices against barrier level for down and up barriers, with the vanilla price and Monte Carlo error bars.

## Limitations

- The core pricer assumes constant volatility, rate, and dividend yield. Only the volatility surface module models a smile or term structure of volatility.
- Greeks are analytic for European options only. American Greeks are not provided.
- Monte Carlo covers European payoffs only in the core. Barrier options are European and single-barrier.
- The Ho-Lee model allows negative rates and uses one normal volatility. It does not fit a volatility skew or mean reversion.
- The volatility surface interpolates linearly and does not enforce no-arbitrage conditions. A production surface needs a parameterisation such as SVI and calendar and butterfly checks.
- Barrier prices are Monte Carlo estimates with a discrete-monitoring correction. For a production check, compare them with the Reiner-Rubinstein closed form.
- Callable bonds assume the issuer calls whenever calling is optimal, with no call costs or make-whole provisions.
- Implied volatility uses bisection, which is robust but slower than Newton's method.
- The binomial tree is O(N²) in time and uses `std::pow` per node, so very large step counts are slow.

## Possible extensions

- Asian and lookback payoffs in Monte Carlo.
- Local or stochastic volatility models (Heston, SABR).
- Newton-Raphson implied volatility with a bisection fallback.
- Hull-White or Black-Derman-Toy short-rate models for fixed income.
- Python bindings via pybind11.
