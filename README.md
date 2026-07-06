# Summative: Mortgage Valuation Model

## (a) Constructing Interest Rate Trees

The University Building Society (UBS) seeks to price a 10-year fixed rate mortgage product for UK university students and employees. Doing so requires a model of how interest rates may evolve over the life of the mortgage. This is essential for two reasons.

First, the fair mortgage rate must reflect the time value of money across the full distribution of possible future interest rate outcomes, rather than being based solely on the current term structure. Second, the mortgage contains an embedded prepayment option: if market rates decline and refinancing becomes optimal, the borrower may repay early. Valuing this option requires a forward-looking model of short rate dynamics.

To capture this uncertainty, two binomial interest rate tree models are constructed: the Ho-Lee (HL) model and the Black-Derman-Toy (BDT) model. In both cases, the short rate, defined as the continuously compounded interest rate applying over a single period, evolves through up or down movements at each time step under the risk-neutral measure. By calibrating the trees to observed market prices of government bonds, the models are forced to fit the current term structure exactly, thereby providing an internally consistent basis for pricing mortgage cash flows and the optionality of the contract.

Both models are implemented using the following assumptions:

- Semi-annual time steps: $\Delta = 0.5$ years
- Continuously compounded short rates
- A 10-year horizon: $T = 10$ years
- Number of steps: $N = 20$
- Risk-neutral probabilities: $q^* = 0.5$, up and down moves are equally likely

These assumptions balance tractability and realism. Semi-annual discretisation is sufficiently fine to capture the medium-term evolution of rates while remaining computationally efficient for tree-based valuation.

The Bank of England publishes daily estimates of the nominal government liability spot curve, derived from the observed UK Gilt prices using the Svensson parametric model. As it is derived from actively traded government securities, the curve reflects real market pricing and incorporates current expectations of inflation, monetary policy and risk premia. While the curve is model based rather than directly observable at all maturities and may therefore introduce some smoothing assumptions, it remains the standard benchmark for risk-free discounting in GBP and is appropriate for calibration purposes in interest rate modelling.

Spot rates are extracted at 20 semi-annual maturities from 0.5 to 10 years. These rates are already quoted as continuously compounded annual rates, which is the convention required by both the HL and BDT models.

The full set of spot rates on 30th January 2026 is shown below:

| Maturity (Yrs) | Spot Rate (%) | Maturity (Yrs) | Spot Rate (%) |
|:---:|:---:|:---:|:---:|
| 0.5 | 3.48 | 5.5 | 4.04 |
| 1.0 | 3.55 | 6.0 | 4.11 |
| 1.5 | 3.59 | 6.5 | 4.18 |
| 2.0 | 3.63 | 7.0 | 4.25 |
| 2.5 | 3.67 | 7.5 | 4.31 |
| 3.0 | 3.72 | 8.0 | 4.38 |
| 3.5 | 3.78 | 8.5 | 4.45 |
| 4.0 | 3.84 | 9.0 | 4.51 |
| 4.5 | 3.90 | 9.5 | 4.57 |
| 5.0 | 3.97 | 10.0 | 4.63 |

The yield curve is upward sloping throughout, the spread between the 6-month and 10-year spot rates is approximately 115 basis points, rising from 3.48% to 4.63%. This positive slope implies that the market prices in higher future rates over longer horizons, a feature that the calibrated trees must be able to reproduce.

From each spot rate, the corresponding zero-coupon bond price per £100 nominal is calculated as:

$$P(k) = 100 \times e^{-r(k)T(k)}$$

Where $r(k)$ is the continuously compounded spot rate and $T(k) = k\Delta$ is the maturity in years.

These 20 bond prices form the calibration targets for the interest rate trees. Exact matching is crucial, as if the model cannot replicate the current term structure, it cannot be relied upon for pricing related interest rate products consistently.

Both models require a volatility parameter, $\sigma$, which controls the scale of short-rate movements at each step in the tree. The choice of volatility is a central modelling input as it directly affects the dispersion of future rate paths and, therefore, the value of the embedded prepayment option. A lower volatility would understate refinancing risk, while a higher volatility would exaggerate it.

Rather than estimating volatility from historical data, $\sigma$ is taken from the GBP swaption implied volatility surface. This choice is preferable for two reasons. First, implied volatility is forward-looking, it summarises the market's current assessment of future rate uncertainty, which is the required quantity for pricing a new mortgage contract. Second, historical volatility estimates based on data contained between 2009-2021 are heavily influenced by the prolonged near-zero rate environment post the Great Financial Crisis and would likely understate the uncertainty embedded in the post-tightening environment.

A swaption is an option to enter an interest rate swap at a predetermined fixed rate. Its market price depends on expected future rate variability. By inverting the pricing model used to value the swaption, similar to the process to extract implied volatility from equity options using Black-Scholes, we obtain the volatility parameter consistent with the observed market price.

We obtained the sigma values from Bloomberg GBP RFR BVOL Cube (Default) via SWPM, using the GBP OIS Swaption Volatility Surface on 30th January 2026. A 6-month expiry is chosen to align with the semi-annual tree step ($\Delta = 0.5$), ensuring consistency between the model's discretisation and the volatility input. The 10-year tenor matches the maturity of the mortgage. Using the same market date as the yield curve also ensures that the term structure and volatility assumptions are drawn from a common market snapshot.

Because the HL and BDT models use different rate dynamics, normal and lognormal respectively, they require different volatility values.

### Ho-Lee sigma (normal volatility)

The Ho-Lee model is based on additive, Gaussian shocks to the short rate. The relevant volatility measure is therefore normal volatility, quoted in basis points per year, which captures the absolute variability of rates.

The Bloomberg normal volatility reading: 68.02 bp/year

$$\sigma_{HL} = \frac{68.02}{10{,}000} = 0.006802 \ (0.6802\%)$$

### BDT Sigma (Black/lognormal volatility)

The BDT model is specified in log-rate, meaning shocks are multiplicative in levels. The appropriate volatility measures are therefore Black volatility, which captures proportional rate uncertainty.

The Bloomberg Black volatility reading: 16.23%

$$\sigma_{BDT} = 0.1623 \ (16.23\%)$$

One limitation is that only a single at-the-money volatility is used from one point on the swaption surface. In reality, implied volatility varies across both strikes and expiries. As a result, the calibration does not capture volatility smile or term-structure effects. Nonetheless, this simplification is appropriate for one-factor binomial tree models and gives a forward-looking rate of volatility compared to that of historical estimates.

### Model Descriptions

The Ho-Lee model assumes that the short rate evolves under the risk-neutral measure according to:

$$r(i+1,j) = r(i,j) + \theta_i \Delta + \sigma_{HL}\sqrt{\Delta}$$

$$r(i+1,j+1) = r(i,j) + \theta_i \Delta - \sigma_{HL}\sqrt{\Delta}$$

Where $i$ denotes the time step and $j$ indexes the node within that step, with $j = 0$ corresponding to the highest-rate node.

The key features of the Ho-Lee model are that volatility is constant across both time and states, implying that uncertainty does not depend on the level of interest rates. Secondly, the drift term $\theta_i$ is allowed to vary through time and is calibrated sequentially to fit the observed yield curve. The model also implies a normal distribution of future rates, which introduces symmetry but also permits negative outcomes. Finally, the tree recombines, meaning that an up move followed by a down move reaches the same node as a down move followed by an up move, keeping the tree computationally manageable.

The main weakness with the Ho-Lee model is that the normality assumption becomes increasingly unrealistic over long horizons, particularly when rates approach the lower nodes. Negative interest rates are highly unlikely in reality.

The BDT model is specified in log-rate space, with up and down nodes represented by the following equations respectively:

$$z(i+1,j) = z(i,j) + \theta_i \Delta + \sigma_{BDT}\sqrt{\Delta}$$

$$z(i+1,j+1) = z(i,j) + \theta_i \Delta - \sigma_{BDT}\sqrt{\Delta}$$

The BDT model therefore guarantees strictly positive rates at every node, which is a key practical advantage over the Ho-Lee model. In addition, the implied distribution of future rates is lognormal rather than normal, producing a right-skewed distribution with a heavy upper tail. This means that high-rate scenarios become more extreme as the horizon lengthens, while rates remain bounded away from large negative values.

The model's strength is that it is economically more realistic treatment of positivity and proportional rate movements. However, the lognormal structure of the tree can also generate very large upward rate scenarios over long horizons, which may appear extreme even if they arise with low probability under the risk-neutral measure.

Both models are calibrated using the same general algorithm, adapted appropriately to normal dynamics in Ho-Lee and lognormal dynamics in BDT.

The initial short rate is set equal to the 6-month spot rate implied by the first zero-coupon bond price.

$$r(0,0) = \frac{\log(P(1)/100)}{\Delta} = 3.48\%$$

This ensures that the tree begins at the observed current short rate.

Let $q(i,j)$ denote the Arrow-Debreu state price, defined as the risk-neutral present value today of receiving £1 if node $(i,j)$ is reached. The initial state price is:

$$q(0,0) = 1$$

How the state prices evolve according to:

$$q(i+1,j) = 0.5\,q(i,j)e^{-r(i,j)\Delta} + 0.5\,q(i,j-1)e^{-r(i,j-1)\Delta}$$

The model price of a zero-coupon bond maturing at step $k$ is then given by:

$$P(k) = 100 \sum_j q(k,j)e^{-r(k,j)\Delta}$$

For each time step, the drift parameter $\theta_i$ is chosen such that the model exactly matches the market price of the bond maturing at step $i+2$. Once $\theta_i$ is obtained, the rate tree is updated for step $i+1$, the corresponding state prices are computed, and the next bond price can be matched.

The final column at step $N$ is filled using the drift from step $N-1$. Although this means the terminal step is not calibrated independently in the same way as earlier steps, the calibration checks confirm that the resulting bond prices still match the market data exactly.

### Results

```
Table: The Risk Neutral Ho-Lee Interest Rate Tree
BoE Nominal Spot Curve, 30 Jan 2026 | sigma = 0.6802% | delta = 0.5 | N = 20
================================================================================================================================
 Time T   0.0    0.5    1.0    1.5    2.0    2.5    3.0    3.5    4.0    4.5    5.0    5.5    6.0    6.5    7.0    7.5    8.0    8.5    9.0    9.5   10.0
Period i    0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20
θi(×100)         0.2812 0.1035 0.1658 0.1681 0.2904 0.3527 0.2550 0.2573 0.4597 0.3020 0.3043 0.3066 0.3089 0.0112 0.5935 0.3158 -0.0418 0.2805 0.2828

 j=0    3.48   4.10   4.63   5.20   5.76   6.39   7.05   7.66   8.26   8.98   9.61  10.24  10.87  11.51  12.00  12.77  13.41  13.87  14.49  15.12  15.74
 j=1           3.14   3.67   4.24   4.80   5.43   6.08   6.69   7.30   8.01   8.65   9.28   9.91  10.55  11.04  11.81  12.45  12.91  13.53  14.16  14.78
 j=2                  2.71   3.27   3.84   4.47   5.12   5.73   6.34   7.05   7.68   8.32   8.95   9.59  10.07  10.85  11.49  11.95  12.57  13.19  13.82
 j=3                         2.31   2.88   3.50   4.16   4.77   5.38   6.09   6.72   7.35   7.99   8.62   9.11   9.89  10.53  10.99  11.61  12.23  12.85
 j=4                                1.92   2.54   3.20   3.81   4.42   5.13   5.76   6.39   7.03   7.66   8.15   8.93   9.57  10.03  10.65  11.27  11.89
 j=5                                       1.58   2.24   2.85   3.46   4.17   4.80   5.43   6.07   6.70   7.19   7.96   8.60   9.06   9.69  10.31  10.93
 j=6                                              1.27   1.88   2.49   3.20   3.84   4.47   5.10   5.74   6.23   7.00   7.64   8.10   8.72   9.35   9.97
 j=7                                                     0.92   1.53   2.24   2.87   3.51   4.14   4.78   5.26   6.04   6.68   7.14   7.76   8.38   9.01
 j=8                                                            0.57   1.28   1.91   2.55   3.18   3.81   4.30   5.08   5.72   6.18   6.80   7.42   8.04
 j=9                                                                   0.32   0.95   1.58   2.22   2.85   3.34   4.12   4.76   5.22   5.84   6.46   7.08
 j=10                                                                        -0.01   0.62   1.26   1.89   2.38   3.16   3.79   4.25   4.88   5.50   6.12
 j=11                                                                              -0.34  -0.03   0.93   1.42   2.19   2.83   3.29   3.91   4.54   5.16
 j=12                                                                                     -0.67  -0.03   0.45   1.23   1.87   2.33   2.95   3.57   4.20
 j=13                                                                                            -0.99  -0.51   0.27   0.91   1.37   1.99   2.61   3.23
 j=14                                                                                                   -1.47  -0.69  -0.05   0.41   1.03   1.65   2.27
 j=15                                                                                                          -1.65  -1.02  -0.56   0.07   0.69   1.31
 j=16                                                                                                                 -1.98  -1.52  -0.90  -0.27   0.35
 j=17                                                                                                                        -2.48  -1.86  -1.24  -0.61
 j=18                                                                                                                               -2.82  -2.20  -1.58
 j=19                                                                                                                                      -3.16  -2.54
 j=20                                                                                                                                             -3.50
```

The calibrated drift parameters are predominately positive, which is consistent with the upward-sloping yield curve. The model therefore requires a generally positive drift adjustment to reconcile current short rates with higher long-term discount rates. The small negative value at step 17 reflects a small flattening in the curve rather than a reversal in its overall shape. More broadly, the variation in $\theta_i$ across time shows the subtle curvature of the term structure.

The spread between upper and lower nodes widens approximately linearly over time, reflecting the constant normal volatility assumption. By step 20, the spread equals 19.24 percentage points. Negative rates emerge from step 10 onwards in the lower part of the tree, highlighting the main weakness of the Ho-Lee model for long dated mortgage valuation.

```
Table: The Risk Neutral BDT Interest Rate Tree
BoE Nominal Spot Curve, 30 Jan 2026 | sigma = 16.2300% | delta = 0.5 | N = 20
================================================================================================================================
 Time T   0.0    0.5    1.0    1.5    2.0    2.5    3.0    3.5    4.0    4.5    5.0    5.5    6.0    6.5    7.0    7.5    8.0    8.5    9.0    9.5   10.0
Period i    0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20
θi(×100)         6.5978 1.5012 3.1213 3.0840 6.1011 7.3721 4.7700 4.6826 9.0097 5.2918 5.2048 5.1321 5.0727 -0.7308 10.3494 4.9743 -1.4675 4.3584 4.3800

 j=0    3.48   4.03   4.56   5.19   5.92   6.84   7.96   9.14  10.50  12.32  14.19  16.33  18.79  21.62  24.16  28.54  32.82  36.54  41.88  48.02  55.05
 j=1           3.21   3.62   4.13   4.70   5.44   6.33   7.27   8.35   9.79  11.28  12.98  14.94  17.19  19.21  22.69  26.09  29.04  33.29  38.17  43.76
 j=2                  2.88   3.28   3.74   4.32   5.03   5.78   6.63   7.78   8.96  10.32  11.88  13.66  15.27  18.03  20.74  23.09  26.47  30.34  34.78
 j=3                         2.61   2.97   3.44   4.00   4.59   5.27   6.19   7.13   8.20   9.44  10.86  12.14  14.33  16.48  18.35  21.04  24.12  27.65
 j=4                                2.36   2.73   3.18   3.65   4.19   4.92   5.66   6.52   7.50   8.63   9.65  11.40  13.10  14.59  16.72  19.17  21.98
 j=5                                       2.17   2.53   2.90   3.33   3.91   4.50   5.18   5.96   6.86   7.67   9.06  10.42  11.60  13.29  15.24  17.47
 j=6                                              2.01   2.31   2.65   3.11   3.58   4.12   4.74   5.45   6.10   7.20   8.28   9.22  10.57  12.11  13.89
 j=7                                                     1.83   2.11   2.47   2.85   3.28   3.77   4.34   4.85   5.72   6.58   7.33   8.40   9.63  11.04
 j=8                                                            1.67   1.96   2.26   2.60   3.00   3.45   3.85   4.55   5.23   5.82   6.68   7.65   8.78
 j=9                                                                   1.56   1.80   2.07   2.38   2.74   3.06   3.62   4.16   4.63   5.31   6.08   6.98
 j=10                                                                        1.43   1.65   1.89   2.18   2.43   2.87   3.31   3.68   4.22   4.84   5.55
 j=11                                                                              1.31   1.50   1.73   1.93   2.29   2.63   2.93   3.35   3.84   4.41
 j=12                                                                                     1.20   1.38   1.54   1.82   2.09   2.33   2.67   3.06   3.50
 j=13                                                                                            1.09   1.22   1.44   1.66   1.85   2.12   2.43   2.79
 j=14                                                                                                   0.97   1.15   1.32   1.47   1.68   1.93   2.21
 j=15                                                                                                          0.91   1.05   1.17   1.34   1.54   1.76
 j=16                                                                                                                 0.83   0.93   1.06   1.22   1.40
 j=17                                                                                                                        0.74   0.85   0.97   1.11
 j=18                                                                                                                               0.67   0.85   0.88
 j=19                                                                                                                                      0.61   0.70
 j=20                                                                                                                                             0.56
```

The BDT drift parameters are much larger in magnitude than those in Ho-Lee, which reflects how the calibration is performed in log-rate. They also exhibit greater variation over time, indicating that more substantial adjustments are required to match the observed yield curve. The negative drift values at steps 14 and 17 show local term structure curvature rather than a broad downward shift in the level of rates.

All rates remain positive, with the minimum rate being 0.56% and the maximum reaching 55.05%. The spread is nearly three times the corresponding spread in Ho-Lee. This is because shocks are multiplicative, rate dispersion grows non-linearly as time rises. The BDT tree therefore produces a much more skewed and dispersed distribution of future rates.

## (b) Valuing the Mortgage Contract

To determine the par mortgage rate, $c^*$, that UBS should charge on its 10-year fixed rate mortgage such that the contract is priced at par (£100) at origination. The par rate represents the equilibrium rate that balances lender profitability with borrower demand.

The mortgage is modelled as a fixed rate, fully amortising loan with an American style prepayment option held by the borrower. Borrowers are assumed to behave optimally, refinancing whenever it is financially beneficial to do so and default risk is ignored. The par mortgage rate is determined separately under the Ho-Lee and BDT interest rate trees calibrated above, allowing for a direct comparison of how different interest rate dynamics affect mortgage pricing.

To value the mortgage, the outstanding principal balance at each time step must first be determined, as this represents the strike price of the borrower's prepayment option. Each payment is broken down into an interest component and a principal component, with the outstanding balance declining accordingly.

Using backward induction on the calibrated interest rate trees, we are able to value the contracts. Firstly, the mortgage is valued assuming no prepayment, treating it as a standard fixed rate bond, representing the value of the mortgage cash flows in the absence of borrower optionality.

Second, the borrower's prepayment right is modelled as an American call option. At each node, the borrower compares the value of immediate exercise with the continuation value. The exercise value is given by finding the maximum of the value of the no prepayment minus the outstanding balance, and zero. While the continuation value is the discounted expected future option value. The option value is therefore the maximum of these two quantities at each node. Prepayment occurs when the present value of remaining payments exceeds the outstanding balance, making refinancing advantageous.

The value of the mortgage to UBS is obtained by subtracting the option value from the no-prepayment value.

The par mortgage rate is determined by solving for the rate that ensures the present value is equal to £100. Brent's method is used to find the root of the equation. The root finding is calculated on the full mortgage value including the option. Ignoring the option would lead to an underestimation of the required coupon rate.

| Model | Par Rate (c*) | Payment | Value No Prepayment | Option Value | Full Mortgage Value |
|:---|:---:|:---:|:---:|:---:|:---:|
| Ho-Lee | 4.678% | £6.33 | £102.42 | £2.42 | £100.00 |
| BDT | 4.726% | £6.35 | £102.66 | £2.66 | £100.00 |

Both models produce a full mortgage value of £100, confirming that the root-finding process has correctly identified the par mortgage rate. The no-prepayment value exceeds par in both models, reflecting that the mortgage coupon is set above the current short rate. This value decreases in high-rate nodes due to heavier discounting and increases in low-rate nodes, demonstrating the inverse relationship between bond prices and interest rates.

The prepayment option is most valuable in low-rate environments, where refinancing is optimal and declines both in high-rate states and toward maturity. At nodes where prepayment is optimal, the mortgage value is capped at the outstanding balance, reflecting the borrower's ability to refinance. This creates a price ceiling and shows the negative convexity of mortgage contracts.

The BDT model produces a higher mortgage rate (4.726%) than the Ho-Lee model (4.6782%). This difference is because of the higher value of the prepayment option under BDT, £2.66 compared to £2.42 under Ho-Lee. Since the mortgage must be valued at par, a higher option value requires a higher coupon rate.

This result can be explained by the different interest rate movements of the two models. The BDT model generates greater rate dispersion and a more realistic distribution of future interest rates, increasing the likelihood of refinancing opportunities. In contrast, the Ho-Lee model, produces a narrower distribution and allows negative rates, which are less economically plausible. As a result, the prepayment option is more valuable under BDT, leading to a higher required mortgage rate.

The par mortgage rate can be interpreted as consisting of two key components. The first is compensation for the time value of money, reflecting the risk-free yield curve. The second is compensation for the borrower's embedded prepayment option, which is a cost to the lender. This option value is effectively spread over the life of the mortgage through a higher coupon rate.

The resulting par rate between 4.68%-4.73% is above both the current short rate and the 10-year spot rate, reflecting the upward sloping yield curve and the economic cost of borrower flexibility. The difference between the Ho-Lee and BDT results highlights that mortgage pricing is sensitive not only to the current term structure, but also to the assumptions about the distribution of future interest rates.

## (c) Constructing MBSs

Constructing securities backed by these mortgages, including a pass-through (PT) security, a principal-only (PO) strip, and an interest-only (IO) strip. The pass-through security entitles investors to all cash flows generated by the group of mortgages, whereas the PO strip receives only principal repayments and the IO strip receives only interest payments. All three securities have the same 10-year maturity as the underlying mortgages and are valued over 20 semi-annual periods.

To recover the cost of securitisation and retain an economic margin, UBS offers the pass-through securities at a coupon rate 50 basis points below the par mortgage rate derived earlier.

The pass-through rate for Ho-Lee is 4.1782% and for BDT it is 4.2264%.

The three securities represent a complete decomposition of the mortgage pool's cash flows. At each payment date, the total cash flow from the mortgage consists of interest and scheduled principal repayment. The pass-through security is equal to the sum of the PO and IO strips at every node.

The three securities are valued by backward induction on the calibrated interest rate trees, following the same general framework as in our mortgage valuation. The borrower's prepayment decision is inherited directly from the mortgage valuation, whenever the mortgage value equals the outstanding balance, prepayment is optimal, and the mortgage is terminated. MBS investors do not control this decision; they simply receive the cash flows that result from the borrower's exercise behaviour.

At each node there are two possible cases. If prepayment occurs, the borrower repays the outstanding balance in full. In that case, both the pass-through and principal-only investors would receive the remaining principal balance, while the interest-only investors receive nothing because the mortgage is paid off and no further interest is paid. The IO strip is highly exposed to prepayment risk, since falling rates accelerate refinancing and thereby shortening the stream of future interest payments.

If prepayment does not occur, the value of each security is computed as the discounted expected value of its future cash flows plus the relevant period cash flow. Thus, the PT security receives both pass-through interest and scheduled principal, the PO receives scheduled principal only and the IO receives pass-through interest only. The terminal condition for all three securities is zero at maturity, since the mortgages are fully paid off then.

In addition to fair values, spot rate durations are computed to measure the interest rate sensitivity of each security. Duration is calculated as a one-step finite difference using the first period up and down nodes. This approximates the percentage change in price per unit in the short rate.

| Security | Ho-Lee Price | BDT Price | Difference | HL Duration | BDT Duration |
|:---|:---:|:---:|:---:|:---:|:---:|
| PT | £98.99 | £99.03 | £0.0346 | 2.9001 | 3.3448 |
| PO | £90.67 | £90.90 | £0.0229 | 17.7373 | 19.9514 |
| IO | £8.32 | £8.13 | (£0.1943) | -158.8042 | -182.4394 |

The pass-through security is valued at £98.99 under Ho-Lee and £99.03 under BDT, in both cases slightly below par. This discount is the direct consequence of the 50-basis point securitisation spread retained by UBS. As the MBS investor receive a coupon below the fair mortgage rate, the present value of the pass-through cash flows is correspondingly below £100. It reflects the fact that investors are not receiving the full coupon paid by borrowers.

The PT durations of 2.90 years under Ho-Lee and 3.34 years under BDT indicate moderate interest rate sensitivity. These durations are lower than would be expected for an otherwise similar non-callable bond, because prepayment shortens effective maturity when rates fall. The embedded prepayment option therefore dampens the duration of the pass-through security.

The PO strip is valued at £90.67 under Ho-Lee and £90.90 under BDT, representing a substantial discount from par. This is because the PO investor receives no interest payments and must wait to recover principal over time. The discounted price therefore reflects the time value of money applied to a delayed stream of principal repayments.

Unlike the PT, the PO benefits from faster prepayment. When rates fall and borrowers refinance, the outstanding principal is returned sooner, which increases the present value of the PO cash flows. As a result, the PO has a large positive duration of 17.74 years under Ho-Lee and 19.95 under BDT. In this sense, the PO behaves like a highly rate-sensitive claim on principal acceleration.

This feature is visible directly in the tree. Under low-rate states, where prepayment is more likely, the PO value rises sharply because principal is returned earlier. In the Ho-Lee tree at step $i = 4$, for example, the PO is worth £63.24 at the highest-rate node but £83.56 at the lowest-rate node, with the latter equal to the outstanding balance, confirming that immediate prepayment occurs.

The IO strip is valued at £8.32 under Ho-Lee and £8.13 under BDT. Its value is relatively small because it receives only interest payments, which decline as the mortgage amortises, and those payments cease entirely if the borrower prepays. The IO is therefore the security most exposed to prepayment risk.

This exposure is why it has a strongly negative duration, -158.80 years under Ho-Lee and -182.44 years under BDT, meaning that its value rises when rates rise and falls when rates fall. The reason is because lower rates increase the likelihood of refinancing, which shortens the future interest stream and reduces the value of the IO despite the benefit of lower discounting. Conversely, when rates rise, prepayment becomes less likely, allowing the IO to continue receiving interest payments for longer.

This effect is clearly visible in the valuation tree. Under Ho-Lee for example, the IO has positive value at high-rate nodes but falls to zero at sufficiently low-rate nodes because prepayment occurs immediately and interest payments stop. This extreme sensitivity explains why IO strips may be used as hedging instruments against the positive duration of other fixed income positions.

Relative to Ho-Lee, the BDT model produces a slightly higher value for the PT and PO securities and a slightly lower value for the IO strip. This is because the BDT model generates wider rate dispersion than Ho-Lee. This increases the frequency and value of low-rate scenarios in which borrowers refinance. Earlier prepayment benefits PO investors, since principal is returned sooner, but harms IO investors, since the stream of interest payments is cut short. The PT combines both effects, so the increase in PO value and decrease in IO value partly offset one another, leaving only a modest net difference in PT value.

The BDT model also produces larger durations for all three securities. The pass-through is more rate-sensitive, the PO has even higher positive duration and the IO becomes even more negatively sensitive. This reflects the larger convexity and stronger rate sensitivity due to the lognormal nature of BDT.

The 50-basis point retained spread generates a modest but reliable stream of securitisation income while transferring most of the interest rate and prepayment risk to MBS investors. The precise present value of that spread differs slightly across the two models because the mortgage and pass-through rates differ.

The three securities appeal to different classes of investors because they isolate different components of mortgage cash flow risk. The pass-through security is most suitable for investors such as pension funds and insurers seeking relatively stable cash flows and moderate duration exposure. It provides both interest and principal repayments.

The PO strip is attractive to investors who expect rates to decline, since it benefits from accelerated prepayment and earlier return of principal. In effect, it is a leveraged exposure to lower rates and faster refinancing. Its discount price and high positive duration make it particularly attractive to investors seeking long duration exposure without purchasing a very long-dated bond.

The IO strip is most attractive either as a hedge or to investors who expect rates to remain high or rise further. As its value increases when prepayment slows, the IO has strongly negative duration and can offset the positive duration of other mortgage related or fixed income positions. This makes it a useful tool for managing prepayment and interest rate risk in portfolios.

A key limitation of all three valuations is their dependence on the prepayment model. Prepayment is assumed to be fully rational and interest-rate driven. In practice, borrowers prepay for many additional reasons and often fail to refinance even when it is financially optimal. This means that actual MBS performance may differ from the values generated by our models.

## (d) Interest Rate Trees with MC

For this question, we used the Ho-Lee and BDT trees calculated in part A. To create a forward simulation for the Monte Carlo, we created a 100,000 by 20 matrix of random variables between 0 and 1, where each row is a different path simulation of the tree, and each column represents a period in the tree. When starting the simulation, if the first value in the row is less than 0.5, the rate moves up, if it is greater than 0.5, the rate moves down. This is repeated for each 20 periods in the tree (each 20 columns of the matrix), until one simulated path from the start node at 0 to the end node at 20 is found. This process is then repeated 100,000 times for each row in the matrix, creating 100,000 simulated paths. All paths simulated are then stored and applied to both the Ho-Lee and BDT trees found in part A, this ensures that both models are directly comparable since the paths taken are identical, with the only difference coming from the different interest rates at each node in the tree. The mean rate is the average rate at T=10 of all 100,000 simulations. The histograms below show the rates at year 10, against the probability density function found from the frequency of the simulated paths reaching the specified rate.

Although both models found a similar mean interest rate, they have significantly different distributions, as shown below. The Ho-Lee model has a normal distribution, with a tighter range and higher density around the mean. This is because the Ho-Lee adds a fixed shock at each step, therefore the rate-spread is bounded and symmetric. The BDT on the other hand, is skewed rightwards, this is because BDT adds shocks to the log rate, so the rate-spread space grows asymmetrically through compounding. Overall, the BDT is lognormal and does not allow for negative rates, unlike the Ho-Lee model. The BDT model also has a larger range than the Ho-Lee model and lower density around the mean. Both models have relatively unrealistic interest rates at the tails, with the Ho-Lee having a rate of -3.5%, and the BDT having a rate of 55%, however the mass at these regions is minimal on the histogram, this is because the probability of taking 20 consecutive up or down movements is very small. The difference between the mean rates, of 6.13% for Ho-Lee and 6.33% for BDT can be explained by the BDT having more mass on higher rates than Ho-Lee and the Ho-Lee having some mass in the negative rates. Furthermore, the BDT being log normal leads to the mean always being above the median. Since both simulations are discrete and not continuous, they show distinctive spikes rather than a smooth distribution.

> **Figure (d):** Histograms of short rates at year 10 — Ho-Lee Model (Normal / Additive), Mean = 6.13%; BDT Model (Lognormal / Multiplicative), Mean = 6.33%.

## (e) Valuing the Mortgage with MC

Similar to part D, we created a 100,000 by 20 matrix of random variables between 0 and 1 to create a new set of simulated paths, keeping independence from part D. Starting at i=0 for each path, the one period discount factor is calculated, with the first payment computed from the mortgage rate found in part B, and discounted back one period to i=0. From i=1 onwards the continuation value is compared to against the outstanding balance. To find the continuation value at each period we use the mortgage value with prepayment tree found in part B. The formula for the continuation value is:

$$hold_{i,j} = e^{-r_{i,j}\cdot\Delta}\left(\tfrac{1}{2}V[i+1][j] + \tfrac{1}{2}V[i+1][j+1] + pmt\right)$$

Where $e^{-r_{i,j}\cdot\Delta}$ is the discount factor for the next period, and $V[i+1][j]$ is the value in the upstate of the mortgage and $V[i+1][j+1]$ is the value in the downstate, each with equal probability of occurring, and $pmt$ is the fixed semi-annual payment made by the borrower, all values are found from part B. The value $hold_{i,j}$ is compared against the value of the outstanding principal at time i found in part B. If $hold_{i,j}$ is less than the outstanding balance, then continuation is cheaper and UBS will receive $pmt$ at the end of period i, this is then discounted using the cumulative discount factor multiplied by the discount factor for period i to find the present value at i=0, and the path moves to the next node. If $hold_{i,j}$ is greater or equal to the outstanding principal, then prepayment is optimal and the mortgage is prepaid at the start of period i and UBS will receive the outstanding principal, this is then discounted to time i=0 using cumulative discount factor. This process is repeated for each node in the path until either the mortgage is prepaid early or matures at the end of the term at i=20. For each row in the 100,000 by 20 matrix, the present value of the future payments is calculated using the process above then stored. We then calculate the mean and the confidence interval for all present values of the 100,000 paths simulated.

We find the point estimate as the mean of all the present values of each 100,000 paths simulated, which is the mortgage price. We find that the point estimate is £100.002985, with a 95% confidence interval of [99.987977, 100.017993]. This value is consistent with part B since the mortgage value of £100.0000 is within the 95% confidence interval found using the Monte Carlo simulation. Although there is a slight difference in value between the backwards induction method of part B and the simulations we used here (0.002985 difference), this does not invalidate the result since part B uses an exact probability weighted average, whereas the Monte Carlo uses random sampling. By the law of large numbers, using 100,000 simulations significantly reduces the gap in price converging on the true expectation, however with a finite set of simulations there is still some small sampling error.

## (f) Prepayment Modelling

### Modelling Mortgage Prepayment: The Baseline Framework

Having established the cash flow structure of a fixed-rate mortgage, we now explore how the possibility of early repayment can be incorporated into its price. Pricing the mortgage requires accounting for the fact that the borrower may choose to repay the outstanding principal early, ending the lender's claim to future cash flows. This prepayment option is central to how mortgages are priced.

### Backward induction on the interest rate tree

The mortgage is priced on a binomial interest rate tree calibrated in the earlier questions. If prepayment occurs, the pass-through security returns the outstanding principal $L_i$, since the borrower retires the remaining balance in full, with all future payments ceasing. If prepayment does not occur, the security's value at that node is the risk-neutral discounted expectation of its value at the two previous nodes, plus the scheduled cash flow arriving at step $i+1$:

$$P_{i,j}^{PT} = e^{-r_{i,j}\Delta}\left(\tfrac{1}{2}P_{i+1,j}^{PT} + \tfrac{1}{2}P_{i+1,j+1}^{PT} + CF(i+1)\right)$$

$CF(i+1)$ denotes the cash flow received at the next step: this is the sum of the interest payment and scheduled principal repayment. Starting from the terminal nodes, where the remaining balance is zero and the security is worth nothing, we work backwards through the entire tree to recover the fair value at $t=0$.

The formula above treats prepayment as a binary event at each node: either the mortgage is entirely repaid, or it is not. In reality, prepayment is better modelled as a partial reduction in the pool's outstanding balance at every period. This is where the prepayment rate enters the cash flow calculation directly. In our model with semi-annual steps, we utilise the Semi-Annual Mortality rate (SaM). This is the fraction of the remaining outstanding balance expected to prepay within each half-year period. The outstanding balance entering period $i+1$ is reduced by both the scheduled amortisation, and also by the unscheduled prepayment:

$$L_{i+1} = L_i \cdot (1 - \text{SaM}) - \text{scheduled principal}_{i+1}$$

This means $CF(i+1)$ is a function of SaM: a higher prepayment rate shrinks the pool balance faster, reducing future interest payments. This is because interest accrues on a smaller outstanding balance while the return of the principal is accelerated. The SaM is derived from the annual CPR via:

$$\text{SaM} = 1 - \sqrt{1 - \text{CPR}}$$

### Issues with fixed CPR

In the baseline model, like the PSA convention, CPR is a fixed, deterministic schedule that is not affected by the interest rate environment. This means the backward induction formula above applies the same prepayment-adjusted cash flows in every node at a given time step; the current level of $r_{i,j}$ has no bearing on how aggressively the pool prepays. Evidently this is a strong and unrealistic assumption. Borrowers are more likely to prepay when their existing coupon rate is significantly above the prevailing market rate, and far less likely to do so when rates have risen. A fixed CPR cannot capture this, nor can it reflect the accumulated history of the rate path, the age of the mortgage, or the degree to which a pool has already been depleted of its most rate-sensitive borrowers. It is precisely this limitation that motivates the extended model developed in the following sections.

### Modelling CPR

We shall now outline how CPR will be modelled to better reflect the prepayment behaviour observed in reality. We introduce a three-factor model to account for the: Refinancing Incentives, Seasoning, and Burnout effects.

#### Refinancing Incentives

The refinancing incentive captures the economic logic that a borrower holding a fixed-rate mortgage will want to prepay and replace it with a cheaper loan whenever current market rates fall meaningfully below their contractual coupon. From the borrower's perspective, the mortgage is a liability. When rates fall, that liability becomes expensive relative to what the market now offers. Rationally, the borrower should extinguish it and refinance at the lower rate, capturing the interest saving over the remaining life of the loan. The incentive is not just a preference, but an option-like payoff: the borrower holds an embedded call on their own debt, and the deeper that option moves in-the-money, the more probable prepayment is.

(Berger, et al., 2021) model the refinancing incentive as a rate gap (the difference between the borrower's contractual coupon and the prevailing market rate on a comparable new loan). Using a linked borrower-loan panel covering approximately 180 million US mortgage records from 2005 to 2017, they uncover striking results: Prepayment hazard is low and essentially flat for loans with negative rate gaps (borrowers with no financial incentive to move rarely prepay). As this gap turns positive and widens, the hazard rises sharply. It then plateaus at around 2 percent per month for rate gaps above approximately 100 basis points. This "step-like" shape is the empirical signature of the incentive effect: prepayment is state-dependent, responding systematically and non-linearly to rate incentives.

**Specification:**

In our model, we capture the refinancing incentive through a rate-spread term entering the conditional prepayment rate. Letting $c$ denote the mortgage coupon, and $r_t$ the prevailing refinancing rate at time $t$.

$$\beta_1 \cdot incentive = \beta_1 \cdot (c - r_t)$$

We set $\beta_1 > 0$, so that a wider positive spread (a deeper in-the-money refinancing option) raises the CPR, just as (Berger, et al., 2021) identify.

(Richard & Roll, 1989) argue that the theoretically correct measure is the ratio $C/R$, rather than the spread $c - r_t$, since the borrower's decision is fundamentally an NPV comparison (how many times more expensive is the old liability relative to refinancing), rather than an additive one. However, (Richard & Roll, 1989) note, $C/R$, and the spread are nearly collinear for mid-life mortgages in moderate rate environments. With the spread specification being standard across the empirical literature, we retain $\beta_1 \cdot (c - r_t)$ as a well-precedented formulation.

#### Seasoning

As a mortgage ages, its prepayment rate tends to rise from near zero at origination toward a higher steady-state level before eventually plateauing. Newly issued mortgages are held by a borrower who has just gone through the effort and cost of taking out a loan, making immediate prepayment unlikely regardless of the rate environment. Over time, as the borrower settles into the property, and becomes more financially aware of their options, the probability of prepayment naturally increases. Seasoning therefore reflects the simple reality that prepayment propensity is not constant from day one, it builds as the mortgage matures.

(Richard & Roll, 1989) show that prepayment rates rise from near zero at origination, reaching a stable level only as the mortgage matures. Crucially, they note that the speed of this ramp depends on the coupon relative to current refinancing rates. They show that a premium pool (C/R = 1.2) seasons fully in around 30 months, while a current-coupon pool takes nearly five years and a discount pool can take nine or more. Borrowers with a deep incentive to move are also quick to season; those with little incentive drag the process out, confirming seasoning as a genuine, rate-dependent empirical phenomenon.

**Specification:**

Following the PSA convention, the seasoning factor ramps linearly from zero to one over the first 30 months and is flat thereafter. Our model is implemented on a semi-annual tree, so each time step, $i$, corresponds to $t = 0.5i$ years, or equivalently $6i$ months. The 30-month PSA ramp therefore spans exactly 5 years, which is 10 semi-annual steps. The seasoning contribution to CPR is:

$$\beta_2 \cdot seasoning(i) = \beta_2 \cdot \min\left(\frac{i}{10}, 1\right)$$

where $i$ is the step number on the tree. The denominator of 10 arises directly from this time-step conversion: 30 months / 6 months per step = 10 steps. For $i \geq 10$, the factor equals 1, meaning the pool is fully seasoned and age no longer adds to CPR through this channel.

#### Burnout

Burnout refers to the tendency of a mortgage pool's prepayment rate to slow over time, even when the refinancing incentive remains strong. Burnout arises because of borrower heterogeneity: mortgage pools are not made up of identical borrowers. When rates fall and the refinancing option moves in-the-money, the most rate-sensitive borrowers (those with the lowest transaction costs, and greatest financial awareness) prepay first. What remains is increasingly concentrated in borrowers who have not acted. If rates stay low or fall further, this residual group prepays at a slower rate than the original pool would have because the most reactive movers have left the pool.

(Richard & Roll, 1989) first coined the term "premium burnout", documenting that pools which had spent extended periods deep in-the-money showed attenuated prepayment responses to subsequent rate declines. However, they treat burnout as an exogenously given calibrated multiplier without a structural explanation as to why it arose. (Hall, 2000) provides this explanation, framing burnout as a problem of unobserved heterogeneity. He argues that borrowers differ in unobservable ways (risk preferences, financial sophistication, etc). As the most prepayment-prone borrowers exit the pool, the surviving population shift systematically toward the slow end. Hall's gamma mixture model improves from previous approaches, integrating a continuous distribution over unobserved borrower types directly into the likelihood, producing burnout endogenously. He finds that the gamma mixture consistently outperforms the previous logit specifications across rate cycles, particularly after rates stabilise following a decline.

**Specification:**

We capture burnout through the pool's cumulative in-the-money exposure since origination. At each semi-annual step $i$, we sum all prior periods in which the refinancing spread was positive, weighted by the step size:

$$\beta_3 \cdot burnout_i = \beta_3 \cdot \sum_{\tau=0}^{i-1} \max(c - r_t, 0) \cdot \Delta$$

where $\Delta = 0.5$, $c$ is the pool coupon, and $r_t$ is the refinancing rate at step $t$. Since $\beta_3 < 0$, this term reduces CPR. The economic logic follows that: the larger the cumulative in-the-money exposure, the more rate-sensitive borrowers will have already left the pool along that particular rate path. This leaves a surviving population increasingly concentrated in slow prepayers. A pool that has spent many periods with rates well below the coupon will have a large burnout term, and correspondingly depressed future prepayment – even if current rates remain attractive.

This specification is path-dependent: two nodes at the same point on the interest rate tree may carry very different burnout values depending on the rate history that led there. As such, we would expect the BDT-calibrated tree to mis-price the burnout effect. Instead, the Monte Carlo simulation should "correctly" price the effect as each simulated path accumulates its own burnout history independently.

### The Logistic CPR Model

We model CPR using a logistic transformation of a linear index $Z_i^s$:

$$\text{CPR}_i^s = \frac{C_{\max}}{1 + e^{-Z_i^s}}$$

where $Z_i^s$ is the linear index aggregating all three prepayment drivers:

$$Z_i^s = \alpha + \beta_1 \cdot \text{Incentive}_i^s + \beta_2 \cdot \text{Seasoning}_i + \beta_3 \cdot \text{Burnout}_i^s$$

There are two strong justifications for this functional form, both grounded in the empirical literature.

**Boundedness.** CPR is a conditional prepayment rate and must lie in $[0, C_{\max}]$, where $C_{\max} \leq 1$. It also cannot be negative, nor exceed its ceiling. The logistic transformation enforces these bounds automatically regardless of the values taken by the incentive, seasoning, and burnout terms.

The empirical prepayment hazard documented in the literature is observed to be non-linear in the refinancing incentive. (Berger, et al., 2021), estimating a non-parametric hazard across 180 million loans, find a pronounced step-like shape: CPR is low and flat for negative rate gaps, rises sharply as the gap turns positive, then plateaus at around 2% per month for large positive gaps. This is the shape a logistic function produces, making it a natural approximation to the empirically observed hazard. (Richard & Roll, 1989) similarly document a highly non-linear, convex relationship between CPR and C/R. Similarly, (Schwartz & Torous, 1989) adopt a log-logistic baseline hazard, noting that a purely linear specification cannot capture the way prepayment accelerates once refinancing becomes sufficiently attractive. (Hall, 2000) also works within a logit framework as his baseline because the conditional probability of prepayment is naturally bounded between zero and one.

Notably, the right-tail asymptote of the logistic function, where CPR flattens and ceases to respond to further increases in $Z_i^s$, captures a key empirical observation. It reflects the idea that even with a very strong refinancing incentive, prepayment rates plateau rather than continue rising, consistent with (Berger, et al., 2021), who find an empirical plateau at around 2% monthly. This highlights the burnout phenomenon: the pool's most rate-sensitive borrowers have already exited, so additional levels of incentive yield diminishing returns. The left-tail asymptote exhibits similar logic for the seasoning effect: newly issued mortgages have a CPR near zero not because the incentive is zero, but because the pool has not yet had time for any borrowers to act. Therefore, the logistic form does not simply impose convenient mathematical bounds; its shape actively reflects the economic forces of our model.

### Calibrating the Model

**1. $C_{max}$:**

The first step is to set our CPR ceiling: this is the maximum annualised prepayment rate that they logistic model is asymptotical towards. Using quarterly Mortgage Lenders & Administrators Return data, dating back to 2007 (Bank of England, 2026), we constructed annualised CPR from MLAR Table 1.21 using the stock-flow identity.

$$\text{Redemptions}_t = \text{Stock}_{t-1} + \text{Gross Advances}_t - \text{Stock}_t$$

$$\text{CPR}_t = 1 - \left(1 - \frac{\text{Redemptions}_t}{\text{Stock}_{t-1}}\right)^4$$

A peak CPR of 29.6% was observed in 2007 Q3; as such, we set $C_{max} = 0.3$.

**2. $C_{baseline}$:**

Next, we set out baseline CPR – our empirical anchor. This is the CPR that is observed for a fully seasoned, at-the-money pool with no burnout history. As we are looking to value a 2026 mortgage, we calibrate this value with the Bank of England MLAR data from 2022 onwards, of which the average CPR is 14.52%; hence we adopt $C_{baseline} = 0.145$.

**3. $\tilde{\alpha}$:**

$\tilde{\alpha}$ is the effective intercept for a fully seasoned, at-the-money, no-burnout loan:

$$\tilde{\alpha} = \alpha + \beta_2 \ ; \quad \beta_1 \cdot incentive = \beta_3 \cdot burnout = 0$$

$$C_{baseline} = \frac{C_{max}}{1 + e^{-\tilde{\alpha}}}$$

Plugging in our CPR parameters & solving yields $\tilde{\alpha} = -0.0666$.

**4. $\beta_2$ (seasoning)**

This parameter dictates how much CPR increases as a loan ages from origination to full seasoning. (Richard & Roll, 1989) document that a par coupon pool ($C/R = 1.0$) takes approximately 60 months to season fully, with prepayment rates near zero at origination rising to their plateau level by month 60. They show that at a constant refinancing rate, a newly issued par pool prepays at approximately 1–2% CPR while the same pool at full seasoning prepays at approximately 12–15% CPR. Taking the midpoint of this empirical range gives a seasoned-to-newborn ratio of approximately 12x, which we adopt as our calibration target.

At zero incentive and zero burnout, the index reduces to $Z = \alpha + \beta_2 \cdot seasoning_i$, giving:

$$\text{CPR}_{seasoned} = \frac{C_{max}}{1 + e^{-(\alpha + \beta_2)}} = \frac{C_{max}}{1 + e^{-\tilde{\alpha}}}$$

$$\text{CPR}_{newborn} = \frac{C_{max}}{1 + e^{-\alpha}} = \frac{C_{max}}{1 + e^{-(\tilde{\alpha} - \beta_2)}}$$

Imposing the 12x condition & solving, we yield $\beta_2 = 3.10$; Therefore, $\alpha = -3.17$.

**5. $\beta_1$ (incentive)**

This is the sensitivity of CPR to the refinancing incentive spread. The target for $\beta_1$ is drawn from (Richard & Roll, 1989). Their paper reports actual average CPR for seasoned pools across a range of coupon rates and prevailing refinancing rate environments, allowing the pure incentive effect to be isolated by comparing prepayment rates for the same pool as market rates move relative to the contractual coupon.

To calibrate $\beta_1$, we focus on pools where the coupon rate is approximately equal to the prevailing refinancing. For such pools, a 100bp decline in the prevailing refinancing rate, a 100bp move into the money, is associated with an increase in actual average CPR from approximately 7.6% to 13.7%, a multiplicative increase of 1.8x. This relationship reflects the non-linear S-shaped response of prepayment to refinancing incentive. It is important to note that these observations are drawn from seasoned pools, meaning the seasoning effect is already fully embedded in the reported CPR figures. We therefore evaluate the 1.8x condition at our model's fully seasoned, zero burnout reference point.

$$\text{CPR}_{+100bps} = 1.8 \times \text{CPR}_{baseline} = 1.8 \times 14.5\% = 26.1\%$$

Now solving for $\beta_1$:

$$\frac{C_{max}}{1 + e^{-(\tilde{\alpha} - 0.01\beta_1)}} = 0.261; \quad \beta_1 = 197$$

**6. $\beta_3$ (burnout)**

Two constraints are applied jointly: the first from (Richard & Roll, 1989) pins the speed at which burnout suppresses CPR, and the second from (Hall, 2000) pins the floor below which CPR should not fall regardless of accumulated burnout.

(Richard & Roll, 1989) track prepayment rates for premium pools sustained in an in-the-money environment over time. As rate-sensitive borrowers progressively exit, the pool becomes increasingly concentrated in slow pre-payers, and CPR declines even as the financial incentive remains unchanged. They find that after approximately 40 months of sustained in-the-money exposure, the prepayment rate falls to approximately 50% of its no-burnout level. This implies a burnout multiplier of 0.5, pinning the rate at which accumulated burnout should suppress CPR in our model.

(Hall, 2000) estimates a gamma mixture model that treats burnout as an endogenous consequence of unobservable borrower heterogeneity rather than an explicit proxy. His simulations show that even heavily burned-out pools maintain a prepayment floor of approximately 0.6% SMM per month, equivalent to 7% CPR annually, reflecting irreducible exogenous prepayment from house sales, divorce and job relocation that operates independently of the rate environment. This prevents CPR collapsing to zero at extreme burnout levels, which would materially overvalue the mortgage.

The (Hall, 2000) floor constraint is binding. Setting CPR = 7% at Burnout = 0.10 represents five years of sustained 200bps in-the-money exposure, the most extreme scenario likely to be generated by the BDT paths:

$$\frac{0.3}{1 + e^{-(3.873 - 0.10\beta_3)}} = 0.07; \quad \beta_3 = 51$$

Verifying this beta against (Richard & Roll, 1989), at 7 semi-annual periods of 200bps exposure, Burnout = 0.07 & $\beta_3 = 51$, produces a burnout multiplier of $17.3\%/29.4\% = 0.59$. This is broadly consistent with their 0.50 finding at 40 months. The modest difference reflects evidence from (Hall, 2000) that proxy-based approaches tend to overstate burnout speed.

### Calibrated Model

Having calibrated all parameters to fit the empirical evidence & UK mortgage data, the final CPR model is as follows:

$$\text{CPR}_t^s = \frac{0.3}{1 + e^{-(-3.17 + 197\cdot incentive_t^s + 3.10\cdot seasoning_t^s - 51\cdot burnout_t^s)}}$$

### Pricing the MBS

The Monte Carlo pricing procedure in part (f) is structurally identical to that used in part (e): 100,000 paths are simulated through the BDT interest rate tree by drawing uniform random numbers at each period to determine whether the rate moves up or down, cash flows are discounted along each path using the realised short rates, and the mortgage value is taken as the average across all paths. The single departure from part (e) is that prepayment is now endogenized; the decision to prepay at any node is governed by the three-factor logistic CPR model described above.

At each node i along a simulated path, the model computes a CPR from the prevailing refinancing incentive, the loan's seasoning, and the accumulated burnout, then converts this to a per-period prepayment probability, the SaM. The node therefore has two possible outcomes: the mortgage either prepays, with probability SaM, or it survives to the next period, with probability $1 - \text{SaM}$. To resolve this, the code draws a second uniform random variable U ~ U(0,1). If U < SaM, prepayment occurs, the lender receives the outstanding balance bal[i], and the path terminates. If U ≥ SaM, no prepayment occurs: the lender receives the scheduled coupon and the simulation advances to the next node. Even deeply in-the-money borrowers are not guaranteed to prepay; they simply face a higher threshold SaM.

This construction introduces path dependency, which is what distinguishes part (f) from prior parts. Burnout accumulates at every surviving node, meaning it depends on every rate realisation the path has passed through. Two paths that arrive at the same node (i, j) will carry different burnout levels if their rate histories differ, and will therefore produce different CPRs and SaMs going forward. A node's value depends on the history leading to it. This breaks the recombining structure of the BDT tree entirely and makes backward induction on the tree infeasible, which is why a Monte Carlo simulation is necessary to price this security.

### MBS Values

> **Figure F3 — Part (f): MBS Tranche Prices** (Binary prepayment (part c) vs Logistic CPR (part f)), price in £ per £100:
>
> | Tranche | Part (c) binary | Part (f) logistic CPR |
> |:---|:---:|:---:|
> | Pass-Through (PT) | £99.03 | £98.00 |
> | Principal Only (PO) | £90.90 | £83.34 |
> | Interest Only (IO) | £8.13 | £14.66 |

Here we compare the part (f) valuation to the part (c) valuation. It is noteworthy that the MBS featuring endogenous CPR modelling is worth less than the optimal prepayment MBS. The logistic MBS has a coupon of 4.2552%, representing a 47.12 bp fall relative to the part (c) MBS. This makes sense as (c) is the worst possible case from the lender's perspective, hence they charge a higher rate; borrowers will exercise their right to prepay with 100% efficiency & rationality. In (f), the frictions we have modelled mean that borrowers are sluggish, burnt-out & slow to get started on prepayments. Seasoning suppresses prepayments in the early years, so the lender collects more coupons before the borrower even begins to think about refinancing. Burnout means that after a prolonged low-rate period the pool is exhausted; the people who would have refinanced have already done so, and the remainder are the sticky non-refinancers. So even if rates stay low, prepayment slows. Probabilistic exercise means some borrowers simply never prepay even when it is deeply optimal to do so.

The net effect is that the lender receives more coupons, over a longer horizon, before prepayment terminates the stream compared to the fully rational model. The mortgage is therefore worth more to the lender at any given coupon. Therefore, to price it back at par (£100), you need a lower coupon. Hence the spread is negative.

### Sensitivity Analysis

> **Figure F1 — Part (f): Logistic CPR Model** — Conditional Prepayment Rate as a function of incentive and seasoning (left panel: CPR vs incentive at different Seasoning levels, Burnout = 0; right panel: CPR vs incentive at different Burnout levels, Seasoning = 1.0).

The left panel holds burnout at zero and varies seasoning. It is clear that at Sea = 0 (a brand-new mortgage), the logistic curve sits far to the right. Therefore, you require a very large rate incentive before CPR picks up meaningfully. As seasoning rises to 1.0, the curve shifts left, becoming responsive to rate incentives at progressively smaller rate gaps. All curves converge at the same ceiling (~30% CPR) for large enough incentives. The inflection point of the S-curve is what moves: seasoning effectively shifts the borrower's "readiness to refinance."

The right panel holds seasoning at 1.0 and varies burnout. Here the shift is dramatic and asymmetric: high burnout doesn't just translate the curve, it compresses and flattens it. At BO = 0.40, even 5% in-the-money borrowers are only reaching ~12% CPR. The surviving pool is populated almost entirely by non-refinancers; the S-curve barely lifts off the x-axis. This illustrates exactly why burnout is a first-order effect: it doesn't just reduce CPR uniformly, it selectively removes the most rate-sensitive borrowers from the pool, leaving behind a population that is intrinsically unresponsive.

## References

Bank of England, 2026. *Mortgage Lenders and Administrators Statistics - 2025 Q3.* [Online] Available at: https://www.bankofengland.co.uk/statistics/mortgage-lenders-and-administrators/2025/2025-q3 [Accessed 23 March 2026].

Berger, D., Milbradt, K., Tourre, F. & Vavra, J., 2021. Mortgage Prepayment and Path-Dependent Effects of Monetary Policy. *American Economic Association,* 111(9), pp. 2829-2878.

Hall, A., 2000. Controlling for Burnout in Estimating Mortgage Prepayment Models. *Journal of Housing Economics,* Volume 9, pp. 215-232.

Richard, S. & Roll, R., 1989. Prepayments on fixed-rate mortgage-backed securities. *The Journal of Portfolio Management,* pp. 73-82.

Schwartz, E. & Torous, W., 1989. Prepayment and the Valuation of Mortgage-Backed Securities. *The Journal of Finance,* 44(2), pp. 375-392.

## Appendix

**(d)** Sample paths plotted for Ho-Lee and BDT trees using the first 200 simulations.

> **Figure (Appendix d):** Sample short-rate paths — Ho-Lee Model (mean path rising from 3.48% to ~6.13% over 10 years) and BDT Model (mean path rising from 3.48% to ~6.33% over 10 years, with right-skewed dispersion).

**(e)** Convergence of MC with increasing N.

```
Convergence of MC estimate with number of paths (BDT):
N paths       Mean (£)     Std Error                   95% CI
100          100.366775    0.227360    [ 99.921149, 100.812401]
500          100.217981    0.099997    [100.021988, 100.413975]
1,000        100.004587    0.080984    [ 99.845859, 100.163315]
5,000         99.989990    0.034599    [ 99.922177, 100.057804]
10,000        99.971731    0.024574    [ 99.923567, 100.019895]
50,000        99.995234    0.010934    [ 99.973803, 100.016665]
100,000      100.002985    0.007657    [ 99.987977, 100.017993]
```
