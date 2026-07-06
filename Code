# -*- coding: utf-8 -*-

# =============================================================================
# FM405 SUMMATIVE - FILE 1: INTEREST RATE TREE CALIBRATION
# Parts (a): Ho-Lee and BDT Interest Rate Tree Calibration
# Data: Bank of England Nominal Spot Curve, 30 January 2026
# Note: Run all files on one document since later parts require outputs from
# earlier parts
# =============================================================================

import numpy as np
from scipy.optimize import brentq

# =============================================================================
# MARKET DATA
# =============================================================================

delta = 0.5   # semi-annual steps
N     = 20    # 20 steps = 10 years

# Bank of England nominal spot rates (%), 30 January 2026
spot_rates = {
    0.5:  3.48, 1.0:  3.55, 1.5:  3.59, 2.0:  3.63, 2.5:  3.67,
    3.0:  3.72, 3.5:  3.78, 4.0:  3.84, 4.5:  3.90, 5.0:  3.97,
    5.5:  4.04, 6.0:  4.11, 6.5:  4.18, 7.0:  4.25, 7.5:  4.31,
    8.0:  4.38, 8.5:  4.45, 9.0:  4.51, 9.5:  4.57, 10.0: 4.63
}

# Zero coupon bond prices per Â£100 face value
P = {}
for k, (maturity, rate) in enumerate(spot_rates.items(), start=1):
    P[k] = 100 * np.exp(-(rate / 100) * maturity)

# =============================================================================
# HO-LEE CALIBRATION
# Normal (additive) model: r_{i+1,j} = r_{i,j} + theta_i*delta +/- sigma*sqrt(delta)
# =============================================================================

sigma_HL = 68.02 / 10000    # 68.02 bp/yr normal vol -> 0.006802

r0_HL       = -np.log(P[1] / 100) / delta
r_tree_HL   = np.zeros((N + 1, N + 1))
r_tree_HL[0][0] = r0_HL
q_HL        = np.zeros((N + 1, N + 1))
q_HL[0][0]  = 1.0
theta_HL    = np.zeros(N)

for i in range(N - 1):
    target = P[i + 2]

    def bond_price_given_theta_HL(theta_i, i=i):
        bottom_next = r_tree_HL[i][i] + theta_i * delta - sigma_HL * np.sqrt(delta)
        q_next = np.zeros(i + 2)
        for j in range(i + 2):
            val = 0.0
            if j <= i:
                r_ij = r_tree_HL[i][i] + (i - j) * 2 * sigma_HL * np.sqrt(delta)
                val += 0.5 * q_HL[i][j] * np.exp(-r_ij * delta)
            if j >= 1 and (j - 1) <= i:
                r_ij1 = r_tree_HL[i][i] + (i - (j - 1)) * 2 * sigma_HL * np.sqrt(delta)
                val += 0.5 * q_HL[i][j - 1] * np.exp(-r_ij1 * delta)
            q_next[j] = val
        price = 0.0
        for j in range(i + 2):
            r_next_j = bottom_next + (i + 1 - j) * 2 * sigma_HL * np.sqrt(delta)
            price += q_next[j] * np.exp(-r_next_j * delta) * 100
        return price

    theta_HL[i] = brentq(
        lambda th: bond_price_given_theta_HL(th) - target, -2.0, 5.0
    )

    bottom_next = r_tree_HL[i][i] + theta_HL[i] * delta - sigma_HL * np.sqrt(delta)
    for j in range(i + 2):
        r_tree_HL[i + 1][j] = bottom_next + (i + 1 - j) * 2 * sigma_HL * np.sqrt(delta)

    for j in range(i + 2):
        val = 0.0
        if j <= i:
            r_ij = r_tree_HL[i][i] + (i - j) * 2 * sigma_HL * np.sqrt(delta)
            val += 0.5 * q_HL[i][j] * np.exp(-r_ij * delta)
        if j >= 1 and (j - 1) <= i:
            r_ij1 = r_tree_HL[i][i] + (i - (j - 1)) * 2 * sigma_HL * np.sqrt(delta)
            val += 0.5 * q_HL[i][j - 1] * np.exp(-r_ij1 * delta)
        q_HL[i + 1][j] = val

# Fill final column
bottom_N_HL = r_tree_HL[N - 1][N - 1] + theta_HL[N - 2] * delta - sigma_HL * np.sqrt(delta)
for j in range(N + 1):
    r_tree_HL[N][j] = bottom_N_HL + (N - j) * 2 * sigma_HL * np.sqrt(delta)

# --- Print Ho-Lee tree -------------------------------------------------------
W = 7
print("=" * (8 + W * (N + 1)))
print("Table: The Risk Neutral Ho-Lee Interest Rate Tree")
print(f"BoE Nominal Spot Curve, 30 Jan 2026 | sigma = {sigma_HL*100:.4f}% | delta = {delta} | N = {N}")
print("=" * (8 + W * (N + 1)))

print(f"{'Time T':>8}", end="")
for i in range(N + 1):
    print(f"{i * delta:>{W}.1f}", end="")
print()

print(f"{'Period i':>8}", end="")
for i in range(N + 1):
    print(f"{i:>{W}d}", end="")
print()

print(f"{'Î¸i(Ã100)':>8}", end="")
print(f"{'':>{W}}", end="")
for i in range(N):
    if i <= N - 2:
        print(f"{theta_HL[i] * 100:>{W}.4f}", end="")
    else:
        print(f"{'':>{W}}", end="")
print()
print()

for j in range(N + 1):
    print(f"{'j='+str(j):>8}", end="")
    for i in range(N + 1):
        if j <= i:
            print(f"{r_tree_HL[i][j] * 100:>{W}.2f}", end="")
        else:
            print(f"{'':>{W}}", end="")
    print()

print()
print("=" * 55)
print("Ho-Lee Calibration Check")
print("=" * 55)
print(f"{'k':>4}  {'Maturity':>10}  {'Market':>10}  {'Model':>10}  {'Diff':>10}")
for k in range(1, N + 1):
    maturity = k * delta
    market   = P[k]
    if k < N:
        model = sum(q_HL[k][j] for j in range(k + 1)) * 100
    else:
        model = sum(
            q_HL[N-1][j] * np.exp(-r_tree_HL[N-1][j] * delta)
            for j in range(N)
        ) * 100
    diff = model - market
    flag = "  *last step" if k == N else ""
    print(f"{k:>4}  {maturity:>10.1f}  {market:>10.4f}  {model:>10.4f}  {diff:>10.6f}{flag}")

# =============================================================================
# BDT CALIBRATION
# Lognormal model: log(r_{i+1,j}) = log(r_{i,j}) + theta_i*delta +/- sigma*sqrt(delta)
# =============================================================================

sigma_BDT = 16.23 / 100     # 16.23% Black vol -> 0.1623

r0_BDT      = -np.log(P[1] / 100) / delta
z0_BDT      = np.log(r0_BDT)

z_tree_BDT  = np.zeros((N + 1, N + 1))
r_tree_BDT  = np.zeros((N + 1, N + 1))
z_tree_BDT[0][0] = z0_BDT
r_tree_BDT[0][0] = r0_BDT

q_BDT       = np.zeros((N + 1, N + 1))
q_BDT[0][0] = 1.0
theta_BDT   = np.zeros(N)

for i in range(N - 1):
    target = P[i + 2]

    def bond_price_given_theta_BDT(theta_i, i=i):
        z_bottom_next = z_tree_BDT[i][i] + theta_i * delta - sigma_BDT * np.sqrt(delta)
        q_next = np.zeros(i + 2)
        for j in range(i + 2):
            val = 0.0
            if j <= i:
                r_ij = np.exp(z_tree_BDT[i][i] + (i - j) * 2 * sigma_BDT * np.sqrt(delta))
                val += 0.5 * q_BDT[i][j] * np.exp(-r_ij * delta)
            if j >= 1 and (j - 1) <= i:
                r_ij1 = np.exp(z_tree_BDT[i][i] + (i - (j - 1)) * 2 * sigma_BDT * np.sqrt(delta))
                val += 0.5 * q_BDT[i][j - 1] * np.exp(-r_ij1 * delta)
            q_next[j] = val
        price = 0.0
        for j in range(i + 2):
            z_next_j = z_bottom_next + (i + 1 - j) * 2 * sigma_BDT * np.sqrt(delta)
            r_next_j = np.exp(z_next_j)
            price += q_next[j] * np.exp(-r_next_j * delta) * 100
        return price

    theta_BDT[i] = brentq(
        lambda th: bond_price_given_theta_BDT(th) - target, -5.0, 10.0
    )

    z_bottom_next = z_tree_BDT[i][i] + theta_BDT[i] * delta - sigma_BDT * np.sqrt(delta)
    for j in range(i + 2):
        z_tree_BDT[i + 1][j] = z_bottom_next + (i + 1 - j) * 2 * sigma_BDT * np.sqrt(delta)
        r_tree_BDT[i + 1][j] = np.exp(z_tree_BDT[i + 1][j])

    for j in range(i + 2):
        val = 0.0
        if j <= i:
            r_ij = np.exp(z_tree_BDT[i][i] + (i - j) * 2 * sigma_BDT * np.sqrt(delta))
            val += 0.5 * q_BDT[i][j] * np.exp(-r_ij * delta)
        if j >= 1 and (j - 1) <= i:
            r_ij1 = np.exp(z_tree_BDT[i][i] + (i - (j - 1)) * 2 * sigma_BDT * np.sqrt(delta))
            val += 0.5 * q_BDT[i][j - 1] * np.exp(-r_ij1 * delta)
        q_BDT[i + 1][j] = val

# Fill final column
z_bottom_N_BDT = z_tree_BDT[N - 1][N - 1] + theta_BDT[N - 2] * delta - sigma_BDT * np.sqrt(delta)
for j in range(N + 1):
    z_tree_BDT[N][j] = z_bottom_N_BDT + (N - j) * 2 * sigma_BDT * np.sqrt(delta)
    r_tree_BDT[N][j] = np.exp(z_tree_BDT[N][j])

# --- Print BDT tree ----------------------------------------------------------
print()
print("=" * (8 + W * (N + 1)))
print("Table: The Risk Neutral BDT Interest Rate Tree")
print(f"BoE Nominal Spot Curve, 30 Jan 2026 | sigma = {sigma_BDT*100:.4f}% | delta = {delta} | N = {N}")
print("=" * (8 + W * (N + 1)))

print(f"{'Time T':>8}", end="")
for i in range(N + 1):
    print(f"{i * delta:>{W}.1f}", end="")
print()

print(f"{'Period i':>8}", end="")
for i in range(N + 1):
    print(f"{i:>{W}d}", end="")
print()

print(f"{'Î¸i(Ã100)':>8}", end="")
print(f"{'':>{W}}", end="")
for i in range(N):
    if i <= N - 2:
        print(f"{theta_BDT[i] * 100:>{W}.4f}", end="")
    else:
        print(f"{'':>{W}}", end="")
print()
print()

for j in range(N + 1):
    print(f"{'j='+str(j):>8}", end="")
    for i in range(N + 1):
        if j <= i:
            print(f"{r_tree_BDT[i][j] * 100:>{W}.2f}", end="")
        else:
            print(f"{'':>{W}}", end="")
    print()

print()
print("=" * 55)
print("BDT Calibration Check")
print("=" * 55)
print(f"{'k':>4}  {'Maturity':>10}  {'Market':>10}  {'Model':>10}  {'Diff':>10}")
for k in range(1, N + 1):
    maturity = k * delta
    market   = P[k]
    if k < N:
        model = sum(q_BDT[k][j] for j in range(k + 1)) * 100
    else:
        model = sum(
            q_BDT[N-1][j] * np.exp(-r_tree_BDT[N-1][j] * delta)
            for j in range(N)
        ) * 100
    diff = model - market
    flag = "  *last step" if k == N else ""
    print(f"{k:>4}  {maturity:>10.1f}  {market:>10.4f}  {model:>10.4f}  {diff:>10.6f}{flag}")

print(f"\nMinimum rate in BDT tree: {np.min(r_tree_BDT[r_tree_BDT > 0]) * 100:.4f}%")
print(f"Maximum rate in BDT tree: {np.max(r_tree_BDT) * 100:.4f}%")

# -*- coding: utf-8 -*-


# =============================================================================
# FM405 SUMMATIVE - FILE 2: MORTGAGE VALUATION AND MBS PRICING
# Part (b): Par mortgage rate and mortgage valuation
# Part (c): MBS pricing — Pass-Through, Interest Only, Principal Only
# Part (d): Simulate interest rate paths on binomial tree, plot histograms
# Part (e): Value the mortgage via Monte Carlo and compare to tree
# Note: Run File 1 containing part A first to generate:
#   r_tree_HL, q_HL, r_tree_BDT, q_BDT, delta, N, P
# =============================================================================

import numpy as np
from scipy.optimize import brentq

# =============================================================================
# SHARED HELPER FUNCTIONS
# =============================================================================

def coupon_payment(c, delta, N):
    """
    Semi-annual coupon payment on £100 par mortgage.
    Converts continuously compounded rate c to periodic rate first.
    Matches lecture slide p.16 formula.
    """
    rp  = np.exp(c * delta) - 1
    return 100 * rp / (1 - (1 + rp)**(-N))


def amortisation_schedule(c, delta, N):
    """
    Compute full amortisation schedule.
    Matches Panel B, Table 4 (lecture slide p.21).
    Returns bal, interest, principal arrays and pmt.
    """
    rp        = np.exp(c * delta) - 1
    pmt       = coupon_payment(c, delta, N)
    bal       = np.zeros(N + 1)
    interest  = np.zeros(N)
    principal = np.zeros(N)
    bal[0]    = 100.0
    for i in range(N):
        interest[i]  = bal[i] * rp
        principal[i] = pmt - interest[i]
        bal[i + 1]   = bal[i] - principal[i]
    bal[N] = 0.0
    return bal, interest, principal, pmt


def build_vnp_tree(c, r_tree, delta, N, pmt):
    """
    Build V^np (no-prepayment) tree by backward induction.
    Matches lecture slide p.18, Table 3.
    Terminal condition: V^np[N][j] = 0 (mortgage fully paid).
    """
    Vnp = np.zeros((N + 1, N + 1))
    for i in range(N - 1, -1, -1):
        for j in range(i + 1):
            r         = r_tree[i][j]
            disc      = np.exp(-r * delta)
            Vnp[i][j] = disc * (0.5 * Vnp[i+1][j] + 0.5 * Vnp[i+1][j+1] + pmt)
    return Vnp


def build_option_tree(Vnp, bal, r_tree, delta, N):
    """
    Build prepayment option C tree by backward induction.
    Matches lecture slides pp.20-21, Table 4 Panel A.

    At each node (i,j):
      C_ex   = max(V^np[i][j] - bal[i], 0)
      C_wait = disc * (0.5*C[i+1][j] + 0.5*C[i+1][j+1])
      C[i][j] = max(C_wait, C_ex)

    i=0: cannot prepay at origination
    i=N: option expires worthless
    """
    C = np.zeros((N + 1, N + 1))
    for i in range(N - 1, -1, -1):
        for j in range(i + 1):
            r      = r_tree[i][j]
            disc   = np.exp(-r * delta)
            C_wait = disc * (0.5 * C[i+1][j] + 0.5 * C[i+1][j+1])
            if i == 0:
                C[i][j] = C_wait
            else:
                C_ex    = max(Vnp[i][j] - bal[i], 0.0)
                C[i][j] = max(C_wait, C_ex)
    return C


def find_par_rate(r_tree, delta, N):
    """
    Find par mortgage rate c* such that V[0][0] = 100.
    Uses min(hold, bal) single-pass for root-finding to avoid
    circular dependency, then decomposes into V^np and C for reporting.
    """
    def mortgage_direct(c):
        pmt          = coupon_payment(c, delta, N)
        bal, _, _, _ = amortisation_schedule(c, delta, N)
        V            = np.zeros((N + 1, N + 1))
        for i in range(N - 1, -1, -1):
            for j in range(i + 1):
                r    = r_tree[i][j]
                disc = np.exp(-r * delta)
                hold = disc * (0.5 * V[i+1][j] + 0.5 * V[i+1][j+1] + pmt)
                V[i][j] = hold if i == 0 else min(hold, bal[i])
        return V[0][0]

    return brentq(lambda c: mortgage_direct(c) - 100.0, 0.001, 0.20)


# =============================================================================
# PRINT FUNCTIONS
# =============================================================================

def print_tree_table(tree, N, delta, title, fmt="{:>9.2f}"):
    W = 9
    print()
    print(title)
    print("-" * (6 + W * (N + 1)))
    print(f"{'t =':>6}", end="")
    for i in range(N + 1):
        print(f"{i * delta:>{W}.1f}", end="")
    print()
    print(f"{'i =':>6}", end="")
    for i in range(N + 1):
        print(f"{i:>{W}d}", end="")
    print()
    print(f"{'j =':>6}")
    for j in range(N + 1):
        print(f"{j:>6}", end="")
        for i in range(N + 1):
            if j <= i:
                print(fmt.format(tree[i][j]), end="")
            else:
                print(f"{'':>{W}}", end="")
        print()


def print_amortisation(bal, interest, principal, pmt, delta, N, title):
    W = 10
    print()
    print(title)
    print("-" * (24 + W * (N + 1)))
    print(f"{'t =':>24}", end="")
    for i in range(N + 1):
        print(f"{i * delta:>{W}.1f}", end="")
    print()
    print(f"{'i =':>24}", end="")
    for i in range(N + 1):
        print(f"{i:>{W}d}", end="")
    print()
    print(f"{'Interest Paid':>24}", end="")
    print(f"{'':>{W}}", end="")
    for v in interest:
        print(f"{v:>{W}.4f}", end="")
    print()
    print(f"{'Principal Paid':>24}", end="")
    print(f"{'':>{W}}", end="")
    for v in principal:
        print(f"{v:>{W}.4f}", end="")
    print()
    print(f"{'Outstanding Principal':>24}", end="")
    for v in bal:
        print(f"{v:>{W}.4f}", end="")
    print()


# =============================================================================
# PART (b): PAR MORTGAGE RATE AND VALUATION
# Process:
#   Step 0: Root-find c* using min(hold,bal)
#   Step 1: Amortisation schedule (Panel B, Table 4)
#   Step 2: V^np tree (Table 3)
#   Step 3: Prepayment option C tree (Table 4 Panel A)
#   Step 4: V = V^np - C
#   Step 5: Verify V[0][0] = 100
# =============================================================================

print("=" * 70)
print("PART (b): PAR MORTGAGE RATE AND VALUATION")
print("=" * 70)

results = {}

for label, r_tree in [("Ho-Lee", r_tree_HL), ("BDT", r_tree_BDT)]:

    # Step 0: find par rate
    c = find_par_rate(r_tree, delta, N)

    # Step 1: amortisation schedule
    bal, interest, principal, pmt = amortisation_schedule(c, delta, N)

    # Step 2: V^np tree
    Vnp = build_vnp_tree(c, r_tree, delta, N, pmt)

    # Step 3: prepayment option tree
    C = build_option_tree(Vnp, bal, r_tree, delta, N)

    # Step 4: mortgage value with prepayment
    V = Vnp - C

    results[label] = dict(c=c, pmt=pmt, bal=bal, interest=interest,
                          principal=principal, Vnp=Vnp, C=C, V=V)

    print()
    print("=" * 70)
    print(f"  MODEL: {label}")
    print("=" * 70)
    print(f"  Par mortgage rate c*         : {c*100:.4f}%  (ann., cont. comp.)")
    print(f"  Semi-annual payment per £100 : £{pmt:.4f}")
    print(f"  V^np[0][0]  (no prepayment)  : £{Vnp[0][0]:.6f}")
    print(f"  Option value  C[0][0]        : £{C[0][0]:.6f}")
    print(f"  Mortgage value V[0][0]       : £{V[0][0]:.6f}   (target = £100.00)")

    # Table 3: V^np tree
    print_tree_table(
        Vnp, N, delta,
        title=f"Table 3: Mortgage Value WITHOUT Prepayment — {label} (£ per £100)"
    )

    # Table 4 Panel A: prepayment option tree
    print_tree_table(
        C, N, delta,
        title=f"Table 4 Panel A: Prepayment Option Tree — {label} (£ per £100)"
    )

    # Mortgage value WITH prepayment
    print_tree_table(
        V, N, delta,
        title=f"Table: Mortgage Value WITH Prepayment — {label} (£ per £100)"
    )

    # Table 4 Panel B: amortisation schedule
    print_amortisation(
        bal, interest, principal, pmt, delta, N,
        title=f"Table 4 Panel B: Amortisation Schedule — {label}"
    )

    print()
    print(f"  Consistency: V^np - C = {Vnp[0][0]:.6f} - {C[0][0]:.6f}"
          f" = {Vnp[0][0] - C[0][0]:.6f}  (= V[0][0] = {V[0][0]:.6f})")


# =============================================================================
# PART (c): MBS PRICING
# Securities: Pass-Through (PT), Interest Only (IO), Principal Only (PO)
# Pass-through rate = par mortgage rate - 50bp for each model
#
# Prepayment decision: driven by mortgage tree from part (b)
# The homeowner prepays based on their mortgage liability at c*
# MBS investors passively receive resulting cash flows
# =============================================================================

print()
print("=" * 70)
print("PART (c): MBS PRICING")
print("=" * 70)

# Extract results from part (b)
c_HL    = results['Ho-Lee']['c']
c_BDT   = results['BDT']['c']
bal_HL  = results['Ho-Lee']['bal']
bal_BDT = results['BDT']['bal']
V_HL    = results['Ho-Lee']['V']
V_BDT   = results['BDT']['V']

# Pass-through rates: 50bp below par mortgage rate
c_PT_HL  = c_HL  - 0.0050
c_PT_BDT = c_BDT - 0.0050

print(f"\n  Ho-Lee  par mortgage rate : {c_HL*100:.4f}%")
print(f"  Ho-Lee  pass-through rate : {c_PT_HL*100:.4f}%  (= par - 50bp)")
print(f"\n  BDT     par mortgage rate : {c_BDT*100:.4f}%")
print(f"  BDT     pass-through rate : {c_PT_BDT*100:.4f}%  (= par - 50bp)")


def build_mbs_trees(c_mort, c_pt, r_tree, bal, V_mort, delta, N):
    """
    Build PT, IO and PO trees following lecture slides pp.24-31.

    Prepayment condition: re-derived from the hold value at each node.
    A node prepays when the present value of continuing (hold) >= outstanding
    balance, which is the optimal exercise condition from part (b).
    This is more robust than checking V_mort == bal due to floating point.

    If PREPAYMENT at node (i,j):
      PT = bal[i]   (outstanding principal returned to PT investor)
      PO = bal[i]   (all principal returned to PO investor)
      IO = 0        (interest payments stop — IO gets nothing)

    If NO PREPAYMENT at node (i,j):
      PT: CF = PT interest (at c_pt) + scheduled principal
      PO: CF = scheduled principal only
      IO: CF = PT interest only (at c_pt)
    """
    # Amortisation at mortgage rate c_mort
    rp_mort         = np.exp(c_mort * delta) - 1
    pmt_mort        = coupon_payment(c_mort, delta, N)
    principal_sched = np.array([pmt_mort - bal[i] * rp_mort for i in range(N)])

    # PT interest payments at pass-through rate c_pt
    rp_pt       = np.exp(c_pt * delta) - 1
    interest_pt = np.array([bal[i] * rp_pt for i in range(N)])

    # Cash flows per period
    CF_PT = interest_pt + principal_sched   # PT: interest + principal
    CF_PO = principal_sched.copy()          # PO: principal only
    CF_IO = interest_pt.copy()              # IO: interest only

    # Detect prepayment nodes by recomputing hold at each node
    # Prepayment occurs when hold >= bal[i] (borrower exercises refinancing option)
    prepay = np.zeros((N + 1, N + 1), dtype=bool)
    for i in range(1, N):
        for j in range(i + 1):
            hold = np.exp(-r_tree[i][j] * delta) * (
                0.5 * V_mort[i+1][j] + 0.5 * V_mort[i+1][j+1] + pmt_mort
            )
            if hold >= bal[i] - 1e-8:
                prepay[i][j] = True

    # Build PT, PO, IO trees by backward induction
    PT = np.zeros((N + 1, N + 1))
    PO = np.zeros((N + 1, N + 1))
    IO = np.zeros((N + 1, N + 1))

    for i in range(N - 1, -1, -1):
        for j in range(i + 1):
            r    = r_tree[i][j]
            disc = np.exp(-r * delta)

            if prepay[i][j]:
                PT[i][j] = bal[i]
                PO[i][j] = bal[i]
                IO[i][j] = 0.0
            else:
                PT[i][j] = disc * (0.5 * PT[i+1][j] + 0.5 * PT[i+1][j+1] + CF_PT[i])
                PO[i][j] = disc * (0.5 * PO[i+1][j] + 0.5 * PO[i+1][j+1] + CF_PO[i])
                IO[i][j] = disc * (0.5 * IO[i+1][j] + 0.5 * IO[i+1][j+1] + CF_IO[i])

    return PT, PO, IO, CF_PT, CF_PO, CF_IO, interest_pt, principal_sched


def spot_rate_duration(tree, r_tree):
    """
    Spot rate duration using one-step finite difference from root.
    Matches lecture slide p.31 formula:
      D ≈ -1/P0 * (P[1,u] - P[1,d]) / (r[1,u] - r[1,d])
    """
    P0  = tree[0][0]
    P1u = tree[1][0]
    P1d = tree[1][1]
    r1u = r_tree[1][0]
    r1d = r_tree[1][1]
    if abs(P0) < 1e-10 or abs(r1u - r1d) < 1e-10:
        return float('nan')
    return -1/P0 * (P1u - P1d) / (r1u - r1d)


def print_cf_schedule(int_pt, prin, bal, delta, N, title):
    W = 10
    print()
    print(title)
    print("-" * (26 + W * (N + 1)))
    print(f"{'t =':>26}", end="")
    for i in range(N + 1):
        print(f"{i * delta:>{W}.1f}", end="")
    print()
    print(f"{'i =':>26}", end="")
    for i in range(N + 1):
        print(f"{i:>{W}d}", end="")
    print()
    print(f"{'PT Interest (pass-through)':>26}", end="")
    print(f"{'':>{W}}", end="")
    for v in int_pt:
        print(f"{v:>{W}.4f}", end="")
    print()
    print(f"{'Scheduled Principal':>26}", end="")
    print(f"{'':>{W}}", end="")
    for v in prin:
        print(f"{v:>{W}.4f}", end="")
    print()
    print(f"{'Outstanding Principal':>26}", end="")
    for v in bal:
        print(f"{v:>{W}.4f}", end="")
    print()


# Run for both models
mbs_results = {}

for label, c_mort, c_pt, r_tree, bal, V_mort in [
    ("Ho-Lee", c_HL,  c_PT_HL,  r_tree_HL,  bal_HL,  V_HL),
    ("BDT",    c_BDT, c_PT_BDT, r_tree_BDT, bal_BDT, V_BDT),
]:
    PT, PO, IO, CF_PT, CF_PO, CF_IO, int_pt, prin = build_mbs_trees(
        c_mort, c_pt, r_tree, bal, V_mort, delta, N
    )

    D_PT = spot_rate_duration(PT, r_tree)
    D_PO = spot_rate_duration(PO, r_tree)
    D_IO = spot_rate_duration(IO, r_tree)

    mbs_results[label] = dict(
        PT=PT, PO=PO, IO=IO,
        int_pt=int_pt, prin=prin,
        D_PT=D_PT, D_PO=D_PO, D_IO=D_IO,
        c_mort=c_mort, c_pt=c_pt
    )

    print()
    print("=" * 70)
    print(f"  PART (c) — {label} MODEL")
    print("=" * 70)
    print(f"  Mortgage rate c*          : {c_mort*100:.4f}%")
    print(f"  Pass-through rate c_PT    : {c_pt*100:.4f}%  (= c* - 50bp)")
    print()
    print(f"  {'Security':<10}  {'Price (£)':>12}  {'Duration (yrs)':>16}")
    print(f"  {'-'*10}  {'-'*12}  {'-'*16}")
    print(f"  {'PT':<10}  {PT[0][0]:>12.6f}  {D_PT:>16.4f}")
    print(f"  {'PO':<10}  {PO[0][0]:>12.6f}  {D_PO:>16.4f}")
    print(f"  {'IO':<10}  {IO[0][0]:>12.6f}  {D_IO:>16.4f}")
    print()
    print(f"  Consistency: PO + IO = {PO[0][0]:.6f} + {IO[0][0]:.6f}"
          f" = {PO[0][0]+IO[0][0]:.6f}  (should = PT = {PT[0][0]:.6f})")

    # Table 6 Panel A: PT tree
    print_tree_table(
        PT, N, delta,
        title=f"Table 6 Panel A: Pass-Through Value Tree — {label} (£ per £100)"
    )

    # Table 7 Panel A: PO tree
    print_tree_table(
        PO, N, delta,
        title=f"Table 7 Panel A: Principal Only (PO) Value Tree — {label} (£ per £100)"
    )

    # Table 7 Panel B: IO tree
    print_tree_table(
        IO, N, delta,
        title=f"Table 7 Panel B: Interest Only (IO) Value Tree — {label} (£ per £100)"
    )

    # Cash flow schedule
    print_cf_schedule(
        int_pt, prin,
        bal_HL if label == "Ho-Lee" else bal_BDT,
        delta, N,
        title=f"Panel B: Cash Flow Schedule — {label}"
    )

# =============================================================================
# COMPARISON TABLE: BDT vs Ho-Lee
# =============================================================================
print()
print("=" * 70)
print("  COMPARISON: BDT vs Ho-Lee")
print("=" * 70)
print(f"  {'Security':<10}  {'Ho-Lee Price':>14}  {'BDT Price':>12}  "
      f"{'Difference':>12}  {'HL Duration':>12}  {'BDT Duration':>13}")
print(f"  {'-'*10}  {'-'*14}  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*13}")

for sec in ['PT', 'PO', 'IO']:
    hl_p  = mbs_results['Ho-Lee'][sec][0][0]
    bdt_p = mbs_results['BDT'][sec][0][0]
    hl_d  = mbs_results['Ho-Lee'][f'D_{sec}']
    bdt_d = mbs_results['BDT'][f'D_{sec}']
    diff  = bdt_p - hl_p
    print(f"  {sec:<10}  {hl_p:>14.6f}  {bdt_p:>12.6f}  "
          f"{diff:>12.6f}  {hl_d:>12.4f}  {bdt_d:>13.4f}")

# =============================================================================
# FM405 SUMMATIVE - FILE 3: MONTE CARLO SIMULATION
# Part (d): Simulate interest rate paths on binomial tree, plot histograms
# Part (e): Value the mortgage via Monte Carlo and compare to tree
#
# METHODOLOGY:
#   MC simulates paths DIRECTLY on the calibrated binomial tree.
#   At each period i, draw RAND() ~ Uniform[0,1]:
#     If RAND() < 0.5  -> up move   -> node j unchanged   (r_{i+1, j})
#     If RAND() >= 0.5 -> down move -> node j increases 1 (r_{i+1, j+1})
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# PARAMETERS
# =============================================================================

N_sim = 100000
seed  = 42
rng   = np.random.default_rng(seed)

c_HL    = results['Ho-Lee']['c']
c_BDT   = results['BDT']['c']
bal_HL  = results['Ho-Lee']['bal']
bal_BDT = results['BDT']['bal']

# =============================================================================
# PART (d): SIMULATE INTEREST RATE PATHS ON BINOMIAL TREE
# =============================================================================
#   Draw RAND() ~ U[0,1] at each step.
#   RAND() < 0.5  -> up move   -> j unchanged   -> r_{i+1, j}
#   RAND() >= 0.5 -> down move -> j increases 1 -> r_{i+1, j+1}
#
# Rates are read directly from the calibrated tree arrays.
# =============================================================================

print("=" * 70)
print("PART (d): MONTE CARLO INTEREST RATE SIMULATION ON BINOMIAL TREE")
print("=" * 70)
print()
print("  Methodology: simulate paths directly on calibrated binomial tree")
print("  RAND() < 0.5  -> up move (j stays);  RAND() >= 0.5 -> down move (j+1)")

# Draw uniform random variables for part (d) visualisation only
# A separate RNG stream is used for part (e) to keep the two analyses
# statistically independent (Issue D-1)
U_d        = rng.uniform(size=(N_sim, N))
up_moves_d = U_d < 0.5   # True = up, False = down


def simulate_tree_paths(r_tree, up_moves, N_sim, N):
    """
    Simulate N_sim paths directly on the binomial tree.
    At each step i, track node index j and read rate from r_tree[i][j].
    Returns paths array (N_sim, N+1) of actual tree node rates.
    """
    paths    = np.zeros((N_sim, N + 1))
    j_paths  = np.zeros(N_sim, dtype=int)
    paths[:, 0] = r_tree[0][0]

    for i in range(N):
        j_paths = np.where(up_moves[:, i], j_paths, j_paths + 1)
        paths[:, i + 1] = r_tree[i + 1][j_paths]

    return paths


paths_HL  = simulate_tree_paths(r_tree_HL,  up_moves_d, N_sim, N)
paths_BDT = simulate_tree_paths(r_tree_BDT, up_moves_d, N_sim, N)

# --- Terminal distributions ---
r_T_HL  = paths_HL[:, N]  * 100
r_T_BDT = paths_BDT[:, N] * 100

print(f"\n  Ho-Lee terminal rates at T = {N*delta:.0f} years:")
print(f"    Mean        : {r_T_HL.mean():.4f}%")
print(f"    Std dev     : {r_T_HL.std():.4f}%")
print(f"    Min         : {r_T_HL.min():.4f}%")
print(f"    Max         : {r_T_HL.max():.4f}%")
print(f"    % negative  : {(r_T_HL < 0).mean()*100:.2f}%")

print(f"\n  BDT terminal rates at T = {N*delta:.0f} years:")
print(f"    Mean        : {r_T_BDT.mean():.4f}%")
print(f"    Std dev     : {r_T_BDT.std():.4f}%")
print(f"    Min         : {r_T_BDT.min():.4f}%")
print(f"    Max         : {r_T_BDT.max():.4f}%")
print(f"    % negative  : {(r_T_BDT < 0).mean()*100:.2f}%")

# --- Verify MC matches tree: ZCB prices (lecture notes slide 70) ---
print(f"\n  Verification: MC zero-coupon bond prices vs market")
print(f"  {'Maturity':>10}  {'Market':>10}  {'MC (HL)':>10}  {'MC (BDT)':>10}")
for k in [2, 4, 6, 8, 10, 14, 18, 20]:
    if k > N:
        break
    market = P[k]
    mc_hl  = np.mean(np.exp(-np.sum(paths_HL[:,  :k] * delta, axis=1))) * 100
    mc_bdt = np.mean(np.exp(-np.sum(paths_BDT[:, :k] * delta, axis=1))) * 100
    print(f"  {k*delta:>10.1f}  {market:>10.4f}  {mc_hl:>10.4f}  {mc_bdt:>10.4f}")

# --- Figure 1: Histograms ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle(
    f"Figure 1: Distribution of Short Rate at T = {N*delta:.0f} Years\n"
    f"N = {N_sim:,} Simulations on Binomial Tree | BoE Spot Curve, 30 Jan 2026",
    fontsize=12
)

axes[0].hist(r_T_HL, bins=100, color='steelblue', edgecolor='none',
             alpha=0.8, density=True)
axes[0].axvline(r_T_HL.mean(), color='navy', linestyle='--', linewidth=1.5,
                label=f'Mean = {r_T_HL.mean():.2f}%')
axes[0].axvline(0, color='red', linestyle=':', linewidth=1.2, label='r = 0')
axes[0].set_title('Ho-Lee Model (Normal / Additive)')
axes[0].set_xlabel('Short Rate (%)')
axes[0].set_ylabel('Density')
axes[0].legend()

axes[1].hist(r_T_BDT, bins=100, color='darkorange', edgecolor='none',
             alpha=0.8, density=True)
axes[1].axvline(r_T_BDT.mean(), color='darkred', linestyle='--', linewidth=1.5,
                label=f'Mean = {r_T_BDT.mean():.2f}%')
axes[1].set_title('BDT Model (Lognormal / Multiplicative)')
axes[1].set_xlabel('Short Rate (%)')
axes[1].set_ylabel('Density')
axes[1].legend()

plt.tight_layout()
plt.savefig('figure1_terminal_rate_histograms.png', dpi=150, bbox_inches='tight')
plt.show()
print("\n  Figure 1 saved.")

# --- Figure 2: Sample paths ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle(
    "Figure 2: Sample Interest Rate Paths (first 200 simulations)\n"
    "Paths traverse actual binomial tree nodes | BoE Spot Curve, 30 Jan 2026",
    fontsize=12
)

t_grid = np.arange(N + 1) * delta

for k in range(200):
    axes[0].plot(t_grid, paths_HL[k]  * 100, color='steelblue',
                 alpha=0.05, linewidth=0.5)
axes[0].plot(t_grid, paths_HL.mean(axis=0) * 100, color='navy',
             linewidth=2, label='Mean path')
axes[0].axhline(0, color='red', linestyle=':', linewidth=1)
axes[0].set_title('Ho-Lee Model')
axes[0].set_xlabel('Years')
axes[0].set_ylabel('Short Rate (%)')
axes[0].legend()

for k in range(200):
    axes[1].plot(t_grid, paths_BDT[k] * 100, color='darkorange',
                 alpha=0.05, linewidth=0.5)
axes[1].plot(t_grid, paths_BDT.mean(axis=0) * 100, color='darkred',
             linewidth=2, label='Mean path')
axes[1].set_title('BDT Model')
axes[1].set_xlabel('Years')
axes[1].set_ylabel('Short Rate (%)')
axes[1].legend()

plt.tight_layout()
plt.savefig('figure2_sample_paths.png', dpi=150, bbox_inches='tight')
plt.show()
print("  Figure 2 saved.")


# =============================================================================
# PART (e): MORTGAGE VALUATION VIA MONTE CARLO (BDT ONLY)
# =============================================================================
# Prepayment: since paths visit actual tree nodes, we read V[i][j] directly
#   from the mortgage value tree computed in part (b).
#   Prepay if V[i][j] <= bal[i] — exactly the condition from the tree.
#   This makes the MC and tree completely consistent on the same probability
#   space. Any remaining difference is pure sampling noise, measured by the CI.
# =============================================================================

print()
print("=" * 70)
print("PART (e): MORTGAGE VALUATION VIA MONTE CARLO (BDT)")
print("=" * 70)
print()
print("  Model: BDT")
print("  Prepayment: hold >= bal[i] recomputed from tree V values")
print("  Discounting: path-specific accumulated discount factors")
print("  Note: independent RNG stream from part (d) for statistical separation")

# Fresh RNG for part (e) — independent from part (d) draws (Issue D-1)
rng_e      = np.random.default_rng(seed + 1)
U_e        = rng_e.uniform(size=(N_sim, N))
up_moves   = U_e < 0.5


def mc_mortgage_tree(r_tree, V_mort, bal, c, delta, N, N_sim, up_moves):
    """
    Value the mortgage along simulated binomial tree paths.

    Cash flow timing (matching tree backward induction exactly):
      - Tree: V[i][j] = disc_i * (0.5*V[i+1,j] + 0.5*V[i+1,j+1] + pmt)
      - pmt is a cash flow received at the END of period i, discounted
        back by disc_i from period i+1 to period i, then further back
        to t=0 by cum_disc. Total discount: cum_disc * disc_i * pmt.
      - On prepayment at period i: V[i][j] = bal[i], which is a present
        value AT period i (not at i+1). The tree discounts it by disc_{i-1}
        (one step back to period i-1). In MC this corresponds to cum_disc,
        which is the product of disc factors for periods 0..i-1.
        Total discount: cum_disc * bal[i].
      - The two cases use different discount factors because they represent
        cash flows at different points in time: pmt arrives at end of period i
        while bal[i] is already expressed as a value at the start of period i.
        This is NOT an inconsistency — it correctly mirrors the tree.

    Prepayment range: checked for all 1 <= i <= N-1 (Issue E-2 fix).
    At i=N-1: V[N][j]=0 for all j, so hold = disc*(0.5*0 + 0.5*0 + pmt) = disc*pmt.
    """
    pmt      = coupon_payment(c, delta, N)
    pv_paths = np.zeros(N_sim)

    for sim in range(N_sim):
        j        = 0
        cum_disc = 1.0
        pv       = 0.0

        for i in range(N):
            r    = r_tree[i][j]
            disc = np.exp(-r * delta)

            # Prepayment check: all interior periods 1 <= i <= N-1
            # At i=N-1: V[N][j]=0, hold = disc*pmt
            if i > 0:
                j_next = j      # up move node
                j_dn   = j + 1  # down move node (capped at N for safety)
                v_up   = V_mort[i+1][j_next] if i+1 <= N else 0.0
                v_dn   = V_mort[i+1][j_dn]   if i+1 <= N else 0.0
                hold   = disc * (0.5 * v_up + 0.5 * v_dn + pmt)
                if hold >= bal[i] - 1e-8:
                    # bal[i] is value AT period i -> discount by cum_disc only
                    pv += cum_disc * bal[i]
                    break

            # No prepayment: pmt paid at END of period i
            pv += cum_disc * disc * pmt
            cum_disc *= disc

            # Advance node: up = j stays, down = j+1
            if not up_moves[sim, i]:
                j += 1

        pv_paths[sim] = pv

    return pv_paths


# --- BDT ---
print("\n  Running Monte Carlo valuation")
pv_BDT = mc_mortgage_tree(
    r_tree_BDT, results['BDT']['V'],
    bal_BDT, c_BDT, delta, N, N_sim, up_moves
)

mc_mean_BDT = pv_BDT.mean()
mc_se_BDT   = pv_BDT.std() / np.sqrt(N_sim)
ci_lo_BDT   = mc_mean_BDT - 1.96 * mc_se_BDT
ci_hi_BDT   = mc_mean_BDT + 1.96 * mc_se_BDT
tree_BDT    = results['BDT']['V'][0][0]

print(f"\n  BDT Results:")
print(f"  {'N simulations':<30}: {N_sim:,}")
print(f"  {'Point estimate':<30}: £{mc_mean_BDT:.6f}")
print(f"  {'Standard error':<30}: £{mc_se_BDT:.6f}")
print(f"  {'95% CI':<30}: [£{ci_lo_BDT:.6f}, £{ci_hi_BDT:.6f}]")
print(f"  {'Tree value from part (b)':<30}: £{tree_BDT:.6f}")
print(f"  {'Difference (MC - Tree)':<30}: £{mc_mean_BDT - tree_BDT:.6f}")
in_ci = ci_lo_BDT <= tree_BDT <= ci_hi_BDT
print(f"  {'Tree value in 95% CI':<30}: {'YES' if in_ci else 'NO'}")

# --- Convergence table ---
print()
print("  Convergence of MC estimate with number of paths (BDT):")
print(f"  {'N paths':<12}  {'Mean (£)':>10}  {'Std Error':>12}  {'95% CI':>28}")
for cp in [100, 500, 1000, 5000, 10000, 50000, 100000]:
    if cp > N_sim:
        break
    m  = pv_BDT[:cp].mean()
    se = pv_BDT[:cp].std() / np.sqrt(cp)
    print(f"  {cp:<12,}  {m:>10.6f}  {se:>12.6f}"
          f"  [{m-1.96*se:>10.6f}, {m+1.96*se:>10.6f}]")

# -*- coding: utf-8 -*-


# =============================================================================
# FM405 SUMMATIVE - FILE 1: INTEREST RATE TREE CALIBRATION
# Parts (a): Ho-Lee and BDT Interest Rate Tree Calibration
# Data: Bank of England Nominal Spot Curve, 30 January 2026
# Note: Run all files on one document since later parts require outputs from
# earlier parts
# =============================================================================

import numpy as np
from scipy.optimize import brentq

# =============================================================================
# MARKET DATA
# =============================================================================

delta = 0.5   # semi-annual steps
N     = 20    # 20 steps = 10 years

# Bank of England nominal spot rates (%), 30 January 2026
spot_rates = {
    0.5:  3.48, 1.0:  3.55, 1.5:  3.59, 2.0:  3.63, 2.5:  3.67,
    3.0:  3.72, 3.5:  3.78, 4.0:  3.84, 4.5:  3.90, 5.0:  3.97,
    5.5:  4.04, 6.0:  4.11, 6.5:  4.18, 7.0:  4.25, 7.5:  4.31,
    8.0:  4.38, 8.5:  4.45, 9.0:  4.51, 9.5:  4.57, 10.0: 4.63
}

# Zero coupon bond prices per Â£100 face value
P = {}
for k, (maturity, rate) in enumerate(spot_rates.items(), start=1):
    P[k] = 100 * np.exp(-(rate / 100) * maturity)

# =============================================================================
# HO-LEE CALIBRATION
# Normal (additive) model: r_{i+1,j} = r_{i,j} + theta_i*delta +/- sigma*sqrt(delta)
# =============================================================================

sigma_HL = 68.02 / 10000    # 68.02 bp/yr normal vol -> 0.006802

r0_HL       = -np.log(P[1] / 100) / delta
r_tree_HL   = np.zeros((N + 1, N + 1))
r_tree_HL[0][0] = r0_HL
q_HL        = np.zeros((N + 1, N + 1))
q_HL[0][0]  = 1.0
theta_HL    = np.zeros(N)

for i in range(N - 1):
    target = P[i + 2]

    def bond_price_given_theta_HL(theta_i, i=i):
        bottom_next = r_tree_HL[i][i] + theta_i * delta - sigma_HL * np.sqrt(delta)
        q_next = np.zeros(i + 2)
        for j in range(i + 2):
            val = 0.0
            if j <= i:
                r_ij = r_tree_HL[i][i] + (i - j) * 2 * sigma_HL * np.sqrt(delta)
                val += 0.5 * q_HL[i][j] * np.exp(-r_ij * delta)
            if j >= 1 and (j - 1) <= i:
                r_ij1 = r_tree_HL[i][i] + (i - (j - 1)) * 2 * sigma_HL * np.sqrt(delta)
                val += 0.5 * q_HL[i][j - 1] * np.exp(-r_ij1 * delta)
            q_next[j] = val
        price = 0.0
        for j in range(i + 2):
            r_next_j = bottom_next + (i + 1 - j) * 2 * sigma_HL * np.sqrt(delta)
            price += q_next[j] * np.exp(-r_next_j * delta) * 100
        return price

    theta_HL[i] = brentq(
        lambda th: bond_price_given_theta_HL(th) - target, -2.0, 5.0
    )

    bottom_next = r_tree_HL[i][i] + theta_HL[i] * delta - sigma_HL * np.sqrt(delta)
    for j in range(i + 2):
        r_tree_HL[i + 1][j] = bottom_next + (i + 1 - j) * 2 * sigma_HL * np.sqrt(delta)

    for j in range(i + 2):
        val = 0.0
        if j <= i:
            r_ij = r_tree_HL[i][i] + (i - j) * 2 * sigma_HL * np.sqrt(delta)
            val += 0.5 * q_HL[i][j] * np.exp(-r_ij * delta)
        if j >= 1 and (j - 1) <= i:
            r_ij1 = r_tree_HL[i][i] + (i - (j - 1)) * 2 * sigma_HL * np.sqrt(delta)
            val += 0.5 * q_HL[i][j - 1] * np.exp(-r_ij1 * delta)
        q_HL[i + 1][j] = val

# Fill final column
bottom_N_HL = r_tree_HL[N - 1][N - 1] + theta_HL[N - 2] * delta - sigma_HL * np.sqrt(delta)
for j in range(N + 1):
    r_tree_HL[N][j] = bottom_N_HL + (N - j) * 2 * sigma_HL * np.sqrt(delta)

# --- Print Ho-Lee tree -------------------------------------------------------
W = 7
print("=" * (8 + W * (N + 1)))
print("Table: The Risk Neutral Ho-Lee Interest Rate Tree")
print(f"BoE Nominal Spot Curve, 30 Jan 2026 | sigma = {sigma_HL*100:.4f}% | delta = {delta} | N = {N}")
print("=" * (8 + W * (N + 1)))

print(f"{'Time T':>8}", end="")
for i in range(N + 1):
    print(f"{i * delta:>{W}.1f}", end="")
print()

print(f"{'Period i':>8}", end="")
for i in range(N + 1):
    print(f"{i:>{W}d}", end="")
print()

print(f"{'Î¸i(Ã100)':>8}", end="")
print(f"{'':>{W}}", end="")
for i in range(N):
    if i <= N - 2:
        print(f"{theta_HL[i] * 100:>{W}.4f}", end="")
    else:
        print(f"{'':>{W}}", end="")
print()
print()

for j in range(N + 1):
    print(f"{'j='+str(j):>8}", end="")
    for i in range(N + 1):
        if j <= i:
            print(f"{r_tree_HL[i][j] * 100:>{W}.2f}", end="")
        else:
            print(f"{'':>{W}}", end="")
    print()

print()
print("=" * 55)
print("Ho-Lee Calibration Check")
print("=" * 55)
print(f"{'k':>4}  {'Maturity':>10}  {'Market':>10}  {'Model':>10}  {'Diff':>10}")
for k in range(1, N + 1):
    maturity = k * delta
    market   = P[k]
    if k < N:
        model = sum(q_HL[k][j] for j in range(k + 1)) * 100
    else:
        model = sum(
            q_HL[N-1][j] * np.exp(-r_tree_HL[N-1][j] * delta)
            for j in range(N)
        ) * 100
    diff = model - market
    flag = "  *last step" if k == N else ""
    print(f"{k:>4}  {maturity:>10.1f}  {market:>10.4f}  {model:>10.4f}  {diff:>10.6f}{flag}")

# =============================================================================
# BDT CALIBRATION
# Lognormal model: log(r_{i+1,j}) = log(r_{i,j}) + theta_i*delta +/- sigma*sqrt(delta)
# =============================================================================

sigma_BDT = 16.23 / 100     # 16.23% Black vol -> 0.1623

r0_BDT      = -np.log(P[1] / 100) / delta
z0_BDT      = np.log(r0_BDT)

z_tree_BDT  = np.zeros((N + 1, N + 1))
r_tree_BDT  = np.zeros((N + 1, N + 1))
z_tree_BDT[0][0] = z0_BDT
r_tree_BDT[0][0] = r0_BDT

q_BDT       = np.zeros((N + 1, N + 1))
q_BDT[0][0] = 1.0
theta_BDT   = np.zeros(N)

for i in range(N - 1):
    target = P[i + 2]

    def bond_price_given_theta_BDT(theta_i, i=i):
        z_bottom_next = z_tree_BDT[i][i] + theta_i * delta - sigma_BDT * np.sqrt(delta)
        q_next = np.zeros(i + 2)
        for j in range(i + 2):
            val = 0.0
            if j <= i:
                r_ij = np.exp(z_tree_BDT[i][i] + (i - j) * 2 * sigma_BDT * np.sqrt(delta))
                val += 0.5 * q_BDT[i][j] * np.exp(-r_ij * delta)
            if j >= 1 and (j - 1) <= i:
                r_ij1 = np.exp(z_tree_BDT[i][i] + (i - (j - 1)) * 2 * sigma_BDT * np.sqrt(delta))
                val += 0.5 * q_BDT[i][j - 1] * np.exp(-r_ij1 * delta)
            q_next[j] = val
        price = 0.0
        for j in range(i + 2):
            z_next_j = z_bottom_next + (i + 1 - j) * 2 * sigma_BDT * np.sqrt(delta)
            r_next_j = np.exp(z_next_j)
            price += q_next[j] * np.exp(-r_next_j * delta) * 100
        return price

    theta_BDT[i] = brentq(
        lambda th: bond_price_given_theta_BDT(th) - target, -5.0, 10.0
    )

    z_bottom_next = z_tree_BDT[i][i] + theta_BDT[i] * delta - sigma_BDT * np.sqrt(delta)
    for j in range(i + 2):
        z_tree_BDT[i + 1][j] = z_bottom_next + (i + 1 - j) * 2 * sigma_BDT * np.sqrt(delta)
        r_tree_BDT[i + 1][j] = np.exp(z_tree_BDT[i + 1][j])

    for j in range(i + 2):
        val = 0.0
        if j <= i:
            r_ij = np.exp(z_tree_BDT[i][i] + (i - j) * 2 * sigma_BDT * np.sqrt(delta))
            val += 0.5 * q_BDT[i][j] * np.exp(-r_ij * delta)
        if j >= 1 and (j - 1) <= i:
            r_ij1 = np.exp(z_tree_BDT[i][i] + (i - (j - 1)) * 2 * sigma_BDT * np.sqrt(delta))
            val += 0.5 * q_BDT[i][j - 1] * np.exp(-r_ij1 * delta)
        q_BDT[i + 1][j] = val

# Fill final column
z_bottom_N_BDT = z_tree_BDT[N - 1][N - 1] + theta_BDT[N - 2] * delta - sigma_BDT * np.sqrt(delta)
for j in range(N + 1):
    z_tree_BDT[N][j] = z_bottom_N_BDT + (N - j) * 2 * sigma_BDT * np.sqrt(delta)
    r_tree_BDT[N][j] = np.exp(z_tree_BDT[N][j])

# --- Print BDT tree ----------------------------------------------------------
print()
print("=" * (8 + W * (N + 1)))
print("Table: The Risk Neutral BDT Interest Rate Tree")
print(f"BoE Nominal Spot Curve, 30 Jan 2026 | sigma = {sigma_BDT*100:.4f}% | delta = {delta} | N = {N}")
print("=" * (8 + W * (N + 1)))

print(f"{'Time T':>8}", end="")
for i in range(N + 1):
    print(f"{i * delta:>{W}.1f}", end="")
print()

print(f"{'Period i':>8}", end="")
for i in range(N + 1):
    print(f"{i:>{W}d}", end="")
print()

print(f"{'Î¸i(Ã100)':>8}", end="")
print(f"{'':>{W}}", end="")
for i in range(N):
    if i <= N - 2:
        print(f"{theta_BDT[i] * 100:>{W}.4f}", end="")
    else:
        print(f"{'':>{W}}", end="")
print()
print()

for j in range(N + 1):
    print(f"{'j='+str(j):>8}", end="")
    for i in range(N + 1):
        if j <= i:
            print(f"{r_tree_BDT[i][j] * 100:>{W}.2f}", end="")
        else:
            print(f"{'':>{W}}", end="")
    print()

print()
print("=" * 55)
print("BDT Calibration Check")
print("=" * 55)
print(f"{'k':>4}  {'Maturity':>10}  {'Market':>10}  {'Model':>10}  {'Diff':>10}")
for k in range(1, N + 1):
    maturity = k * delta
    market   = P[k]
    if k < N:
        model = sum(q_BDT[k][j] for j in range(k + 1)) * 100
    else:
        model = sum(
            q_BDT[N-1][j] * np.exp(-r_tree_BDT[N-1][j] * delta)
            for j in range(N)
        ) * 100
    diff = model - market
    flag = "  *last step" if k == N else ""
    print(f"{k:>4}  {maturity:>10.1f}  {market:>10.4f}  {model:>10.4f}  {diff:>10.6f}{flag}")

print(f"\nMinimum rate in BDT tree: {np.min(r_tree_BDT[r_tree_BDT > 0]) * 100:.4f}%")
print(f"Maximum rate in BDT tree: {np.max(r_tree_BDT) * 100:.4f}%")
