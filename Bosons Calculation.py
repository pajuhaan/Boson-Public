#!/usr/bin/env python3
"""
Relator Z4 mass triad (W, Z, H)
==============================

This script computes, in a fully deterministic way:

1) The closed Z4 parent scale M4 from the Relator Z_N ladder (N=4).
2) The "minimal" (main-text) mass triad (MW, MZ, MH).
3) The "refined" (appendix) overlap-corrected triad (MW, MZ, MH) and Z-pole proxies.
4) A clean comparison table vs PDG reference values, including percent errors.

All numerical inputs are declared up-front.
No fitted parameters are introduced in either model.

Run:
    python relator_z4_triads.py
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Optional

# -----------------------------
# 0) Numerical inputs (EDIT HERE)
# -----------------------------

# Fundamental constants used in the ladder (CODATA-style inputs; set to the values you use in the paper)
ALPHA_FS = 1.0 / 137.035999084  # fine-structure constant alpha
ME_GEV = 0.00051099895000       # electron rest energy m_e c^2 in GeV

# PDG reference values (EDIT to your preferred PDG edition/rounding)
PDG_MASSES_GEV: Dict[str, float] = {
    "W": 80.3692,
    "Z": 91.1880,
    "H": 125.20,
}

# Optional PDG electroweak targets (if you want error columns for them)
# If you do not want these comparisons, set them to None.
PDG_SIN2_THETA_EFF: Optional[float] = None  # e.g. 0.23153 (update as needed)
PDG_GV_L: Optional[float] = None            # e.g. -0.03783 (update as needed)
PDG_GA_L: Optional[float] = None            # e.g. -0.50111 (update as needed)
PDG_A_L: Optional[float] = None             # e.g. 0.147  (update as needed)

# -----------------------------
# 1) Helper utilities
# -----------------------------

def rel_error_percent(pred: float, ref: float) -> float:
    """Relative error in percent; delta = 100*(pred/ref - 1)."""
    return 100.0 * (pred / ref - 1.0)

def fmt(x: float, nd: int = 10) -> str:
    """Consistent floating formatting."""
    return f"{x:.{nd}f}"

@dataclass(frozen=True)
class Triad:
    MW: float
    MZ: float
    MH: float

# -----------------------------
# 2) Relator Z4 ladder (closed)
# -----------------------------

def compute_Z4_parent_scale_M4(alpha_fs: float, me_gev: float) -> Dict[str, float]:
    """
    Computes the closed N=4 rung M4 from:

    K_parallel^(1)(1/4) = 409/(105*pi^2),
    K_parallel^(1)(1/2) = 0 (PV-null half-twist),
    S4 = 6 K_parallel^(1)(1/4),
    P_J^(4) = (N-1)/(3N) = 1/4,
    K_parallel^(4) = 6/pi + P_J^(4)*S4,
    sigma_A^(4) = 3 / K_parallel^(4),
    sigma_{A,e} = pi/2,
    M4 = m_e * exp( (sigma_{A,e} - sigma_A^(4)) / (4*alpha) ).
    """
    pi = math.pi

    # PV kernel values (closed)
    K1_a14 = 409.0 / (105.0 * pi**2)  # K_parallel^(1)(1/4)
    # K1_a12 = 0.0  # K_parallel^(1)(1/2), not needed explicitly

    # N=4
    N = 4.0
    P_J_4 = (N - 1.0) / (3.0 * N)      # = 1/4
    S4 = 6.0 * K1_a14
    K4 = 6.0 / pi + P_J_4 * S4

    sigmaA4 = 3.0 / K4
    sigmaA_e = pi / 2.0

    M4 = me_gev * math.exp((sigmaA_e - sigmaA4) / (4.0 * alpha_fs))

    return {
        "K1_a14": K1_a14,
        "S4": S4,
        "P_J_4": P_J_4,
        "K4": K4,
        "sigmaA4": sigmaA4,
        "sigmaA_e": sigmaA_e,
        "M4": M4,
    }

# -----------------------------
# 3) Minimal Z4 triad (main text)
# -----------------------------

def triad_minimal(M4: float) -> Dict[str, float]:
    """
    Minimal parameter-free splitting:
      s_R^2 = (2/3)*(1/3) = 2/9,
      c_R = sqrt(1 - s_R^2) = sqrt(7)/3,
      M4 = MW + MZ,
      MW/MZ = c_R,
      MH = M4/sqrt(2).
    """
    sR2 = 2.0 / 9.0
    cR = math.sqrt(1.0 - sR2)  # sqrt(7)/3

    MZ = M4 / (1.0 + cR)
    MW = cR * M4 / (1.0 + cR)
    MH = M4 / math.sqrt(2.0)

    return {
        "sR2": sR2,
        "cR": cR,
        "MW": MW,
        "MZ": MZ,
        "MH": MH,
    }

# -----------------------------
# 4) Refined Z4 triad (overlap-corrected appendix)
# -----------------------------

def triad_refined(M4: float, cR: float) -> Dict[str, float]:
    """
    Overlap-corrected version:
      r = MW/MZ = c_R,
      q_WW = 1/4,
      q_WZ = r^2/(1+r^2)^2 (=63/256 when r^2=7/9),
      Δ4 = (M4/(2π^2)) (q_WW^2 + 2 q_WZ^2),
      corrected additivity: MW + MZ = M4 - Δ4,
      solve with MW = r MZ,
      Higgs: λ0 = 1/8, q_eff^2 = (1/3)(q_WW^2 + 2 q_WZ^2),
             λ_R = λ0 (1 + (3/4) q_eff^2),
             MH = 2 sqrt(λ_R) M4.
    """
    pi = math.pi
    r = cR

    q_WW = 1.0 / 4.0
    q_WW2 = q_WW**2

    q_WZ = (r**2) / ((1.0 + r**2)**2)
    q_WZ2 = q_WZ**2

    Delta4 = (M4 / (2.0 * pi**2)) * (q_WW2 + 2.0 * q_WZ2)

    MZ = (M4 - Delta4) / (1.0 + r)
    MW = r * (M4 - Delta4) / (1.0 + r)

    q_eff2 = (1.0 / 3.0) * (q_WW2 + 2.0 * q_WZ2)
    lam0 = 1.0 / 8.0
    lamR = lam0 * (1.0 + (3.0 / 4.0) * q_eff2)
    MH = 2.0 * math.sqrt(lamR) * M4

    return {
        "r": r,
        "q_WW": q_WW,
        "q_WZ": q_WZ,
        "Delta4": Delta4,
        "q_eff2": q_eff2,
        "lambda0": lam0,
        "lambdaR": lamR,
        "MW": MW,
        "MZ": MZ,
        "MH": MH,
    }

# -----------------------------
# 5) Refined Z-pole proxies (s_eff^2, rho_R, gV/gA, A_l)
# -----------------------------

def z_pole_proxies() -> Dict[str, float]:
    """
    Implements your appendix definitions:
      K_parallel^(1)(0) = 2/pi,
      K_parallel^(1)(1/4) = 409/(105*pi^2),
      δ4 = (1/2) K(1/4)/K(0),
      q0 = 1/4,
      s_eff^2 = 1/(4 + δ4 + 2 q0^4),
      ρ_R = 1 + (δ4/4) q0^2,
      gA = sqrt(ρ_R) T3 with T3=-1/2,
      gV = sqrt(ρ_R)(T3 - 2 Q s_eff^2) with Q=-1,
      A_l = 2 gV gA / (gV^2 + gA^2).
    """
    pi = math.pi
    K10 = 2.0 / pi
    K114 = 409.0 / (105.0 * pi**2)
    delta4 = 0.5 * (K114 / K10)

    q0 = 1.0 / 4.0
    s_eff2 = 1.0 / (4.0 + delta4 + 2.0 * (q0**4))
    rhoR = 1.0 + (delta4 / 4.0) * (q0**2)

    T3 = -0.5
    Q = -1.0
    gA = math.sqrt(rhoR) * T3
    gV = math.sqrt(rhoR) * (T3 - 2.0 * Q * s_eff2)  # = sqrt(rhoR)*(-1/2 + 2 s_eff2)

    A_l = (2.0 * gV * gA) / (gV**2 + gA**2)

    return {
        "delta4": delta4,
        "s_eff2": s_eff2,
        "rhoR": rhoR,
        "gV_l": gV,
        "gA_l": gA,
        "A_l": A_l,
    }

# -----------------------------
# 6) Main execution; tables
# -----------------------------

def main() -> None:
    # --- Ladder -> M4
    z4 = compute_Z4_parent_scale_M4(ALPHA_FS, ME_GEV)
    M4 = z4["M4"]

    # --- Minimal model
    minres = triad_minimal(M4)
    # --- Refined model
    ref = triad_refined(M4, minres["cR"])
    # --- Z-pole proxies
    zp = z_pole_proxies()

    # --- Build a comparison table for masses
    rows = []
    for name, pdg_val in PDG_MASSES_GEV.items():
        pred_min = minres[f"M{name}"] if f"M{name}" in minres else minres[name]
        pred_ref = ref[f"M{name}"] if f"M{name}" in ref else ref[name]

        rows.append({
            "Particle": name,
            "PDG [GeV]": pdg_val,
            "Minimal pred [GeV]": pred_min,
            "Minimal δ [%]": rel_error_percent(pred_min, pdg_val),
            "Refined pred [GeV]": pred_ref,
            "Refined δ [%]": rel_error_percent(pred_ref, pdg_val),
        })

    # Print headline numbers
    print("\n=== Relator Z4 ladder (closed) ===")
    print(f"alpha_fs              = {fmt(ALPHA_FS, 15)}")
    print(f"m_e c^2 [GeV]         = {fmt(ME_GEV, 15)}")
    print(f"K_parallel^(1)(1/4)   = {fmt(z4['K1_a14'], 15)}")
    print(f"S4                    = {fmt(z4['S4'], 15)}")
    print(f"K_parallel^(4)        = {fmt(z4['K4'], 15)}")
    print(f"sigma_A^(4)           = {fmt(z4['sigmaA4'], 15)}")
    print(f"sigma_{'{'}A,e{'}'}          = {fmt(z4['sigmaA_e'], 15)}")
    print(f"M4 [GeV]              = {fmt(M4, 12)}")

    print("\n=== Minimal splitting (main text) ===")
    print(f"s_R^2                 = {fmt(minres['sR2'], 12)}")
    print(f"c_R                   = {fmt(minres['cR'], 12)}")
    print(f"MW, MZ, MH [GeV]       = {fmt(minres['MW'], 9)}, {fmt(minres['MZ'], 9)}, {fmt(minres['MH'], 9)}")

    print("\n=== Refined splitting (overlap-corrected appendix) ===")
    print(f"Delta4 [GeV]          = {fmt(ref['Delta4'], 12)}")
    print(f"MW, MZ, MH [GeV]       = {fmt(ref['MW'], 9)}, {fmt(ref['MZ'], 9)}, {fmt(ref['MH'], 9)}")
    print(f"lambda_R              = {fmt(ref['lambdaR'], 12)}")

    print("\n=== Z-pole proxy outputs (appendix definitions) ===")
    print(f"delta4                = {fmt(zp['delta4'], 12)}")
    print(f"s_eff^2               = {fmt(zp['s_eff2'], 12)}")
    print(f"rho_R                 = {fmt(zp['rhoR'], 12)}")
    print(f"gA_l                  = {fmt(zp['gA_l'], 12)}")
    print(f"gV_l                  = {fmt(zp['gV_l'], 12)}")
    print(f"A_l                   = {fmt(zp['A_l'], 12)}")

    # On-shell mixing proxy from the mass ratio (same in both models because MW/MZ=r=c_R is imposed)
    s_os2 = 1.0 - (minres["MW"] / minres["MZ"])**2
    print("\n=== On-shell proxy ===")
    print(f"s_os^2 = 1 - (MW/MZ)^2 = {fmt(s_os2, 12)}  (should equal s_R^2 by construction)")

    # Print the mass table
    try:
        import pandas as pd  # type: ignore
        df = pd.DataFrame(rows)
        # nicer ordering
        df = df[["Particle", "PDG [GeV]", "Minimal pred [GeV]", "Minimal δ [%]",
                 "Refined pred [GeV]", "Refined δ [%]"]]
        print("\n=== Mass comparison table ===")
        print(df.to_string(index=False, justify="center", float_format=lambda x: f"{x: .6f}"))
    except Exception:
        # Fallback if pandas is not available
        print("\n=== Mass comparison table (fallback) ===")
        header = f"{'Particle':>8} | {'PDG':>10} | {'Min pred':>12} | {'Min δ%':>8} | {'Ref pred':>12} | {'Ref δ%':>8}"
        print(header)
        print("-" * len(header))
        for r in rows:
            print(
                f"{r['Particle']:>8} | "
                f"{r['PDG [GeV]']:>10.6f} | "
                f"{r['Minimal pred [GeV]']:>12.6f} | "
                f"{r['Minimal δ [%]']:>8.4f} | "
                f"{r['Refined pred [GeV]']:>12.6f} | "
                f"{r['Refined δ [%]']:>8.4f}"
            )

    # Optional EW target comparison
    if any(x is not None for x in [PDG_SIN2_THETA_EFF, PDG_GV_L, PDG_GA_L, PDG_A_L]):
        print("\n=== Optional EW-coupling comparisons (fill PDG_* constants at top) ===")
        if PDG_SIN2_THETA_EFF is not None:
            print(f"sin^2(theta_eff): pred={fmt(zp['s_eff2'], 9)}  PDG={fmt(PDG_SIN2_THETA_EFF, 9)}"
                  f"  δ[%]={fmt(rel_error_percent(zp['s_eff2'], PDG_SIN2_THETA_EFF), 6)}")
        if PDG_GV_L is not None:
            print(f"gV_l:            pred={fmt(zp['gV_l'], 9)}  PDG={fmt(PDG_GV_L, 9)}"
                  f"  δ[%]={fmt(rel_error_percent(zp['gV_l'], PDG_GV_L), 6)}")
        if PDG_GA_L is not None:
            print(f"gA_l:            pred={fmt(zp['gA_l'], 9)}  PDG={fmt(PDG_GA_L, 9)}"
                  f"  δ[%]={fmt(rel_error_percent(zp['gA_l'], PDG_GA_L), 6)}")
        if PDG_A_L is not None:
            print(f"A_l:             pred={fmt(zp['A_l'], 9)}  PDG={fmt(PDG_A_L, 9)}"
                  f"  δ[%]={fmt(rel_error_percent(zp['A_l'], PDG_A_L), 6)}")

if __name__ == "__main__":
    main()