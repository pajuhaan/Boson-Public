# Boson-Public — Relator $\mathbb{Z}_N$ Bouquet Ladder

[![GitHub](https://img.shields.io/badge/GitHub-pajuhaan%2FBoson--Public-181717?logo=github)](https://github.com/pajuhaan/Boson-Public)
[![Stars](https://img.shields.io/github/stars/pajuhaan/Boson-Public?style=social)](https://github.com/pajuhaan/Boson-Public/stargazers)
[![Last Commit](https://img.shields.io/github/last-commit/pajuhaan/Boson-Public)](https://github.com/pajuhaan/Boson-Public/commits/main)

Public repository: https://github.com/pajuhaan/Boson-Public

This repository contains a small, reproducible Python reference implementation of the **Relator $\mathbb{Z}_N$ bouquet scan** used to compute

1) the closed ladder “rung” masses $m_N$ from an electron anchor $(m_e,\alpha)$,  
2) the **$\mathbb{Z}_4$ electroweak parent scale** $M_4$ and a **parameter-free mass triad** $(M_W,M_Z,M_H)$ in a minimal baseline,  
3) a **refined $\mathbb{Z}_4$ variant** with closed overlap corrections and a coupling-level normalization (Appendix-style model),  
4) the **$\mathbb{Z}_5$ first ultra-heavy rung** $M_5$ and its **minimal splitting** into two new charged vector modes at the $\mathcal{O}(10\!\-\!20)\,\mathrm{TeV}$ scale.

I keep the code modular and explicit so the numerical chain is readable and auditable.

---

## What is computed

### 1) Locked $\mathbb{Z}_N$ bouquet data

For an $N$-loop phase-locked bouquet, the discrete twisted sectors are indexed by $k=1,\dots,N-1$ with twist offsets

$$
a_k \;:=\; \min\!\left(\frac{k}{N},\,1-\frac{k}{N}\right)\in\Bigl[0,\tfrac12\Bigr].
$$

### 2) Angular reduction constant and dilation exponent

Define

$$
S_N \;:=\; 3\sum_{k=1}^{N-1}K_{\parallel}^{(1)}(a_k),
\qquad
P_J^{(N)} \;:=\;\frac{N-1}{3N},
$$

$$
K_{\parallel}^{(N)} \;:=\;\frac{6}{\pi}+P_J^{(N)}S_N,
\qquad
\sigma_A^{(N)} \;:=\;\frac{3}{K_{\parallel}^{(N)}}.
$$

Here $K_{\parallel}^{(1)}(a)$ is the one-sector PV-reduction kernel (closed evaluations are used at the special offsets needed for $N=3,4,5$).

### 3) Relator mass ladder (electron-anchored)

The ladder mass at rung $N$ is computed from

$$
\frac{m_N}{m_e}
\;=\;
\exp\!\left(\frac{\sigma_{A,e}-\sigma_A^{(N)}}{4\alpha}\right),
\qquad
\sigma_{A,e}\equiv\frac{\pi}{2}.
$$

The implementation supports SI-consistent inputs and the standard particle-physics convention where masses are reported as rest-energy scales in GeV.

---

## Model variants implemented

### A) Minimal baseline ($\mathbb{Z}_4$ triad)

From $M_4=m_4$, a purely geometric mixing proxy is introduced for $\mathbb{Z}_4$

$$
s_{\mathrm R}^2=\frac{2}{9},
\qquad
c_{\mathrm R}=\sqrt{1-s_{\mathrm R}^2}=\frac{\sqrt7}{3}.
$$

A minimal on-shell identification is then used

$$
\frac{M_W}{M_Z}=c_{\mathrm R},
\qquad
M_4=M_W+M_Z,
$$

which yields $(M_W,M_Z)$ in closed form, and a minimal radial-mode estimate gives $M_H$.

### B) Refined $\mathbb{Z}_4$ overlap-corrected model (Appendix-style)

A closed overlap correction $\Delta_4$ is computed from a Gaussian-collar overlap invariant $q_{ij}$ and inserted into a corrected additivity

$$
M_W+M_Z=M_4-\Delta_4,
$$

with the same ratio constraint $M_W/M_Z=c_{\mathrm R}$.  
This variant also adds a minimal coupling-level normalization giving $s_{\mathrm{eff}}^2$, $\rho_{\mathrm R}$, and the leptonic couplings $\bar g_V^\ell,\bar g_A^\ell$ and $A_\ell$.

### C) $\mathbb{Z}_5$ ultra-heavy rung and new charged sector (prediction)

The code computes $M_5=m_5$ from the same ladder relation and performs a minimal two-mode split controlled by closed PV weights

- two complex-conjugate pairs $(k=1,4)$ and $(k=2,3)$,  
- two massive complex charged vector modes $\mathcal{V}_{5,1}^{\pm}$ and $\mathcal{V}_{5,2}^{\pm}$ with masses $(M_{5,1},M_{5,2})$,  
- an optional radial scalar $\mathcal{S}_5$ (baseline estimate).

The $\mathbb{Z}_5$ sector is presented as a leading-order prediction; widths, lifetimes, and inter-topology mixing are not modeled here.

---

## Installation

```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell
# .venv\Scripts\Activate.ps1

pip install -r requirements.txt
```

---

## Usage

Run the full computation and print tables

```bash
python scripts/run_all.py
```

Expected output includes

- rung masses $m_N$ for $N=3,4,5$,  
- $\mathbb{Z}_4$ minimal baseline $(M_W,M_Z,M_H)$ and deviations vs PDG,  
- $\mathbb{Z}_4$ refined overlap-corrected masses and coupling proxies,  
- $\mathbb{Z}_5$ spectrum summary $(M_5, M_{5,1}, M_{5,2})$ and optional $\mathcal{S}_5$.

---

## Notes and limitations

- This is a leading-order implementation of the Relator bouquet scan and the minimal splitting logic.  
- It does not compute lifetimes, total widths, radiative corrections, detector-level observables, or detailed micro-structural selection effects.  
- Inter-topology transitions (schematically $\mathcal{B}_5\to\mathcal{B}_4$) are treated as suppressed in the minimal baseline and are not used to fix branching fractions.  
- The point of this repo is reproducibility of the closed-form pipeline, not a precision electroweak fit.

---

## Reference numerical checkpoints (sanity checks)

The repository is expected to reproduce the following checkpoints (up to rounding and constant precision)

- $\mathbb{Z}_4$ parent rung: $M_4 \approx 173.2276\ \mathrm{GeV}$  
- Minimal $\mathbb{Z}_4$ triad (baseline)  
  - $M_W \approx 81.1791\ \mathrm{GeV}$  
  - $M_Z \approx 92.0485\ \mathrm{GeV}$  
  - $M_H \approx 122.4904\ \mathrm{GeV}$  
- Refined $\mathbb{Z}_4$ overlap-corrected  
  - $M_W \approx 80.42394\ \mathrm{GeV}$  
  - $M_Z \approx 91.19218\ \mathrm{GeV}$  
  - $M_H \approx 125.27037\ \mathrm{GeV}$  
- $\mathbb{Z}_5$ rung and split (prediction)  
  - $M_5 \approx 32.9261\ \mathrm{TeV}$  
  - $M_{5,1}\approx 22.1554\ \mathrm{TeV}$  
  - $M_{5,2}\approx 10.7707\ \mathrm{TeV}$  
  - optional $\mathcal{S}_5 \sim 20.8\ \mathrm{TeV}$ (baseline)

---

## How to cite


