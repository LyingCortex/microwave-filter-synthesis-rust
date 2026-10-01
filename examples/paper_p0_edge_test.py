#!/usr/bin/env python3
"""DECISIVE p0-edge test.

The printed Section-V design is self-consistent (verified: unitary poly identity
D*D* ~= |N|2+|p0t Pt|2+|p0r Pr|2 to 1e-3 at high order; S-product (aN+bD)/2 = Pt*Pr).
The code's 2x2 (eq.13) is *derived from* |S11|^2 = eps at a point. So the CORRECT
edges are exactly the points where the printed design's antenna |S11| = 0.1 (20 dB).

This script:
  1. loads the printed Np, D, Pt, Pr, p0t, p0r, n0, n2, b0;
  2. finds where |S11(s)| = 0.1 along the j-axis (the RL spec is pinned);
  3. checks whether the code's edges (tx_edge=-j, rx_edge=j(Δ+1)) coincide;
  4. runs the code's exact 2x2 at the printed state and compares to |p0t|2, |p0r|2.
"""
import numpy as np
J = 1j

# printed (from paper_check_printed2.py)
Np = np.array([0.76+0.65j, 0.033-0.039j, 2.166+1.854j, 0.09-0.105j, 2.31+1.976j,
               0.0886-0.103j, 1.14+0.977j, 0.0391-0.04561j, 0.265+0.226j,
               0.0076-0.0089j, 0.026+0.022j, 5.5e-4-6.48e-4j, 9.29e-4+7.95e-4j,
               1.055e-5-1.23e-5j, 8.87e-6+7.59e-6j], complex)
D = np.array([1, 1.77-0.051j, 4.417-0.104j, 5.3-0.242j, 6.724-0.3119j,
              5.61-0.376j, 4.425-0.326j, 2.51-0.24j, 1.25-0.14j, 0.44-0.065j,
              0.128-0.02248j, 0.024-5.7e-3j, 3.27e-3-9.5e-4j, 2.26e-4-9.21e-5j,
              1.123e-5-4.25e-6j], complex)
Pt = np.array([1, 0.455+3.9j, -6.027+1.53j, -1.97-4.6j, 1.833-1.21j,
               0.359+0.351j, -0.0275+0.046j, -0.0018-0.000579j], complex)
Pr = np.array([1, 0.43-4.29j, -7.337-1.6j, -2.31+6.37j, 2.95+1.6j,
               0.565-0.7j, -0.076-0.091j, -0.005+0.0027j], complex)
p0t = 4.717e-4 + 1.743e-4j
p0r = 4.336e-4 + 1.6e-4j
n_j, b0 = 1.47, -0.171
n2 = n_j**2
n0 = (1-J*n2*b0)/(1+J*n2*b0)

def ep(c, z): return np.polyval(c, z)

print(f"n0={n0:.6f}  |p0t|2={abs(p0t)**2:.6e}  |p0r|2={abs(p0r)**2:.6e}")

# roots of Pt, Pr -> where do the channel polynomials live (domain probe)
rt, rr = np.roots(Pt), np.roots(Pr)
print("Pt roots imag:", np.round(np.sort(rt.imag), 3))
print("Pr roots imag:", np.round(np.sort(rr.imag), 3))

# |S11(s)| = |Np(s)|/|D(s)|  (since Np = n0*N_true, Np/D = S11)
ys = np.linspace(-4.5, 2.0, 2000)
s11 = np.abs([ep(Np, J*y) for y in ys]) / np.abs([ep(D, J*y) for y in ys])
# find zeros of (|S11| - 0.1)
target = 0.1
d = s11 - target
signs = np.sign(d)
crossings = []
for i in range(len(ys)-1):
    if signs[i]==0: crossings.append(ys[i]); continue
    if signs[i]*signs[i+1] < 0:
        # linear interp
        crossings.append(ys[i] - d[i]*(ys[i+1]-ys[i])/(d[i+1]-d[i]))
print("\n|S11|=0.1 (20 dB) crossings at s=j*y, y =", np.round(crossings, 4))

# code's edges. normalized_separation for the code's geom:
# delta_center = rx_center - tx_center. code tx=15.15-15.35 center 15.25,
# rx=14.9-15.1 center 15.0. delta=-0.25GHz, avg_bw=0.2GHz
delta = 15.0-15.25; avg_bw = 0.2
Delta = 2.0*delta/avg_bw
print(f"\ncode normalized_separation Delta={Delta:.4f}")
print(f"code tx_edge  = -j  -> y = -1.0")
print(f"code rx_edge  = j(Delta+1) -> y = {Delta+1:.4f}")

# run the code's exact 2x2 at the printed state (common domain = printed domain)
# Pt_common, Pr_common = printed Pt, Pr (already in common domain per S-product check)
def n2_2(N, Pt, Pr, e1, e2, eps1, eps2):
    N1=abs(ep(N,e1))**2; Pt1=abs(ep(Pt,e1))**2; Pr1=abs(ep(Pr,e1))**2
    N2=abs(ep(N,e2))**2; Pt2=abs(ep(Pt,e2))**2; Pr2=abs(ep(Pr,e2))**2
    b1=N1*(1-eps1)/eps1; b2=N2*(1-eps2)/eps2
    det=Pt1*Pr2-Pt2*Pr1
    p0t2=(b1*Pr2-b2*Pr1)/det
    p0r2=(Pt1*b2-Pt2*b1)/det
    return p0t2,p0r2,det

eps = 10**(-20/10)
print("\n-- code's edges: tx_edge=-j (TX RL), rx_edge=j(Delta+1) (RX RL) --")
# NOTE: which channel's RL at which edge? code: row1=tx edge w/ eps_tx, row2=rx edge w/ eps_rx
# but printed Pt=TX? Let's try both assignments.
for (e_tx, e_rx, label) in [ (J*-1.0, J*(Delta+1), "code outer?"),
                             (J*(Delta+1), J*-1.0, "swap") ]:
    p0t2,p0r2,det = n2_2(Np/abs(Np[0]), Pt, Pr, e_tx, e_rx, eps, eps)
    print(f"  {label}: det={det:.3e} p0t2={p0t2:.4e} p0r2={p0r2:.4e}")

print("\n-- test: does 2x2 at the *found* |S11|=0.1 crossings reproduce |p0t|2,|p0r|2? --")
if len(crossings)>=2:
    for (a,b) in [(crossings[0],crossings[1]),(crossings[-2],crossings[-1]),(crossings[0],crossings[-1])]:
        p0t2,p0r2,det = n2_2(Np/abs(Np[0]), Pt, Pr, J*a, J*b, eps, eps)
        print(f"  y=[{a:.3f},{b:.3f}] det={det:.3e} p0t2={p0t2:.4e} p0r2={p0r2:.4e}")
print(f"  target: |p0t|2={abs(p0t)**2:.4e} |p0r|2={abs(p0r)**2:.4e}")
