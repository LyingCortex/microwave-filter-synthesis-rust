#!/usr/bin/env python3
"""Probe the paper's common domain from the printed polynomials.
Compute |S11|, |S21_TX|=|p0t Pt/D|, |S21_RX|=|p0r Pr/D| vs s=j*y and
identify where each channel's passband actually sits."""
import numpy as np
J = 1j
Np = np.array([0.76+0.65j,0.033-0.039j,2.166+1.854j,0.09-0.105j,2.31+1.976j,
              0.0886-0.103j,1.14+0.977j,0.0391-0.04561j,0.265+0.226j,
              0.0076-0.0089j,0.026+0.022j,5.5e-4-6.48e-4j,9.29e-4+7.95e-4j,
              1.055e-5-1.23e-5j,8.87e-6+7.59e-6j],complex)
D = np.array([1,1.77-0.051j,4.417-0.104j,5.3-0.242j,6.724-0.3119j,5.61-0.376j,
             4.425-0.326j,2.51-0.24j,1.25-0.14j,0.44-0.065j,0.128-0.02248j,
             0.024-5.7e-3j,3.27e-3-9.5e-4j,2.26e-4-9.21e-5j,1.123e-5-4.25e-6j],complex)
Pt = np.array([1,0.455+3.9j,-6.027+1.53j,-1.97-4.6j,1.833-1.21j,0.359+0.351j,
              -0.0275+0.046j,-0.0018-0.000579j],complex)
Pr = np.array([1,0.43-4.29j,-7.337-1.6j,-2.31+6.37j,2.95+1.6j,0.565-0.7j,
              -0.076-0.091j,-0.005+0.0027j],complex)
p0t=4.717e-4+1.743e-4j; p0r=4.336e-4+1.6e-4j
n2=1.47**2; b0=-0.171; n0=(1-J*n2*b0)/(1+J*n2*b0)
def ep(c,z): return np.polyval(c,z)
ys=np.linspace(-2.0,2.0,401)
print(" y      |S11|    |S21_TX|   |S21_RX|")
for y in np.arange(-1.4,1.41,0.2):
    s=J*y
    S11=abs(ep(Np,s))/abs(ep(D,s))
    s21t=abs(p0t*ep(Pt,s)/ep(D,s))
    s21r=abs(p0r*ep(Pr,s)/ep(D,s))
    print(f"{y:+.1f}  {S11:.4f}   {s21t:.4f}     {s21r:.4f}")
