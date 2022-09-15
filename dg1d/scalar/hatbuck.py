import numpy as np
from bucklev import *
from scipy import optimize

xmin, xmax = -1.0, 1.0
def hatbuck(x):
    f = np.empty_like(x)
    for i,xx in enumerate(x):
        if xx < -0.5 or xx > 0.0:
            f[i] = 0.0
        else:
            f[i] = 1.0
    return f

# Finding the exact solution

u_s  = np.sqrt(a_buck/(1.0+a_buck)) # (f(u_s)-f(0))/(u_s-0)=f'(u_s)
u_ss = 1.0-1.0/np.sqrt(1.0+a_buck) # (f(u_ss)-f(1))/(u_ss-1)=f'(u_ss)

# Inverse of f' for computing rarefactions
# Inverse of f' restricted to [u_buck,1], an interval that contains [u_s,1]
def inv_f_s(v):
    # Inverse of f' at v equals root of this polynomial in [0.5,1]
    def p(u):
        value = v*(u**2+a_buck*(1.0-u)**2)**2-2.0*a_buck*u*(1.0-u)
        return value
    output = optimize.brentq(p, u_s, 1.0) # Gives root of polynomial
    return output

# Inverse of f' restricted to [0,u_buck], an interval that contains [0,u_ss]
def inv_f_ss(v):
    # Inverse of f' at v equals root of this polynomial in [0.5,1]
    def p(u):
        value = v*(u**2+a_buck*(1.0-u)**2)**2-2.0*a_buck*u*(1.0-u)
        return value
    output = optimize.brentq(p, 0.0, u_buck) # Gives root of polynomial
    return output

# Only works for time until the rarefactions intersect
def initial_condition(x, t=0.0):
    y = np.empty_like(x)
    for i, xx in enumerate(x):
        if xx <= -0.5:
            y[i] =  0.0
        elif -0.5 < xx <= -0.5+fprime(u_ss)*t:
            y[i] = inv_f_ss((xx+0.5)/t)
        elif -0.5+fprime(u_ss)*t < xx <= 0.0:
            y[i] = 1.0
        elif 0.0 < xx <= fprime(u_s)*t:
            y[i] = inv_f_s(xx/t)
        else:
            y[i] = 0.0
    return y

