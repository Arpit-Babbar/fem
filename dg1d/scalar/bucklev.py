'''
Buckley-Leverett model
'''
import numpy as np

# convex in [0,u_bucklev] and concave in [u_bucklev,1]
# u_bucklev depends on a_bucklev
a_buck = 0.25
u_buck = 0.287141

# f = u^2/(u^2 + a(1-u)^2)
def flux(x,u):
    a = a_buck
    return  u**2/(u**2 + a*(1.0-u)**2)

def fprime(u):
    a = a_buck
    L = u**2 + a*(1.0-u)**2
    return 2.0*a*u*(1.0-u) / L**2

rhs = rhst = rhstt = rhsttt = rhstttt = lambda x,t : 0.0

def max_speed_bucklev(ul, ur):
    umin = max(min(ul,ur), 0.0)
    umax = min(max(ul,ur), 1.0)
    if umin > u_buck or umax < u_buck:
        return max(fprime(umin),fprime(umax))
    else:
        return fprime(u_buck)

# Rusanov flux
def rusanov(x, ul, ur):
    a = max_speed_bucklev(ul,ur)
    fl = flux(x, ul)
    fr = flux(x, ur)
    return 0.5*(fl + fr) - 0.5*a*(ur - ul)

def upwind(x, ul, ur):
    return flux(x, ul)

# Max speed based on cell average values
def max_speed(u):
    smax = fprime(u_buck)
    smax = max(smax, np.abs(fprime(u)).max())
    return smax