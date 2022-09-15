import numpy as np

xmin, xmax = -1.0, 1.0

def initial_condition(x, t= 0.0):
    return np.sin(4*np.pi*x)


