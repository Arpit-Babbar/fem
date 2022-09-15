import numpy as np

xmin, xmax = 0.0, 1.0

def initial_condition(x, t= 0.0):
    return 1 + np.sin(2*np.pi*x)


