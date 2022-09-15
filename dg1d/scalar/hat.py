import numpy as np

xmin, xmax = 0.0, 1.0

def initial_condition(x, t= 0.0):
    f = np.empty_like(x)
    for i,xx in enumerate(x):
        if xx < 0.25 or xx > 0.75:
            f[i] = 0.0
        else:
            f[i] = 1.0
    return f
