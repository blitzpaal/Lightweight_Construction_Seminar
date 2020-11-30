from scipy.optimize import NonlinearConstraint, Bounds, differential_evolution, rosen
import numpy as np
from Shaft_dimensioning import calculate_shaft_strength

def symetric(x):
    return x[0] - x[-1]

def balanced(x):
    return np.sum(x)

nlc = NonlinearConstraint(balanced, 0.0, 0.0)
# specify limits using a `Bounds` object.   
bounds = Bounds([-90., -90., -90., -90., -90., -90., -90., -90., -90.], [90., 90., 90., 90., 90., 90., 90., 90., 90.])

#bounds = Bounds([-90., -90., -90., -90.], [90., 90., 90., 90.], )

#result = differential_evolution(CLT_Shaft, bounds, constraints=(nlc))
result = differential_evolution(calculate_shaft_strength, bounds)
print(result.x, result.fun)