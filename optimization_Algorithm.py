from scipy.optimize import NonlinearConstraint, Bounds, differential_evolution, rosen
import numpy as np
from CLT_Test import CLT_Shaft

# specify limits using a `Bounds` object.
bounds = Bounds([0., 0., 0., 0.], [90., 90., 90., 90.])

result = differential_evolution(CLT_Shaft, bounds)
print(result.x, result.fun)