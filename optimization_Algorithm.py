from scipy.optimize import NonlinearConstraint, Bounds, differential_evolution, shgo, dual_annealing, minimize
import numpy as np
from Shaft_dimensioning import calculate_shaft_strength, compose_stack
from CLT_calculation import calc_Q_0, CLT_ABD

# Material data
t_ply = 0.125  # ply thickness in mm
E_11 = 126000  # Longitudinal tensile modulus in MPa
E_22 = 9000  # Transverse tensile modulus in MPa
G_12 = 4600  # Shear modulus in MPa
v_12 = 0.3  # Poisson’s ratio 1

# Stiffness matrix of UD-Layer
Q_0 = calc_Q_0(E_11, E_22, G_12, v_12)

def balanced(stack_angle):
    if symetric == True:
        stack_angle = np.concatenate((stack_angle, np.flip(stack_angle)))
    
    stack = compose_stack(stack_angle, t_ply)
    stack = compose_stack(stack_angle, t_ply)

    ABD = CLT_ABD(stack, Q_0)

    return ABD[0:2,2]

# Specify limits using a `Bounds` object
bounds = Bounds([-90., -90., -90., -90., -90.], [90., 90., 90., 90., 90.])

# Constraints for balanced laminate (two options)
balanced_laminate = NonlinearConstraint(balanced, 0.0, 0.0)
balanced = False

# Constraint for symetric laminate
symetric = True
"""
# Global optimization
glob_result = differential_evolution(calculate_shaft_strength, bounds=bounds, args=(balanced, symetric))
#glob_result = differential_evolution(calculate_shaft_strength, bounds=bounds, args=(balanced, symetric), constraints=(balanced_laminate))
print(glob_result.x, glob_result.fun)
"""
# Local optimization
x0 = np.array([81.98268888, -78.81994856,  47.28577036,  39.81491606,  -9.35838018])

#loc_result = minimize(calculate_shaft_strength, x0, tol=1e-6, bounds=bounds, args=(balanced, symetric))
loc_result = minimize(calculate_shaft_strength, x0, tol=1e-6, bounds=bounds, args=(balanced, symetric), constraints=(balanced_laminate))

print(loc_result.x, loc_result.fun)
