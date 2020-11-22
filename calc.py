import numpy as np
import math


t = 0.125   # in mm
E11z = 126000   # in MPa
E22z = 9000
E11d = 126000
E22d = 9000
G12 = 4600
v12 = 0.3
v21 = 0.021
rho = 1570  # in kg/m^3
Rm1z = 1250     # in MPa
Rm2z = 50
Rm1d = 1260
Rm2d = 210
Rm12 = 84


# Stiffness matrix for UD layer: 0deg
Q_0 = np.array([[E11z/(1-v12*v21), (v21*E11z)/(1-v12*v21), 0],
                [(v12*E22z)/(1-v12*v21), E22z/(1-v12*v21), 0],
                [0, 0, G12]])


def stiffness_matrix(a_deg):
    a_deg = 45  # layer orientation in degrees
    a = a_deg * (math.pi / 180)  # layer orientation in rad

    # transformation matrices t, (t_s is actually the already inverted t_s)
    t_s = np.array([[math.cos(a) ** 2, math.sin(a) ** 2, -2 * math.sin(a) * math.cos(a)],
                    [math.sin(a) ** 2, math.cos(a) ** 2, 2 * math.sin(a) * math.cos(a)],
                    [math.sin(a) * math.cos(a), -math.sin(a) * math.cos(a), math.cos(a) ** 2 - math.sin(a) ** 2]])

    t_e = np.array([[math.cos(a) ** 2, math.sin(a) ** 2, math.sin(a) * math.cos(a)],
                    [math.sin(a) ** 2, math.cos(a) ** 2, -math.sin(a) * math.cos(a)],
                    [-2 * math.sin(a) * math.cos(a), 2 * math.sin(a) * math.cos(a),
                     math.cos(a) ** 2 - math.sin(a) ** 2]])

    return np.matmul(np.matmul(t_s, Q_0), t_e)


def abd(x, y):
#for i in range(0, 3):
    a = y*2 + x + y*2
    b = y*-3 + y*3
    d = y*5.2 + x*0.08 + y*5.2
    return np.array([[a, b],
                     [b, d]])
