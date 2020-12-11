import numpy as np
import math

from CLT_calculation import calc_Q_0, CLT_ABD, Shell_Engineering_Constants, CLT_Stress
from Rotational_speed_Buckling import twisted_critical_speed, bending_critial_speed, torsion_buckling

# Geometry and Load data
l = 500  # shaft length in mm
l_f = 35 # flange length in mm
l_s = 500-2*l_f # free length betwwen flanges in mm
d_i = 45  # inner diameter in mm
T_n = 400*1000  # net torque in Nmm
RF = 2  # reserve factor
T = T_n*RF  # torque with reserve factor in Nmm
N_n = 5000/60  # net rotational speed in 1/s
N = N_n*RF  # rotational speed with reserve factor in 1/s

# Material data
t_ply = 0.125  # ply thickness in mm
E_11 = 126000  # Longitudinal tensile modulus in MPa
E_22 = 9000  # Transverse tensile modulus in MPa
G_12 = 4600  # Shear modulus in MPa
v_12 = 0.3  # Poisson’s ratio 1
rho = 1570 * 10**-12  # Density in Ns^2/mm^4
R_m1Z = 1250  # Longitudinal tensile strength in MPa
R_m2Z = 50  # Transverse tensile strength in MPa
R_m1D = 1260  # Longitudinal compressive strength in MPa
R_m2D = 210  # Transverse compressive strength in MPa
R_m12 = 84  # Shear strength in MPa

# Puck slope parameters
p_tp_ten = 0.35
p_tp_com = 0.3
p_tt_ten = 0.25
p_tt_com = 0.25

# Laminate data
stack = None
F = None
Q_0 = None
ABD = None

def stress_vector(T, stack, d_i):
    # since there is only torque n_x, n_y, m_x, m_y and m_xy = 0
    F = np.zeros(6)
    r_i = d_i / 2
    t = np.sum(stack[:,1])
    r_o = r_i + t
    A_m = np.pi*((r_o + r_i) / 2) ** 2
    tau_xy = T/(2*A_m*t)
    n_xy = tau_xy * t
    F[2] = n_xy

    return F

def compose_stack(stack_angle, t_ply):
    stack = np.zeros((stack_angle.shape[0],3))
    stack[:,0] = stack_angle
    stack[:,1] = t_ply
    stack[:,0] = np.deg2rad(stack[:,0]) # transform ply angle from degree to radian

    z_ref = np.sum(stack[:,1]) / 2 # z-position of reference plane, relative to bottom of stack
    for i in range(stack.shape[0]):
        if i==0:
            stack[i,2] = -z_ref+stack[i,1]
        elif i>0:
            stack[i,2] = stack[i-1,2]+stack[i,1]

    return stack

def calculate_shaft_strength(stack_angle, balanced, symetric):
    if balanced == True:
        stack_angle_bal = np.zeros(stack_angle.shape[0]*2)
        i=0
        for angle in stack_angle:
            stack_angle_bal[i] = angle
            stack_angle_bal[i+1] = -angle
            i = i+2
        stack_angle = stack_angle_bal
    
    if symetric == True:
        stack_angle = np.concatenate((stack_angle, np.flip(stack_angle)))
    
    stack = compose_stack(stack_angle, t_ply)

    F = stress_vector(T, stack, d_i)

    # Stiffness matrix of UD-Layer
    Q_0 = calc_Q_0(E_11, E_22, G_12, v_12)

    ABD = CLT_ABD(stack, Q_0)

    t = np.sum(stack[:,1])

    E_Ax, E_Ay, G_Axy, v_Axy = Shell_Engineering_Constants(ABD, t)

    f_E_FF, f_E_IFF = CLT_Stress(stack, Q_0, ABD, F, R_m1Z, R_m2Z, R_m1D, R_m2D, R_m12, p_tp_ten, p_tp_com, p_tt_ten, p_tt_com)

    f_E_max = np.max((f_E_FF,f_E_IFF))

    N_crit_twist = twisted_critical_speed(d_i, t, l_s, rho, G_Axy)
    N_crit_bend = bending_critial_speed(d_i, t, l_s, rho, E_Ax)
    RF_N = N / min(N_crit_twist, N_crit_bend)

    T_crit = torsion_buckling(d_i, t, l_s, E_Ax, E_Ay)
    RF_T = T / T_crit

    return max(f_E_max, RF_N, RF_T)

balanced = False
symetric = True
print(calculate_shaft_strength(np.array([ 79.0,  85.0, -68.0,  23.0, -18.0]), balanced, symetric))