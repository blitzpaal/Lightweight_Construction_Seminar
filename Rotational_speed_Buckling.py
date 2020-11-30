import numpy as np
import math

# Acc. to Roloff/Matek S.408
def twisted_critical_speed(d_i, t, l_s, rho, G_Axy):
    r_i = d_i / 2
    r_o = r_i + t
    d_o = 2*r_o

    I_p = np.pi/2 * (r_o**4 - r_i**4)
    J = 1/8 * (d_o**2 + d_i**2) * rho * 2 * np.pi * (r_o+r_i)/2 * t * l_s

    phi = 180/np.pi * l_s/(G_Axy*I_p)

    N_crit = 72.3 * (1/(phi*J))**0.5

    return N_crit

# Acc. to Haberhauer/Bodenstein S.307
def bending_critial_speed(d_i, t, l_s, rho, E_Ax):
    r_i = d_i / 2
    r_o = r_i + t

    A = np.pi * (r_o**2-r_i**2)
    I_b = np.pi * (r_o**4-r_i**4) / 4

    omega_crit = np.pi**2/l_s**2 * ((E_Ax*I_b)/(rho*A))**0.5

    N_crit = omega_crit / (2*np.pi)

    return N_crit

"""
    I_p = np.pi/2 * (r_o**4 - r_i**4)
    I_A = 2*I_p

    m = rho * l_s * np.pi*(r_o**2-r_i**2)

    f_max = (rho * m * l_s**3) / (48 * E_Ax * I_A)

    N_crit = 946 * (1/f_max)**0.5
    """

# Acc. to "The Use of Continuous Fiber Composites in Driveshafts", Duane V. Byerly, 1996
def torsion_buckling(d_i, t, l_s, E_Ax, E_Ay):
    #r_i = d_i / 2
    #r_o = r_i + t
    #r_m = (r_o+r_i)/2

    #T_crit = 2 * np.pi * r_m**2 * t * 0.272 * (E_Ax*E_Ay**3)**(1/4) * (t/r_m)**(3/2)
    
    d_o = d_i + 2*t
    d_m = (d_o+d_i)/2
    
    T_crit = 1.854/(l_s)**0.5 * E_Ax**0.375 * E_Ay**0.625 * t**2.25 * d_m**1.25

    return T_crit