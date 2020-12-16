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

def bending_critial_speed(d_i, t, l_s, rho, E_Ax):
    r_i = d_i / 2
    r_o = r_i + t

    A = np.pi * (r_o**2-r_i**2)
    I_b = np.pi * (r_o**4-r_i**4) / 4

    # Acc. to Haberhauer/Bodenstein S.307
    #omega_crit = np.pi**2/l_s**2 * ((E_Ax*I_b)/(rho*A))**0.5
    #N_crit = omega_crit / (2*np.pi)

    d_m = d_i + t

    # Acc. to Dickhut
    N_crit = 0.5*np.pi/8**0.5 * d_m/l_s**2 * (E_Ax/rho)**0.5

    return N_crit

def torsion_buckling(d_i, t, l_s, E_Ax, E_Ay, D_yy):
    r_i = d_i / 2
    r_o = r_i + t
    r_m = (r_o+r_i)/2
    
    d_o = d_i + 2*t
    d_m = (d_o+d_i)/2
    
    # Acc. to "The Use of Continuous Fiber Composites in Driveshafts", Duane V. Byerly, 1996
    # T_crit = 1.854/(l_s)**0.5 * E_Ax**0.375 * E_Ay**0.625 * t**2.25 * d_m**1.25

    # Acc. to "Handbook of structural stability", Column Research Committee of Japan, Tokyo, Corona Publishing, 1971
    # T_crit = 2 * np.pi * r_m**2 * t * 0.272 * (E_Ax*E_Ay**3)**(1/4) * (t/r_m)**(3/2)

    k_l = 0.760
    k_s = 1.030 # fixed bearing (0.925 for floating bearing)

    # Acc. to Dickhut / "Instability of orthotropic cylindrical shells under combined torsion and hydrostatic pressure", Simitses, 1967
    T_crit = k_l * k_s * (2*np.pi**3)/12**(3/8) * (r_m**(5/4)*t**(3/8))/l_s**0.5 * E_Ax**(3/8) * D_yy**(5/8)

    return T_crit