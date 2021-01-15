import numpy as np
import math

def calc_T_e(ply):
    T_e = np.array([[np.cos(ply[0])**2, np.sin(ply[0])**2, np.sin(ply[0])*np.cos(ply[0])],
    [np.sin(ply[0])**2, np.cos(ply[0])**2, -np.sin(ply[0])*np.cos(ply[0])],
    [-2*np.sin(ply[0])*np.cos(ply[0]), 2*np.sin(ply[0])*np.cos(ply[0]), np.cos(ply[0])**2-np.sin(ply[0])**2]])

    return T_e

def calc_T_s(ply):
    T_s = np.array([[np.cos(ply[0])**2, np.sin(ply[0])**2, 2*np.sin(ply[0])*np.cos(ply[0])],
    [np.sin(ply[0])**2, np.cos(ply[0])**2, -2*np.sin(ply[0])*np.cos(ply[0])],
    [-np.sin(ply[0])*np.cos(ply[0]), np.sin(ply[0])*np.cos(ply[0]), np.cos(ply[0])**2-np.sin(ply[0])**2]])

    return T_s

def calc_Q_0(E_11, E_22, G_12, v_12):
    v_21 = v_12*E_22/E_11  # Poisson’s ratio 2
    Q_0 = np.array([[E_11/(1-v_12*v_21),v_21*E_11/(1-v_12*v_21),0],
    [v_12*E_22/(1-v_12*v_21),E_22/(1-v_12*v_21),0],
    [0,0,G_12]])

    return Q_0

def CLT_ABD(stack, Q_0):
    # Transform Laminate stiffness matrix of each layer
    Q_trans = np.empty((Q_0.shape[0],Q_0.shape[1],stack.shape[0]))

    for i in range(stack.shape[0]):
        T_e = calc_T_e(stack[i])
        T_s = calc_T_s(stack[i])

        Q_trans[:,:,i] = np.linalg.inv(T_s) @ Q_0 @ T_e

    # ABD Matrix calculation
    A = np.sum(Q_trans*stack[:,1],2)
    B = np.sum(Q_trans*stack[:,1]*(stack[:,2]-stack[:,1]/2),2)
    D = np.sum(Q_trans*(stack[:,1]**3/12+stack[:,1]*(stack[:,2]-stack[:,1]/2)**2),2)
    #print(A)
    #print(B)
    #print(D)

    # Assemble ABD Matrix
    ABD = np.vstack((np.hstack((A, B)), np.hstack((B, D))))
    #print(ABD)

    return ABD

def Shell_Engineering_Constants(ABD, t):
    ABD_inv = np.linalg.inv(ABD)

    E_Ax = 1 / (ABD_inv[0,0]*t)
    E_Ay = 1 / (ABD_inv[1,1]*t)
    G_Axy = 1 / (ABD_inv[2,2]*t)
    v_Axy = -ABD_inv[0,1] / ABD_inv[0,0]
    v_Ayx = -ABD_inv[0,1] / ABD_inv[1,1]

    return E_Ax, E_Ay, G_Axy, v_Axy, v_Ayx

def Plate_Engineering_Constants(ABD, t):
    ABD_inv = np.linalg.inv(ABD)

    E_Dx = 12 / (ABD_inv[3,3]*t**3)
    E_Dy = 12 / (ABD_inv[4,4]*t**3)
    G_Dxy = 12 / (ABD_inv[5,5]*t**3)
    v_Dxy = -ABD_inv[3,4] / ABD_inv[3,3]
    v_Dyx = -ABD_inv[3,4] / ABD_inv[4,4]

    return E_Dx, E_Dy, G_Dxy, v_Dxy, v_Dyx

def CLT_Stress_Puck(stack, Q_0, ABD, F, R_m1Z, R_m2Z, R_m1D, R_m2D, R_m12, p_tp_ten, p_tp_com, p_tt_ten, p_tt_com):
    # Strain in global coordinate system
    eps_kap_0 = np.linalg.inv(ABD) @ F
    eps_0 = eps_kap_0[:3]
    kap_0 = eps_kap_0[3:]

    # Global and local strain and local stress in each layer
    eps_xy = np.zeros((stack.shape[0],eps_0.shape[0]))
    eps_12 = np.zeros((stack.shape[0],eps_0.shape[0]))
    sig_12 = np.zeros((stack.shape[0],eps_0.shape[0]))

    # Stress exposure for each layer with Puck criterion
    f_E_FF = np.zeros(stack.shape[0])
    f_E_IFF = np.zeros(stack.shape[0])

    for i in range(stack.shape[0]):
        eps_xy[i,:] = eps_0 + stack[i,2] * kap_0

        T_e = calc_T_e(stack[i])

        eps_12[i,:] = T_e @ eps_xy[i,:]

        sig_12[i,:] = Q_0 @ eps_12[i,:]

        # Puck fibre failure criterion
        if sig_12[i,0] >= 0:
            f_E_FF[i] = np.abs(sig_12[i,0]) / R_m1Z
        elif sig_12[i,0] < 0:
            f_E_FF[i] = np.abs(sig_12[i,0]) / R_m1D

        if sig_12[i,1] >= 0: # Mode A
            f_E_IFF[i] = ((1-p_tp_ten*R_m2Z/R_m12)**2 * (sig_12[i,1]/R_m2Z)**2 + (sig_12[i,2]/R_m12)**2)**0.5 + p_tp_ten*sig_12[i,1]/R_m12
        else:
            R_tt_A = R_m12/(2*p_tp_com) * ((1+2*p_tp_com*R_m2D/R_m12)**0.5 - 1)
            tau_21c = R_m12 * (1+2*p_tt_com)**0.5
            if sig_12[i,1] < 0 and 0 <= np.abs(sig_12[i,1]/sig_12[i,2]) <= np.abs(R_tt_A/tau_21c): # Mode B
                f_E_IFF[i] = ((sig_12[i,2]/R_m12)**2 + (p_tp_com/R_m12*sig_12[i,1])**2)**0.5 + p_tp_com*sig_12[i,1]/R_m12
            elif sig_12[i,1] < 0 and 0 <= np.abs(sig_12[i,2]/sig_12[i,1]) <= np.abs(tau_21c/R_tt_A): # Mode C
                f_E_IFF[i] = ((sig_12[i,2] / (2*(1+p_tt_com)*R_m12))**2 + (sig_12[i,1]/R_m2D)**2) * R_m2D/(-sig_12[i,1])

        #print("Puck FF: " + str(f_E_FF[i]) + " Puck IFF: " + str(f_E_IFF[i]))           
        
    #print(sig_12)

    return f_E_FF, f_E_IFF

def CLT_Stress_TW(stack, Q_0, ABD, F, R_m1Z, R_m2Z, R_m1D, R_m2D, R_m12):
    # Strain in global coordinate system
    eps_kap_0 = np.linalg.inv(ABD) @ F
    eps_0 = eps_kap_0[:3]
    kap_0 = eps_kap_0[3:]

    # Global and local strain and local stress in each layer
    eps_xy = np.zeros((stack.shape[0],eps_0.shape[0]))
    eps_12 = np.zeros((stack.shape[0],eps_0.shape[0]))
    sig_12 = np.zeros((stack.shape[0],eps_0.shape[0]))

    # Stress exposure for each layer with Tsai-Wu criterion
    f_E_TW = np.zeros(stack.shape[0])

    for i in range(stack.shape[0]):
        eps_xy[i,:] = eps_0 + stack[i,2] * kap_0

        T_e = calc_T_e(stack[i])

        eps_12[i,:] = T_e @ eps_xy[i,:]

        sig_12[i,:] = Q_0 @ eps_12[i,:]

        # Tsai-Wu failure criterion
        f_E_TW[i] = sig_12[i,0]**2/(R_m1Z*R_m1D) - sig_12[i,0]*sig_12[i,1]/(R_m1Z*R_m1D*R_m2Z*R_m2D)**0.5 + sig_12[i,1]**2/(R_m2Z*R_m2D) + sig_12[i,2]**2/R_m12**2 + sig_12[i,0]*(1/R_m1Z-1/R_m1D) + sig_12[i,1]*(1/R_m2Z-1/R_m2D)

        #print("Tsai-Wu: " + str(f_E_TW[i]))       
        
    #print(sig_12)

    return f_E_TW