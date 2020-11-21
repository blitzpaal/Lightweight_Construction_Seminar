# Program to calculate strength of rotor shaft
# Author: Paul Straehle

import numpy as np
import matplotlib.pyplot as plt

# Geometry and Load data
l = 500 # length in mm
d_i = 45 # inner diameter in mm
T_n = 400*1000 # net torque in Nmm
RF = 2 # reserve factor
T = T_n*RF # torque with reserve factor in mm
N = 5000*np.pi/60 # roational speed in rad/s

# Material data
t_ply = 0.125 # ply thickness in mm
E_11 = 126000 # Longitudinal tensile modulus in MPa
E_22 = 9000 # Transverse tensile modulus in MPa
G_12 = 4600 # Shear modulus in MPa
v_12 = 0.3 # Poisson’s ratio 1
v_21 = v_12*E_22/E_11 # Poisson’s ratio 2
rho = 1570 * 10**-9 # Density in kg/mm^3
R_m1Z = 1250 # Longitudinal tensile strength in MPa
R_m2Z = 50 # Transverse tensile strength in MPa
R_m1D = 1260 # Longitudinal compressive strength in MPa
R_m2D = 210 # Transverse compressive strength in MPa
R_m12 = 84 # Shear strength in MPa

# Laminate stack
# [angle, thickness, height] Laminate buildup from bottom to top with orientation (in degree), thickness (in mm) and z-position (in mm, calculated later) of each layer
stack = np.array([[0,t_ply,0],[45,t_ply,0],[45,t_ply,0],[90,t_ply,0],[-45,t_ply,0],[-45,t_ply,0],
[-45,t_ply,0],[-45,t_ply,0],[90,t_ply,0],[45,t_ply,0],[45,t_ply,0],[0,t_ply,0]],float)
stack[:,0] = np.deg2rad(stack[:,0]) # transform ply angle from degree to radian

# Stiffness matrix of UD-Layer
Q_0 = np.array([[E_11/(1-v_12*v_21),v_21*E_11/(1-v_12*v_21),0],
[v_12*E_22/(1-v_12*v_21),E_22/(1-v_12*v_21),0],
[0,0,G_12]])

# Transforme Laminate stiffness matrix of each layer
Q_trans = np.empty((Q_0.shape[0],Q_0.shape[1],stack.shape[0]))

z_ref = np.sum(stack[:,1]) / 2 # z-position of reference plane, starting from bottom
for i in range(stack.shape[0]):
    if i==0:
        stack[i,2] = -z_ref+stack[i,1]
    elif i>0:
        stack[i,2] = stack[i-1,2]+stack[i,1]

    T_s = np.array([[np.cos(stack[i,0])**2, np.sin(stack[i,0])**2, 2*np.sin(stack[i,0])*np.cos(stack[i,0])],
    [np.sin(stack[i,0])**2, np.cos(stack[i,0])**2, -2*np.sin(stack[i,0])*np.cos(stack[i,0])],
    [-np.sin(stack[i,0])*np.cos(stack[i,0]), np.sin(stack[i,0])*np.cos(stack[i,0]), np.cos(stack[i,0])**2-np.sin(stack[i,0])**2]])

    T_e = np.array([[np.cos(stack[i,0])**2, np.sin(stack[i,0])**2, np.sin(stack[i,0])*np.cos(stack[i,0])],
    [np.sin(stack[i,0])**2, np.cos(stack[i,0])**2, -np.sin(stack[i,0])*np.cos(stack[i,0])],
    [-2*np.sin(stack[i,0])*np.cos(stack[i,0]), 2*np.sin(stack[i,0])*np.cos(stack[i,0]), np.cos(stack[i,0])**2-np.sin(stack[i,0])**2]])

    Q_trans[:,:,i] = np.linalg.inv(T_s).dot(Q_0).dot(T_e)

# ABD Matrix calculation
A = np.sum(Q_trans*stack[:,1],2)
B = np.sum(Q_trans*stack[:,1]*(stack[:,2]-stack[:,1]/2),2)
D = np.sum(Q_trans*(stack[:,1]**3/12+stack[:,1]*(stack[:,2]-stack[:,1]/2)**2),2)