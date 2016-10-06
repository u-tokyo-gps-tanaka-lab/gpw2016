# -*- coding: utf-8 -*-
# dc_simulation.py
# Katsuki Ohto

import math
import numpy as np

import dc

# simulation constant
FRICTION_RINK = 12.009216 # ( = g * mu )

DX_V_R = 0.001386965639075
DY_V_R = 0.041588442394742

#DR_V_R = np.hypot(DX_V_R, DY_V_R)
DR_V_R = 0.041611563471

PHI = math.atan2(DY_V_R, DX_V_R)

ALPHA = 0.03333564

#B = (math.cos(SPIRAL_ALPHA) ** 2) / (1 - (math.cos(SPIRAL_ALPHA) ** 2))
B = 29.986811440344486

#A = DX_V_R * (dc.calc_r(dc.XY_TEE, dc.XY_THROW) / DR_V_R) / (math.cos(PHI) * math.exp(SPIRAL_B * PHI))
A = 3.45628911574e-19

MAT_ALPHA = [np.array([[math.cos(ALPHA), -math.sin(ALPHA)],
                       [math.sin(ALPHA), math.cos(ALPHA)]]),
             np.array([[math.cos(-ALPHA), -math.sin(-ALPHA)],
                       [math.sin(-ALPHA), math.cos(-ALPHA)]])]

# curve function
def calc_r_by_theta(theta):
    return A * math.exp(B * theta)

def calc_theta_by_r(r):
    return math.log(r / A) / B

def calc_r_by_v(v):
    return DR_V_R * v * v

def calc_v_by_r(r):
    return math.sqrt(r / DR_V_R)

def calc_v_by_theta(theta):
    return calc_v_by_r(calc_r_by_theta(theta))

def calc_theta_by_v(v):
    return calc_theta_by_r(calc_r_by_v(v))

def calc_xy_by_vxy(vxy, spin):
    rate = DR_V_R * calc_r(vxy)
    return rate * np.matmul(MAT_ALPHA[spin], vxy)

def calc_goal(xy, vxy, w):
    dxy = calc_xy_by_vxy(vxy, int(w < 0))
    return xy + dxy

def calc_v_by_t(t):
    return t * FRICTION_RINK

def calc_t_by_v(v):
    return v / FRICTION_RINK

def calc_vtheta_by_theta(theta):
    return -(theta + ALPHA)

def calc_vtheta_by_v(v):
    return calc_vtheta_by_theta(calc_theta_by_v(v))

def calc_vxy_by_theta(theta):
    v = calc_v_by_theta(theta)
    vtheta = calc_vtheta_by_theta(theta)
    return np.array([], dtype = np.float32)

def calc_next_theta_by_theta(theta):
    return theta

asleep_distance2_table = np.empty(2 ** 12, dtype = np.float32)
awake_distance2_table = np.empty(2 ** 12, dtype = np.float32)


if __name__ == '__main__':

    # test
    print(DR_V_R)
    
    a = ALPHA
    print(a)
    print(math.cos(a) / math.sqrt(1 - (math.cos(a) ** 2)))
    b = B
    print(math.acos(b / math.sqrt(b * b + 1)))

    print(PHI)

    print(A)

    print(calc_r_by_theta(PHI))

    print(calc_theta_by_r(dc.calc_r(dc.XY_TEE, dc.XY_THROW)))

    print(calc_v_by_theta(1.53))

    print(calc_vtheta_by_v(31.0))





