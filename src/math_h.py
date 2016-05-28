import math
import numpy as np
import scipy.constants as sci
import scipy.special as sp
import scipy.optimize as so


class sc:  # Scientific constants
    m_e = sci.electron_mass
    e_0 = sci.e
    hbar = sci.hbar
    eV = sci.electron_volt
    meV = 10e-3*sci.electron_volt
    pi = sci.pi


def C_x(a, b, angle):
    return (a ** 2) * (math.cos(angle) ** 2) + (b ** 2) * (math.sin(angle) ** 2)


def C_y(a, b, angle):
    return (a ** 2) * (math.sin(angle) ** 2) + (b ** 2) * (math.cos(angle) ** 2)


def C_xy(a, b, angle):
    return (a ** 2 - b ** 2) * (math.sin(2 * angle))


def realignAng(M_xx, M_yy, a, b, angle):

    C_X = C_x(a, b, angle)
    C_Y = C_y(a, b, angle)
    C_XY = C_xy(a, b, angle)
    return 0.5 * math.atan((M_xx * M_yy * C_XY) /
                           ((M_xx ** 2) * C_X - (M_yy ** 2) * C_Y))


def transformMaxtrix(M_xx, M_yy, a, b, angle):

    phi = realignAng(M_xx, M_yy, a, b, angle)

    realign_mat = np.zeros((2, 2))
    realign_mat[0][0] = realign_mat[1][1] = math.cos(phi)
    realign_mat[1][0] = math.sin(phi); realign_mat[0][1] = -1.0 * math.sin(phi)

    mass_mat = np.zeros((2, 2))
    mass_mat[0][0] = M_xx; mass_mat[1][1] = M_yy

    angle_mat = np.zeros((2, 2))
    angle_mat[0][0] = realign_mat[1][1] = math.cos(angle)
    angle_mat[1][0] = -1.0 * math.sin(angle)
    angle_mat[0][1] =  math.sin(angle)

    t_matrix = np.linalg.multi_dot([realign_mat, mass_mat, angle_mat])
    return t_matrix # return transform matrix


def major_eff(a, b, M_xx, M_yy, angle):

    phi = realignAng(M_xx, M_yy, a, b, angle)

    C_X = C_x(a, b, angle)
    C_Y = C_y(a, b, angle)
    C_XY = C_xy(a, b, angle)
    A = (a * b * M_xx * M_yy) / ((M_xx ** 2) * C_X * (math.sin(phi) ** 2) +
        (M_yy ** 2) * C_Y * (math.cos(phi) ** 2) - C_XY * math.sin(2 * phi))
    return A


def minor_eff(a, b, M_xx, M_yy):

    phi = realignAng(M_xx, M_yy, a, b, angle)

    C_X = C_x(a, b, angle)
    C_Y = C_y(a, b, angle)
    C_XY = C_xy(a, b, angle)
    B = (a * b * M_xx * M_yy) / ((M_xx ** 2) * C_X * (math.cos(phi) ** 2) +
        (M_yy ** 2) * C_Y * (math.sin(phi) ** 2) + C_XY * math.sin(2 * phi))
    return B


def Xi_boundary(A_eff, B_eff):
    return math.atan(B_eff/A_eff)











