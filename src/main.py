from src.math_h import *

# -------------------------------------------------------------------------- #
#
#
#  2D Schrodinger equation solver for elliptic boundary
#
#
# -------------------------------------------------------------------------- #


# input valuable
# ellipses x^2/a^2 + y^2/b^2 = 1 for inner boundary, x^2/c^2 + y^2/d^2 = 1 for outer boundary

majorAxisLength = 5e-9
a_b_ratio = 1.1
angle = (sc.pi/180)*60

a_major = majorAxisLength
b_minor = majorAxisLength / a_b_ratio
oxide_thickness = 2e-9

c_major = a_major + oxide_thickness
d_minor = b_minor + oxide_thickness

m_xx = 0.910 # Relative electron mass for xx direction
m_yy = 0.160 # Relative electron mass for yy direction

e_xx = m_xx * sc.m_e  # effective electron mass for x-axis (SI unit)
e_yy = m_yy * sc.m_e  # effective electron mass for y-axis (SI unit)

M_xx = math.sqrt(m_xx)  # Magnification coefficient of x-axis
M_yy = math.sqrt(m_yy)  # Magnification coefficient of y-axis


def length_conv(majorAxisLength, a_b_ratio, m_xx, m_yy, angle):

    a_major = majorAxisLength
    b_minor = majorAxisLength / a_b_ratio

    M_xx = math.sqrt(m_xx)  # Magnification coefficient of x-axis
    M_yy = math.sqrt(m_yy)  # Magnification coefficient of y-axis

    A = major_eff(a, b, M_xx, M_yy)
    B = minor_eff(a, b, M_xx, M_yy)

    return [A, B]




def solve_elliptic(a, b, m_xx, m_yy, angle, t_ox):

    def solve_boundary():

         se_conf = sp.mathieu_even_coef(m, q)
         ce_conf = sp.mathieu_odd_coef(m, q)

         def se_bnd_cond(m, q, z):
             return (sp.mathieu_modcem1(m, q, z)[1] + sp.mathieu_modcem2(m, q, z)[1]
             /sp.mathieu_modcem1(m, q, z)[0] + sp.mathieu_modcem2(m, q, z)[0])

         def ce_bnd_solve(m, q, z):
             return (sp.mathieu_modsem1(m, q, z)[1] + sp.mathieu_modsem2(m, q, z)[1]
             /sp.mathieu_modsem1(m, q, z)[0] + sp.mathieu_modsem2(m, q, z)[0])


