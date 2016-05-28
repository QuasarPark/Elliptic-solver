from src.math_h import *

I_s = 1.44e-17

V_T = sci.Boltzmann*300/sci.e

I = I_s*(math.exp(0.5/V_T)-1)

print(V_T)
print("Answer of Problem 1 is", I)

def func1(x):
    return x**2 - 0.75*x + 0.01

ans1 = so.newton(func1, x0=1.0, maxiter=10000, tol=1e-20)
print(ans1)