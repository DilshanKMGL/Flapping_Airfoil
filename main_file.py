import function_file as ff
import sympy as sp
import numpy as np
from scipy.optimize import newton
from scipy.optimize import fsolve

airfoil = 'NACA4412'
N, r, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane = ff.read_data(airfoil)

v_func, v_derive = ff.get_v_function(Gkn, r, center_circle)
tolerance = 1e-10
iteration = 50

z = sp.symbols('z', real=False)
func = v_func - trailing_edge_z
derive_func = v_derive
v_value = ff.newton(func, derive_func, center_circle + r, tolerance, iteration)
u_value = ff.get_u_value(v_value, center_circle, r)
print(v_value)
print(u_value)

####
