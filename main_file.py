import function_file as ff
import sympy as sp
import time

start = time.time()

# ------ airfoil data
airfoil = 'NACA4412'
N, r, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane = ff.read_data(airfoil)
# ------ free stream velocity
free_velocity = 20.0
free_aoa = 2.0
free_aoa = sp.rad(free_aoa).evalf()
# ------ plunging parameters
pl_amplitude = 2
pl_frequency = 3
# ------ time step
time_step = 0.02
current_time = 0.00
iteration = 50
# ------ new vortex
distance = 0.05
angle = 0
# ------ data store
circulation = []
le_vortex_strength = []
vortex_z = []
vortex_v = []
vortex_u = []
# derivatives

v_func, v_derive = ff.get_v_function(Gkn, r, center_circle)

func = v_func - trailing_edge_z
derive_func = v_derive
'''
tolerance = 1e-8
iteration = 50
v_value = ff.newton(func, derive_func, center_circle + r, tolerance, iteration)
u_value = ff.get_u_value(v_value, center_circle, r)
'''

####
