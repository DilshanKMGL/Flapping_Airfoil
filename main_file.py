import function_file as ff
import sympy as sp
import time
import cmath

start = time.time()

# ------ airfoil data
airfoil = 'NACA4412'
N, r, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane = ff.read_data(airfoil)
# ------ free stream velocity
free_velocity = 20.0
free_aoa = 0.0
free_aoa = sp.rad(free_aoa).evalf()
# ------ plunging parameters
pl_amplitude = 2
pl_frequency = 3
# ------ time step
time_step = 0.02
current_time = 0.00
iteration = 50
# ------ new vortex
distance = 0.005
angle = 0
# ------ data store
circulation = []
te_vortex_strength = []
te_vortex_z = []
te_vortex_u = []
trailing_edge_pre_pos = trailing_edge_z
# ------ derivatives
v, u = sp.symbols('v, u', real=False)
z_fun, v_fun = ff.get_v_function(Gkn, r, center_circle)
first_derivative = 1 / sp.diff(v_fun, u)
second_derivative = sp.diff(1 / sp.diff(v_fun, u), u) / sp.diff(v_fun, u)
# ------ trailing edge
#te_v_value_init = ff.newton(z_fun-trailing_edge_z, sp.diff(z_fun, v), center_circle + r, 1e-8, 50)
#te_u_value_init = ff.get_u_value(te_v_value_init, center_circle, r)

# ---- iteration code

fun_1 = ff.plunging(pl_amplitude, pl_frequency, current_time)
func = cmath.polar(fun_1[1] + free_velocity)
velocity = func[0]
aoa = free_aoa + func[1]

search_point = center_circle + r
# te_value is a list [z value, u value]
te_value = ff.get_trailing_edge(z_fun, trailing_edge_pre_pos, search_point, center_circle, r, fun_1[0])
trailing_edge_pre_pos = te_value[0]
new_vortex_z, new_vortex_u = ff.new_vortex_position(z_fun, trailing_edge_z, search_point, angle, distance, center_circle, r)
te_vortex_u.append(new_vortex_u)
te_vortex_z.append(new_vortex_z)

# create function to make complex potential function
vortex_function = ff.create_circulation(center_circle).subs(v, (u * r)/center_circle).evalf()
le_vortex_sum = sum(te_vortex_strength)
vortex_function += ff.create_vortex(te_vortex_u[-1], le_vortex_sum)
vortex_function += ff.get_freestream(velocity, aoa)
temp_func = [ff.get_vortex(te_vortex_u[vortex], te_vortex_strength[vortex])
                 for vortex in range(len(te_vortex_strength))]
vortex_function += sum(temp_func)
new_circulation = ff.calculate_circulation(vortex_function, te_value[1])
circulation.append(sp.re(new_circulation))
te_vortex_strength.append(sp.re(-le_vortex_sum - new_circulation))
print(circulation)
print(te_vortex_strength)





'''
v_func, v_derive = ff.get_v_function(Gkn, r, center_circle)
func = v_func - trailing_edge_z
derive_func = v_derive
tolerance = 1e-8
iteration = 50
v_value = ff.newton(func, derive_func, center_circle + r, 1e-8, 50)
u_value = ff.get_u_value(v_value, center_circle, r)
'''