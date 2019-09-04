import numpy as np
import sympy as sp
import cmath


def read_data(heading):
    heading = str(heading) + '_data.txt'
    file1 = open(heading, 'r')
    line = file1.readlines()

    N = int(line[1])
    r = float(line[3])
    center_circle = complex(line[5].replace('i', 'j'))
    trailing_edge_z = complex(line[7].replace('i', 'j'))
    Gkn = [complex(index.replace('i', 'j')) for index in line[9][:len(line[9]) - 3].split(' ')]
    z_plane = np.asarray([complex(index.replace('i', 'j')) for index in line[11][:len(line[11]) - 3].split(' ')])
    v_plane = np.asarray([complex(index.replace('i', 'j')) for index in line[13][:len(line[13]) - 3].split(' ')])
    u_plane = np.subtract(v_plane, center_circle) / r

    return N, r, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane


def get_v_function(Gkn, radius, center):
    n = np.arange(1, len(Gkn) + 1)
    v = sp.symbols('v', real=False)
    z = v + sum(Gkn * (radius / (v - center)) ** n)
    z_dev = sp.diff(z, v)
    return z, z_dev


def newton(f, Df, x0, epsilon, max_iter):
    """Approximate solution of f(x)=0 by Newton's method.

    Parameters
    ----------
    f : function
        Function for which we are searching for a solution f(x)=0.
    Df : function
        Derivative of f(x).
    x0 : number
        Initial guess for a solution f(x)=0.
    epsilon : number
        Stopping criteria is abs(f(x)) < epsilon.
    max_iter : integer
        Maximum number of iterations of Newton's method.

    Returns
    -------
    xn : number
        Implement Newton's method: compute the linear approximation
        of f(x) at xn and find x intercept by the formula
            x = xn - f(xn)/Df(xn)
        Continue until abs(f(xn)) < epsilon and return xn.
        If Df(xn) == 0, return None. If the number of iterations
        exceeds max_iter, then return None.

    Examples
    --------
    #>>> f = lambda x: x**2 - x - 1
    #>>> Df = lambda x: 2*x - 1
    #>>> newton(f,Df,1,1e-8,10)
    Found solution after 5 iterations.
    1.618033988749989
    """

    v = sp.symbols('v', real=False)
    xn = x0
    for n in range(0, max_iter):
        fxn = f.subs([(v, xn)]).evalf()  # f(xn)
        if abs(fxn) < epsilon:
            # print('Found solution after', n, 'iterations.')
            return xn
        Dfxn = Df.subs([(v, xn)]).evalf()  # Df(xn)
        if Dfxn == 0:
            # print('Zero derivative. No solution found.')
            return None
        xn = (xn - fxn / Dfxn).evalf()
    print('Exceeded maximum iterations. No solution found.')
    return None


def get_u_value(v_value, center, r):
    u_value = (v_value - center) / r
    return u_value


def get_v_value(u_value, center, r):
    v_value = u_value * r + center
    return v_value
