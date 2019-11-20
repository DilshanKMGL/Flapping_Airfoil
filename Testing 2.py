import numpy as np
from matplotlib import pyplot as plt

te_vortex_u = np.array([5,10,15,20])
le_vortex_u = np.array([2,4,6])
te_vortex_strength = np.array([1,2,3,4])
le_vortex_strength = np.array([3,6,9])

le_ss = np.tile(le_vortex_strength, (len(te_vortex_strength), 1))  # triling edge vortex strength
le_uu_conj = np.conjugate(np.tile(le_vortex_u, (len(te_vortex_u), 1)))  # conjugate for calculation
le_uu = np.tile(le_vortex_u, (len(te_vortex_u), 1))
te_u = np.tile(te_vortex_u, (len(le_vortex_u), 1)).transpose()

print(le_ss)
print(te_u)
print(le_uu)
print(sum(le_ss*(le_u-le_uu)))