import numpy as np
import matplotlib.pyplot as plt

from map_tools import *
from map_cycle import *
from map_run import *

plt.rcParams["font.size"]=10
plt.rcParams["lines.markersize"]=3
plt.rcParams["figure.figsize"]=[12, 6]
plt.rcParams["figure.dpi"]=200

### create a Run object with a run number, a range of cycles, the current data format, and ask it to look for rings and to calibrate.
remnant_field_map = Run(648, np.arange(0, 90), dtformat='v3swap', runtype='rings', calib=True)

# /!\ note on maps directory: Run looks for maps in "../maps_bin/"
# /!\ note on calibration: takes by default the first rh0=0, z=0 ring it encounters as a calibration ring for B_rho and B_phi
# and the second one for B_z (it assumes that the fluxgate has been rotated by 90°)
# if no second calibration ring then it arbitrarily calibrates Bz with the first ring, which suppresses m=0 modes

# Cycle object can be extracted from a run by [rho, z] index:
some_ring = remnant_field_map.cycles[70, -39]

# plot example for some ring
fig, ax = some_ring.simplePlot('phi', 'Brho')
plt.show()

# extract Glms up to order l=3 and associated normalized errors
lmax=3
G, G_err = remnant_field_map.getGlm(lmax, source='rho', norm=True)
# /!\ note on 'source': uses by default the z probe for m=0 modes, then 'source' tells which probes to use for other modes.
# G[l][lmax + 1 + m] is G_lm
print("G_10 = {:.2e} +/- {:.2e}".format(G[1][lmax+1+0], G_err[1][lmax+1+0]))

# plot Glm bars of 0.1 width for l=1 and all m
remnant_field_map.plotGm(l=1, w=0.1)
plt.show()
