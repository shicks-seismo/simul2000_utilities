#!/usr/bin/env python

import numpy as np
import scipy.interpolate as spint
import pandas as pd

# Start of parameters to define
NEW_X_NODES = [-740.0, -240.0, -210.0, -180.0, -150.0, -120.0, -90.0, -65.0,
               -40.0, -15.0, 10.0, 35.0, 65.0, 95.0, 125.0, 155.0, 185.0,
               225.0, 300.0, 740.0]
NEW_Y_NODES = [-740.0, -100.0, 100.0, 740.0]
NEW_Z_NODES = [0.0]
MOD_IN = "vel.xyz"
# End of parameters to define

# Read in xzv file and then sort
in_num_mod = pd.read_csv(MOD_IN, usecols=[0, 1, 2],
                         delim_whitespace=True, header=None)
in_num_mod = in_num_mod.sort_values(by=[2, 0])

# Define grid points
x = np.unique(np.array(list(in_num_mod[0])))
z = np.unique(np.array(list(in_num_mod[1])))

v = np.array(list(in_num_mod[2]))
v = np.reshape(v, (len(z), len(x)))
# Interpolate velocity
vel_interp = spint.interp2d(x=x, y=z, z=v,
                            kind="linear")

# Define node points for output
nodes_x_out = NEW_X_NODES
nodes_y_out = NEW_Y_NODES
nodes_z_out = NEW_Z_NODES

# Define node points for averaging and prepare output MOD
MOD_OUT = open("MOD_ave", "w")

# Write first line
MOD_OUT.write(" 1.0{:3g} {:2g} {:2g}  2\n".format(
    len(nodes_x_out), len(nodes_y_out), len(nodes_z_out)))

# Write x nodes
for n, node in enumerate(nodes_x_out):
    if n == len(nodes_x_out) - 1:
        MOD_OUT.write("{:6.1f}\n".format(node))
    else:
        MOD_OUT.write("{:6.1f} ".format(node))

# Write y nodes
for n, node in enumerate(nodes_y_out):
    if n == len(nodes_y_out) - 1:
        MOD_OUT.write("{:6.1f}\n".format(node))
    else:
        MOD_OUT.write("{:6.1f} ".format(node))

# Write z nodes
for n, node in enumerate(nodes_z_out):
    if n == len(nodes_z_out) - 1:
        MOD_OUT.write("{:6.1f}\n".format(node))
    else:
        MOD_OUT.write("{:6.1f} ".format(node))

# Write zero lines
MOD_OUT.write(" {:2g} {:2g} {:2g}\n".format(0, 0, 0))
MOD_OUT.write(" {:2g} {:2g} {:2g}\n".format(0, 0, 0))

# Write velocities to file
for nz, z in enumerate(nodes_z_out):
    for ny, y in enumerate(nodes_y_out):
        for nx, x in enumerate(nodes_x_out):
            if nx == len(nodes_x_out) - 1:
                MOD_OUT.write("{:5.2f}\n".format(vel_interp([x, y, z])[0]))
            else:
                MOD_OUT.write("{:5.2f}".format(vel_interp([x, y, z])[0]))

MOD_OUT.close()
