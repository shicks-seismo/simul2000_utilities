#!/usr/bin/env python


def mk_nodes_file(velomod, mod):
    """
    Make nodes-data file read by res2spread
    """

    w = open("nodes-data", "w")
    for n, l in enumerate(open(velomod)):
        if n < 4:
            w.write(l)
    for n, l in enumerate(open(mod)):
        if n >= 4:
            w.write(l)
    w.close()


def add_bound_nodes_resindices(in_file, out_file):
    """
    Add in boundary nodes to res2spread -R resolution row-grid indices output
    """

    from simul2000_grid_utils import grid_vals_2_list

    # Add in bounding node points
    nodes_x = [[float(l.rstrip()[i:i+6]) for i in range(
        0, len(l.rstrip()), 6)]
        for n, l in enumerate(open(in_file)) if n == 1][0]
    nodes_y = [[float(l.rstrip()[i:i+6]) for i in range(
        0, len(l.rstrip()), 6)]
        for n, l in enumerate(open(in_file)) if n == 2][0]
    nodes_z = [[float(l.rstrip()[i:i+6]) for i in range(
        0, len(l.rstrip()), 6)]
        for n, l in enumerate(open(in_file)) if n == 3][0]
    nodes_x[0:0] = [-740.0]
    nodes_x.append(740.0)
    nodes_y[0:0] = [-740.0]
    nodes_y.append(740.0)
    nodes_z[0:0] = [-740.0]
    nodes_z.append(740.0)

    # Now make file
    out = open(out_file, "w")
    out.write(" 1.0 {:2g} {:2g} {:2g}  2\n".format(
        len(nodes_x), len(nodes_y), len(nodes_z)))

    for x, x_node in enumerate(nodes_x):
        if x < len(nodes_x)-1:
            out.write("{:6.1f}".format(x_node))
        else:
            out.write("{:6.1f}\n".format(x_node))
    for y, y_node in enumerate(nodes_y):
        if y < len(nodes_y)-1:
            out.write("{:6.1f}".format(y_node))
        else:
            out.write("{:6.1f}\n".format(y_node))
    for z, z_node in enumerate(nodes_z):
        if z < len(nodes_z)-1:
            out.write("{:6.1f}".format(z_node))
        else:
            out.write("{:6.1f}\n".format(z_node))

    out.write("  0  0  0\n")

    # Get resolution values
    res_vals = grid_vals_2_list(file_in="res.tmp", point_width=5, min_line=5)

    # Now write values
    i = 0
    for z, z_node in enumerate(nodes_z):
        for y, y_node in enumerate(nodes_y):
            for x, x_node in enumerate(nodes_x):
                if (x == 0 or x == len(nodes_x)-1 or y == 0 or
                        y == len(nodes_y)-1 or z == 0 or
                        z == len(nodes_z)-1):
                    resnew = 0.00
                else:
                    resnew = res_vals[i]
                    i = i + 1

                if x == len(nodes_x)-1:
                    out.write("{:5.2f}\n".format(resnew))
                else:
                    out.write("{:5.2f}".format(resnew))
