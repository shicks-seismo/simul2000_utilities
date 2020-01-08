#!/usr/bin/env python


def get_node_points(mod):
    """
    Get list of nodes (x, y, z) from MOD / velomod file
    """

    nodes_x_in = [[float(i) for i in l.split()]
                  for n, l in enumerate(open(mod)) if n == 1][0]
    nodes_y_in = [[float(i) for i in l.split()]
                  for n, l in enumerate(open(mod)) if n == 2][0]
    nodes_z_in = [[float(i) for i in l.split()]
                  for n, l in enumerate(open(mod)) if n == 3][0]

    return(nodes_x_in, nodes_y_in, nodes_z_in)


def grid_vals_2_list(file_in, point_width, min_line, ny, nz, section):
    """
    Convert values on MOD / velomod grid to single list
    point_width = fixed width of each value
    Min_line = First line to start on (starting at 0)
    ny = Number of y nodes
    nz = Number of z nodes
    section - which block to get - e.g. section=1=> vp; section=2=>vp/vs
    """
    if section == 1:
        vals = [float(l.rstrip()[i:i+point_width])
                for n, l in enumerate(open(file_in))
                if n >= min_line and n < (ny*nz + min_line)
                for i in range(0, len(l.rstrip()), point_width)]
    elif section == 2:
        vals = [float(l.rstrip()[i:i+point_width])
                for n, l in enumerate(open(file_in))
                if n >= (ny*nz + min_line) and n < (2*ny*nz + min_line)
                for i in range(0, len(l.rstrip()), point_width)]

    return vals
