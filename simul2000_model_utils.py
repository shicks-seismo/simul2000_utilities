#!/usr/bin/env python
""""
Routines for dealing with simul2000 model files
"""


def get_grid_orig_rot(stns_file):
    """
    Find the geographic origin of a simul model

    Parameters
    ----------
    stns_file : Path and filename of STNS file with rotation parameters on
                first line (Unit 02 file)

    Returns
    -------
    lat : Latitude in decimal degrees of grid origin
    lon : Longitude in decimal degrees of grid origin
    rot : Bearing (in degrees clockwise) of grid rotation
    """
    f = open(stns_file, "r")
    for n, l in enumerate(f):
        if n == 0:
            hemi_EW = l[12]
            if hemi_EW == " ":
                lon_mult = -1
            else:
                lon_mult = 1
            hemi_NS = l[2]
            if hemi_NS == "S":
                lat_mult = -1
            else:
                lat_mult = 1
            lat = lat_mult * (int(l.split()[0]) + float(l.split()[1])/60)
            lon = lon_mult * (int(l.split()[2]) + float(l.split()[3])/60)
            rot = 360 - float(l.split()[4])
            return(lat, lon, rot)


def convert_between_diff_grids(stns_file_in, mod_file_in, mod_file_out,
                               stns_file_out, meth="2D"):
    """
    Convert between 3-D simul2000 models with different origins and grids.
    Generates output MOD (for simul input) and VELOMOD (for plotting) files

    Parameters
    ----------
    stns_file_in : Path and filename of STNS file with rotation parameters for
                   input grid on first line (Unit 02 file)
    mod_file_in : Path and filename of input MOD file
    stns_file_out : Path and filename of STNS file with rotation parameters for
                    output grid on first line (Unit 02 file)
    mod_file_out : Path and filename of output MOD file
                   (containing grid to convert to)
    """
    from scipy.interpolate import LinearNDInterpolator
    orig_lat_in, orig_lon_in, rot_in = get_grid_orig_rot(stns_file_in)
    _, _, _, nodes_in, _, vals_in = model2list(
            mod_file_in, stns_file=stns_file_in, out_coords="GEOG")
    v_interp = LinearNDInterpolator(nodes_in, vals_in)
    z_nodes_out, y_nodes_out, x_nodes_out, nodes_out, _, _ = model2list(
            mod_file_out, stns_file=stns_file_out, out_coords="GEOG")
    orig_lat_out, orig_lon_out, rot_out = get_grid_orig_rot(stns_file_out)
    outfile_MOD = "converted_MOD"
    outfile_VELOMOD = "converted_VELOMOD"
    write_header(z_nodes_out, y_nodes_out, x_nodes_out,
                 outfile_MOD, outfile_VELOMOD)
    w = open(outfile_MOD, "a")
    v = open(outfile_VELOMOD, "a")
    # Write the interpolated values
    i = 0
    for nz, z in enumerate(z_nodes_out):
        for ny, y in enumerate(y_nodes_out):
            for nx, x in enumerate(x_nodes_out):
                if nx == len(x_nodes_out) - 1:
                    w.write("{:5.2f}\n".format(float(v_interp(nodes_out[i]))))
                    v.write("{:5.2f}\n".format(float(v_interp(nodes_out[i]))))
                else:
                    w.write("{:5.2f}".format(float(v_interp(nodes_out[i]))))
                    v.write("{:5.2f}".format(float(v_interp(nodes_out[i]))))
                i = i+1


def write_header(nodes_z_out, nodes_y_out, nodes_x_out,
                 outfile_MOD, outfile_VELOMOD):
    """
    Write header linest to output MOD (for input) and VELOMOD (for plotting)
    """
    w = open(outfile_MOD, "w")
    v = open(outfile_VELOMOD, "w")
    # Write top line
    w.write(" 1.0 {:2g} {:2g} {:2g}  2\n".format(
            len(nodes_x_out), len(nodes_y_out), len(nodes_z_out)))
    v.write(" 1.0 {:2g} {:2g} {:2g}  2\n".format(
            len(nodes_x_out), len(nodes_y_out), len(nodes_z_out)))
    # Write x nodes
    for n, node in enumerate(nodes_x_out):
        if n == len(nodes_x_out) - 1:
            w.write("{:6.1f}\n".format(node))
            v.write("{:6.1f}\n".format(node))
        else:
            w.write("{:6.1f} ".format(node))
            v.write("{:6.1f}".format(node))

    # Write y nodes
    for n, node in enumerate(nodes_y_out):
        if n == len(nodes_y_out) - 1:
            w.write("{:6.1f}\n".format(node))
            v.write("{:6.1f}\n".format(node))
        else:
            w.write("{:6.1f} ".format(node))
            v.write("{:6.1f}".format(node))

    # Write z nodes
    for n, node in enumerate(nodes_z_out):
        if n == len(nodes_z_out) - 1:
            w.write("{:6.1f}\n".format(node))
            v.write("{:6.1f}\n".format(node))
        else:
            w.write("{:6.1f} ".format(node))
            v.write("{:6.1f}".format(node))
    w.write("  0  0  0\n  0  0  0\n")
    v.write("  0  0  0\n  0  0  0\n")


def xy_2_latlon(z_nodes, y_nodes, x_nodes, orig_lat, orig_lon, rot):
    """
    Grid grid point to geographic latitude and longitude

    Parameters
    ----------
    z_nodes : List of floats containing node positions in z-direction (km)
    y_nodes : List of floats containing node positions in y-direction (km)
    x_nodes : List of floats containing node positions in x-direction (km)
    orig_lat : Latitude in decimal degrees of grid origin
    orig_lon : Longitude in decimal degrees of grid origin
    rot : Bearing (in degrees clockwise) of grid rotation

    Returns
    -------
    points_out : List of lists containing floats of longitude, latitude, depth
                 Longitude and latitude values are in decimal degrees format

    """
    import geopy
    import geopy.distance
    start = geopy.Point(orig_lat, orig_lon)
    points_out = []
    for z in z_nodes:
        for y in y_nodes:
            dy = geopy.distance.distance(kilometers=y)
            shifted_dy = dy.destination(point=start, bearing=rot)
            for x in x_nodes:
                rot_x = rot - 90
                dx = geopy.distance.distance(kilometers=x)
                shifted_dx = dx.destination(point=shifted_dy, bearing=rot_x)
                points_out.append([shifted_dx.longitude, shifted_dx.latitude,
                                  z])
    return(points_out)


def model2list(infile, stns_file="", quan_type="V", phase="P",
               out_coords="GEOG", outfile=""):
    """
    Get list of nodes (x, y, z) and velocities from simul model file.
    Automatically detects whether input model is of MOD (Unit03) or
    VELOMOD (Unit23) type.

    Parameters
    ----------
    infile : file path and name of input simul model file to read
    stnsfile: STNS file with grid origin (optional - for out_coords="GEOG")
    quan_type : "V" = velocity; "Q" = attenuation
    vel_type : "P" = vp; "S" = vs; "PS" = vp/vs ratio
    out_coords : grid-system to output node points with model values
                 "GEOG" (geographic i.e. lat, lon or "CART" i.e. cartesian)

    Returns
    ------
    x_nodes : list of node points in x-direction relative to origin (in km)
    y_nodes : list of node points in y-direction relative to origin (in km)
    z_nodes : list of node points in z-direction relative to origin (in km)
    """

    quan_type = quan_type.upper()
    phase = phase.upper()

    mod_type = "MOD"
    try:
        [[float(i) for i in l.split()]
         for n, l in enumerate(open(infile)) if n == 1]
    except ValueError:
        mod_type = "VELOMOD"
    if mod_type == "MOD" and phase == "S":
        raise Exception("MOD format does not contain Vs values")

    # Read in header information
    if mod_type == "MOD":
        x_nodes = [[float(i) for i in l.split()]
                   for n, l in enumerate(open(infile)) if n == 1][0]
        y_nodes = [[float(i) for i in l.split()]
                   for n, l in enumerate(open(infile)) if n == 2][0]
        z_nodes = [[float(i) for i in l.split()]
                   for n, l in enumerate(open(infile)) if n == 3][0]
    elif mod_type == "VELOMOD":
        x_nodes = [float(l.rstrip()[i:i+6])
                   for n, l in enumerate(open(infile)) if n == 1
                   for i in range(0, len(l.rstrip()), 6)]
        y_nodes = [float(l.rstrip()[i:i+6])
                   for n, l in enumerate(open(infile)) if n == 2
                   for i in range(0, len(l.rstrip()), 6)]
        z_nodes = [float(l.rstrip()[i:i+6])
                   for n, l in enumerate(open(infile)) if n == 3
                   for i in range(0, len(l.rstrip()), 6)]
    # Now read in the model values
    l_start = [n for n, li in enumerate(open(infile))
               if li == "  0  0  0\n"][-1]+1  # Find model starting line
    # Define where to start reading from and to
    if phase == "P" or quan_type == "Q":
        seg = 1
    elif ((phase == "S" and mod_type == "VELOMOD")
          or (phase == "PS" and mod_type == "MOD")):
        n_start = l_start
        n_end = l_start + 1*len(y_nodes)*len(x_nodes)
        seg = 2
    elif phase == "PS" and mod_type == "VELOMOD":
        seg = 3
    n_start = l_start + (seg-1) * len(y_nodes)*len(z_nodes)
    n_end = l_start + seg * len(y_nodes) * len(z_nodes)
    # Read in model values
    vals = [float(j) for n, l in enumerate(open(infile))
            if n_start <= n < n_end for j in l.split()]
    # Get grid origin and rotation if converted to geographic coordinates
    if out_coords == "GEOG":
        orig_lat, orig_lon, rot = get_grid_orig_rot(stns_file)
        geog_grid = xy_2_latlon(z_nodes, y_nodes, x_nodes,
                                orig_lat, orig_lon, rot)
    # Set-up output file if required
    if outfile != "":
        if out_coords == "GEOG":
            w = open(outfile+".latlong", "w")
        else:
            w = open(outfile+".xyz", "w")
    # Now find corresponding grid points and write to file if required
    i = 0
    nodes = []
    nodes_vals = []
    vals_all = []
    for z in z_nodes:
        for y in y_nodes:
            for x in x_nodes:
                if out_coords != "GEOG":
                    nodes_vals.append([x, y, z, vals[i]])
                    nodes.append([x, y, z])
                    vals_all.append(vals[i])
                    if outfile:
                        w.write("{:6.1f} {:6.1f} {:6.1f} {:6.2f}\n".format(
                                x, y, z, vals[i]))
                else:
                    nodes_vals.append([geog_grid[i][0], geog_grid[i][1],
                                      geog_grid[i][2], vals[i]])
                    nodes.append([geog_grid[i][0], geog_grid[i][1],
                                  geog_grid[i][2]])
                    vals_all.append(vals[i])
                    if outfile:
                        w.write("{:8.4f} {:8.4f} {:5.1f} {:6.2f}\n".
                                format(geog_grid[i][0], geog_grid[i][1],
                                       geog_grid[i][2], vals[i]))
                i += 1
    return(x_nodes, y_nodes, z_nodes, nodes, nodes_vals, vals)
