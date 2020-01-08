#!/usr/bin/env python


def add_tstar_randomnoise(f27in, noise_stdev):

    """
    Add random noise to synthetic t* measurements
    """

    import numpy as np

    # Prepare files
    w = open("f27_plusnoise", "w")
    f = open(f27in, "r")

    # Loop over each line in input file
    for l in f:

        # Write header line
        if l[0] == "1" or l[5:6] == "0":
            w.write(l)

        # Get measurements
        else:
            dataline = [l.rstrip()[i:i+18]
                        for i in range(0, len(l.rstrip()), 18)]
            for n, obs in enumerate(dataline):
                station = obs[0:4]
                pha_wgt = obs[5:7]
                t_star = float(obs[7:])

                # Add noise
                t_star_noise = t_star + np.random.normal(scale=noise_stdev)

                # Now write new tstar value

                write_format = "{:4s} {:2s} {:10.4f}".format(station, pha_wgt,
                                                             t_star_noise)
                if n == len(dataline)-1:
                    w.write("{:}\n".format(write_format))
                else:
                    w.write("{:}".format(write_format))


def add_arrivaltime_randomnoise(f27in, noise_stdev):

    """
    Add random noise to synthetic arrival times based on pick uncertainty
    noise_stdev = [pick_uncert_Q0, pick_uncert_Q1, pick_uncert_Q2,
                   pick_uncert_Q3, pick_uncert_Q4]
    """

    import numpy as np

    # Prepare files
    w = open("f27_plusnoise", "w")
    f = open(f27in, "r")

    # Loop over each line in input file
    for l in f:

        # Write header line and separator line
        if l[0] == "1" or l[0:1] == "0":
            w.write(l)

        # Get measurements
        else:
            dataline = [l.rstrip()[i:i+14]
                        for i in range(0, len(l.rstrip()), 14)]
            for n, obs in enumerate(dataline):
                station = obs[0:4]
                pha_wgt = obs[5:8]
                wgt = int(obs[7:8])
                arrival_time = float(obs[8:])

                # Add noise based on arrival time weight
                arrival_time_noise = (arrival_time + np.random.normal(
                    scale=noise_stdev[wgt]))

                # Now write new arrival time
                write_format = "{:4s} {:3s}{:6.2f}".format(station, pha_wgt,
                                                            arrival_time_noise)
                if n == len(dataline)-1:
                    w.write("{:}\n".format(write_format))
                else:
                    w.write("{:}".format(write_format))
