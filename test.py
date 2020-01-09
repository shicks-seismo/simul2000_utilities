#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 13:42:26 2020

@author: sph1r17
"""

from simul2000_grid_utils import convert_between_diff_grids
convert_between_diff_grids("STNS_in", "mod.med15_15p.out", "MOD_Lidong_nofix",
                           "STNS", meth="2D")