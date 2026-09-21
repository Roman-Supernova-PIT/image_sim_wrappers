#!/usr/bin/env python
#
# Created Mar 23 2026
# Wrapper for Megan's imagelib toolkit to make png image
# from asdf file
# ===========================================================================

import os, sys


import logging, time, datetime, yaml, glob
import argparse
#import asdf

#from astropy.coordinates import SkyCoord
#from astropy.time import Time
#from astropy.table import Table, vstack, join, MaskedColumn
#from astropy import units as u
#from astropy.visualization import simple_norm
#import copy
#import galsim
#import importlib
#import matplotlib.pyplot as plt
#import numpy as np

import imagelib


# ===========
def get_args():

    parser_inputs = argparse.ArgumentParser()
    
    msg = "Required: wildcard (or name) for asdf images to process"
    parser_inputs.add_argument("wildcard_asdf_name", help=msg, default=None)

    args = parser_inputs.parse_args()

    return args
    # end get_args


# ==================================
if __name__ == "__main__":

    args     = get_args()

    asdf_list = glob.glob(args.wildcard_asdf_name)
    #print(f" Process {asdf_list}")
    for asdf_file in asdf_list:
        png_file = asdf_file.split(".")[0]+".png"
        print(f" Create figure {png_file}")
        imagelib.mkfigure(asdf_file, plotname = png_file) 
    
    # === END: ===
