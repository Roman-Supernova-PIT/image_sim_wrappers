#!/usr/bin/env python

import os, sys, glob, argparse, logging, time, datetime
import numpy as np


DEFAULT_WILDCARD_L1_FILES = "SNPIT*L1.asdf"

# ============================
def setup_logging():
    logging.basicConfig(level=logging.INFO,
                        format="[%(levelname)6s ] %(message)s")     
    return

def get_args():
    parser_inputs = argparse.ArgumentParser()
    
    msg = "name (or wildcard) of input L1 asdf file (default={DEFAULT_WILDCARD_L1_FILES})"
    parser_inputs.add_argument("--l1_asdf_file", help = msg, type=str,
                               default=DEFAULT_WILDCARD_L1_FILES)

    msg = "split job (0 to N-1) and total number of split jobs (N)"
    parser_inputs.add_argument("--jobsplit", "-j", help = msg, type=int, default=None, nargs="+");

    msg = "print list of L1 images to process, but do not process"
    parser_inputs.add_argument("--check", "-c", help=msg, action="store_true")    
    
    args = parser_inputs.parse_args()

    return args

def get_L1_asdf_noL2(args):
    l1_asdf_file = args.l1_asdf_file
    
    L1_asdf_list_all = sorted( glob.glob(l1_asdf_file) )    
    n_L1 = len(L1_asdf_list_all)
    logging.info(f"Found {n_L1} files with {L1_asdf_list_all}")

    if args.jobsplit:
        
        ijob = args.jobsplit[0]
        njob = args.jobsplit[1]
        L1_asdf_list_all = [val for i, val in enumerate(L1_asdf_list_all) if i % njob == ijob]
        n_L1 = len(L1_asdf_list_all)
        logging.info(f"APPLY SPLIT: ijob = {ijob} of {njob}")
        logging.info(f"AFTER SPLIT: Found {n_L1} files with {L1_asdf_list_all}")  

    logging.info(f"N_EXPECT_ASDF:  {n_L1}")
    
    L1_asdf_list = []
    for L1_asdf_file in L1_asdf_list_all:
        L2_asdf_file, L2_png_file = get_L2_filenames(L1_asdf_file)
        exist_L2 = os.path.exists(L2_asdf_file)
        if exist_L2 :
            logging.info(f" {L1_asdf_file} -> L2 exists")
        else:
            logging.info(f" {L1_asdf_file} -> L2 does NOT exist")
            L1_asdf_list.append(L1_asdf_file)
            
    return L1_asdf_list

def get_L2_filenames(L1_asdf_file):
    # return name of L2 asdf file and L2 png file
    L2_asdf_file = L1_asdf_file.replace('L1','L2')
    L2_png_file  = L2_asdf_file.replace('asdf','png')
    return L2_asdf_file, L2_png_file

def print_flashy_header(msg):

    logging.info(f"")
    logging.info(f"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    logging.info(f"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    logging.info(f"@@@@@@@ {msg}")
    logging.info(f"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    logging.info(f"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")    
    logging.info(f"")

    return

def run_romancal(args, L1_asdf_file):

    L2_asdf_file, L2_png_file = get_L2_filenames(L1_asdf_file)

    msg_run = f"Run romancal on {L1_asdf_file}"

    if args.check:
        logging.info(f"  {msg_run}")
        return

    t_start = time.time()
    # - - - - - - 
    # print flashy header to help visually find start of verbose L2 logging
    print_flashy_header(msg_run)
    
    result = ExposurePipeline.call(L1_asdf_file)
    
    print_flashy_header(f"Export L2 product to {L2_asdf_file}")
    result.to_asdf(L2_asdf_file)
    
    logging.info(f" Create figure {L2_png_file}")

    imagelib.mkfigure(L2_asdf_file, plotname = L2_png_file)  

    t_end = time.time()
    t_sim = t_end - t_start # seconds                                                                   
    logging.info(f"CPU({L2_asdf_file}):  {t_sim:.0f} seconds ")
    
    return

# ==================================
if __name__ == "__main__":

    setup_logging()
    args  = get_args()

    # fetch list of L1 asdf files that don't have associated L2
    L1_asdf_noL2_list = get_L1_asdf_noL2(args)
    logging.info('')
    
    if not args.check:
        logging.info(f"Import romancal and imagelib ...")
        from romancal.pipeline import ExposurePipeline
        import imagelib
    
    for L1_asdf_file in L1_asdf_noL2_list:
        run_romancal(args, L1_asdf_file)
    
    logging.info('\nDone.')

    # === END: ===
