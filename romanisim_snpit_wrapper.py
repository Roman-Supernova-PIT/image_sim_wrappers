#!/usr/bin/env python
#
# for help desjard@stsci.edu  (Tyler)
#
# Galaxy catalog columns from notebook:
# gal cols = ['ra', 'dec', 'type', 'n', 'half_light_radius', 'pa', 'ba',
#              'F062', 'F087', 'F106', 'F129', 'F146', 'F158', 'F184', 'F213']
#
# MA table:
#  https://roman-docs.stsci.edu/roman-instruments/the-wide-field-instrument/observing-with-the-wfi/wfi-multiaccum-ma-tables
#
# Questions;
#   - is lensing applied ?
#   - is Galactic extinction applied ?
#   - how to retrieve subset of objects overlaid on image ?
#   - slurm-like system to distribute jobs?  CPU limit per user?
#   - is there a dumb-fast option to quickly test infrastructure ?
#
#
# Jul 10 2026: begin minor changes to embed in bigger submit-pipeline;
#            e.g., add several os.path.expandvars
#
# ===========================================================================

import os, sys
from contextlib import redirect_stdout

import logging, time, datetime, yaml, glob, subprocess, gzip
import argparse
import asdf

from   astropy.coordinates import SkyCoord
from   astropy.time import Time
from   astropy.table import Table, vstack, join, MaskedColumn
from   astropy import units as u
from   astropy.visualization import simple_norm
from   astropy.table import unique, join

import copy
import galsim
import importlib
#import matplotlib.pyplot as plt
import numpy as np




# ============================================================================
# define hard-coded inputs ... perhaps later will read from input config file


WRITE_TRUTH_TABLE_ONLY = False  # if true, skip sim after writing truth table

ROTATE_GALAXIES_90DEG = True  # Apr 25 2026

BAND_LIST_SNANASIM = [ 'R062-R','Z087-Z','Y106-Y','J129-J','H158-H','F184-F','K213-K' ]
BAND_LIST_HOSTLIB  = [ 'R_obs', 'Z_obs', 'Y_obs', 'J_obs', 'H_obs', 'F_obs', 'K_obs' ]
BAND_LIST_SOC      = [ 'F062',  'F087',  'F106',  'F129',  'F158',  'F184',  'F213' ]
BAND_LIST_INPUT    = [ b[0] for b in BAND_LIST_HOSTLIB ]  # single-char represenation for user input

# note that a_Sersic and b_Sersic are replaced with a0,b0 or a1,b1
HOSTLIB_COLUMNS = [ 'GALID', 'RA_GAL', 'DEC_GAL', 'a_rot',
                    'n_Sersic', 'a_Sersic', 'b_Sersic', 'w1_Sersic', 'MAGNIFICATION' ] + \
                    BAND_LIST_HOSTLIB

SOC_COLUMNS     = [ 'id',    'ra',     'dec',     'pa',
                    'n',     'a',      'b',     'w1', 'MAGNIFICATION' ] + \
                    BAND_LIST_SOC


PREFIX_OUTFILE   = 'SNPIT_'
PREFIX_TRUTH     = 'TRUTH_'

LABEL_STAR      = "STAR"
LABEL_GALAXY    = "GALAXY"
LABEL_TRANSIENT = "TRANSIENT"
# - - - - -
# keys in config input file
KEY_POINTING_FILE      = "POINTING_FILE"

KEY_STAR_CATALOG       = "STAR_CATALOG"
KEY_GALAXY_CATALOG     = "GALAXY_CATALOG"
KEY_TRANSIENT_CATALOG  = "TRANSIENT_CATALOG"

KEY_SNANA_TRANSIENT_PATH = "SNANA_TRANSIENT_PATH"
KEY_MJD_TEMPLATE         = "MJD_TEMPLATE"
KEY_SNANA_HOSTLIB_FILE   = "SNANA_HOSTLIB_FILE"

STAR_CATNAME_GAIA = "GAIA"
STAR_CATNAME_SYN  = "SYN"

KEY_RA_CEN  = "RA_CEN"
KEY_DEC_CEN = "DEC_CEN"
KEY_RADIUS  = "RADIUS"
KEY_ROLL    = "ROLL"

# - - - - - -
STRING_L1 = "L1"
STRING_L2 = "L2"

CAL_DICT = { STRING_L1: 'uncal',   STRING_L2: 'cal' }

TYPE_SERSIC = "SER"   # for galaxies
TYPE_PSF    = "PSF"   # for point sources (transients and stars

# ============================
def setup_logging():
    logging.basicConfig(level=logging.INFO,
                        format="[%(levelname)6s %(message)s")     

    #logging.getLogger("matplotlib").setLevel(logging.ERROR)
    #logging.getLogger("seaborn").setLevel(logging.ERROR)
    return

# ===========
def get_args():
    parser_inputs = argparse.ArgumentParser()

    msg = "print config HELP"
    parser_inputs.add_argument("--HELP", "-H", help=msg, action="store_true")
    
    msg = "Required: name of sim config input file"
    parser_inputs.add_argument("input_config_file", help=msg, nargs="?", default=None)

    msg = "Required: image level (L1 or L2)"
    parser_inputs.add_argument("--level", "-l",
                               help=msg, default=None, type=str )
    
    msg = "Required: SCA number (1-18)"
    parser_inputs.add_argument("--scanum", "-s",
                               help=msg, default=None, type=int )

    msg = f"Required: BAND {' '.join(BAND_LIST_INPUT)}"
    parser_inputs.add_argument("--band_snana", "-b",
                               help=msg, default=None, type=str )    
    
    msg = "optional: MJD shift (w.r.t transient sim) for image sim "
    parser_inputs.add_argument("--mjd_shift", "-m", help=msg, default=0, type=float )

    msg = f"optional jobsplit by mjd; default = 0 1. Examples: 0 4 or 1 4 or 2 4 or 3 4"
    parser_inputs.add_argument("--jobsplit_mjd",
                               help=msg, default=None, type=int, nargs="+" )

    msg = f"optional jobsplit by texpose with Texpose -> Texpose/Nsplit"
    parser_inputs.add_argument("--jobsplit_texpose",
                               help=msg, default=None, type=int, nargs="+" )    
    
    msg = "snid to dump astropy table rows (and enable --check mode to skip doing sim)"
    parser_inputs.add_argument("--snid_dump", help=msg, default=None, type=int )

    msg = "galid to dump astropy table rows (and enable --check mode to skip doing sim)"
    parser_inputs.add_argument("--galid_dump", help=msg, default=None, type=int )

    msg = "mjd to dump entire astropy table (and enable --check mode to skip doing sim)"
    parser_inputs.add_argument("--mjd_dump", help=msg, default=0.0, type=float )
    
    msg = "name of start-dump star csv table file"
    parser_inputs.add_argument("--star_dump", help=msg, default=None, type=str )

    # - - - -
    msg = "Check loops over band/MJD.SCA but do not submit them"
    parser_inputs.add_argument("--check", "-c", help=msg, action="store_true")
    
    msg = "Run quick_sim_test with prescaled sources and no transients"
    parser_inputs.add_argument("--quick_sim_test", "-q", help=msg, action="store_true")

    
    args = parser_inputs.parse_args()

    # load band_soc corresponding to single-char band_snana
    args.band_soc = None
    for b_soc, b_snana in zip(BAND_LIST_SOC, BAND_LIST_INPUT):
        if b_snana == args.band_snana:
            args.band_soc = b_soc

    if args.snid_dump or args.galid_dump or args.mjd_dump>1.0 :
        args.check = True

    args.jobsplit = None
    if args.jobsplit_mjd:     args.jobsplit = args.jobsplit_mjd
    if args.jobsplit_texpose: args.jobsplit = args.jobsplit_texpose

    if args.HELP:
        print_config_HELP()
        
    return args
    # end get_args

def print_config_HELP():

    msg = f"""
    # =================================================================
    #      config file HELP for romanisim_snpit_wrapper.py
    # =================================================================
    
    # this block is used by sbatch_prep_romanisim+romancal.py to prepare slurm jobs
    # for romanisim_snpit_wrapper.py. 
    # Z*4 -> allocate x4 more cores for Z
    # F/4 -> divide each F-band exposure into 4 sub-exposures
    
    SBATCH_PREP:
      SCA_RANGE: 1 18       # range of SCA numbers to process
      BANDS:     R  Z*4 Y J H F/4  
      NODES:     mem-med   
      WALLTIME:  "4:00:00"
  
    # ---------------------------------------------------
    # keys for romanisim_snpit_wrapper.py

    # define MA number for exposure times
    MA_TEXPOSE_MAP_FILE:  $INPUT_CATALOG_DIR/MA_TEXPOSE_MAP.DAT

    # note FoV is about 0.8 x 0.4 ~ .281 deg^2
    SKY_REGION:
      RA_CEN:   242.5   # degrees
      DEC_CEN:   54.4
      RADIUS:     1.0  
      ROLL:      10.0  # degress

    # IMGNUM in SNANA sim is matched to pointing file to get RA,DEC,ROLL
    POINTING_FILE:  $INPUT_CATALOG_DIR/IMSIM_for_JPTEST.POINTING.gz

    # recommend reading star cat from pre-made csv file
    STAR_CATALOG:
    - $INPUT_CATALOG_DIR/STARS_GAIA_ROMAN_NORTH_radius1.csv.gz
    - $INPUT_CATALOG_DIR/STARS_SYN_NORTH_PS10.csv.gz
    #- GAIA   # catalog query is slow, so use this only with --star_dump <cat_file>
    #- SYN    # idem

    # galaxy catalog is the HOSTLIB used for SNANA sim
    GALAXY_CATALOG:
      SNANA_HOSTLIB_FILE:   $INPUT_CATALOG_DIR/IMSIM_for_JPTEST.HOSTLIB.gz

    # artificially reduce number of galaxies by requiring last GALID digits;
    #GALID_LAST_DIGITS:   174  81  # require last GALID digits to be 174 or 81

    # read transients from SNANA sim data folder
    TRANSIENT_CATALOG:
      SNANA_TRANSIENT_PATH:  $INPUT_CATALOG_DIR/snana_sim_sundial_pilot+core

    
    # select MJD ranges to process
    # MJD_RANGES:
    - 61890 61980   # CORE 
    - 61499 61531   # PILOT

    # NMJD_PROCESS: 1  # use this for quick/debug outputs

    """

    print(f"{msg}")

    sys.exit("\n Done.")
    return

def read_config(args):
    path_expand = os.path.expandvars(args.input_config_file)
    logging.info(f"Reading YAML input from {path_expand}")
    with open(path_expand) as f:
        return yaml.safe_load(f.read())

# ---------------------------
def print_banner(banner, level=1):
    logging.info('')

    if level == 1:
        logging.info('# ================================================================')
    elif level == 2:
        logging.info('# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
    elif level == 9:
        logging.info('# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ')
        logging.info('# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ')        
        
    logging.info(f"{banner}")
    return

def get_sky_region(config):
    SKY_REGION = config['SKY_REGION']
    ra     = SKY_REGION[KEY_RA_CEN]
    dec    = SKY_REGION[KEY_DEC_CEN]
    radius = SKY_REGION[KEY_RADIUS]  # Search radius in degrees

    if KEY_ROLL in SKY_REGION :
        roll = SKY_REGION[KEY_ROLL]
    else:
        roll = 0.0
    
    
    return ra, dec, radius, roll

def get_mjd_ranges(config):

    mjdmin_list = []
    mjdmax_list = []

    MJD_RANGES = config.setdefault('MJD_RANGES', ['0 999999'] )
    for mjd_range in config['MJD_RANGES']:
        mjd_list = mjd_range.split()
        mjdmin_list.append( float(mjd_list[0]) )
        mjdmax_list.append( float(mjd_list[1]) )

    return mjdmin_list, mjdmax_list

def get_star_cat_csv(csv_file):
    
    print_banner(f"Read star catalog from {csv_file}")
    star_cat = Table.read(csv_file, format='ascii.csv', comment='#' )

    n_star = len(star_cat)
    logging.info(f"Number of stars in csv_file: {n_star} ")
    logging.info(f"Catalog columns: {star_cat.colnames} ")
    
    return star_cat 

def get_star_cat_GAIA(args, config):
    
    print_banner("Read star catalog from GAIA")

    from   astroquery.gaia import Gaia
    from   romanisim import gaia, bandpass, catalog
    
    # get Gaia star cat; read from DB or from csv file (depends on size)
    ra, dec, radius, roll = get_sky_region(config)

    t0 = time.time()
    query = f'SELECT * FROM gaiadr3.gaia_source WHERE distance({ra}, {dec}, ra, dec) < {radius}'

    logging.info(f"Query: {query}")
    job    = Gaia.launch_job_async(query)
    result = job.get_results()
    t1 = time.time()
    t_query = t1 - t0 # seconds
    logging.info(f"CPU(GAIA-QUERY):  {t_query:.0f} seconds ")
    
    logging.info('Filter the Gaia results for stars and exclude bright stars')
    result = result[result['classprob_dsc_combmod_star'] >= 0.7]
    result = result[result['phot_g_mean_mag'] > 16.5]

    # Set the observation time
    obs_time = '2026-10-31T00:00:00'

    logging.info('Make the Roman I-Sim formatted catalog')
    gaia_catalog = gaia.gaia2romanisimcat(result, Time(obs_time),
                                          fluxfields=set(bandpass.galsim2roman_bandpass.values()) )

    logging.info('Reject anything with missing fluxes or positions')
    names = [f for f in gaia_catalog.dtype.names if f[0] == 'F']
    names += ['ra', 'dec']

    bad = np.zeros(len(gaia_catalog), dtype='bool')
    for b in names:
        bad = ~np.isfinite(gaia_catalog[b])
        if hasattr(gaia_catalog[b], 'mask'):
            bad |= gaia_catalog[b].mask
        gaia_catalog = gaia_catalog[~bad]


    # add id column defined as 100,000 * ra * dec
    new_id = np.abs(100000 * gaia_catalog['ra'] * gaia_catalog['dec']).astype(int) 
    gaia_catalog['id'] = new_id
    new_order = ['id'] + [name for name in gaia_catalog.colnames if name != 'id']
    gaia_catalog = gaia_catalog[new_order]


    n_star = len(gaia_catalog)
    logging.info(f"Number of Gaia stars in catalog: {n_star} ")
    logging.info(f"Gaia catalog columns: {gaia_catalog.colnames} ")

    print(f"{gaia_catalog[0:5]}" )
    sys.stdout.flush()

    
    if args.star_dump :
        csv_file = args.star_dump
        logging.info(f"\t write {csv_file}")
        gaia_catalog.write(csv_file, format='ascii.csv', overwrite=True)

    return gaia_catalog
    # end get_star_cat_GAIA


def get_star_cat_SYN(args, config):
    
    ra, dec, radius, roll = get_sky_region(config)
    band_soc        = args.band_soc

    nstar_per_degsq = 30000
    nstar_per_degsq *= 20 # xxx temp for LMC

    n_star    = int((3.14*radius*radius) * nstar_per_degsq) 
    faint_mag = 24.
    seed      = 45345

    print_banner(f"Make SYNthetic star catalog: n_star={n_star}   m < {faint_mag}")
    
    from romanisim  import gaia, catalog

    star_cat = catalog.make_stars(SkyCoord(ra, dec, unit='deg'), n_star,
                                  radius=radius, index=5/3., faintmag=faint_mag, 
                                  truncation_radius=None,
                                  bandpasses = BAND_LIST_SOC, rng=None, seed=seed)

    # add id column defined as 100,000 * ra * dec
    new_id = np.abs(100000 * star_cat['ra'] * star_cat['dec']).astype(int) 
    star_cat['id'] = new_id
    new_order = ['id'] + [name for name in star_cat.colnames if name != 'id']
    star_cat = star_cat[new_order]
    
    if args.star_dump:
        csv_file = args.star_dump
        logging.info(f"\t write {csv_file}")
        star_cat.write(csv_file, format='ascii.csv', overwrite=True)
    
    print(f"\n{star_cat[0:5]}")
    sys.stdout.flush()
    
    return star_cat

def get_galaxy_cat_HOSTLIB(args, config, hostlib_file, iser):

    # input iser = 0 or 1 = Sersic profile (e.g., disk & bulge)
    
    print_banner(f"Read galaxy catalog with Sersic component {iser}")

    hostlib_col  = HOSTLIB_COLUMNS
    ra, dec, radius, roll = get_sky_region(config)
    
    # replace a_Sersic and b_Sersic ..
    sersic_replace_list = [ 'n', 'a', 'b' ]
    for s in sersic_replace_list:
        #logging.info(f"\t replace {s} with {s}{iser} ")
        s_old = f"{s}_Sersic"
        s_new = f"{s}{iser}_Sersic"
        hostlib_col  = [ s_new if item == s_old else item for item in hostlib_col]

    logging.info(f"HOSTLIB file to read:    {hostlib_file}")
    logging.info(f"HOSTLIB columns to read: {hostlib_col} ")
    
    t_all = Table.read(hostlib_file, format='ascii', include_names = hostlib_col )
    ntot_gal = len(t_all)

    t_select = apply_radius_cut(config, t_all, 'RA_GAL', 'DEC_GAL')

    # - - - - - -
    ncut_gal = len(t_select)
    logging.info(f"Number of galaxies (all -> radiusCut): {ntot_gal} -> {ncut_gal}")
    logging.info(f"Rename columns for query:")
    
    for colname_h, colname_soc in zip(hostlib_col, SOC_COLUMNS):
        if colname_h in t_select.colnames:  # note that w0_Sersic or w1_sersic could be missing
            t_select.rename_column(colname_h, colname_soc)

    # for iser = 0, create w0_Sersic column = 1 - w1
    if iser == 0:
        t_select['w0'] = 1.0 - t_select['w1']
        t_select['w0'] = [float(f"{x:.4f}") for x in t_select['w0']] # avoid spurious extra digits

    var_w = f"w{iser}"  # w0 or w1, used to select fraction of flux for this Sersic component
    logging.info(f"Convert ABmag columns to maggie flux = 10**(-0.4*ABmag) * SersicFraction")
    for mag in BAND_LIST_SOC:  # here is a mag
        if not reject_band(args,mag):
            logging.info(f"\t convert {mag} ABmag to maggie flux")
            t_select[mag] = t_select[var_w] * 10.0**(-0.4*t_select[mag])  # convert to maggie

    
    logging.info(f"Add column for half_light_radius = a and  ba = b/a:")
    t_select['half_light_radius'] = t_select['a']
    t_select['ba'] = t_select['b'] / t_select['a']
    t_select['ba'] = [float(f"{x:.4f}") for x in t_select['ba']] # avoid spurious extra digits
        
    logging.info(f"Add column for type = {TYPE_SERSIC} ")
    t_select['type'] = TYPE_SERSIC

    if ROTATE_GALAXIES_90DEG:  # Apr 25 2026
        logging.info("Rotate galaxy pointing angle (pa) by 90 degrees")  
        t_select['pa'] += 90.0
        mask = t_select['pa'] > 180   
        t_select['pa'][mask] -= 180   # avoid angle > 180

    logging.info(f"Table columns for query: {t_select.colnames}")

    # check pre-scale option using GALID pattern match
    t_final = pattern_match_galid(config,t_select)
    
    print(f"{t_final[0:5]}")
    sys.stdout.flush()
    
    #sys.exit(f"\n xxx DEBUG EXIT xxx")
    
    return t_final  # end get_galaxy_cat_HOSTLIB

def pattern_match_galid(config,t_in):

    t_out = t_in
    #KEY = 'PATTERN_MATCH_GALID'
    KEY = 'GALID_LAST_DIGITS'
    if KEY not in config: return t_out

    # add GALID column that is string (id is int)
    t_in['GALID'] = t_in['id'].astype(str)
    
    t_list = []
    pattern_list = config[KEY].split()
    for pattern in pattern_list:
        lenp  = len(pattern)
        mask  = [s[-lenp:] == pattern for s in t_in['GALID']]
        t_list.append(t_in[mask])

    t_out = vstack(t_list)
    t_out.remove_column('GALID')
            
    n_in  = len(t_in)
    n_out = len(t_out)
    logging.info(f" Require GALID end with {pattern_list}: n_row = {n_in} --> {n_out}")

    #sys.exit(f"\n xxx DEBUG STOP xxx ")
    return t_out

def apply_mjd_cut(config, t, COLNAME_MJD ):

    mjdmin_list, mjdmax_list = get_mjd_ranges(config)
    mask_sum = None
    for mjdmin, mjdmax in zip(mjdmin_list,mjdmax_list):
        mask0 = t[COLNAME_MJD] >= (mjdmin-0.0001)
        mask1 = t[COLNAME_MJD] <= (mjdmax+0.0001)
        mask  = mask0 & mask1
        #sys.exit(f"\n xxx mask = {mask[0:20]}")
        if mask_sum is None:
            mask_sum = mask
        else:
            mask_sum = mask_sum | mask


    t_cut = t[mask_sum]
    return t_cut

def apply_radius_cut(config, t, COLNAME_RA, COLNAME_DEC):

    # return table with radius cut around central RA,DEC
    
    # apply radius cut 
    ra, dec, radius, roll = get_sky_region(config)
    
    t_coords = t[ COLNAME_RA, COLNAME_DEC ]
    center = SkyCoord(ra= ra*u.deg, dec = dec*u.deg, frame='icrs')

    #  Create a SkyCoord object for the events
    event_coords = SkyCoord(ra=t_coords[COLNAME_RA]*u.deg,
                            dec=t_coords[COLNAME_DEC]*u.deg, frame='icrs')

    #  Calculate the separation (on-sky angular distance)
    sep = center.separation(event_coords)

    # 5. Apply the cut
    mask = sep < radius * u.deg
    t_select = t[mask]

    return t_select

def prefix_asdf(mjd, sca, band, level ):
    # compute visit id = int(mjd * 1E4)
    visit_id = int(mjd * 1.0E4 + 0.5)
    prefix   = f"{PREFIX_OUTFILE}VISIT{visit_id}_WFI{sca:02d}_{band}_{level}"
    return prefix

def run_sim(args, full_cat, info_dict ):

    sca           = info_dict['sca']
    band_soc      = info_dict['band_soc']
    mjd_transient = info_dict['mjd']
    ma_num        = info_dict['ma_num']
    ra            = info_dict['ra']
    dec           = info_dict['dec']
    roll          = info_dict['roll']

    # run the sim !!!

    mjd_image = mjd_transient + args.mjd_shift
    t_mjd     = Time(mjd_image, format='mjd')
    obs_date  = t_mjd.isot
    seed      = 7  # Galsim random number generator seed for reproducibility
    
    str_level = args.level     
    level     = int(str_level[-1])
    cal_level = CAL_DICT[str_level]
    
    prefix       = prefix_asdf(mjd_image, sca, band_soc, str_level )
    prefix_truth = prefix.replace(PREFIX_OUTFILE, PREFIX_TRUTH)
    
    fnam_asdf  = f'{prefix}.asdf'  
    fnam_png   = f'{prefix}.png'
    fnam_log   = f'{prefix}.log'  
    fnam_truth = f'{prefix_truth}.dat.gz'
    
    # skip if file exists to enable auto-contunuation by re-submitting
    if os.path.exists(fnam_asdf):
        print_banner(f"Skip romanIsim for already existing {fnam_asdf}", level = 2);
        return fnam_asdf
    
    print_banner(f"Run the sim for {fnam_asdf}:", level=2)
    logging.info(f"\t mjd={mjd_image} -> {obs_date}")
    
    if args.check :
        logging.info(f"\t SKIP romanIsim for --check mode.")
        return fnam_asdf
    
    # xxx mark delete ra, dec, radius, roll = get_sky_region(config)
    
    logging.info('\t parse args_sim')
    # Set other arguments for use in Roman I-Sim. The code expects a specific format for these,
    # so this is a little complicated looking.
    parser = argparse.ArgumentParser()

    parser.set_defaults(usecrds=True, psftype="stpsf", level=level, filename=fnam_asdf,
                        drop_extra_dq=True, sca=sca, 
                        bandpass=band_soc, pretend_spectral=None)

    args_sim = parser.parse_args([])

    logging.info('\t Set reference files to None for CRDS')
    for k in parameters.reference_data:
        parameters.reference_data[k] = None

    logging.info('\t Set Galsim RNG object')
    rng = galsim.UniformDeviate(seed)

    logging.info(f'\t Set default persistance information')
    persist = persistence.Persistence()

    logging.info(f'\t Set metadata')
    metadata = ris.set_metadata(date=obs_date, bandpass=band_soc, sca=sca,
                                ma_table_number=ma_num, usecrds=True)
    
    logging.info(f'\t Update the WCS info ')
    wcs.fill_in_parameters(metadata, SkyCoord(ra, dec, unit='deg', frame='icrs'),
                           boresight=False, pa_aper=roll )

    # load goofy obserbation_id
    obs_id = int(1.0e10*mjd_image) + abs(int(ra*dec)) + int(roll)
    metadata['observation'] = {}
    metadata['observation']['observation_id'] = str(obs_id)

    # get list of ids overlaid (code from Thushara, waiting to be ingested by NASA)
    write_truth_table(fnam_truth, args, mjd_image, metadata, full_cat)

    if WRITE_TRUTH_TABLE_ONLY:
        logging.info(f"\t SKIP romanIsim for WRITE_TRUTH_TABLE_ONLY")
        return fnam_asdf
    
    # - - - - - - -  -
    #sys.exit(f"\n xxx metadata = \n{metadata}")
    
    t_start = time.time()
    logging.info(f'\t run the sim')
    #import pdb; pdb.set_trace() # DEBUGGER
    
    ris.simulate_image_file(args_sim, metadata, full_cat, rng, persist)
        
    t_end = time.time()

    t_sim = t_end - t_start # seconds
    logging.info(f"CPU({prefix}):  {t_sim:.0f} seconds ")

    
    logging.info(f" Create figure {fnam_png}")
    imagelib.mkfigure(fnam_asdf, plotname = fnam_png)  # crashes; no mkfigure attribute?
    
    return fnam_asdf


def write_truth_table(fnam_truth, args, mjd, metadata, full_cat):
    
    # Created Aug 4 2026
    # cut-and-paste code from Thushara
    #
    # MAYBE TO-DO: sum galaxy sersic components to get total mag

    logging.info(f"Write TRUTH table to {fnam_truth}")
    
    band_soc = args.band_soc  # this is table column flux in maggies
    scanum   = args.scanum
    
    imwcs   = wcs.get_wcs(metadata, usecrds=True)
    image   = galsim.ImageF(roman.n_pix, roman.n_pix, wcs=imwcs, xmin=0, ymin=0)
    trimcat = imcode.trim_objlist(full_cat, image)
    xpos, ypos = image.wcs._xy(np.radians(trimcat['ra']), np.radians(trimcat['dec']))
    ignore_distant_sources = 10
    t_keep         = imcode.in_bounds(xpos, ypos, image.bounds, ignore_distant_sources)
    t_kept_catalog = trimcat[t_keep]
    
    # fetch and store x,y detector locations of "KEPT" objects
    x_kept, y_kept = image.wcs._xy(np.radians(t_kept_catalog['ra']), np.radians(t_kept_catalog['dec']))
    t_kept_catalog['x_det'] = x_kept
    t_kept_catalog['y_det'] = y_kept
    
    #nx = len(x_kept);  ny=len(y_kept)
    #print(f"\n xxx nx/ny = {nx}/{ny}  \n xxx x_kept = {x_kept} \n xxx y_kept = {y_kept} \n")
    # - - - -
    
    # select columns to write out to csv file
    col_noband_list = ['id', 'ra', 'dec', 'x_det', 'y_det', 'n', 'label' ]
    col_save_list   = col_noband_list + [ band_soc] 
    t_kept_cols     = t_kept_catalog[ col_save_list ]

    mask_galaxy    = t_kept_cols['label'] == LABEL_GALAXY
    mask_star      = t_kept_cols['label'] == LABEL_STAR
    mask_transient = t_kept_cols['label'] == LABEL_TRANSIENT

    t_kept_gal     = t_kept_cols[mask_galaxy]
    if len(t_kept_gal) > 0:
        grouped        = t_kept_gal.group_by('id')
        t_kept_gal_sum = grouped['id', band_soc].groups.aggregate(np.sum)
        t_other_cols   = unique(t_kept_gal[col_noband_list], keys='id')
        t_kept_gal_merged = join(t_kept_gal_sum, t_other_cols, keys='id')
        t_kept_gal_merged['n'] = -1  # sersic index has no meaning for flux sum over profiles
        n_gal         = len(t_kept_gal_merged)         # count 1 per galaxy id
        n_gal_sersic  = len(t_kept_cols[mask_galaxy])  # count both sersic profiles
    else:
        n_gal = 0; n_gal_sersic = 0
        
    n_star        = len(t_kept_cols[mask_star])
    n_tra         = len(t_kept_cols[mask_transient])
    
    logging.info(f"\t TRUTH Nobj(transient, star)       = {n_tra} {n_star} ")
    logging.info(f"\t TRUTH Nobj(galaxy_sersic, galaxy) = {n_gal_sersic}  {n_gal}")
    
    # for original SERSIC galaxy rows, rename GALAXY -> GALAXY_SERSIC;
    # GALAXY will be used for flux-sum over sersic profiles.
    t_kept_cols['label']              = t_kept_cols['label'].astype("<U20")

    if n_gal > 0:
        t_kept_cols['label'][mask_galaxy] = LABEL_GALAXY + "_SERSIC"
        t_kept_final = vstack( [ t_kept_cols, t_kept_gal_merged ] )
    else:
        t_kept_final = t_kept_cols

    #sys.exit(f"\n xxx t_kept_gal_sum = \n{t_kept_gal_sum} \n xxx other_cols = \n{t_other_cols}")
    
    # define columns to add
    colname_mag    = 'mag'
    colname_maggie = 'flux_maggie'
    colname_mjd    = 'mjd'
    colname_band   = 'band'
    colname_sca    = 'sca'

    t_kept_final[colname_maggie]  = t_kept_final[band_soc]
    t_kept_final[colname_mag]     = -2.5*np.log10(t_kept_final[colname_maggie])  # restore mag
    t_kept_final[colname_mjd]     = mjd
    t_kept_final[colname_band]    = band_soc
    t_kept_final[colname_sca]     = scanum
    del t_kept_final[band_soc]  # remove redundant column that is now called maggie

    # set formats
    t_kept_final[colname_maggie].format = '%.4e'
    t_kept_final[colname_mag].format    = '%.4f'
    t_kept_final['ra'].format           = '%.6f'
    t_kept_final['dec'].format          = '%.6f'
    t_kept_final['x_det'].format        = '%08.3f'
    t_kept_final['y_det'].format        = '%08.3f'        

    #logging.info(f"Write truth table to {fnam_truth}")
    #t_kept_final.write(fnam_truth, format='ascii.csv', delimiter=',', overwrite=True)

    with gzip.open(fnam_truth, 'wt', encoding='utf-8') as f:
        t_kept_final.write(f, format='ascii.csv')
        
    return
    
def get_transient_cat_SNANA(transient_path, genversion, args, config):

    logging.info('')
    logging.info(f" @@@@@ Prepare catalog for {genversion} @@@@@ ")
    
    # return full transient catalog with all models and all MJDs
    sim_path   = f"{transient_path}/{genversion}"

    head_list  = sorted(glob.glob(f"{sim_path}/*_HEAD.FITS*"))
    phot_list  = sorted(glob.glob(f"{sim_path}/*_PHOT.FITS*"))

    t_dict = {}
    n_file = 0
    for head_file, phot_file in zip(head_list,phot_list):
        n_file += 1
        cat_name = genversion + f"-{n_file:02d}"
        t_dict[cat_name] = \
            get_transient_cat_SNANA_file(genversion, head_file, phot_file, args, config)

    logging.info('')
    logging.info(f"  Stack {n_file} sub-catalogs:")
    t_cat_all = cat_vstack(t_dict)

    print_col_list = ['MJD', 'ra', 'dec', args.band_soc, 'TEXPOSE']
    #print_col_list += [ 'IMGNUM' ]
    #print(f"\n{t_cat_all[print_col_list][0:5]} \n")

    print(f"\n{t_cat_all[0:5]} \n")
    sys.stdout.flush()

    return t_cat_all

def get_transient_cat_SNANA_file(genversion, head_file, phot_file, args, config):

    base_name = os.path.basename(head_file).split('_HEAD')[0]
    logging.info(f"  Process {genversion}/{base_name}")
    
    # read HEAD file
    head_col_list = [ 'SNID', 'RA', 'DEC', 'NOBS', 'PTROBS_MIN', 'PTROBS_MAX' ]
    t_head_all = Table.read(head_file, format='fits' )
    nrow_head_orig = len(t_head_all)

    
    t_head_all = apply_radius_cut(config, t_head_all, 'RA', 'DEC')
    nrow_head_cut = len(t_head_all)
    logging.info(f"\t Number of transients (all -> RadiusCut) = {nrow_head_orig} -> {nrow_head_cut}")
    del t_head_all.meta['EXTNAME']  # remove junk colummn to avoid inner join confusion later
    
    t_head         = t_head_all[head_col_list]
    t_head['SNID'] = t_head['SNID'].astype(int)
    t_coords       = t_head[ 'SNID', 'RA', 'DEC']
    
    # read PHOT file
    phot_col_list  = [ 'MJD', 'BAND',  'SIM_MAGOBS', 'TEXPOSE' ]
    if config['t_pointing'] :
        phot_col_list += [ 'IMGNUM', 'DETNUM' ]  # IMGNUM in SNANA is POINTING id

    t_phot_all     = Table.read(phot_file, format='fits' )
    nrow_phot_orig = len(t_phot_all)
        
    del t_phot_all.meta['EXTNAME']
    t_phot     = t_phot_all[phot_col_list]
    # do NOT remove pad rows, otherwise PTROBS_MIN/MAX won't work
    
    del t_head_all, t_phot_all

    # construct SNID columm for phot table
    nrow_head = len(t_head)
    nrow_phot = len(t_phot)
    logging.info(f"\t nrow[head,phot] = {nrow_head} , {nrow_phot}")
    snid_phot_list = np.zeros(nrow_phot, dtype='int')
    for snid, ptrobs_min, ptrobs_max in \
            zip( t_head['SNID'], t_head['PTROBS_MIN'], t_head['PTROBS_MAX'] ):
        snid_phot_list[ptrobs_min-1:ptrobs_max] = int(snid)

    
    t_phot['SNID'] = snid_phot_list  
    t_phot.sort([ 'SNID', 'MJD', 'BAND' ])  # useful for visual debug

    # apply MJD cut after glueing SNID 
    t_phot        = apply_mjd_cut(config, t_phot, 'MJD' )
    nrow_phot_cut = len(t_phot)
    logging.info(f"\t Number of MJD + PAD rows (all -> MJDcut) = " \
                 f"{nrow_phot_orig} -> {nrow_phot_cut}" )

    if nrow_phot_cut == 0:
        sys.exit(f"ERROR: no photometry rows selected for transients.")
        
    # now remove pad rows with MJD = -777
    mask_to_keep = t_phot['MJD'] > 0.0
    t_phot       = t_phot[mask_to_keep]
    t_phot['BAND'] = [x.strip() for x in t_phot['BAND']]  # remove pad spaces
    
    # join HEAD and PHOT tables to append ra,dec to PHOT table
    logging.info(f"\t join coords to LCphot table ... ")
    t_phot_merge = join(t_phot, t_coords, keys = 'SNID', join_type='inner')
    t_phot_merge.sort([ 'SNID', 'MJD', 'BAND' ]) # for visual debug 
    nrow_merge = len(t_phot_merge )

    names_orig = ['RA', 'DEC']
    names_soc  = ['ra', 'dec']
    for old, new in zip(names_orig, names_soc):
        t_phot_merge.rename_column(old, new)

    print_debug = False
    if print_debug:
        print(f"\n xxx head  table = \n{t_head[:5]}")    
        print(f"\n xxx phot  table = \n{t_phot[495:505]}")
        print(f"\n xxx merge table = \n{t_phot_merge[495:505]}")
        sys.stdout.flush()
        
    # - - - - -
    # prepare t_cat for each band, then add them
    t_cat = {}

    logging.info(f"\t convert SIM_MAGOBS to maggie flux ... ")
    for band_snana, band_soc in zip(BAND_LIST_SNANASIM, BAND_LIST_SOC):
        mask_band   = t_phot_merge['BAND'] == band_snana
        t_phot_band = t_phot_merge[mask_band] #  rows with selected band_snana

        if len(t_phot_band) > 0 and not reject_band(args,band_soc): 
            t_phot_band[band_soc] = 10.0**(-0.4*t_phot_band['SIM_MAGOBS'])
            t_phot_band.remove_column('BAND')
            t_phot_band.remove_column('SIM_MAGOBS')
            #cat_name = base_name + '-' + band_soc
            cat_name = genversion + '-' + band_soc
            t_cat[cat_name] = t_phot_band

    if len(t_cat) == 0:
        sys.exit(f"\n ERROR: no {args.band_snana} photometry to stack in get_transient_cat_SNANA_file")
        
    t_cat_all = cat_vstack(t_cat) # stack all bands
    t_cat_all['type'] = TYPE_PSF

    missing_col_list =  [ 'n', 'half_light_radius', 'ba' ] 
    for col in missing_col_list:
        t_cat_all[col] = -1

    t_cat_all.rename_column("SNID", 'id')


    t_cat_all = split_texpose(t_cat_all, args)
    #sys.exit(f"\n xxx t_cat_all = \n{t_cat_all[0:20]}" )
                        
    return t_cat_all
    # end get_transient_cat_SNANA

def get_mjd_template_cat(mjd_template_list, texpose, args, config):

    # for each MJD template, put zero-flux source in middle of FP

    logging.info("")
    logging.info(f"  Prepare MJD_TEMPLATES for {mjd_template_list}")
    
    band_soc              = args.band_soc
    ra, dec, radius, roll = get_sky_region(config)
    
    n_mjd     = len(mjd_template_list)
    ra_list   = [ ra  ] * n_mjd
    dec_list  = [ dec ] * n_mjd
    flux_list = [ 1.0E-30 ] * n_mjd

    columns_list = [ mjd_template_list, ra_list, dec_list, flux_list ]

    cat = Table(columns_list, names=('MJD', 'ra', 'dec', band_soc))
    cat['TEXPOSE']             = texpose
    cat['type']                = TYPE_PSF
    cat['id']                  =  1
    cat['n']                   = -1
    cat['half_light_radius']   = -1
    cat['ba']                  = -1

    cat = split_texpose(cat, args)
    
    print(f"{cat}")
    sys.stdout.flush()
    
    return cat
# end get_mjd_template_cat

def split_texpose(t_cat, args):

    # if we are splitting by Texpose, then each obs is expanded into
    # n_split obs, each with Texpose/n_split exposure time.
    # This is designed for F814 band that is expected to be split
    # into 4 obs per visit.
    #
    # Input t_cat is astropy table with SNANA transients
    #
    # Beware: this algorithm is slow; takes 30 seconds to sploit 15,000 rows to 60,000 rows
    
    if not args.jobsplit_texpose: return t_cat

    t0 = time.time()
    
    n_split = args.jobsplit[1]  # Divide Texpose -> Texpose / n_splot

    logging.info(f"")
    logging.info(f"\t Split TEXPOSE -> TEXPOSE / {n_split}")
    
    # manually construct rows for new table
    nrow_orig  = len(t_cat)
    colnames   = t_cat.colnames
    split_table_list = []
    new_mjd_list     = np.array([])
    new_texpose_list = np.array([])
    
    for row in t_cat:
        texpose_orig  = row['TEXPOSE']   # note mjd_orig is centered on TEXPOSE
        #texpose_orig  = 1000.0   # xxx REMOVE
        texpose_split = texpose_orig/ n_split

        mjd_off0_sec      = -texpose_orig/2.0 + texpose_split/2.0
        mjd_off1_sec      = +texpose_orig/2.0 + texpose_split/2.0
        mjd_offsets_sec   = np.arange(mjd_off0_sec, mjd_off1_sec, texpose_split)
        
        split_mjd_list      = row['MJD'] + mjd_offsets_sec/86400.0
        split_texpose_list  = np.array( [texpose_split]*n_split )
        
        #sys.exit(f"\n xxx split_mjd_list = {split_mjd_list} ")
        new_mjd_list     = np.append(new_mjd_list,     split_mjd_list)
        new_texpose_list = np.append(new_texpose_list, split_texpose_list)
        
    # - - - -
    #sys.exit(f"\n xxx new_mjd_list = \n{new_mjd_list[0:12]} \n" \
    #         f" xxx new_texpose_list = {new_texpose_list[0:12]} \n")

    # vstack n_split duplicate tables
    repeated_indices = np.repeat(np.arange(nrow_orig), n_split)
    t_cat_new          = t_cat[repeated_indices]

    # replace MJD and TEXPOSE columns of new table
    t_cat_new['MJD']     = new_mjd_list
    t_cat_new['TEXPOSE'] = new_texpose_list
    
    #print(f"\n xxx split_texpose: tcat_new = \n{t_cat_new[0:12]}")
    
    nrow_split = len(t_cat_new) 
    ratio      = nrow_split / nrow_orig
    if ratio != n_split :
        sys.exit(f"\n ERROR: nrow_split / nrow_orig = {ratio}, but n_split = {n_split}")

    t_split = time.time() - t0
    logging.info(f"\t Done splitting {nrow_orig} rows into {nrow_split} rows.")
    logging.info(f"\t CPU(split_texpose):  {t_split:.0f} seconds ")

    return t_cat_new

# end split_texpose

def cat_vstack(cat_dict, mjd=None):

    # stack all catalogs in cat_dict.
    # If mjd is not None, pick out these rows if MJD column exists (for transients only)

    cat_list = []
    for cat_name, cat in cat_dict.items():

        if mjd is not None and 'MJD' in cat.colnames:
            cat  = cat[ cat['MJD'] == mjd ] 
            is_transient = True
        else:
            is_transient = False   # star or galaxy            

        n_row = len(cat)
        if n_row > 0:
            cat_list.append(cat)
            logging.info(f"\t add '{cat_name}' to vstack with {n_row} rows")
            if is_transient :
                n = min(n_row,3)
                id_list = cat['id'][:n].tolist()
                logging.info(f"\t\t (first few transient IDs: {id_list})")

    full_cat = vstack( cat_list )
    return full_cat


def get_unique_mjd_by_band(args, config, cat_dict, tr_list):

    # return dictionary where each key is band, and item is sorted list of unique MJDs
    # e.g, unique_mjd_dict['F062'] = [ 54000, 54010, 54020 ...]
    # tr_list is list of transient dictionary keys so that stars and galaxies are ignored.

    print_banner(f"Find unique MJDs for each band for {tr_list} ", level = 2 )

    nmjd_process = config.setdefault('NMJD_PROCESS',9999)

    if nmjd_process < 9999:
        logging.info(f"Truncate image processing : NMJD_PROCESS = {nmjd_process}")
    
    unique_mjd_dict = {}
    for b in BAND_LIST_SOC:

        if reject_band(args,b): continue
        
        unique_mjd_list    = []
        n_mjd_tot = 0
        for tr in tr_list:  # loop over transient models and survey components
            t_temp   = cat_dict[tr]
            if b in t_temp.colnames:
                #t_temp   = t_temp[ ~t_temp[b].mask ]  # require valid flux                               
                mjd_list = t_temp['MJD'].tolist()[0:nmjd_process]
                n_mjd_tot += len(mjd_list)
                unique_mjd_list += sorted(list(set(mjd_list)))

        # get final unique list among all transient models for this band
        unique_mjd_list  = sorted(list(set(unique_mjd_list)))
        n_mjd_unique = len(unique_mjd_list)
        logging.info(f"  {b}: found {n_mjd_unique:4d} unique MJDs among {n_mjd_tot:8d} observations.")
        logging.info(f"\t unique_mjd_list: {unique_mjd_list}")

        # Jul 13 2026: check job-split option
        if args.jobsplit:
            ijob = args.jobsplit[0]
            njob = args.jobsplit[1]
            unique_mjd_list = [val for i, val in enumerate(unique_mjd_list) if i % njob == ijob]
            
            logging.info(f"\t APPLY SPLIT: ijob = {ijob} of {njob}")
            logging.info(f"\t AFTER SPLIT unique_mjd_list: {unique_mjd_list} ")

        unique_mjd_dict[b] = unique_mjd_list            
        
    #sys.exit(f"\n xxx DEBUG DUMP from unique")
    return unique_mjd_dict
# end get_unique_mjd_by_band

def get_unique_mjd_pointing(args, config):

    # return list of MJDs in pointing file.
    # This util is called if there are no transients

    t_pointing = config.setdefault('t_pointing',None)

    if not t_pointing:
        sys.exit(f"\n ERROR: with no transients, must have {KEY_POINTING_FILE} in config.")

    unique_mjd_dict = {}
    unique_texpose_dict = {}
    for b_soc, b_snana in zip(BAND_LIST_SOC, BAND_LIST_INPUT):
        if b_snana == args.band_snana:        
            t_match      = t_pointing[t_pointing['FILTER'] == b_snana ]
            mjd_list     = t_match['MJD'].tolist()
            texpose_list = t_match['TEXPOSE'].tolist()
            
            unique_mjd_dict[b_soc]     = mjd_list
            unique_texpose_dict[b_soc] = texpose_list
    
    return unique_mjd_dict, unique_texpose_dict

def expand_cat_by_mjd(args, config, cat_dict):

    # since there are no transients to define MJD and IMGNUM,
    # use pointing file to define these quantities. For each
    # row in cat_dict astropy tables, expand into row for each MJD
    # in pointing file.

    t_pointing = config.setdefault('t_pointing',None)
    if not t_pointing: return cat_dict
    
    t_match      = t_pointing[t_pointing['FILTER'] == args.band_snana ]
    mjd_list     = t_match['MJD'].tolist()
    texpose_list = t_match['TEXPOSE'].tolist()
    imgnum_list  = t_match['IMGNUM'].tolist()
    n_mjd = len(mjd_list)

    cat_expand_dict = {}
    for key, cat_orig in cat_dict.items():
        key_base = os.path.basename(key)
        logging.info(f" no transients -> expand {key_base} cat with {n_mjd} MJDs")
        cat_expand = cat_orig[np.repeat(np.arange(len(cat_orig)), n_mjd)]
        cat_expand["MJD"]     = np.tile(mjd_list,     len(cat_orig))
        cat_expand["TEXPOSE"] = np.tile(texpose_list, len(cat_orig))
        cat_expand["IMGNUM"]  = np.tile(imgnum_list,  len(cat_orig))
        cat_expand["DETNUM"]  = -9
        cat_expand_dict[key] = cat_expand
        #sys.exit(f"\n xxx cat_orig = \n{cat_dict}\n\n xxx cat_expand = \n{cat_expand} ")
    
    return cat_expand_dict

def print_asdf_expect(args, unique_mjd_dict):

    logging.info(f"")
    bands = list(unique_mjd_dict.keys())
    for band in bands:
        # Process only one band input by user, but maybe later allow few bands
        if reject_band(args, band) : continue
        n_expect = len(unique_mjd_dict[band])
        logging.info(f"N_EXPECT_ASDF: {n_expect}")
        for mjd in unique_mjd_dict[band]:
            prefix = prefix_asdf(mjd, args.scanum, band, args.level)
            logging.info(f"EXPECT_ASDF:  {prefix}.asdf")

    logging.info(f"")            
    return
    
def reject_band(args,b):
    if b == args.band_soc:  return False
    return True

def parse_texpoe_ma_map(args, config):
    
    ma_dict = {}
    KEY = 'MA_TEXPOSE_MAP_FILE'
    if KEY not in config:
        sys.exit(f"\n FATAL ERROR: missing {KEY} key in {args.input_config_file}")

    ma_file = os.path.expandvars(config[KEY])
    print_banner(f"Read map of T_expose vs. MA number from {ma_file}")
    
    with open(ma_file) as f:
        ma_dict = yaml.safe_load(f.read())['MA_DICT']

    # - - -
    #ma = get_ma_number(87.0,  ma_dict, True)
    #ma = get_ma_number(433.0, ma_dict, True)
    #ma = get_ma_number(748, ma_dict, True)
    #sys.exit(f"\n xxx DEBUG STOP \n xxx ma_dict = \n{ma_dict}")
    
    return ma_dict


def parse_pointing_file(args, config):

    df = None
    pointing_file = config.setdefault(KEY_POINTING_FILE,None)
    if pointing_file is None: return None

    pointing_file = os.path.expandvars(pointing_file)
    t_pointing    = Table.read(pointing_file, format='ascii.csv', comment='#' )

    n_row = len(t_pointing)
    print_banner(f"Read {n_row} pointings from {pointing_file}" )

    return t_pointing

def get_pointing_info(args, config, cat, check_pointing_dict, vbose):

    # return ra, dec, roll
    t_pointing = config['t_pointing']

    # if IMGNUM and DETNUM are in catalog, this means SNANA transients
    # are include, so do check that SNANA and POINTING file agree
    DO_CHECK_SNANA = 'IMGNUM' in cat.colnames  and 'DETNUM' in cat.colnames

    if args.jobsplit_texpose:
        # visit is split into multiple exposure, so cannot do accurate MJD check
        TOL_MJD = 0.2  # days
    else:
        # prepare for precise MJD tolerance check
        TOL_MJD = 0.0001  # days
    
    if t_pointing :  # from pointing file
        imgnum_list = cat['IMGNUM'].tolist()
        detnum_list = cat['DETNUM'].tolist()
        imgnum_snana  = next((x for x in imgnum_list if isinstance(x,int)), None)
        detnum_snana  = next((x for x in detnum_list if isinstance(x,int)), None)

        t_match       = t_pointing[t_pointing['IMGNUM'] == imgnum_snana ]
        ra            = t_match['RA_WFI_CENTER'].value[0]
        dec           = t_match['DEC_WFI_CENTER'].value[0]
        roll          = t_match['ROLL'].value[0]
        band_ptg      = t_match['FILTER'].value[0]
        mjd_ptg       = t_match['MJD'].value[0]
        
        band_check      = check_pointing_dict['band']
        mjd_check       = check_pointing_dict['mjd']
        pass_check_band = (band_check == band_ptg)
        pass_check_mjd  = abs(mjd_check-mjd_ptg) < TOL_MJD
        
        string_list     = [ 'BAND',     'MJD'  ]
        val_found_list  = [ band_ptg,   mjd_ptg ]
        val_check_list  = [ band_check, mjd_check]
        pass_check_list = [ pass_check_band, pass_check_mjd ]
        if vbose:
            logging.info(f"\t input IMGNUM={imgnum_snana} -> ra / dec / roll = {ra} / {dec} / {roll}")

        nerr = 0
        for val_found, val_check, string, pass_check in \
            zip(val_found_list,val_check_list, string_list, pass_check_list ):
            if not pass_check :
                val_dif = val_found-val_check
                logging.info(f"ERROR for {string}: expect {val_check:.4f} but found {val_found:.4f} in POINTING file\n\t diff(MJD) = {val_dif:4f} ")
                nerr += 1

        if nerr == 0 :
            logging.info("\t Consistent SNANA sim and POINTING file for MJD & BAND")
        else:
            sys.exit(f"\n ABORT on value checks between SNANA and POINTING file")
                     
    else:
        ra, dec, radius, roll = get_sky_region(config)  # from config file
    
    return ra, dec, roll


def get_texpose_transient(cat):
    texpose_list = cat['TEXPOSE'].tolist()
    texpose      = next((x for x in texpose_list if isinstance(x, float)), None)    
    return texpose

def get_ma_number(texpose, ma_dict, vbose):

    texpose_list = list(ma_dict.keys())
    i_texpose = int(texpose)
    closest_texpose = min(texpose_list, key=lambda x: abs(x - i_texpose))
    ma_num = ma_dict[closest_texpose]

    if vbose:
        logging.info(f"\t input Texpose={texpose} -> " \
                     f"map Texpose={closest_texpose} and ma_num={ma_num}")
        
    return ma_num
        
def quick_sim_test(cat_dict, config):

    print_banner(f"QUICK SIM TEST WITH GAIA AND GALAXIES ONLY (NO TRANSIENTS)")
    
    full_cat  = cat_vstack(cat_dict)

    n_subset = 100  # keep this many rows
    random_indices = np.random.choice(len(full_cat), n_subset, replace=False)
    subset_cat     = full_cat[random_indices]

    mjd       = 54685
    sca       = 1          # 1-18
    ma_num    = 1002
    ra, dec, radius, roll = get_sky_region(config)    
    for band in BAND_LIST_SOC:

        info_dict = { 'sca':sca,  'band_soc':band_soc, 'mjd':mjd, 'ma_num':ma_num, 
                      'ra':ra, 'dec':dec, 'roll':roll   }
        
        dummy_file = run_sim(args, subset_cat, info_dict )   
        sys.exit('\n xxx DEBUG STOP after processing {band} xxx\n')

    return

def do_galid_dump(args,cat_dict):

    galid_dump     = args.galid_dump
    full_cat       = cat_vstack(cat_dict)
    galid_dump_cat =   full_cat[ full_cat['id'] == galid_dump ]

    logging.info(f" Astropy Table dump for galid = {galid_dump}:")

    col_list = [ 'id', 'ra', 'dec', 'a', 'b', 'n', 'w0', 'w1', 'half_light_radius', 'pa', 'ba',  args.band_soc ]
    print(f"{galid_dump_cat[col_list]} \n")
    sys.exit("Bye Bye.")

    return


def do_dump_full_cat_to_disk(args, mjd, full_cat):

    imjd = int(10000.0*mjd+0.5)
    dump_file = f"full_cat_{args.band_snana}_{imjd}.csv"

    logging.info(f"")
    logging.info(f" Dump entire astropy table to {dump_file}")

    new_order = ['id'] + [name for name in full_cat.colnames if name != 'id']
    t_reordered = full_cat[new_order]
    
    t_reordered.write(dump_file, format='ascii.csv', overwrite=True)
    
    #sys.exit("\n BYE.")
    return

def import_romanisim_tools():
    logging.info("")
    logging.info("!!!!!! import romanisim tools !!!!!!!!!! ")

    global romanisim, galsim, imagelib
    
    logging.info(f"\t import roman from galsim")
    from galsim      import roman   # 8.04.2026  for truth table

    logging.info(f"\t import romanisim")    
    from romanisim   import image as imcode  # 8.04.2026 for truth table
    from romanisim   import ris_make_utils as ris
    from romanisim   import gaia, catalog, log, persistence

    logging.info(f"\t import romanisim.models")    
    from romanisim.models import bandpass, wcs, parameters

    logging.info(f"\t import inject_sources_into_l2")        
    from romanisim.image  import inject_sources_into_l2

    logging.info(f"\t import imagelib")        
    import imagelib

    logging.info("   Done with import romanisim tools")
    return

    
# ==================================
if __name__ == "__main__":

    setup_logging()
    logging.info(f"Begin {' '.join(sys.argv)}")

    args     = get_args()
    config   = read_config(args)    
    cat_dict = {}

    if not args.star_dump:
        mjdmin_list, mjdmax_list = get_mjd_ranges(config)
        config['MA_DICT']        = parse_texpoe_ma_map(args, config)
        config['t_pointing']     = parse_pointing_file(args, config)

    # - - - - -
    STAR_CATALOG      = config.setdefault(KEY_STAR_CATALOG,[])
    GALAXY_CATALOG    = config.setdefault(KEY_GALAXY_CATALOG,{})
    TRANSIENT_CATALOG = config.setdefault(KEY_TRANSIENT_CATALOG,{})
    HAS_TRANSIENTS    = len(TRANSIENT_CATALOG) > 0
    
    # - - - - - - - 
    if len(STAR_CATALOG) > 0:

        # check already prepared csv files
        csv_list = [k for k in STAR_CATALOG if ".csv" in str(k)]
        for csv in csv_list:
            csv = os.path.expandvars(csv)
            key = 'star_' + csv
            cat_dict[key] = get_star_cat_csv(csv)
            cat_dict[key]['label'] = LABEL_STAR
            
        # check GAIA DR3 catalog
        if STAR_CATNAME_GAIA in STAR_CATALOG:
            key = 'star_' + STAR_CATNAME_GAIA
            cat_dict[key]  = get_star_cat_GAIA(args,config)
            cat_dict[key]['label'] = LABEL_STAR
            
        # check synthetic star cat
        if STAR_CATNAME_SYN in STAR_CATALOG:
            key = 'star_' + STAR_CATNAME_SYN
            cat_dict[key]  = get_star_cat_SYN(args,config)   
            cat_dict[key]['label'] = LABEL_STAR

        if args.star_dump: sys.exit('\n Bye bye after star dump')
    # - - - -

    for key, item in GALAXY_CATALOG.items() :
        if key == KEY_SNANA_HOSTLIB_FILE:
            gal_cat = os.path.expandvars(item)
            cat_dict['HOSTLIB_SERSIC0'] = get_galaxy_cat_HOSTLIB(args, config, gal_cat, 0)
            cat_dict['HOSTLIB_SERSIC1'] = get_galaxy_cat_HOSTLIB(args, config, gal_cat, 1)
            cat_dict['HOSTLIB_SERSIC0']['label'] = LABEL_GALAXY
            cat_dict['HOSTLIB_SERSIC1']['label'] = LABEL_GALAXY        

            if args.galid_dump:
                do_galid_dump(args,cat_dict)

        # ... check other source of galaxy catalog ???


    # - - - - - -
    if args.quick_sim_test:
        quick_sim_test(cat_dict, config)  # GAIA and Galaxies only; no transients    

    # - - - - -
    print_banner('# ============== TRANSIENTS ================')
    cat_transient_dict = {}
    unique_mjd_dict   = { args.band_soc : [ mjdmin_list[0] ] }  # default is one epoch if no transients
    texpose = -9.0
    
    for key, item in TRANSIENT_CATALOG.items() :
        if key == KEY_SNANA_TRANSIENT_PATH:
            transient_path  = os.path.expandvars(item)
            genversion_list = []
            item_list       = sorted(glob.glob(f"{transient_path}/*"))  # sim folders + junk files
            
            for item in item_list :
                genv = os.path.basename(item)
                if os.path.isdir(f"{transient_path}/{genv}"):
                    genversion_list.append(genv)
                    temp_dict = get_transient_cat_SNANA(transient_path, genv, args, config)
                    cat_dict[genv] = temp_dict
                    cat_dict[genv]['label'] = LABEL_TRANSIENT
                    if texpose < 0:
                        texpose = get_texpose_transient(cat_dict[genv])  # for optional templates below
                                
            print_banner('# ============== DONE with SNANA TRANSIENTS ================')
            
        if key == KEY_MJD_TEMPLATE :
            # mjds for templates; no transients here
            mjd_template_list = [float(x) for x in item.split() ]
            temp_dict = get_mjd_template_cat(mjd_template_list, texpose, args, config)
            cat_dict['MJD_TEMPLATES'] = temp_dict
            genversion_list += ['MJD_TEMPLATES']
            print_banner('# ============== DONE with MJD_TEMPLATEs ================')

    # - - - -

    if not HAS_TRANSIENTS and config['t_pointing']:
        cat_dict = expand_cat_by_mjd(args, config, cat_dict)
        
    if HAS_TRANSIENTS and len(cat_dict) > 0 :
        unique_mjd_dict = get_unique_mjd_by_band(args, config, cat_dict, genversion_list)
    else:
        unique_mjd_dict, texpose_pointing_dict = \
            get_unique_mjd_pointing(args, config)  # rely on pointing file if no transients
    
    
    # print list of expected asdf files
    print_asdf_expect(args, unique_mjd_dict)
    
    # - - - - - - -  -
    # loop over bands and unique MJDs and run sim for each
    sca = args.scanum
    n_sim_tot = 0
    snid_dump = args.snid_dump
    if snid_dump:  snid_dump_cat_dict = {}
    
    # SMDC
    if not args.check:
        logging.info("")
        logging.info("!!!!!! import romanisim tools !!!!!!!!!! ")
        from galsim      import roman   # 8.04.2026  for truth table
        from romanisim   import image as imcode  # 8.04.2026 for truth table
        from romanisim   import ris_make_utils as ris
        from romanisim   import gaia, catalog, log, persistence
        from romanisim.models import bandpass, wcs, parameters
        from romanisim.image  import inject_sources_into_l2
        import imagelib
        logging.info("   Done with import romanisim tools")

    
    for band_soc, band_snana in zip(BAND_LIST_SOC, BAND_LIST_SNANASIM):
        b_snana      =  band_snana[-1]  # single char band; e.g, 'R', 'Z' etc ...

        # Process only one band input by user, but maybe later allow few bands
        if reject_band(args, band_soc) : continue

        n_sim_expect = len(unique_mjd_dict[band_soc])
        n_sim = 0
        for imjd, mjd in enumerate(unique_mjd_dict[band_soc]):
            n_sim += 1
            n_sim_tot += 1
            
            banner = \
                f"{n_sim} of {n_sim_expect} for {band_soc}-{b_snana}: " \
                f"Stack catalogs for    MJD={mjd:.4f}"
            print_banner(banner, level=9)
            
            full_cat    = cat_vstack(cat_dict, mjd = mjd)

            if HAS_TRANSIENTS:
                texpose  = get_texpose_transient(full_cat)
            else:
                texpose  = texpose_pointing_dict[band_soc][imjd]
            
            ma_num      = get_ma_number(texpose, config['MA_DICT'], True)
            check_pointing_dict = { 'band':b_snana, 'mjd': mjd }
            ra,dec,roll = get_pointing_info(args, config, full_cat, check_pointing_dict, True)

            info_dict = { 'sca':sca,  'band_soc':band_soc, 'mjd':mjd, 'ma_num':ma_num, 
                          'ra':ra, 'dec':dec, 'roll':roll   }
            filename_asdf = run_sim(args, full_cat, info_dict)

            #if abs(args.mjd_dump-mjd) < .0001 :
            if args.mjd_dump > 1.0 :
                do_dump_full_cat_to_disk(args, mjd, full_cat)
            
            if snid_dump:
                snid_dump_cat =   full_cat[ full_cat['SNID'] == snid_dump ]
                if len(snid_dump_cat) > 0:
                    snid_dump_cat_dict[f"{snid_dump}-{mjd}"] = snid_dump_cat
            
    # - - - - -
    if snid_dump:
        snid_dump_cat = cat_vstack(snid_dump_cat_dict, None)
        logging.info(f"")
        logging.info(f" Astropy Table dump for snid = {snid_dump}:")
        print(f"{snid_dump_cat}\n")
        
    logging.info('')
    logging.info(f"Done with {n_sim_tot} SCA image sims")
    
    # === END: ===
