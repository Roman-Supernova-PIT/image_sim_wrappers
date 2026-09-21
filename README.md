# image_sim_wrappers

## Quick-start guide with slurm (beware: tested only on SMDC)

- Simulate catalog-level transients with SNANA. The SIMLIB/cadence library used
  by SNANA should be synced with a POINTING file containing\
  &nbsp;&nbsp;&nbsp;&nbsp;    RA_WFI_CEN, DEC_WFI_CEN, ROLL
  for each visit that is identified by an image number (IMGNUM). This pointing info
  is passed to romanisim for each MJD/IMGNUM. If there is no POINTING file, you can
  hard wire a fixed RA,DEC,ROLL in the wrapper-config file (see below).

   To simulate only stars and/or galaxies, a manually created POINTING file is
   required to define a list of MJDs, since there are no transients to define
   the MJDs.
   
   SNANA is not yet working on SMDC, so need to run this on another cluster
    
- Copy the SNANA HOSTLIB and output sim data folder to SMDC

- Construct input config for the roman code wrappers using help from  
   **romanisim_snpit_wrapper.py --HELP**

- Prepare slurm jobs with\
    **sbatch_prep_romanisim+romancal.py --config_file sim_science.config --prep**\
        or\
    **sbatch_prep_romanisim+romancal.py --c sim_science.config -p**\

- Launch romanisim jobs with  **./RUN1_ALL_SIM.sh**;
   monitor progress in STATUS_SIM.DAT, which includes WALLTIME and failure stats.

   Don't panic if there is no STATUS file for a while; it won't appear until at
   least one slurm task has finished.
      
- launch romancal with **./RUN2_ALL_CAL.sh**; it will wait for RUN1_ALL_SIM.sh to finish.
   Monitor status with STATUS_CAL.DAT
   
- When all romancal jobs are done, make grand summary with **./RUN3_SUMMARY.sh**

- If all looks good, clean up some of the mess by creating BACKUP*tar files
   using **./RUN4_CLEAN.sh **

## Wrapper Descriptions:

- **romanisim_snpit_wrapper.py**\
   Read galaxies and transients from SNANA sim data folder, and
   read GAIA+Synthetic stars from CSV file; this input is translated
   into required astropy format for romanisim. The MJDs in the
   transient sim determine the simulated MJDs, with optional
   input MJD_RANGE cuts. Output includes
    - L1 asdf file for each MJD in transient file
    - TRUTH table of overlaid objects in each L1 file;\
         id, ra, dec, x_det, y_det, n, label, flux_maggie, mag, mjd, band, sca
	 
  For help on the input config file,\
      **romanisim_snpit_wrapper.py --HELP**
   
- **romancal_snpit_wrapper.py**\
  read L1 file (or wildcard for list) and run romancal to produce
  L2 file for each L1. Beware that there is currently a sublte WCS
  problem in the L2 output, and a hacky fix is needed for campari
  and phrosty.


- **sbatch_prep_romanisim+romancal.py**\
  Read config file for romanisim_snpit_wrapper.py, and prepare slurm jobs
  (SMDC) for both romanisim_snpit_wrapper.py & romancal_snpit_wrapper.py.
  Also produces simple bash scripts to submit all of the jobs; see
  RUN*sh above in quic-start guide.
  The sbatch instructions are optionally included in the SBATCH_PREP
  block of the config input for romanisim_snpit_wrapper.py.

  This process includes monitor process that updates\
  &nbsp;&nbsp;&nbsp;&nbsp;    STATUS_SIM.DAT\
  &nbsp;&nbsp;&nbsp;&nbsp;    STATUS_CAL.DAT\
  to show progress of slurm jobs.

  For more help, see SBATCH_PREP block in  **romanisim_snpit_wrapper.py --HELP**
  

