#!/usr/bin/env python

# Created Jul 10 2026
# Read bands and SCA range from sim config, and prepare
# sbatch and sh file for each task. ALso produce
# RUN_ALL_SIM.sh and RUN_ALL_CAL.sh to launch everthing.
#
# How to use the RUN* outputs of this script:
#
# ./RUN1_SIM.sh
# ./RUN2_CAL.sh &    [it will wait until SIM is done so you can sleep]
#
# when SIM and CAL are all done (see STATUS_SIM.DAT and STATUS_CAL.DAT)
#  ./RUN3_SUMMARY.sh
#  ./RUN4_CLEAN.sh
#
import os, sys, glob, argparse, logging, time, datetime, yaml
import subprocess


NTASK_UPD_STDOUT = 50
CONFIG_SBATCH_PREP = "SBATCH_PREP"

SIF_FILE     = "/data/snpit/roman-snpit-env-cpu.sif"

TASKNAME_SIM = "SIM"
TASKNAME_CAL = "CAL"
TASKNAME_LIST = [ TASKNAME_SIM, TASKNAME_CAL ]
TASKNAME_SUMMARY = 'SUMMARY'
TASKNAME_CLEAN   = 'CLEAN'

PREFIX_SBATCH_DICT = {
    TASKNAME_SIM : "run_sim",
    TASKNAME_CAL : "run_cal"
}

CODE_WRAPPER_DICT = {
    TASKNAME_SIM : f"$HOME/bin/romanisim_snpit_wrapper.py",
    TASKNAME_CAL : f"$HOME/bin/romancal_snpit_wrapper.py",
}

RUNNUM_ORDER_DICT = {
    TASKNAME_SIM     : 1,
    TASKNAME_CAL     : 2,
    TASKNAME_SUMMARY : 3,
    TASKNAME_CLEAN   : 4
}

#SUBMIT_ORDER_DICT = {
#    TASKNAME_SIM : 1,
#    TASKNAME_CAL : 2
#}

# define map between SNANA single-char band and SOC band name
BAND_LIST_DICT  = {  
    'R': 'F062',  'Z': 'F087',  'Y': 'F106',  'J':'F129',  'H': 'F158',  'F': 'F184',  'K': 'F213'
}


BACKUP_PREFIX = 'BACKUP'

# ============================
# ============================
def setup_logging():
    logging.basicConfig(level=logging.INFO,
                        format="[%(levelname)6s %(filename)11.11s ] %(message)s")  
    return

def get_args():
    parser_inputs = argparse.ArgumentParser()
    
    msg = "name of input config file"
    parser_inputs.add_argument("--config_file", "-c", help=msg, type=str, default=None)

    # - - - actions: prep, clean, status
    msg = "action: prepare run scripts for sbatch"
    parser_inputs.add_argument("--prep", "-p", help=msg, action="store_true")

    msg = "action: get status update for SIM or CAL"
    parser_inputs.add_argument("--status", "-s", help=msg, type=str, default=None)

    msg = "action: launch RUN2_CAL when SIM is done (combine with --status)"
    parser_inputs.add_argument("--run_cal_when_done", help=msg, action="store_true")
    
    msg = "action: make final summary table needed for data base ingestion (after romancal finishes)"
    parser_inputs.add_argument("--summary", help=msg, action="store_true")

    msg = "action: cleanup by making backup tar files (after romancal finishes)"
    parser_inputs.add_argument("--clean", help=msg, action="store_true")

    msg = "action: print info for each run* file processed"
    parser_inputs.add_argument("--verbose", help=msg, action="store_true")          
    
    args = parser_inputs.parse_args()

    narg_missing = 0
    requred_arg_list = ['config_file'] 
    for arg_name in requred_arg_list:
        value = getattr(args, arg_name, None)
        if value is None:
            logging.info(f"ERROR: missing required arg --{arg_name}.")
            narg_missing += 1

    if narg_missing:
        sys.exit(f"\n\t ABORT on {narg_missing} missing required args.")
        
    return args


def compress_files(flag, dir_name, wildcard, name_backup, wildcard_keep ):

    # name of tar file is BACKUP_{name_backup}.tar
    # Inputs
    #  flag > 0 -> compress
    #  flag < 0 -> uncompress
    #  dir_name -> cd to this directory
    #  wildcard -> include these files in tar file
    #  name_backup -> tar file name is BACKUP_{name_backup}.tar
    #        if name_backup has .tar extension, the use this name with BACKUP prefix.
    #  wildcard_keep -> do NOT remove these files
    #

    if '.tar' in name_backup:
        tar_file   = name_backup
    else:
        tar_file   = f"{BACKUP_PREFIX}_{name_backup}.tar"

    targz_file = f"{tar_file}.gz"
    cddir      = f"cd {dir_name}"

   # be careful if wildcard string is too short; don't want rm *
    if len(wildcard) < 3 :
        msgerr = []
        msgerr = f"wildcard = '{wildcard}' is dangerously short string"
        msgerr = f"that could result in removing too much."
        msgerr = f"Provide longer wildcard string."
        log_assert(False,msgerr)

    if flag > 0 :
        cmd_tar  = f"tar -cf {tar_file} {wildcard} "
        cmd_gzip = f"gzip {tar_file}"

        if len(wildcard_keep) == 0 :
            cmd_rm   = f"rm {wildcard}"
        else:
            # remove all wildcard files EXCEPT for wildcard_keep
            cmd_rm = f"find {wildcard} ! -name '{wildcard_keep}' " + \
                     "-type f -exec rm {} +"

        cmd_all  = f"{cddir} ; {cmd_tar} ; {cmd_gzip} ; {cmd_rm} "
    else:
        cmd_unpack = f"tar -xzf {targz_file}"
        cmd_rm     = f"rm {targz_file}"
        cmd_all    = f"{cddir} ; {cmd_unpack} ; {cmd_rm} "

    os.system(cmd_all)

    # end compress_files

def get_log_file_name(prefix):

    log_file = prefix + '.log'  # default

    # if log_file exists, change the name to preserve previous log files
    # e.g., if prefix.log exists, return prefix.log-1.
    # e.g., if prefix.log and prefix.log-1 exists, return prefix.log-2 ...etc ...
    
    log_list = sorted(glob.glob(f"{log_file}*"))
    n_log = len(log_list)
    if n_log > 0:
        log_file = f"{log_file}-{n_log}"
        
    return log_file


def prep_task(TASK, b, sca, SBATCH_PREP, args, i_split, i_task):

    # Created .sbatch and .sh file for task specified by:
    #   TASK        = SIM or CAL
    #   b           = single-char snana band name ; e.g 'J' or 'Z'
    #   sca         = sca number (1-18)
    #   SBATCH_PREP = sbatch info from user config
    #   i_split     = split task for this band (0 to Nsplit-1)
    #   i_task      = task number (0 to NTASK-1) used to determine which NODE
    
    config_file      = args.config_file
    n_split          = SBATCH_PREP['BAND_NSPLIT_DICT'][b]
    
    PREFIX      = PREFIX_SBATCH_DICT[TASK]
    CODE        = CODE_WRAPPER_DICT[TASK]
    
    b_soc       = BAND_LIST_DICT[b]        # e.g.,  b(SNANA)=R -> b_soc = F062
    sca2d       = f"{sca:02d}"
    
    prefix      = f"{PREFIX}_{sca2d}-{b}"
    if n_split > 1:
        prefix += f"_split{i_split}"
    
    bash_file   = f"{prefix}.sh"
    sbatch_file = f"{prefix}.sbatch"
    done_file   = f"{prefix}.done"
    log_file    = get_log_file_name(prefix)
    
    n_node      = SBATCH_PREP['N_NODE']
    i_node      = i_task % n_node

    node        = SBATCH_PREP['NODE_LIST'][i_node] # round-robin node for each task
    walltime    = SBATCH_PREP['WALLTIME']

    
    # start with sbatch script
    if args.verbose:
        logging.info(f" Prepare {sbatch_file} ")
    else:
        if i_task % NTASK_UPD_STDOUT == 0 and i_task > 0 and TASK == TASKNAME_SIM:
            logging.info(f" Prepare i_task = {i_task}")
    
    with open(sbatch_file,"wt") as s:
        s.write(f"#!/bin/bash\n")
        s.write(f"#SBATCH -p {node} \n")
        s.write(f"#SBATCH --time {walltime} \n")
        s.write(f"#SBATCH --output {log_file} \n")
        s.write(f"\n")
        s.write(f"singularity run \\\n")
        s.write(f"\t{SIF_FILE} \\\n")
        s.write(f"\t/bin/bash  \\\n")
        s.write(f"\t{bash_file} \n")

    # next write the bash script that call romanisim wrapper

    with open(bash_file,"wt") as f:
        f.write(f"source $HOME/.bashrc \n")
        f.write(f"source $HOME/new_venv/bin/activate\n")
        f.write(f"\n")
        
        if TASK == TASKNAME_SIM:
            arg_list = [ f"{config_file}", "--level L1", f"--scanum {sca}", f"--band_snana {b}" ]
            
            keyname_jobsplit   = SBATCH_PREP['BAND_SPLIT_KEYNAME'][b]
        else:
            arg_list = [ f"--l1_asdf_file 'SNPIT*WFI{sca2d}_{b_soc}_L1.asdf'" ] # fragile alert
            keyname_jobsplit   = '--jobsplit'  # for romancal, split does not depend on MJD vs texpose

        if n_split > 1:
            arg_split = f"{keyname_jobsplit} {i_split} {n_split} "
            arg_list += [ arg_split ]
            
        f.write(f"{CODE} \\\n")
        for arg in arg_list:
            f.write(f"\t{arg} \\\n")
        
        f.write(f"\n")

        # touch done file
        f.write(f"touch {done_file}\n")

        # update STATUS
        cmd_status = SBATCH_PREP['cmd_status_dict'][TASK]
        f.write(f"{cmd_status}\n")
        
        f.write(f"\n")
        
    return sbatch_file, node
# end prep_task

def prep_task_driver(args, SBATCH_PREP):

    band_list = SBATCH_PREP['BAND_LIST']  # list of single char filter names
    sca_min   = SBATCH_PREP['SCA_MIN']
    sca_max   = SBATCH_PREP['SCA_MAX']

    sim_sbatch_file_list = []
    cal_sbatch_file_list = []
    i_task = 0

    nodelist_sim_file = nodelist_file_name(TASKNAME_SIM)
    nodelist_cal_file = nodelist_file_name(TASKNAME_CAL)    
    fp_node_sim       = open(nodelist_sim_file, "wt")
    fp_node_cal       = open(nodelist_cal_file, "wt")    

    logging.info(f"")
    logging.info(f"Begin preparing sbatch tasks:")
    for b in band_list:
        n_split = SBATCH_PREP['BAND_NSPLIT_DICT'][b]
        for i_split in range(0,n_split):
            
            for sca in range(sca_min,sca_max+1):            
                sim_sbatch_file, node = \
                    prep_task(TASKNAME_SIM, b, sca, SBATCH_PREP, args, i_split, i_task)
                sim_sbatch_file_list.append(sim_sbatch_file)
                
                cal_sbatch_file, node = \
                    prep_task(TASKNAME_CAL, b, sca, SBATCH_PREP, args, i_split, i_task)
                cal_sbatch_file_list.append(cal_sbatch_file)            

                fp_node_sim.write(f"{sim_sbatch_file}:  {node}\n")
                fp_node_cal.write(f"{cal_sbatch_file}:  {node}\n")                
                i_task += 1

    fp_node_sim.close()
    fp_node_cal.close()    
    # - - - - - - 
    logging.info(f"")
    create_run_all_script(TASKNAME_SIM, sim_sbatch_file_list, SBATCH_PREP)
    create_run_all_script(TASKNAME_CAL, cal_sbatch_file_list, SBATCH_PREP)

    create_postproc_scripts(args, SBATCH_PREP)
    
    set_chmod()

    logging.info(f" Done.")
    return

# end prep_task_driver

def nodelist_file_name(task):
    return f"NODELIST_{task}.DAT"

def parse_config(args):

    config_file = args.config_file
    
    path_expand = os.path.expandvars(config_file)
    logging.info(f"Reading YAML input from {path_expand}")
    with open(path_expand) as f:
        config = yaml.safe_load(f.read())

    SBATCH_PREP = config[CONFIG_SBATCH_PREP]

    # add convenient items based on user input
    SCA_RANGE = SBATCH_PREP['SCA_RANGE']
    SBATCH_PREP['SCA_MIN']   = int(SCA_RANGE.split()[0])
    SBATCH_PREP['SCA_MAX']   = int(SCA_RANGE.split()[1])

    SBATCH_PREP = parse_config_bands(SBATCH_PREP)

    SBATCH_PREP = parse_config_nodes(SBATCH_PREP)
    
    #sys.exit(f"\n xxx SBATCH_PREP = \n{SBATCH_PREP}")
    
    # store commands to monitor status
    jobname = sys.argv[0]
    conf_arg   = f"--config_file {config_file}"  # common args for SIM and CAL tasks
    cmd_status_dict = {
        TASKNAME_SIM : f"{jobname} {conf_arg} --status {TASKNAME_SIM}", #  ??  --run_cal_when_done",
        TASKNAME_CAL : f"{jobname} {conf_arg} --status {TASKNAME_CAL}"
    }
    SBATCH_PREP['cmd_status_dict'] = cmd_status_dict
    
    return SBATCH_PREP

def get_prefix_sbatch(task, sca, b, i_split):
    PREFIX    = PREFIX_SBATCH_DICT[task]
    prefix   = f"{PREFIX}_{sca:02d}-{b}"

    return prefix

def parse_config_nodes(SBATCH_PREP):

    SBATCH_PREP['NODE_LIST'] = SBATCH_PREP['NODES'].split()
    SBATCH_PREP['N_NODE']    = len(SBATCH_PREP['NODE_LIST'])

    return SBATCH_PREP

def parse_config_bands(SBATCH_PREP):

    BANDS          = SBATCH_PREP['BANDS']
    band_list_full = BANDS.split()  # includes extra chars; e.g. Z*2 or F/4
    SBATCH_PREP['BAND_LIST']      = [ b[0] for b in band_list_full ]  # keep single char only
    SBATCH_PREP['BAND_NSPLIT_DICT']   = {}
    SBATCH_PREP['BAND_SPLIT_KEYNAME'] = {}  
    
    # check each band for split; e.g. Z*2 -> split Z band into 2 tasks
    for band_full in band_list_full:

        njobsplit = 1  # default is one job per band
        b         = band_full[0]  # first char is single-char representation
        SBATCH_PREP['BAND_SPLIT_KEYNAME'][b] = None
            
        if '*' in band_full:
            njobsplit = int(band_full[-1])  # assumes number is single digit
            SBATCH_PREP['BAND_SPLIT_KEYNAME'][b] = '--jobsplit_mjd'
            
        if '/' in band_full:
            njobsplit = int(band_full[-1])  # assumes number is single digit
            SBATCH_PREP['BAND_SPLIT_KEYNAME'][b] = '--jobsplit_texpose'
            
        SBATCH_PREP['BAND_NSPLIT_DICT'][b] = njobsplit

    #sys.exit(f"\n xxx SBATCH_PREP = \n{SBATCH_PREP}")
    return SBATCH_PREP

def set_chmod():
    logging.info(f" Set chmod +x run* RUN* ")
    os.system(f"chmod +x run*")
    os.system(f"chmod +x RUN*")    
    return

def create_postproc_scripts(args, SBATCH_PREP):

    # here PPTASK refers to postproc task (summary or clean);
    # NOT sim or cal image task.
    
    PPTASK_LIST = [ TASKNAME_SUMMARY, TASKNAME_CLEAN ]
    
    for pptask in PPTASK_LIST:
        ORDER  = RUNNUM_ORDER_DICT[pptask]
        script = f'RUN{ORDER}_{pptask}.sh'

        str_image_tasks = ' & '.join(TASKNAME_LIST)
        logging.info(f" ./{script} for final {pptask} after {str_image_tasks} have finished")
    
        with open(script,"wt") as s:
            key = '--' + pptask.lower()
            s.write(f"{sys.argv[0]} --config {args.config_file} {key}\n")
        
    return

def get_filename_timestamp(task,when):
    if when == 0:
        return f"TIMESTAMP_START_{task}"
    else:
        return f"TIMESTAMP_DONE_{task}"
# end get_filename_timestamp

    
def create_run_all_script(task_name, sbatch_file_list, SBATCH_PREP):

    ORDER       = RUNNUM_ORDER_DICT[task_name]  # 1 for SIM, 2 for CAL

    prefix      = f"RUN{ORDER}_ALL_" + task_name
    script_file = prefix + '.sh'
    start_file  = get_filename_timestamp(task_name,0)
    n_task      = len(sbatch_file_list)
    
    comment = ''
    if task_name == TASKNAME_CAL:
        comment = f'(wait until {TASKNAME_SIM} tasks finish)'
        
    logging.info(f" ./{script_file} to launch {n_task} sbatch jobs for {task_name}    {comment}")

    with open(script_file, "wt") as s:

        # for CAL task, what for SIM to finish
        if  task_name == TASKNAME_CAL :
            sim_done_file = get_filename_timestamp(TASKNAME_SIM,1)
            s.write(f"while [ ! -f {sim_done_file} ]; do sleep 60; done\n\n")
            
        s.write(f"touch {start_file}\n\n")
        
        for sbatch_file in sbatch_file_list:
            s.write(f"sbatch {sbatch_file}\n")

        # leave message to monitor status
        cmd_status = SBATCH_PREP['cmd_status_dict'][task_name]
        s.write(f"\n")
        s.write(f"echo # ---------------------------------------------------- \n")        
        s.write(f"echo Monitor status with command \n")
        s.write(f"echo {cmd_status} \n")
    return
# end create_run_all_script

def get_status_driver(args, SBATCH_PREP):

    TASK        = args.status.upper()

    if TASK not in list(PREFIX_SBATCH_DICT.keys() ):
        sys.exit(f"Invalid task name = {TASK} for status. \nValid task names are {TASKNAME_LIST}") 
        
    STATUS_FILE = 'STATUS_' + TASK + '.DAT'
    logging.info(f"Start preparing {STATUS_FILE}")

    # read node for each sbatch file
    NODELIST_FILE = nodelist_file_name(TASK)
    with open(NODELIST_FILE) as f:
        nodelist_map = yaml.safe_load(f.read())
        
    STATUS_WAIT = "WAIT"
    STATUS_RUN  = "RUN"
    STATUS_DONE = "DONE"
    STATUS_OOM  = "OOM"   # out of memory 
    STATUS_FAIL = "FAIL"  # e.g., Traceback detected

    STORE_PREVIOUS = False
    if STORE_PREVIOUS:
        if os.path.exists(STATUS_FILE) :
            cmd = f"mv {STATUS_FILE} {STATUS_FILE}_PREVIOUS"
            os.system(cmd)

    prefix      = PREFIX_SBATCH_DICT[TASK]
    sbatch_list = sorted(glob.glob(f"{prefix}*.sbatch"))

    status_lines =  []

    status_lines.append(f"TASKNAME             " \
                        f"STATUS    NODE        NIMG_DONE  NIMG_TOT   CPU_TOT(hr)   CPU_PER_IMG(hr)")

    # for each sbatch script, check for corresponding log file
    sum_njob_tot = 0;  sum_njob_done = 0; sum_cpu_tot = 0.0 
    sum_sbatch_fail = 0;  sum_sbatch_oom = 0;

    n_log_expect = len(sbatch_list)
    n_log = 0
    
    for sbatch in sbatch_list:
        PREFIX = sbatch.split('.')[0]

        log_base = f"{PREFIX}.log"
        log_list = sorted(glob.glob(f"{log_base}*"))

        if len(log_list) == 0 : continue
        log_file = log_list[-1]  #  current log is last log in list
        
        exist    = os.path.exists(log_file)
        n_log   += 1
        
        if args.verbose:
            logging.info(f"\t {log_file} exists = {exist}")
        else:
            upd = (n_log > 1 and n_log % NTASK_UPD_STDOUT == 0) or (n_log==n_log_expect)
            if upd:
                prefix_common = PREFIX_SBATCH_DICT[TASK]
                logging.info(f"\t Processed {n_log:4d} of {n_log_expect:4d} {prefix_common}*log files")
        
        
        status = STATUS_WAIT

        node = nodelist_map[sbatch]
        
        if exist:
            njob_tot = 0 ; njob_done = 0; cpu_tot = 0.0

            n_expect   = grep_wrapper('N_EXPECT_ASDF', log_file)            
            cpu_list   = grep_wrapper('CPU(SNPIT',     log_list)
            oom        = grep_wrapper('oom',           log_file)   # check for out of memory
            fail       = grep_wrapper('Traceback',     log_file)   # check for crash
            
            if n_expect:
                njob_tot = int(n_expect[0])

            if cpu_list:
                njob_done = len(cpu_list)
                cpu_tot   = sum( [ float(x) for x in cpu_list ] )
                cpu_tot /= 3600.

            if njob_tot > 0:
                status = STATUS_RUN
                
            if njob_done == njob_tot and njob_tot > 0 :
                status = STATUS_DONE

            if oom:
                status = STATUS_OOM
                sum_sbatch_oom += 1
                
            if fail:
                status = STATUS_FAIL
                sum_sbatch_fail += 1
                
            sum_njob_tot  += njob_tot
            sum_njob_done += njob_done
            sum_cpu_tot   += cpu_tot
            if njob_done > 0:
                cpu_per_img = cpu_tot / njob_done
            else:
                cpu_per_img = 0.0

            line = f"{PREFIX:<20}   {status:<4}  {node:<12}  {njob_done:8d}  {njob_tot:8d}  " \
                f"{cpu_tot:9.1f}   {cpu_per_img:12.2f}"
            status_lines.append(f"{line}")
        
    # - - - - -
    walltime = get_walltime(TASK)
    
    status_lines.append(f"")
    status_lines.append(f"# ----------------------------------------------------------------")    
    status_lines.append(f"SUM_NIMG_TOT:   {sum_njob_tot:4d}    # total number of expected images")
    status_lines.append(f"SUM_NIMG_DONE:  {sum_njob_done:4d}    # total number of produced images")    
    status_lines.append(f"SUM_CPU_TOT:    {sum_cpu_tot:.1f}        # hr")
    status_lines.append(f"WALLTIME:       {walltime:.3f}       # hr")
    status_lines.append(f"N_SBATCH_OOM:   {sum_sbatch_oom:4d}      # out of memory")
    status_lines.append(f"N_SBATCH_FAIL:  {sum_sbatch_fail:4d}      # Traceback detected")
    
    status_all = STATUS_RUN ; comment = 'Some tasks still running or pending'
    if sum_njob_done == sum_njob_tot and sum_njob_tot > 0:
        status_all = STATUS_DONE ; comment = 'All tasks complete with no detectable errors'
    if sum_sbatch_fail > 0 or sum_sbatch_oom > 0:
        status_all = STATUS_FAIL ; comment = 'At least one task has Traceback or OOM'
        
    status_lines.append(f"STATUS_ALL:    {status_all}   # {comment}")
    status_lines.append(f"# ----------------------------------------------------------------")
    
    # - - - - -
    # write all the status lines here as quickly as possible to avoid
    # conflicts with multiple run* tasks 
    with  open(STATUS_FILE,"wt") as s:
        for line in status_lines:
            s.write(f"{line}\n")

    # create grand DONE file for this task when everything is really done
    if status_all == STATUS_DONE:
        cmd_touch_done = "touch " + get_filename_timestamp(TASK,1)
        os.system(cmd_touch_done)  
            
    # check option to launch CAL when SIM is done
    IS_SIM  = (TASK       == TASKNAME_SIM)
    IS_DONE = (status_all == STATUS_DONE)
    if IS_SIM and IS_DONE and args.run_cal_when_done:
        try:
            run_cal_script = glob.glob("RUN*CAL.sh")[0]
            logging.info(f"All {TASK} tasks are done --> launch {run_cal_script}")
            cmd            = f"./{run_cal_script}"
        except:
            pass

    logging.info(f"Done updating {STATUS_FILE}")
    return
# end get_status_driver

def get_walltime(task):

    from pathlib import Path
    from datetime import datetime, timezone

    walltime = 0.0
    start_file = get_filename_timestamp(task,0)

    if not os.path.exists(start_file):
        logging.info(f"")
        logging.info(f"# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        logging.info(f"  WARNING: {start_file} does not exist -> cannot determine wall time.")
        logging.info(f"# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        logging.info(f"")
        return walltime

    stat = Path(start_file).stat()
    created_timestamp = stat.st_ctime
    created_time      = datetime.fromtimestamp(created_timestamp, tz=timezone.utc)
    current_time      = datetime.now(timezone.utc)

    # Calculate elapsed walltime
    walltime = (current_time - created_time).total_seconds()
    walltime /= 3600.0
    
    return walltime
# end get_walltime

def grep_wrapper(key, log_file):

    # regurn list of values after key in log_file.
    grep_list = []
    
    try:
        cmd_list = ['grep', key ]
        if isinstance(log_file, list):
            cmd_list += log_file  # check several logs
        else:
            cmd_list.append(log_file)

        
        grep_output = subprocess.check_output(cmd_list, text=True)
        wd_list = grep_output.split()
        for i, wd in enumerate(wd_list):
            if key in wd:
                val = wd_list[i+1]
                grep_list.append(val)
    except:
        grep_list = None
    
    return grep_list


def summary_driver(args, SBATCH_PREP):

    # make summary table for every L2 image; used for data base ingestion

    #config_file   = args.config_file

    summary_table = "SUMMARY_L2.DAT"

    # get RA and DEC to use as goofy addition to obs_id
    config_file = os.path.expandvars(args.config_file)
    with open(config_file) as f:
        config = yaml.safe_load(f.read())
        #sys.exit(f"\n xxx config = \n{config}")
    RA   = config['SKY_REGION']['RA_CEN']
    DEC  = config['SKY_REGION']['DEC_CEN']
    ROLL = config['SKY_REGION']['ROLL']    
    obs_idoff = abs(int(RA * DEC)) + int(ROLL)
    
    img_list = sorted(glob.glob("SNPIT*L2.asdf"))
    n_img = len(img_list)
    
    with open(summary_table,"wt") as s:
        s.write(f"# obs_id = int(MJD*1E10) + abs(int(RA*DEC)) + int(ROLL)\n")
        s.write(f"# WFI_RA(center)  = {RA}\n")
        s.write(f"# WFI_DEC(center) = {DEC}\n")
        s.write(f"# WFI_ROLL angle  = {ROLL}\n")
        s.write(f"#\n")
        s.write("VARNAMES:  IMGNAME  observation_id   SCA   BAND \n")
        for img in img_list:
            obs_id, scanum, band = extract_info_from_img_filename(img, obs_idoff)
            s.write(f"ROW: {img}  {obs_id}  {scanum}  {band}\n")
                
    
    logging.info(f"See observation summary for {n_img} images in {summary_table}")
    sys.exit()

    return
# end summary_driver

def extract_info_from_img_filename(img_file, obs_idoff):

    # extra info from image file name,
    # Beware that getting info from a filename is generally bad practice.

    visit = img_file.split('VISIT')[1][0:9]
    obs_id = int(visit) * 1000000 + obs_idoff

    scanum = img_file.split('WFI')[1][0:2]

    band = img_file.split('WFI')[1][3:7]
    
    return obs_id, scanum, band

def cleanup_driver(args, SBATCH_PREP):

    logging.info(f" CLEANUP:")
    flag_compress = +1
    dir_name = './'
    wildcard_keep = []
    
    wildcard_list = [
        'run_sim*.sh', 'run_sim*.sbatch', 'run_sim*log*', 'run_sim*done',
        'run_cal*.sh', 'run_cal*.sbatch', 'run_cal*log*', 'run_cal*done',
        'SNPIT*L1.png', 'SNPIT*L2.png', 'SNPIT*cat.parquet'
    ]
    name_backup_list = [
        'sim_sh', 'sim_sbatch', 'sim_log', 'sim_done',
        'cal_sh', 'cal_sbatch', 'cal_log', 'cal_done',
        'SNPIT_L1_png', 'SNPIT_L2_png', 'SNPIT_L1_cat.parquet'
    ]
    

    for wildcard, name_backup in zip(wildcard_list,name_backup_list):
        nf = len(glob.glob(wildcard))
        logging.info(f" Create tarball for {nf:4d}  {wildcard:<20} files")
        compress_files(flag_compress, dir_name, wildcard, name_backup, wildcard_keep )
        
    return

# ==================================
if __name__ == "__main__":

    setup_logging()

    logging.info(f"")
    
    args  = get_args()

    SBATCH_PREP = parse_config(args)    

    if args.prep:
        prep_task_driver(args, SBATCH_PREP)

    elif args.status:
        get_status_driver(args, SBATCH_PREP)

    elif args.summary:
        summary_driver(args, SBATCH_PREP)

    elif args.clean:
        cleanup_driver(args, SBATCH_PREP)             

    else:
        sys.exit(f"\n ERROR: no action specified. Try --prep or --status or --clean")
    # == END: ===
    
