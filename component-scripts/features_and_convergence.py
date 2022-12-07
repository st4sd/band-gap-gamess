#!/usr/bin/env python

# Copyright IBM Inc. 2015, 2019. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
# Author(s):
#   James McDonagh
#   Vassilis Vassiliadis

import glob
import logging
import os

import argparse
import re
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

__title__ = os.path.basename(__file__)
__version__ = 1.0
__copyright__ = "Copyright IBM 2020"


def setup_logger(cwd, loglev="INFO"):
    """
    Make logger setup
    INPUT :
        cwd : the working directory you want the log file to be saved in
    OUTPUT:
        FILE: log file
    """
    # set log level from user
    intloglev = getattr(logging, loglev)
    try:
        intloglev + 1
    except TypeError:
        print("ERROR - cannot convert loglev to numeric value using default of 20 = INFO")
        with open("error_logging.log", "w+") as logerr:
            logerr.write("ERROR - cannot convert loglev to numeric value using default of 20 = INFO")
        intloglev = 20

    # Format routine the message comes from, the leve of the information debug, info, warning, error, critical
    # writes all levels to teh log file Optimizer.log
    logging.raiseExceptions = True
    log = logging.getLogger()
    log.setLevel(intloglev)
    pathlog = os.path.join(cwd, "{}.log".format(__title__.split(".")[0]))

    # File logging handle set up
    filelog = logging.FileHandler("{}".format(pathlog), mode="w")
    filelog.setLevel(intloglev)
    fileformat = logging.Formatter("%(levelname)s - %(name)s - %(message)s")
    filelog.setFormatter(fileformat)

    # Setup handle for screen printing only prints info and above not debugging info
    screen = logging.StreamHandler()
    screen.setLevel(10)
    screenformat = logging.Formatter('%(message)s')
    screen.setFormatter(screenformat)

    # get log instance
    log.addHandler(screen)
    log.addHandler(filelog)

    log.info("The handlers {} logging level {} {}".format(log.handlers, loglev, intloglev))
    log.info('Started {}\n'.format(datetime.now()))

    return log

def plot_convergence(data, headers=None, key="max_gradient", gamess_default_gradient_tolerence=0.0001):
    """
    function to plot convergence metrics
    :param data: numpy array - convergence data with columns matching the headers order
    :param headers: list - headers of the oclumns to say what the data is for
    :param key: string - the header you want to plot to check convegence on
    :param gamess_default_gradient_tolerence: The defult convergence threshold you want to show on the plots
    """

    log = logging.getLogger(__name__)

    # plot 1 - total iteration Vs. energy
    print("Plotting iter against energy")
    fig = plt.figure(figsize=(30, 30))
    try:
        x = np.array(data[:, 1], dtype=np.int)
    except IndexError as err:
        log.warning("Data index error likely the gamess output file is incompleted. convergence plotting can not be completed. {}".format(err))
        return
    y = np.array(data[:, 3], dtype=np.float64)
    plt.plot(x, y, "ro-")
    plt.xlabel("Total iteration count", fontsize=40)
    plt.xticks(fontsize=20)
    plt.ylabel("Energy (Hartrees)", fontsize=40)
    plt.yticks(fontsize=20)
    plt.savefig("total_iter_vs_energy.png")
    plt.close()

    # plot 2 - total iteratiob Vs. max gradient
    print("Plotting iter against max gradient")
    fig = plt.figure(figsize=(20, 20))
    y = np.array(data[:, 4], dtype=np.float64)
    dy = y - np.array([gamess_default_gradient_tolerence]*len(x))
    plt.plot(x, y, "ro-", label="optimization run")
    plt.plot(x, dy, "co-", label="optimization run difference compared to tolerence")
    plt.plot(x, [0.0001] * len(x), "k.-", label="GAMESS default gradient tolerance")
    plt.legend(fontsize=20)
    plt.xlabel("Total iteration count", fontsize=40)
    plt.xticks(fontsize=20)
    plt.ylabel("Max Gradient", fontsize=40)
    plt.yticks(fontsize=20)
    plt.savefig("total_iter_vs_max_gradient.png")
    plt.close()

    # plot 3 - total iteratiob Vs. rms
    print("Plotting iter against rms")
    fig = plt.figure(figsize=(20, 20))
    y = np.array(data[:, 5], dtype=np.float64)
    plt.plot(x, y, "ro-")
    plt.xlabel("Total iteration count", fontsize=40)
    plt.xticks(fontsize=20)
    plt.ylabel("RMS", fontsize=40)
    plt.yticks(fontsize=20)
    plt.savefig("total_iter_vs_rms.png")
    plt.close()

    # plot 4 - total iteration Vs. radius of step
    print("Plotting iter against radius of step")
    fig = plt.figure(figsize=(20, 20))
    y = np.array(data[:, 6], dtype=np.float64)
    plt.plot(x, y, "ro-")
    plt.xlabel("Total iteration count", fontsize=40)
    plt.xticks(fontsize=20)
    plt.ylabel("Radius of step", fontsize=40)
    plt.yticks(fontsize=20)
    plt.savefig("total_iter_vs_radius_of_step.png")
    plt.close()
    
    # plot 5 - total iteration Vs. trust radius
    print("Plotting iter against trust radius")
    fig = plt.figure(figsize=(20, 20))
    y = np.array(data[:, 7], dtype=np.float64)
    plt.plot(x, y, "ro-")
    plt.xlabel("Total iteration count", fontsize=40)
    plt.xticks(fontsize=20)
    plt.ylabel("Trust radius", fontsize=40)
    plt.yticks(fontsize=20)
    plt.savefig("total_iter_vs_trust_radius.png")
    plt.close()

    if headers is not None:

        # cut offs code will slice the min_trajectory to get all point less than these cutoffs and plot them
        cutoffs = [0.01, 0.001,  0.0005, 0.0001]

        # The convergence criteria to visually dipict how close we are
        convergence_criteria = gamess_default_gradient_tolerence

        df = pd.DataFrame(data=data, columns=headers)
        try:
            raw_plot_data = [float(ent) for ent in df[key].values]
        except Exception:
            log.warning("WARNING could not convert data to float some plots have not been completed")
            return

        log.info("Data: {}".format(raw_plot_data))

        trajectories = {}
        min_trajectory = []
        iterations = []
        cutoff_met = []
        cutoff_val = []

        for inx, x in enumerate(raw_plot_data):
            if inx == 0:
                minval = x

            if x < minval:
                minval = x

            min_trajectory.append(minval)
            iterations.append(inx)

        trajectories["min"] = min_trajectory

        for c in sorted(cutoffs):
            print("Cutoff to process: {}".format(c))
            try:
                tmp = [inx for inx, ent in enumerate(raw_plot_data) if ent < c][0]
            except IndexError:
                tmp = np.nan

            cutoff_met.append(tmp)
            cutoff_val.append(c)

        print("Plotting raw infomration")
        fig = plt.figure(figsize=(20,20))
        plt.plot(cutoff_val, cutoff_met, "ro-", label="cutoff Vs. iteration met")
        plt.legend(fontsize=20)
        plt.xlabel("Cutoff", fontsize=40)
        plt.xticks(fontsize=20)
        plt.ylabel("Iteration met", fontsize=40)
        plt.yticks(fontsize=20)
        plt.title("Cutoff Vs. iteration met", fontsize=40)
        plt.grid(True, which="both", ls="-", color='0.7')
        plt.savefig("cutoff_vs_iteration.png")
        plt.close()

        for c in cutoffs:
            trajectories[c] = [ent if ent < c else np.nan for ent in min_trajectory]

        print("Plotting data .....")
        for k, v in trajectories.items():

            print("Plotting information Iteration Vs. {}, cutoff {}".format(key, k))
            # Plot direct data no log
            x = [i for i in range(len(v))]
            fig = plt.figure(figsize=(20, 20))
            plt.plot(x, v, "ro-", label="Iteration Vs. {}".format(key))
            plt.plot(x, [convergence_criteria]*len(x), "ko-", label="convergence criteria")
            plt.legend(fontsize=20)
            plt.xlabel("Iteration", fontsize=40)
            plt.xticks(fontsize=20)
            plt.ylabel("{}".format(key), fontsize=40)
            plt.yticks(fontsize=20)
            plt.grid(True, which="both", ls="-", color='0.7')
            log.info("Value: {}".format(v))
            non_nans = [x for x in v if not np.isnan(x)]
            if len(non_nans) == 0:
                non_nans.append(np.nan)
            plt.title("Iteration Vs. {} (cutoff: {} lowest value: {})".format(key, k, min(non_nans)), fontsize=40)
            plt.savefig("iteration_vs_{}_cutoff_{}.png".format(key, k))
            plt.close()

            print("Plotting information Iteration Vs. log10({}), cutoff {}".format(key, k))
            # Plot log threshold
            x = [i for i in range(len(v))]
            logy = list(map(np.log10, v))
            fig = plt.figure(figsize=(20, 20))
            plt.plot(x, logy, "ro-", label="Iteration Vs. log10({})".format(key))
            plt.plot(x, [np.log10(convergence_criteria)]*len(x), "ko-", label="convergence criteria")
            plt.legend(fontsize=20)
            plt.xlabel("Iteration", fontsize=25)
            plt.xticks(fontsize=20)
            plt.ylabel("log10({})".format(key), fontsize=25)
            plt.yticks(fontsize=20)
            non_nans = [x for x in v if not np.isnan(x)]
            if len(non_nans) == 0:
                min_non_nans = np.nan
                log_non_nans = np.nan
            else:
                log_non_nans = "{:f}".format(np.log10(min(non_nans)))
                min_non_nans = min(non_nans)

            if isinstance(k, str):
                logval = k
            else:
                logval = "{:f}".format(np.log10(k))

            plt.title("Iteration Vs. {} (cutoff: {} (log10:{})  lowest value: {} (log10: {}))".format(key, k, logval, min_non_nans, log_non_nans), fontsize=25)
            plt.grid(True, which="both", ls="-", color='0.7')
            plt.savefig("iteration_vs_log10_{}_cutoff_{}.png".format(key, k))
            plt.close()

            print("Plotting information on a log scale {} Vs. Iteration, cutoff {}".format(key, k))
            # Plot log scale inverted axes
            x = [i for i in range(len(v))]
            fig = plt.figure(figsize=(20, 20))
            plt.plot(v, x, "ro-", label="Iteration Vs. {}".format(key))
            plt.plot([convergence_criteria]*len(x), x, "ko-", label="convergence criteria")
            plt.legend(fontsize=20)
            plt.xlabel("{}".format(key), fontsize=40)
            plt.xticks(fontsize=20)
            plt.xscale("log")
            plt.ylabel("Iteration", fontsize=40)
            plt.yticks(fontsize=20)
            plt.grid(True, which="both", ls="-", color='0.7')
            non_nans = [x for x in v if not np.isnan(x)]
            if len(non_nans) == 0:
                non_nans.append(np.nan)
            plt.title("{} Vs. Iteration (cutoff: {} lowest value: {})".format(key, k, min(non_nans)), fontsize=40)
            plt.savefig("iteration_vs_{}_log_scale_cutoff_{}.png".format(key, k))
            plt.close()

def get_convergence_data(files, line_regex="GRAD. MAX", radius_regex="RADIUS OF STEP TAKEN"):
    """
    A function to get the GAMESS convergence energy data from a GAMESS output file
    :param files: list - GAMESS output files to get convegence data out of 
    :param line_regex: string - The line regex to get the convergence and energetic information from
    :param radius_regex: string - The line regex to get the trust radius from the optimizer from
    """
    log = logging.getLogger(__name__)

    #trust_regex = "CURRENT TRUST RADIUS"
    # grep "GRAD. MAX" out.stdout*
    # out.stdout: NSERCH:   0  E=    -1473.7025045563  GRAD. MAX=  0.0007669  R.M.S.=  0.0002916

    log.info("Using regexs to get data from GAMESS outputfiles\n\tenergy regex: {}\n\ttrust radius reg: {}".format(line_regex,radius_regex))
    data = []
    dfound = False
    headers = ["filename", "total_search_iteration", "iteration", "energy", "max_gradient", "rms", "radius_of_step", "trust_radius"]
    count = 0
    for fin in files:
        with open(fin, "r") as f:
            for line in f:
                if line_regex in line:
                    ents = line.split()
                    d = [fin, count, int(ents[1]), float(ents[3]), float(ents[6]), float(ents[8])]
                    dfound = True
                
                if radius_regex in line and dfound is True:
                    ents = line.split()
                    d.append(float(ents[4]))
                    d.append(float(ents[8]))
                    count = count + 1
                    dfound = False
                    data.append(d)

    log.info("Writing data to file 'convergence_analysis_data.csv'")

    with open("convergence_analysis_data.csv", "w") as fout:
        fout.write("{}\n".format(",".join(headers)))
        for ent in data:
            fout.write("{}\n".format(",".join([str(elm) for elm in ent])))

    data = np.array(data)

    return data, headers

def get_decriptors(gamess_output_file):
    """
    A fuction to get useful energies and quantities out of a GAMESS output file and writes them to a new files named after the data type
    :param gamess_output_file: string - the output GAMESS file name and path  
    """

    log = logging.getLogger(__name__)
    # search for: COORDINATES OF ALL ATOMS ARE (ANGS) + 2 lines for coordinates stop at INTERNUCLEAR DISTANCES (ANGS.)
    #             EQUILIBRIUM GEOMETRY LOCATED Geometry optimzation successful
    #             TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS + 1 for charges and population data stop MULLIKEN SPHERICAL HARMONIC POPULATIONS
    #             (DEBYE) + 1 for dipole moment data dx dy dz D stop END OF PROPERTY EVALUATION

    get_features = False
    geometry = False
    population_and_charges = False
    dipole_moment = False
    bond_orders = False
    valence = False

    coords = []
    popu_charge = []
    dipole = []
    bonds_order = []
    valence_state = []

    geometry_start_re = "COORDINATES OF ALL ATOMS ARE (ANGS)"
    geometry_stop_re = "INTERNUCLEAR DISTANCES (ANGS.)"
    popu_charge_start_re = "TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS"
    popu_charge_stop_re = "MULLIKEN SPHERICAL HARMONIC POPULATIONS"
    dipole_start_re = "(DEBYE)" 
    dipole_stop_re = "END OF PROPERTY EVALUATION" 
    bond_order_start_re = "BOND ORDER AND VALENCE ANALYSIS"
    bond_order_stop_re = " TOTAL       BONDED        FREE" 
    valence_start_re = " ATOM            VALENCE     VALENCE     VALENCE"
    valence_stop_re = " ---------------------" 

    log.info("Attempting to assemble features from {}".format(gamess_output_file))
    with open(gamess_output_file) as fout:
        for inx, line in enumerate(fout):
            if "EQUILIBRIUM GEOMETRY LOCATED" in line:
                log.info("Geometry optimization completed successfully extracting features")
                get_features = True

            elif get_features is True and geometry_start_re in line:
                geometry = True
                geom_start = inx + 2

            elif get_features is True and geometry_stop_re in line:
                geometry = False
            
            elif get_features is True and popu_charge_start_re in line:
                population_and_charges = True
                popu_and_char_start = inx + 1

            elif get_features is True and popu_charge_stop_re in line:
                population_and_charges = False
            
            elif get_features is True and dipole_start_re in line:
                dipole_moment = True
                dipole_start = inx

            elif get_features is True and dipole_stop_re in line:
                dipole_moment = False
            
            elif get_features is True and bond_order_start_re in line:
                bond_orders = True
                bond_order_start = inx + 3

            elif get_features is True and bond_order_stop_re in line:
                bond_orders = False
            
            elif get_features is True and valence_start_re in line:
                valence = True
                valence_start = inx + 1

            elif get_features is True and valence_stop_re in line:
                valence = False


            else:
                pass

            if geometry is True and inx > geom_start:
                coords.append(line.strip())

            if population_and_charges is True and inx >= popu_and_char_start: 
                popu_charge.append(line.strip())

            if dipole_moment is True and inx >= dipole_start:
                dipole.append(line.strip())
            
            if bond_orders is True and inx > bond_order_start:
                bonds_order.append(line.strip())
            
            if valence is True and inx >= valence_start:
                valence_state.append(line.strip())

    if get_features is False:
        log.warning("WARNING Geometry optimzation unsuccessful no descriptors will be output for file {}".format(gamess_output_file))
        return

    log.debug("Coordinates:\n{}\n\nCharges and population analysis:\n{}\n\nDipole moment:\n{}\n\n,Bonds orders:\n{}\n\nValences:\n{}\n\n".format(
                                                                                                                                        coords, 
                                                                                                                                        popu_charge, 
                                                                                                                                        dipole,
                                                                                                                                        bonds_order,
                                                                                                                                        valence_state))
    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    try:
        log.info("Writing coordinate data to coordinates.xyz.\nWriting mass data to mass_data.csv.")
        with open("coordinates.xyz", "w") as fout, open("mass_data.csv", "w") as ffout:
            fout.write("{}\n".format(len(coords)))
            fout.write("Made by {} data from {}\n".format(__title__, gamess_output_file))
            ffout.write("Data from {}\n".format(gamess_output_file))
            ffout.write("atom, masses\n")
            for ent in coords:
                e = ent.split()
                fout.write("{} {} {} {}\n".format(e[0], e[2], e[3], e[4]))
                ffout.write("{},{}\n".format(e[0], e[1]))
    except IndexError as err:
        log.warning("Coordinates and mass features could not be completed due to {}".format(err))

    try:
        log.info("Writing population data to population_analysis.csv.")
        with open("population_analysis.csv", "w") as fout:
            # Keep column numbers consistent
            fout.write("Data_from,{},NA,NA,NA,NA\n".format(gamess_output_file))
            for inx, ent in enumerate(popu_charge):
                if inx == 0:
                    fout.write("ID,{}\n".format(",".join([elm for elm in re.split(r'\s{2,}', ent)])))
                else:
                    fout.write("{}\n".format(",".join([elm for elm in ent.split()])))
    except IndexError as err:
        log.warning("Population features could not be completed due to {}".format(err))

    try:
        log.info("Writing dipole moment data to dipole_moment.csv.")
        with open("dipole_moment.csv", "w") as fout:
            # keep column numbers consistent
            fout.write("Data_from, {},NA,NA\n".format(gamess_output_file))
            for inx, ent in enumerate(dipole):
                if inx == 0:
                    header = "DX,DY,DZ,D_DEBYE"
                    #fout.write("{}\n".format(",".join([elm for elm in re.split(r'\s{2,}', ent)])))
                    fout.write("{}\n".format(header))
                else:
                    fout.write("{}\n".format(",".join([elm for elm in ent.split()])))
    except IndexError as err:
        log.warning("Moment features could not be completed due to {}".format(err))
    
    try:
        log.info("Writing bond order data to bond_ord.csv.")
        with open("bond_ord.csv", "w") as fout:
            # keep column numbers consistent
            fout.write("Data_from, {},NA,NA\n".format(gamess_output_file))
            for inx, ent in enumerate(bonds_order):
                e = ent.split()
                if inx == 0:
                    header = "ATOM,PAIR,DIST,B.ORDER"
                    fout.write("{}\n".format(header))
                # Appears as up to three blocks of 4 entries
                try:
                    fout.write("{},{},{},{}\n".format(e[0].strip(), e[1].strip(), e[2].strip(), e[3].strip()))
                except IndexError as err:
                    log.error("WARNING index error trying to write bonds ord. {}".format(err))

                try:
                    fout.write("{},{},{},{}\n".format(e[4].strip(), e[5].strip(), e[6].strip(), e[7].strip()))
                except IndexError as err:
                    pass
                
                try:
                    fout.write("{},{},{},{}\n".format(e[8].strip(), e[9].strip(), e[10].strip(), e[11].strip()))
                except IndexError as err:
                    pass
    except IndexError as err:
        log.warning("Bond order features could not be completed due to {}".format(err))

    try:
        log.info("Writing valences data to valence.csv.\n----------- Features assembled -----------\n")
        with open("valence.csv", "w") as fout:
            # keep column numbers consistent
            fout.write("Data_from, {},NA,NA,NA\n".format(gamess_output_file))
            for inx, ent in enumerate(valence_state):
                if inx == 0:
                    header = "ID,TOTAL_VALENCE,BONDED_VALENCE,FREE_VALENCE"
                    #fout.write("{}\n".format(",".join([elm for elm in re.split(r'\s{2,}', ent)])))
                    fout.write("{}\n".format(header))
                fout.write("{}\n".format(",".join([elm for elm in ent.split()])))
    except IndexError as err:
        log.warning("Valence features could not be completed due to {}".format(err))

def main():
    """
    """
    
    try:
        usage = "python {} Required Parameters {}\n".format(__title__, " options .....")

        parser = argparse.ArgumentParser(description="Command line binary points script",
                                         usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        group = parser.add_argument_group("Principal arguments")
        group.add_argument("-d", "--directory", action="store", default=os.getcwd(), help="Directory to find all GAMESS output files (out.stdout*)")
        
        group = parser.add_argument_group("Optional arguments")
        group.add_argument("-f", "--file_regex", action="store", default="out.stdout*", help="Regex to find all GAMESS output files in directory")
        group.add_argument("-e", "--energy_regex", action="store", default="GRAD. MAX", help="Regex to energy gradient and rms in GAMESS output files")
        group.add_argument("-t", "--trust_radius_regex", action="store", default="RADIUS OF STEP TAKEN", help="Regex to find trust radii in GAMESS ouput")
        group.add_argument("-g", "--gradient_tolerance", action="store", default=0.0001, help="GAMESS gradient tolerance (OPTTOL)")
        group.add_argument("-n", "--no_plot", action="store_true", default=False, help="If this flag is given no plotting of the convergence data will "
                                                                                            "be done.")
        
        group = parser.add_argument_group("Logging arguments")
        group.add_argument("--loglev", action="store", default="INFO", help="log level")

        op = parser.parse_args()

    except argparse.ArgumentError as ee:
        print("\nERROR - command line arguments are ill defined please check the arguments\n")
        raise ee

    setup_logger(os.getcwd(), op.loglev)
    log = logging.getLogger(__name__)
    log.info("\nOrganiszation : IBM Research UK, Hartree Centre, Daresbury Laboratory\nCreated"
             "       : June 2017\nProgram       : {}\nVersion       : {}\n "
             "--------------------------------\n".format(__title__, __version__))

    log.info("Command line input =\n\t{}".format(op))
    
    log.info("Finding GAMESS ouput files in directory {} using search string {}".format(op.directory, op.file_regex))

    # Change directory to get the full absolute file paths
    sd = os.getcwd()
    os.chdir(op.directory)
    wd = os.getcwd()
    os.chdir(sd)
    
    files = list(filter(os.path.isfile,glob.glob(os.path.join(wd, op.file_regex))))
    files.sort(key=lambda x: os.path.getmtime(x)) 
    #files = sorted(glob.glob(os.path.join(wd, op.file_regex)))
    #if glob.glob(os.path.join(wd, "out.stdout")):
    #    files.append(os.path.join(wd, "out.stdout"))
    log.info("Files found using file match {}:\n{}".format(op.file_regex, files))
    
    # Get the descriptors from the last file i.e. most recently updated
    get_decriptors(files[-1])

    # Use the file list to get the information you want
    data, headers = get_convergence_data(files=files, line_regex=op.energy_regex, radius_regex=op.trust_radius_regex)
    
    if not op.no_plot:
        # plot the convergence data graphs
        plot_convergence(data, headers=headers, gamess_default_gradient_tolerence=op.gradient_tolerance)
    else:
        log.info("No plotting requested for convergence skipping convergence plots")

if __name__ == "__main__":
    main()
