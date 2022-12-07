
# Copyright IBM Inc. 2017, 2018, 2019, 2020. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
# Author(s):
#   James McDonagh
#   Hsiang Han Hsu


import os
import sys
import re
import yaml
import logging
import shutil
import glob
import logging
import re
import numpy as np
import experiment.model.codes
from datetime import datetime

__title__ = os.path.basename(__file__)
__version__ = 1.0
__copyright__ = "Copyright IBM Corp. 2017, 2018, 2019, 2020"


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


def get_consistent_orbs_and_geom(filein, search_start_term_geom="COORDINATES OF SYMMETRY UNIQUE ATOMS",
                                 search_end_term_geom="GENERATED AT", offset_lines=3,
                                 search_start_term_orbitals="$VEC", search_end_term_orbitals=("$GRAD", "POPULA"),
                                 offset_end_lines=-1):
    """
    function to get the complete geometry and obital data for a complete optimization run
    :param filein:
    :return:
    """

    log = logging.getLogger(__name__)

    # Geometry search lines and line offset
    # search_start_term_geom = "COORDINATES OF SYMMETRY UNIQUE ATOMS
    # search_end_term_geom = "GENERATED AT"
    # offset_lines = 3

    # Orital search lines and offset
    # search_start_term_orbitals = "$VEC"
    # search_end_term_orbitals = "$END"

    # offset_end_lines = 1  # minus one line number to get the last molecular orbital

    # get reversed lists
    log.info("Parsing {}".format(filein))
    start_geometry = False
    stop_geometry = False
    start_orbital = False
    stop_orbital = True

    start_geometry = False
    stop_geometry = False
    start_orbital = False
    stop_orbital = False

    start_indexes_geom = []
    end_indexes_geom = []
    start_indexes_orbitals = []
    end_indexes_orbitals = []

    # Loop over the lines and store the geometry and orbitals only if they follow in the expected order
    # This is a bit clumbersome but it seems to work. The loop checks that all lines are found in the correct order
    # or expected order if an unexpected order is found then the loop deletes the most recently added line indexes
    # and attempts to reset
    with open(filein, "r") as fin:
        for inx, l in enumerate(fin):

            if search_start_term_geom in l:  # and stop_orbital is True:
                if start_geometry is False:
                    start_indexes_geom.append(inx + offset_lines)
                    start_geometry = True
                else:
                    log.info("unexpected start geometry found {} {} {} {} {}".format(inx,
                                                                                     start_geometry,
                                                                                     stop_geometry,
                                                                                     start_orbital,
                                                                                     stop_orbital))
                    if start_geometry is True:
                        del start_indexes_geom[-1]

                    if stop_geometry is True:
                        del end_indexes_geom[-1]

                    if start_orbital is True:
                        del start_indexes_orbitals[-1]

                    if stop_orbital is True:
                        del end_indexes_orbitals[-1]

                    start_geometry = False
                    stop_geometry = False
                    start_orbital = False
                    stop_orbital = False

            if search_end_term_geom in l:  # and start_geometry is True:
                if start_geometry is True:
                    end_indexes_geom.append(inx)
                    stop_geometry = True
                else:
                    log.info("unexpected stop geometry found {} {} {} {} {}".format(inx,
                                                                                    start_geometry,
                                                                                    stop_geometry,
                                                                                    start_orbital,
                                                                                    stop_orbital))

                    if start_geometry is True:
                        del start_indexes_geom[-1]

                    if stop_geometry is True:
                        del end_indexes_geom[-1]

                    if start_orbital is True:
                        del start_indexes_orbitals[-1]

                    if stop_orbital is True:
                        del end_indexes_orbitals[-1]

                    start_geometry = False
                    stop_geometry = False
                    start_orbital = False
                    stop_orbital = False

            if search_start_term_orbitals in l:  # and stop_geometry is True:
                if all(x is True for x in [start_geometry, stop_geometry]):
                    start_indexes_orbitals.append(inx)
                    start_orbital = True
                else:
                    log.info("unexpected start orbital found {} {} {} {} {}".format(inx,
                                                                                    start_geometry,
                                                                                    stop_geometry,
                                                                                    start_orbital,
                                                                                    stop_orbital))

                    if start_geometry is True:
                        del start_indexes_geom[-1]

                    if stop_geometry is True:
                        del end_indexes_geom[-1]

                    if start_orbital is True:
                        del start_indexes_orbitals[-1]

                    if stop_orbital is True:
                        del end_indexes_orbitals[-1]

                    start_geometry = False
                    stop_geometry = False
                    start_orbital = False
                    stop_orbital = False

            if any(term in l for term in search_end_term_orbitals):  # and start_orbital is True:
                if all(x is True for x in [start_geometry, stop_geometry, start_orbital]):
                    end_indexes_orbitals.append(inx + offset_end_lines)
                    stop_orbital = True
                else:
                    log.info("unexpected stop orbital found {} {} {} {} {}".format(inx,
                                                                                   start_geometry,
                                                                                   stop_geometry,
                                                                                   start_orbital,
                                                                                   stop_orbital))

                    if start_geometry is True:
                        del start_indexes_geom[-1]

                    if stop_geometry is True:
                        del end_indexes_geom[-1]

                    if start_orbital is True:
                        del start_indexes_orbitals[-1]

                    if stop_orbital is True:
                        del end_indexes_orbitals[-1]

                    start_geometry = False
                    stop_geometry = False
                    start_orbital = False
                    stop_orbital = False

            if all(x is True for x in [start_geometry, stop_geometry, start_orbital, stop_orbital]):
                log.debug("Complete geometry and orbitals found")
                start_geometry = False
                stop_geometry = False
                start_orbital = False
                stop_orbital = False

    start_indexes_geom.reverse()
    end_indexes_geom.reverse()
    start_indexes_orbitals.reverse()
    end_indexes_orbitals.reverse()

    if not all(len(ent) == len(start_indexes_geom) for ent in [start_indexes_geom,
                                                               end_indexes_geom,
                                                               start_indexes_orbitals,
                                                               end_indexes_orbitals]):
        log.warning("WARNING - the lengths of the indexes are different to find start and end of geometry and orbitals")

    log.info("Length start geometry {}, length stop geometry {}, length start orbital {}, "
             "length stop orbital {}".format(len(start_indexes_geom),
                                             len(end_indexes_geom),
                                             len(start_indexes_orbitals),
                                             len(end_indexes_orbitals)))

    log.info("Indexes:\nstart geometry: {}\nstop geometry: {}\nstart orbitals: {}\nstop "
             "orbitals: {}".format(start_indexes_geom,
                                   end_indexes_geom,
                                   start_indexes_orbitals,
                                   end_indexes_orbitals))

    i = 0
    fail = True

    while fail is True:

        if i >= len(start_indexes_geom):
            log.error("No valid geometry and orbital sets can be found")
            return False, False

        start_geom = start_indexes_geom[i]
        stop_geom = end_indexes_geom[i]
        start_orbit = start_indexes_orbitals[i]
        stop_orbit = end_indexes_orbitals[i]
        i = i + 1

        # geometry and orbitals should be separated by 3 lines check for this
        diff = start_orbit - stop_geom
        if diff != 3:
            log.warning("Geometry and orbitals should be separated by 3 lines but are not for "
                        "index {} diff = {}".format(i - 1, diff))
            fail = True

        else:
            log.info("Geometry and orbitals separated by 3 lines as expected")
            fail = False

        # check stop_geom and start_obital are with the range start_geom to stop_orbital
        if start_geom < stop_geom < stop_orbit:
            if start_geom < start_orbit < stop_orbit:
                log.info("Intermediate lines fall in the correct range as expected")
                fail = False
            else:
                log.info("start_orbit line is not in the  correct range")
                fail = True
        else:
            log.info("stop_geom line is not in the  correct range")
            fail = True

        # check stop_orbit is larger than stop_orbital[i]
        if i < len(end_indexes_orbitals):
            if stop_orbit > end_indexes_orbitals[i]:
                log.info("stop_orbit in correct order as expected")
                fail = False
            else:
                log.warning("stop_orbit is greater than the next stop_orbit")
                fail = True
        else:
            log.info("Appears there is no more geometry and orbitals to test against")

    with open(filein) as fin:
        log.info("Successfully found geometry and orbitals\n")
        geometry = [l.strip() for ind, l in enumerate(fin) if start_geom <= ind < stop_geom]
        log.info("Geometry:\n{}\n.....\n{}\nLength {}\n".format(geometry[0], geometry[-1], len(geometry)))
        fin.seek(0)

        orbitals = [l for ind, l in enumerate(fin) if start_orbit <= ind < stop_orbit]
        log.info("Orbitals:\n{}\n.....\n{}\nLength {}\n".format(orbitals[0], orbitals[-1], len(orbitals)))

    log.info("Geometry and Orbitals found? {}".format(not fail))

    return geometry, orbitals


def get_template_data(filename, regex="$DATA"):
    """
    A function to get template data from a template file or input file
    :param filename: The path and filename for a template or GAMESS input file
    """
    log = logging.getLogger(__name__)
    try:
        with open(filename, "r") as fin:
            raw_data = [line.strip() for line in fin]
    except FileNotFoundError:
        log.error("ERROR - file not found for template or input file {}".format(filename))
        return False

    # look for $DATA in the raw_data if found assume the file is an input file
    # if not found assume it is a template file

    data_index = None
    for inx, ent in enumerate(raw_data):
        if regex in ent:
            data_index = inx

    if data_index is None:
        log.info("\nNo '{}' found assuming this is a template file and using all lines\n".format(regex))
        data = [elm for elm in raw_data]
    else:
        log.info("\n'{}' found assuming this is a input file and using only lines up to {} as method\n".format(regex,
                                                                                                               data_index))
        log.debug(raw_data[0:data_index])
        data = [elm for elm in raw_data[0:data_index] if elm]

    log.debug("Template/input data: {}".format(data))

    if len(data) <= 0:
        log.error("ERROR input or template file {} is empty".format(filename))
        return False

    return data


def count_molecular_orbitals(mol_orbs, gamess_out_file):
    """
    Function to count the number of molecular orbitals in a GAMESS VEC section in a *.dat file
    :param mol_orbs: list of molecular orbitals from a *.dat file
    :param gamess_out_file: the output file to find number of functions from
    for reference https://ask.cyberinfrastructure.org/t/how-can-i-restart-my-gamess-calculation-and-resubmit/259
    for reference http://erg.biophys.msu.ru/wordpress/archives/581
    for reference http://classic.chem.msu.su/cgi-bin/ceilidh.exe/gran/gamess/forum/?C35e9ea938bHW-7937-1296+00.htm
    """
    log = logging.getLogger(__name__)
    log.info(gamess_out_file)

    try:
        previous = mol_orbs[0].split()[0].strip()
        count = 1
    except IndexError as e:
        log.error("ERROR - index error when initizalizing molecular orbital counting")
        # raise e
        return False

    try:
        for inx, line in enumerate(mol_orbs):
            n = line.split()[0].strip()
            if n != previous:
                count = count + 1
                previous = n
    except IndexError as e:
        log.error("ERROR - index error counting orbitals")
        # raise e
        return False

    log.info("Counted molecular orbitals number of molecular orbitals is in dat file {}".format(count))

    additional_search = False
    n_shells = None
    n_basis_functions = None
    n_molecular_orbitals = None
    try:
        with open(gamess_out_file, "r") as fin:
            for line in fin:
                if "TOTAL NUMBER OF BASIS SET SHELLS" in line:
                    n_shells = line.split("=")[1].strip()
                elif "TOTAL NUMBER OF BASIS FUNCTIONS" in line and n_basis_functions is None:
                    n_basis_functions = line.split("=")[1].strip()
                elif "NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS" in line and n_basis_functions is None:
                    n_basis_functions = line.split("=")[1].strip()
                    additional_search = True
                elif additional_search is True and "TOTAL NUMBER OF MOS IN VARIATION SPACE" in line:
                    additional_search = False
                    n_molecular_orbitals = line.split("=")[1].strip()

    except FileNotFoundError:
        log.error("GAMESS output file {} not found".format(gamess_out_file))
        return False
    except IndexError:
        log.error("Index error parsing basis and orbital infomration")
        return False

    log.info("GAMESS calculation had {} shells, {} basis functions, {} Molecular orbitals in the \
            variational space".format(n_shells, n_basis_functions, n_molecular_orbitals))

    # TODO: These should match but seem to frequently not in DFT need to understand why.
    if count != n_basis_functions:
        log.warning("required {} basis functions counted {} in $VEC".format(n_basis_functions, count))

    # NOTE: based on http://erg.biophys.msu.ru/wordpress/archives/581 GAMESS need only read the basis functions
    # of the number of shells to restart. For simplicity you can copy all function and us norb to tell GAMESS
    # to only read the shells it needs.

    return n_shells


def write_file(filename, geometry=None, template_data=None, molecular_orbitals=None, norbs=0, restart=None):
    """
    A function to write a new file with the last geometry and moleclar orbitals used
    example file at https://www2.chemistry.msu.edu/courses/cem888/gamess-instr.html
    :param filename: The filename to write to
    :param geometry: Last geometry tried by the optimization
    :param template_data: template input file data
    """
    log = logging.getLogger(__name__)
    try:
        date_and_time = datetime.now()
        timestamp = date_and_time.strftime("%m_%d_%Y-%H_%M_%S")
        with open(filename, "r") as fin, open("{}_{}".format(filename, timestamp), "w") as fout:
            original_input = [line for line in fin]
            for line in original_input:
                fout.write("{}".format(line))
    except IndexError:
        log.warning("WARNING - an error occured when trying to make a backup of the old input file. \
        The file will be overwritten")

    try:
        paths = os.getcwd().split(os.path.sep)
        name = [ent for ent in paths if ".instance" in ent][0]
    except IndexError:
        name = None

    with open(filename, "w") as fout:
        edit = False
        for line in template_data:
            if "GUESS" in line:
                log.info("guess line found '{}'".format(line))
                fout.write(" $GUESS GUESS=MOREAD NORB={} $END\n".format(norbs))

            elif "CONTRL" in line or edit is True:
                log.info("control line found '{}'".format(line))
                # commonly CONTRL is over two lines if so write the line as
                # is and make edit True so that we apply the change to the next line
                if "END" in line:
                    edit = False
                    if "QMTTOL" not in line:
                        fout.write(" {} QMTTOL=1.0E-05 $END\n".format(line[:-4]))
                    else:
                        fout.write(" {}\n".format(line.strip()))
                else:
                    fout.write(" {}\n".format(line.strip()))
                    edit = True

            elif "STATPT" in line and "DXMAX" not in line:
                log.info("stat line found '{}'".format(line))
                fout.write(" {} DXMAX=0.05 $END\n".format(line[:-4]))

            else:
                fout.write(" {}\n".format(line.strip()))

        fout.write("\n")
        fout.write(" $DATA\n")
        fout.write("Restart {} {} instance {}\n".format(restart, timestamp, name))
        fout.write(" C1\n")
        for line in geometry:
            fout.write(" {}\n".format(line.strip()))
        fout.write(" $END\n")
        for line in molecular_orbitals:
            fout.write("{}".format(line))
        fout.write(" $END\n")


def get_vec_data(filename):
    """
    A function to get the molecular orbitals for the vec section of the input from a *.dat file
    :param filename: the *.dat filename
    """
    log = logging.getLogger(__name__)
    search_start_term = "VEC"
    search_end_term = "END"
    # search_end_term = " POPULATION ANALYSIS"
    offset_end_lines = 1  # minus one line number to get the last molecular orbital

    with open(filename, "r") as fin:
        try:
            start = max([inx for inx, l in enumerate(fin) if search_start_term in l])
            log.info("Start {}".format(start))
        except ValueError:
            log.error("ERROR - value error evaluating max of the line number matching {} in get vec data".format(
                search_start_term))
            return False
        fin.seek(0)

        try:
            end = min([inx for inx, l in enumerate(fin) if search_end_term in l and inx > start])
            log.info("Stop {}".format(end))
        except ValueError:
            log.error("ERROR - value error evaluating max of the line number matching {} in get vec data".format(
                search_end_term))
            return False
        fin.seek(0)

        log.info("Start line {} Stop line {} both found for $VEC data".format(start, end))
        if end <= start:
            log.error("ERROR - $VEC data - start and stop lines are the same or stop is less than start. \
                    Check {}\nstart regex is {}\nstop regex is {}\n.".format(filename, search_start_term,
                                                                             search_end_term))
            return False

        # Note: We want to preserve the line formatting of the entries as it should be correct to use as input
        mol_orbs = [ent for inx, ent in enumerate(fin) if start <= inx < end]

    log.debug("\nMolecular orbital:\n{}".format(mol_orbs))

    if len(mol_orbs) <= 0:
        log.error("ERROR molecular orbitals not found in {}".format(filename))
        return False

    return mol_orbs


def determine_gamess_exit_reason(gamess_out_file, component_type):
    """
    Afunction to attempt to determine the exit reason of GAMESS from the output file
    :param gamess_out_file: the standard GAMESS output file
    :return:
    """
    log = logging.getLogger(__name__)

    run_type = None
    eq_geom = False
    term_norm = False
    ex_grace = False
    stat_point_not_found = False
    ddikick_failed = False
    term_abnorm = False
    orbital_input_error = False
    lin_dep_warning = False

    with open(gamess_out_file, "r") as fin:
        for line in fin:

            # Determine the type of calculation
            if "RUNTYP" in line:
                list_line_contents = line.split()
                try:
                    run_type = [elm.split("=")[1].strip() for elm in list_line_contents if "RUNTYP" in elm][0]
                    # run_type = next((elm.split("=")[1].strip() for elm in list_line_contents if "RUNTYP" in elm))
                except IndexError:
                    log.error("Could not determine run type from GAMESS assuming component name defined it")
                    if re.match(u"^[A-Za-z0-9]*GeometryOptimisation[A-Z,a-z0-9]*", component_type):
                        run_type = "OPTIMIZE"
                    elif re.match(u"^[A-Za-z0-9]*Energy[A-Za-z0-9]*", component_type):
                        run_type = "ENERGY"
                if run_type is None:
                    log.error("Run type unknown setting to optimize")
                    run_type = "OPTIMIZE"

                log.info("Run type = {}".format(run_type))

            # This horrible set of statements defines which key lines have been found and which have not
            if "EQUILIBRIUM GEOMETRY LOCATED" in line and run_type == "OPTIMIZE":
                eq_geom = True

            elif "EXECUTION OF GAMESS TERMINATED NORMALLY" in line:
                term_norm = True

            elif "exited gracefully" in line:
                ex_grace = True

            elif "FAILURE TO LOCATE STATIONARY POINT" in line:
                stat_point_not_found = True

            elif "ddikick.x: Sending kill signal to DDI" in line:
                ddikick_failed = True

            elif "terminated due to error(s)" in line:
                term_abnorm = True

            elif "ERROR: PREMATURE END OF ORBITAL INPUT ENCOUNTERED" in line:
                orbital_input_error = True

            elif "POSSIBLE LINEAR DEPENDENCE PROBLEMS DETECTED, INPUT QMTTOL" in line:
                lin_dep_warning = True

    # Collect fail reason results together
    if run_type == "ENERGY":
        fails = any(ent is True for ent in [ddikick_failed, term_abnorm])
        reasons = "ddick sent kill: {} terminated on error: {}".format(ddikick_failed, term_abnorm)
        log.error("GAMESS has failed for the following reasons {} but can be restarted".format(reasons))
    elif run_type == "OPTIMIZE":
        fails = any(ent is True for ent in [ddikick_failed, term_abnorm, stat_point_not_found])
        reasons = "ddick sent kill: {} terminated on error: {} stationary point located: {} equilibrium point located {}".format(
            ddikick_failed, term_abnorm, stat_point_not_found, eq_geom)
        log.error("GAMESS has failed for the following reasons {} but can be restarted".format(reasons))
    else:
        fails = False

    # using the line found and not found allows us to define
    # the restart should happen or not
    # TODO: This will need adding to to enable better
    if eq_geom and term_norm  and run_type == "OPTIMIZE" and not fails:
        log.info("Calculation determined to be geometry optimization and equilibrium has been reached - SUCCESS")
        return False

    elif term_norm and run_type == "ENERGY" and not fails:
        log.info("Calculation determined to be energy and has completed - SUCCESS")
        return False

    elif fails is True:
        return True

    else:
        log.error("Failed GAMESS run restart will be attempted, but analysis of failure reason cannot be completed "
                  "as the output file is likely incomplete")
        return True


def Restart(workingDirectory, restarts, componentName, log, exitReason, exitCode):
    """
    Flags if a component can be restarted after encountering exitReason
     Also should perform any operations required for the restart to work.
     Examples of common requirements
     - update input files options
     - rename output files of previous run to prevent overwriting (the `restarts` parameter is useful for this)
   Any files created/renamed/moved must be in the workingDir or subdirs
    Parameters:
        workingDirectory: Workflow directory of the component to be restarted
        restarts: The number of times this function has been called for this component
        componentName: The label the workflow engine uses to id this component
        log: A logiging.logger object - use to write output messages
        exitReason: A string. Defines why the components program exited.
                  One of the values in the experiment.model.codes.exitReasons dictionary
        exitCode: A integer. The exit-code returned by the components program
    Returns:
        One of experiment.model.codes.restartContext
        - RestartContextRestartPossible
        - RestartContextHookNotAvailable
        - RestartContextRestartNotRequired
        - RestartContextRestartNotPossible
        - RestartContextHookFailed
        - RestartContextRestartConditionsNotMet
    """
    log.info("Restart:\nworking dir: {}\nNumber of restarts: {}\n"
             "Component name: {}\nLog: {}\nExit reason {}\nExit code {}".format(workingDirectory,
                                                                                restarts, componentName, log,
                                                                                exitReason, exitCode))

    ##############
    restart_threshold = 5  # Will not restart more than this number of times for a single component
    ##############
    delete_old = True  # Deletes all old file except 'excluded_file_types variable entries' from the GAMESS directory
    ##############

    # create a string with out any component number for easier comparison operations
    componentType = re.sub("\d", "", componentName)
    # restart for hit wall time for all simulation compnents

    # variable for the restart threshold
    # TODO: HACK this should be an option to the restart method talk to MJ
    sep = os.path.sep
    split_wd = workingDirectory.split(os.path.sep)
    log.info(split_wd)
    log.info("Default restart threshold {}".format(restart_threshold))
    try:
        # Workdir format is $INSTANCE_DIR_NAME/stages/stage<idx>/<componentname>
        path_instance = os.path.sep.join(split_wd[:-3])
        path_metadata = os.path.join(path_instance, "elaunch.yaml")
        log.info(f"Path to elaunch.yaml file is {path_metadata}")
        with open(path_metadata, "r") as f:
            metadata = yaml.load(f, Loader=yaml.FullLoader)
        variables = metadata.get('variables', {})
        global_vars = variables.get('global', {})
        restart_threshold = int(global_vars.get('restart_threshold', restart_threshold))
    except Exception as err:
        log.warning(f"Unable to extract value of \"restart_threshold\" variable in elaunch.yaml due "
                    f"to {err} - will use value {restart_threshold}")

    if restarts >= restart_threshold:
        log.error("Number of restarts {} is great than or equal to restart threshold {}. \
                This component will now fail".format(restarts, restart_threshold))
        return experiment.model.codes.restartContexts['RestartContextRestartNotPossible']

    # NOTE: the above check must happen first other wise when an experiment is restarted after
    # completing the permitted number the hook will not allow it to continue.
    subdirs = [d for d in os.listdir(workingDirectory) if
               os.path.isdir(os.path.join(workingDirectory, d)) and re.match(r"Run\d*", d)]
    # NOTE: assumes restart_threshold is a fixed number !
    if len(subdirs) >= restarts:
        restarts = 1 + len(subdirs)

    # Check restart method is appropiate i.e. is it a GAMESS component
    retgopt = re.match(u"^[A-Za-z0-9]*GeometryOptimisation[A-Z,a-z0-9]*", componentType)
    retenr = re.match(u"^[A-Za-z0-9]*Energy[A-Z,a-z0-9]*", componentType)

    if retgopt:
        log.info("Component type Geometry optimization")
    elif retenr:
        log.info("Component type Energy")
    else:
        log.info("No restart hook for component")
        return experiment.model.codes.restartContexts['RestartContextHookNotAvailable']

    if retgopt or retenr:

        # make a backup sub-directory
        dirNam = "Run{}".format(restarts)
        path = os.path.join(workingDirectory, dirNam)
        if not os.path.exists(path):
            os.makedirs(os.path.join(workingDirectory, dirNam))
        else:
            log.warning("Path {} already exists, exiting restart and closing experiment".format(path))
            return experiment.model.codes.restartContexts['RestartContextRestartNotPossible']

        # copy all content files to the backup directory
        files = [f for f in os.listdir(workingDirectory) if os.path.isfile(os.path.join(workingDirectory, f))]
        for f in files:
            cpath = os.path.join(workingDirectory, f)
            shutil.copy2(cpath, path)
        log.info("A backup of the {} run is now stored in {}".format(restarts, path))

        # edit input files for restart methods
        gamess_input = glob.glob(os.path.join(workingDirectory, "*.inp"))
        gamess_dat_punch_file = glob.glob(os.path.join(workingDirectory, "*.dat"))
        # TODO: specific name is there a better way?
        gamess_out_file = glob.glob(os.path.join(workingDirectory, "out.stdout"))
        if len(gamess_input) == 0:
            log.error("Restart cannot find the appropriate gamess input (*.inp) file in {}".format(workingDirectory))
            # raise RuntimeError("ERROR - cannot find GAMESS input file *.inp")
            return experiment.model.codes.restartContexts['RestartContextRestartNotPossible']
        elif len(gamess_input) > 1:
            log.error("Restart has found more than one gamess input (*.inp) file in {}\nFiles found {}".format(
                workingDirectory, gamess_input))
            # raise RuntimeError("ERROR- found more than one GAMESS input file *.inp")
            return experiment.model.codes.restartContexts['RestartContextRestartNotPossible']
        else:
            gamess_input = gamess_input[0]

        found_dat = True
        found_out = True

        if len(gamess_dat_punch_file) == 0:
            log.error("Restart cannot find the appropriate gamess punch (*.dat) file in {}".format(workingDirectory))
            # raise RuntimeError("ERROR - cannot find GAMESS punch file *.dat")
            found_dat = False
        elif len(gamess_dat_punch_file) > 1:
            log.error("Restart has found more than one gamess punch (*.dat) file in {}\nFiles found {}".format(
                workingDirectory, gamess_dat_punch_file))
            # raise RuntimeError("ERROR found more than one GAMESS punch file *.dat")
            return experiment.model.codes.restartContexts['RestartContextRestartNotPossible']
        else:
            gamess_dat_punch_file = gamess_dat_punch_file[0]

        if len(gamess_out_file) == 0:
            log.error(
                "Restart cannot find the appropriate gamess out file (out.stdout) file in {}".format(workingDirectory))
            found_out = False
            # raise RuntimeError("ERROR - cannot find GAMESS punch file *.dat")
        elif len(gamess_out_file) > 1:
            log.error("Restart has found more than one gamess out file (out.stdout) file in {}\nFiles found {}".format(
                workingDirectory, gamess_dat_punch_file))
            # raise RuntimeError("ERROR found more than one GAMESS punch file *.dat")
            return experiment.model.codes.restartContexts['RestartContextRestartNotPossible']
        else:
            gamess_out_file = gamess_out_file[0]

        if not found_dat and not found_out:
            # Looks like this never ran just run as normal
            log.warning("Looks like component %s never ran - can rerun without any additional prep" % componentName)
            return experiment.model.codes.restartContexts["RestartContextRestartPossible"]

        # TODO: Add specific error message check and make variations
        # this is harder to do as the documentation is very poor so it is hard to know what to do
        # with specific errors outside of the general cases I have already dealt with.
        # the most common case seems to require changing QMTTOL this is done by default at the
        # momement.
        # determine if restart is needed as this hook is called even on supposedly successful GAMESS components
        restart_needed = determine_gamess_exit_reason(gamess_out_file, componentName)
        restart_context = experiment.model.codes.restartContexts[
            'RestartContextRestartPossible'] if restart_needed is True else experiment.model.codes.restartContexts[
            'RestartContextRestartNotRequired']

        log.info("GAMESS Input: {} GAMESS punch file {} GAMESS out file {}".format(gamess_input, gamess_dat_punch_file,
                                                                                   gamess_out_file))

        if restart_context == experiment.model.codes.restartContexts['RestartContextRestartPossible']:

            log.info("Restart is determined to be required for GAMESS")

            geometry, molecular_orbitals = get_consistent_orbs_and_geom(gamess_dat_punch_file)

            # mol_orbs = get_vec_data(gamess_dat_punch_file)

            if geometry is False:
                return experiment.model.codes.restartContexts['RestartContextRestartNotPossible']

            if molecular_orbitals is False:
                return experiment.model.codes.restartContexts['RestartContextRestartNotPossible']

            template_data = get_template_data(gamess_input)
            if template_data is False:
                return experiment.model.codes.restartContexts['RestartContextRestartNotPossible']

            norbs = count_molecular_orbitals(molecular_orbitals, gamess_out_file)
            if norbs is False:
                return experiment.model.codes.restartContexts['RestartContextRestartNotPossible']

            write_file(gamess_input, geometry=geometry, template_data=template_data,
                       molecular_orbitals=molecular_orbitals, norbs=norbs, restart=restarts)

            if delete_old is True:
                all_fs = [os.path.join(workingDirectory, f) for f in os.listdir(workingDirectory) if
                          os.path.isfile(os.path.join(workingDirectory, f))]
                excluded_file_types = [".inp", ".stdout", ".out", ".log", ".stderr", ".err", ".py"]
                for ent in excluded_file_types:
                    all_fs = [elm for elm in all_fs if ent not in elm]

                for ent in all_fs:
                    os.remove(ent)

            return restart_context

        else:
            log.info("Restart is not required for GAMESS")
            return restart_context


if __name__ == "__main__":
    # For testing
    # Format routine the message comes from, the leve of the information debug, info, warning, error, critical
    # writes all levels to teh log file Optimizer.log
    pathlog = os.path.join(os.getcwd(), 'test.log')
    logging.basicConfig(format='%(message)s', filemode='w', filename=pathlog, level=logging.INFO)

    # Setup handle for screen printing only prints info and above not debugging info
    screen = logging.StreamHandler()
    screen.setLevel(logging.INFO)
    logformat = logging.Formatter('%(levelname)s %(message)s')
    logging.getLogger("").addHandler(screen)
    log = logging.getLogger(__name__)

    log.info('Started {}\n'.format(datetime.now()))
    Restart(os.getcwd(), 4, "GeometryOptimisation", log, "ResourceExhausted", 1)
