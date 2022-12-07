#!/usr/bin/env python

# Copyright IBM Inc. 2015, 2019. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
# Author(s):
#   James McDonagh
#   Hsiang Han Hsu
#   Michael Johnston
#   Vassilis Vassiliadis


import os, re, csv, sys, glob
import numpy as np
import logging
from datetime import datetime

__organisation__ ="IBM Research, 2019"
__title__ = "csv2inp.py"
__copyright__ = "Copyright IBM 2019"

# VV: This causes terminals to clear when grabbing the logs of the pod that's executing this command
# os.system('clear')

# ===== Self-defined functions =====
# num2str
def num2str(num, precision = 0):
    return "%0.*f" % (precision, num)

# Extract text in a string
def catchtxt(string, starttxt, endtxt, errmsg=''):
    len_start = len(starttxt)
    ind1 = string.find(starttxt)
    if ind1 != -1:
        string2 = string[(ind1+len_start):]
        ind2 = string2.find(endtxt)
        outputstr = string2[0:ind2]
    else:
        outputstr = errmsg
    
    return outputstr

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


Ha2eV = 27.21138602  # 1 Ha = 27.21138602 eV
D2au = 0.393430307  # 1 Debye = 0.393430307 a.u.
svinp = 0  # Save NG inp files to "./not_finished" folder
pth = './'  # Default inp files location

if __name__ == "__main__":

    import optparse

    usage = "usage: %prog [options] [output folders]"

    parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)

    parser.add_option("--logLevel", dest="logLevel",
                      help="The level of logging. Default %default",
                      type="int",
                      default=30,
                      metavar="LOGGING")
    parser.add_option("-l", "--labels", dest="labels",
            help="A comma separated list of strings to use as labels for each file passed to this program."
                 " Must be the same number of labels as files. Default: %default ",
                      default="molecule",
                      metavar="LABEL")
    parser.add_option("-s", "--suffix", dest="suffix",
            help="Suffix of the files containing GAMESS standard out: Default % default",
                      default="out",
                      metavar="SUFFIX")
    parser.add_option("-o", "--output", dest="output",
            help="Output file name. Default %default",
                      default="energies.csv",
                      metavar="OUTPUT")

    options, args = parser.parse_args()

    regex_float = '-*\d+\.\d+'
    setup_logger(os.getcwd())
    log = logging.getLogger(__name__)
    log.info("option {} argument {}".format(options, args))

    fcsv = open(options.output, 'w')
    header = ['label', 'completed', 'total-energy', 'homo', 'lumo', 'gap',
              'electric-moments', 'total-time', 'total-time-per-core']

    writer = csv.writer(fcsv, sys.stdout, lineterminator='\n', delimiter=',')
    writer.writerow(header)
    data = []

    lumo_ha = 0
    homo_ha = 0
    E_gap = 0
    data_wct = 0

    # VV: use regex to identify number of ranks and cores, look for line:
    # MPI kickoff will start GAMESS on 1 cores in 1 nodes.

    regex_cores_nodes = re.compile(r"MPI kickoff will start GAMESS on ([0-9]+) cores in ([0-9]+) nodes\.")

    for (f_out,f2) in zip(args,options.labels.split(',')):
        
        g = os.path.join(f_out, "*%s" % options.suffix)
        f_out = glob.glob(g)
        if len(f_out) == 0:
            log.warning('WARNING - No match to %s - skipping this output dir' % g) 
            continue
        else:
            log.info("Files matching regex({}): {}".format(g, f_out))
        
        f_out = f_out[0]

        # ===== Open *.out =====
        try:
            log.info("Extracting information from GAMESS output file: {}".format(f_out))
            fopen = open(f_out)
            txt = fopen.read()
            fopen.close()
        except:
            txt = 'h'*10

        # ===== File name =====
        data_fn = [f2]

        # ===== OK, NG =====
        norm_term = False
        converged = False

        if txt.find('TERMINATED NORMALLY') != -1 or txt.find('gracefully') != -1:
            norm_term = True
        # If not found i.e. == -1 and is an energy calculation
        if txt.find("SCF IS UNCONVERGED") == -1 and txt.find("RUNTYP=ENERGY") != -1:
            converged = True
        # If it is found i.e. equilibrium is located and optimized
        if txt.find("EQUILIBRIUM GEOMETRY LOCATED") != -1 and txt.find("RUNTYP=OPTIMIZE") != -1:
            converged = True

        if converged is False:
            log.error("ERROR - unconverged scf in this calculation")
            ON = 'NG'
        elif norm_term is False:
            log.warning("WARNING - GAMESS completed unusually check the output")
            ON = 'NG'
        else:
            log.info("GAMESS calculation converged and exited normally")
            ON = 'OK'
            
        data_ON = [ON]
        data = []

        if data_ON == ['OK']:
            log.info("Data was found and will now be processed")
            data = data_fn + data_ON

            # ===== core number =====
            # coreno = catchtxt(txt, ' Initiating ', ' compute processes on ')
            match_cores_ranks = regex_cores_nodes.search(txt)
            if match_cores_ranks:
                log.info("found match for cores/ranks{}".format(match_cores_ranks))
                total_cores = 0

                cores, ranks = match_cores_ranks.groups(0)
                cores = int(cores)
                ranks = int(ranks)
                coreno = cores * ranks
            
            # ===== TOTAL ENERGY =====
            start = 'TOTAL ENERGY'
            end = 'TERMINATED NORMALLY'
            ind = np.array([], dtype=int)

            for m in re.finditer(start, txt):
                ind = np.append(ind, m.start(0))

            txt2 = catchtxt(txt[ind[-1]:], start, end)
            data_etot = re.findall(regex_float, txt2)[0]
            data_etot = float(data_etot)*Ha2eV
            data = data + [data_etot]

            # ===== HOMO/LUMO =====
            start = 'NUMBER OF ELECTRONS'
            end = 'ORBITALS ARE OCCUPIED'

            txt2 = catchtxt(txt, start, end)

            # For the case of OPTIMIZE
            if txt.find('MOLECULAR ORBITALS') != -1:
                start_O = 'MOLECULAR ORBITALS'
                end_O = 'gracefully'
                txt_O = catchtxt(txt, start_O, end_O)
            else:
                txt_O = txt
      
            orbhomo = int(re.findall('\d+', txt2)[-1])
            orblumo = orbhomo+1

            n = 5  # there are 5 columns
            d_homo = orbhomo%n
            d_lumo = orblumo%n
            numspace_homo = 10-len(str(orbhomo))+1  # number of spaces between orbitals
            numspace_lumo = 10-len(str(orblumo))+1

            if d_homo == 0:
                d_homo = n

            if d_lumo == 0:
                d_lumo = n

            if txt.find('MOLECULAR ORBITALS') == -1:
                start = 'EIGENVECTORS'
            else:
                start = 'MOLECULAR ORBITALS'

            end = 'gracefully'
            txt3 = catchtxt(txt, start, end)

            # Extract HOMO energy
            start = ' '*numspace_homo + str(orbhomo)
            end = ' '*21 + 'A'

            txt4 = catchtxt(txt3, start, end)
            eng = re.findall(regex_float, txt4)
            homo_ha = eng[d_homo-1]

            # Extract LUMO energy
            start = ' '*numspace_lumo + str(orblumo)
            end = ' '*21 + 'A'

            txt4 = catchtxt(txt3, start, end)
            eng = re.findall(regex_float, txt4)
            lumo_ha = eng[d_lumo-1]

            lumo_eV = float(lumo_ha)*Ha2eV
            homo_eV = float(homo_ha)*Ha2eV
            E_gap = lumo_eV-homo_eV
            data_hl = [num2str(homo_eV,3), num2str(lumo_eV,3), num2str(E_gap,3)]
            data = data + data_hl

            # ===== Dipole Moment (Debye) =====
            ind = txt_O.find('ELECTROSTATIC MOMENTS')
            txt3 = txt_O[ind:]
            ind = txt3.find('(DEBYE)\n')
            txt3 = txt3[ind+len('(DEBYE)\n'):]
            ind = txt3.find('\n')
            txt3 = txt3[:ind]
            dip = re.findall(regex_float, txt3)
            data_dp = dip[3]
            data = data + [data_dp]
            
            # ===== TOTAL WALL CLOCK TIME =====
            start = 'TOTAL WALL CLOCK TIME'
            end = 'SECONDS'

            ind = np.array([], dtype=int)
            for m in re.finditer(start, txt):
                ind = np.append(ind, m.start(0))

            txt2 = catchtxt(txt[ind[-1]:], start, end)
            data_wct = re.findall(regex_float, txt2)[0]
            data_wctco = float(data_wct)/coreno
            
            data = data + [data_wct, num2str(data_wctco,2)]

            print(data)

            writer.writerow(data)
            
        else:
            log.warning("WARNING - No data found result will be skipped")
            data = data_fn + data_ON
            data = data + ['NA']*(len(header)-2)
            writer.writerow(data)
            
    fcsv.close()
    log.info('\nFinished!')
