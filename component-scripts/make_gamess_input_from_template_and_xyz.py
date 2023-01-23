#!/usr/bin/env python

# Copyright IBM Inc. 2015, 2020. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
# Author(s):
#   James McDonagh

# python packages
import argparse
import logging
import os
from datetime import datetime
import logging
import pandas as pd
import re
import xyz_class

#RDKit
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors, rdmolfiles, Descriptors

__version__ = "0.1"
__authors__ = "James L. McDonagh"
__contact__ = "james.mcdonagh@uk.ibm.com"
__title__ = os.path.basename(__file__)
__copyright__ = "Copyright IBM Corp. 2020"

# functions
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

def count_electons(mol, indx=-1):
    """
    This is a function to count the total number of electrons in a molecule and the number of valence electrons
    mol : RDKit molecule object
    indx: int - confomrer index
    """
    log = logging.getLogger(__name__)

    mf = rdMolDescriptors.CalcMolFormula(mol)

    number_of_electrons = 0
    total_charge = 0
    number_unpaired_es = 0

    for atom in mol.GetAtoms():
        atom_n_e = atom.GetAtomicNum()
        atoms_sym = atom.GetSymbol()
        atom_chr = atom.GetFormalCharge()
        # TODO: assuming if the atom has a 'radical electron' it means unpaired accounting for bonding etc
        # https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html
        rade = atom.GetNumRadicalElectrons()

        number_of_electrons = number_of_electrons + atom_n_e - atom_chr
        total_charge = total_charge + atom_chr
        number_unpaired_es = number_unpaired_es + rade

    number_of_valence_electrons = Chem.Descriptors.NumValenceElectrons(mol)
    spin = number_unpaired_es + 1

    if (number_of_electrons % 2) == 0:
        radical = False
    else:
        radical = True
    log.info("Molecular formula: {}".format(mf))
    log.info("Total number of electrons: {}".format(number_of_electrons))
    log.info("Total charge: {}".format(total_charge))
    log.info("From RDKit number of valence electrons {}".format(number_of_valence_electrons))
    log.info("Radical molecule (only valid for organic moleules) {}".format(radical))
    log.info("Spin state is {}\n--------------------------\n".format(spin))

    return [total_charge, spin, number_of_valence_electrons, number_of_electrons, radical]

def xyz(molecule, n_conf=-1):
    """
    Get elements and atomic coordinates
    :molecule: is a molecule object which has had Hydrogens added, has been embeded and has been energy minimized
    :n_conf: the molecule conformation number to make an xyz representation of
    return pandas dataframe of elements and coordinates
    """
    coords = Chem.MolToMolBlock(molecule, confId=n_conf)
    atomic_positions = []
    for ent in coords.split("\n"):
        elms = ent.split()
        try:
            float(elms[0])
            if len(elms) == 16:
                atomic_positions.append([elms[3], float(elms[0]), float(elms[1]), float(elms[2])])
        except ValueError:
            pass
        except IndexError:
            pass

    atomic_positions = pd.DataFrame(atomic_positions, columns=["element", "x", "y", "z"])

    return atomic_positions

def gamess_input_from_template(mol, gamess_xyz_data, name, template, indx=-1, spin=None, charge=None, out_filename="molecule.inp"):
    """
    Function to make a GAMESS input from a template options file
    :param molecule: RDKit mol
    :param gamess_xyz_data: pandas dataframe - element atomic mass x y z
    :param name: molecule name
    :param template: template file to read and use the options from
    :param indx: molecule conformer index
    :return:
    """

    log = logging.getLogger(__name__)

    mf = rdMolDescriptors.CalcMolFormula(mol)
    coordinates = gamess_xyz_data[["elements", "x", "y", "z"]]
    atomic_masses = ["{:.1f}".format(round(float(w), 0)) for w in gamess_xyz_data["masses"]]
    log.info("Atomic mass: {}".format(atomic_masses))
    coordinates.insert(loc=1, column="masses", value=atomic_masses)
    coords_csv = coordinates.to_csv(header=False, index=False, sep=" ")

    with open(template, "r") as fin:
        temp_data = [line.strip() for line in fin if line]

    elec_data = count_electons(mol)
    if spin is not None:
        elec_data[1] = spin.strip()
        log.info("Spin set to user defined value {}".format(elec_data[1]))
    else:
        log.warning("WARNING - setting molecular spin automatically this is not recommended!")

    if charge is not None:
        elec_data[0] = charge.strip()
        log.info("Charge set to user defined value {}".format(elec_data[0]))
    else:
        log.warning("WARNING - setting molecular charge automatically please check.")

    log.info("Output file : {}".format(out_filename))

    with open(out_filename, "w") as fout:
        for line in temp_data:
            log.debug("Line for input files: {}".format(line))
            if "MULT" in line:
                line_split = ["MULT={}".format(elec_data[1]) if "MULT" in ent else ent for ent in line.split()]
                line = " ".join([str(elm) for elm in line_split])
                log.debug("MULT found set to {}".format(elec_data[1]))

            if "ICHARG" in line:
                line_split = ["ICHARG={}".format(elec_data[0]) if "ICHARG" in ent else ent for ent in line.split()]
                line = " ".join([str(elm) for elm in line_split])
                log.debug("ICHARG found set to {}".format(elec_data[0]))

            if str(elec_data[1]) != "1" and "SCFTYP=RHF" in line:
                line_split = ["SCFTYP=ROHF".format(elec_data[0]) if "SCFTYP" in ent else ent for ent in line.split()]
                line = " ".join([str(elm) for elm in line_split])
                log.debug("SCFTYP changes to ROHF due to openshell state")

            fout.write(" {}\n".format(line.strip()))

        #fout.write("\n")
        fout.write(" $DATA\n")
        fout.write("{} {}\n".format(name, re.sub(r"[^\w]", "", mf)))
        fout.write(" C1\n")
        for ent in coords_csv.split("\n"):
            if ent:
                fout.write(" {}\n".format(ent))
        fout.write(" $END\n")

    return out_filename


def main():
    """
    """

    try:
        usage = "python {} {}\n".format(__title__, " options .....")

        parser = argparse.ArgumentParser(description="Command line input options",
                                         usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument("-xp",
                            "--xyz_path",
                            metavar="FILEPATH",
                            action="store",
                            required=True,
                            help="xyz path to read")

        parser.add_argument("-xf",
                            "--xyz_file",
                            metavar="FILEPATH",
                            action="store",
                            required=True,
                            help="xyz file to read")

        parser.add_argument("-g",
                            "--gamess_template",
                            type=str,
                            action="store",
                            required=True,
                            help="Template for a gamess input file. Should contain all keywords and "
                                 "parameters this script will append the atomic positions.")

        parser.add_argument("-sf",
                            "--sdf",
                            type=str,
                            action="store",
                            help="smiles for the molecule",
                            required=True)

        parser.add_argument("-sp",
                            "--sdf_path",
                            type=str,
                            action="store",
                            help="smiles for the molecule",
                            required=True)

        parser.add_argument("-o",
                            "--output_name",
                            type=str,
                            action="store",
                            help="output gamess input file from template",
                            default="molecule.inp")

        parser.add_argument("--spin",
                            type=str,
                            action="store",
                            help="molecules spin (multiplicity)",
                            default=None)

        parser.add_argument("--charge",
                            type=str,
                            action="store",
                            help="molecules charge state",
                            default=None)

        parser.add_argument("--loglev",
                            action="store",
                            default="INFO",
                            help="log level")

        op = parser.parse_args()

    except argparse.ArgumentError as ee:
        print("\nERROR - command line arguments are ill defined please check the arguments\n")
        raise ee

    setup_logger(os.getcwd(), op.loglev)
    log = logging.getLogger(__name__)
    log.info("\nAuthors       : {}\nOrganiszation : IBM Research Europe\nCreated"
             "       : June 2020\nProgram       : {}\nVersion       : {}\n "
             "--------------------------------\n".format(__authors__, __title__, __version__))

    log.info("Command line input =\n\t{}".format(op))

    # These are separated for the ST4SD engine so it is easier to specify
    op.xyz_file = os.path.join(op.xyz_path, op.xyz_file)
    op.sdf = os.path.join(op.sdf_path, op.sdf)

    # Read xyz with custom xyz class
    xyz_reader = xyz_class.XYZParser()
    xyz_reader.parse_xyz_file(op.xyz_file)
    gamess_xyz_data = xyz_reader.get_coordinates_and_symbols_with_atomic_weights()
    log.info(gamess_xyz_data)

    # SDF file is used to pass molecule data so RDKit can be used to get common properties
    mol = Chem.rdmolfiles.MolFromMolFile(op.sdf)

    if mol is not None:
        name = f"{rdMolDescriptors.CalcMolFormula(mol)} {Chem.MolToSmiles(mol)}"
        gamess_input_from_template(mol,
                                   gamess_xyz_data,
                                   name,
                                   op.gamess_template,
                                   indx=-1,
                                   spin=op.spin,
                                   charge=op.charge,
                                   out_filename=op.output_name)

if __name__ == "__main__":
    main()

# mol = (Chem.AddHs(Chem.MolFromSmiles("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O")))
# AllChem.EmbedMolecule(mol)
# with open("{}.sdf".format("glucose"), "w") as fo:
#     fo.write(Chem.MolToMolBlock(mol))