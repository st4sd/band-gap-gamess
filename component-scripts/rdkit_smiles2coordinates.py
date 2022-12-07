#!/usr/bin/env python

# Copyright IBM Inc. 2015, 2019. All Rights Reserved.
# SPDX-License-Identifier: CC-BY-4.0
# Author(s):
#   James McDonagh
# You can find a copy of the Creative Commons Attribution 4.0 License below
# https://creativecommons.org/licenses/by/4.0/legalcode

"""
This script enables molecule 3D structure generation from a SMILES string
https://www.rdkit.org/docs/GettingStartedInPython.html
creative commons sa 4.0 tutorial used to learn rdkit methods which is (C) 2007-2021 by Greg Landrum
https://creativecommons.org/licenses/by-sa/4.0/
"""

# Python packages and utilities
import os
import sys
import re
import argparse
import logging
import csv
import pandas as pd
from datetime import datetime

#RDKit
sys.path.append('/usr/lib/python2.7/dist-packages/')
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions
from rdkit.Chem import rdMolDescriptors

__version__ = "0.1"
__title__ = os.path.basename(__file__)
__copyright__ = "Copyright IBM Corp. 2017, 2018, 2019"


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

def sample_conformations(mol, n_conformers, r_seed=17591, RMS=False):
    """
    Generates n_conformers on the molecule mol
    :mol: RDKit molecule object
    :n_conformers: the number of conformers to generate
    :r_seed: random seedif not given the geometry will differ each time the call is made
    """

    if RMS is False:
        conformer_index = AllChem.EmbedMultipleConfs(mol, numConfs=n_conformers, numThreads=0,randomSeed=r_seed)

    else:
        conformer_index = AllChem.EmbedMultipleConfs(mol, numConfs=n_conformers, numThreads=0,
                                                 randomSeed=r_seed, maxAttempts=10000, pruneRmsThresh=0.25)
    # options to give above to only accept based on RMS thresh and max attempts maxAttempts=10000, pruneRmsThresh=0.25

    Chem.rdMolTransforms.CanonicalizeMol(mol, ignoreHs=False)
    return list(conformer_index)

def energy_minimize_all_confs(mol, max_int=2000):
    """
    energy minimize all conformations of a molecule which are defined in the molecule object
    :mol: molecule object
    :max_int: max number of minimization iterations
    return 1 not all converged OR 0 all converged
    """
    log = logging.getLogger(__name__)

    # Note numThreads=0 used max avaliable threads for processor
    # Minimize all conformations in the molecule object include inter-molecular potentials and terms
    try:
        log.info("Attempting conformer minimzation with MMFF")
        result = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=max_int, numThreads=0,
                                                   ignoreInterfragInteractions=False)
    except Exception as err:
        log.info("Attempting conformer minimzation with UFF as MMFF encountered an err {}".format(err))
        result = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=max_int, numThreads=0,
                                                  ignoreInterfragInteractions=False)


    converged = [ent[0] for ent in result]
    energies = [ent[1] for ent in result]

    if any(ent == 1 for ent in converged):
        log.info("WARNING - At least some conformations have failed the minimization and are not converged")
        log.info("Converged? (0 = success, 1 = failed): {}".format(converged))
        log.info("Minimized energies: {}".format(energies))

    else:
        log.info("Energy minimization of all conformers successful")
        log.info("Minimized energy: {}\n".format(*["conformer {}: Energy {} Hartree\n".format(i, energy) for i,
                                                                                    energy in enumerate(energies)]))

    # Find lowest energy conformer
    min_ener = 0.0
    indx = 0
    for i, (conv, ener) in enumerate(zip(converged, energies)):
        if i == 0:
            min_ener = ener
            indx = i

        if conv == 0 and ener < min_ener:
            min_ener = ener
            indx = i

    # https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html
    log.info("\n".join(["index {} converged (0 true; 1/-1 false) {} energy {}".format(ind, con, en) for ind, (con, en)
                        in enumerate(zip(converged, energies))]))
    log.info("Minimum energy conformer determined to be index '{}' (zero based) with an "
             "energy {}".format(indx, min_ener))

    return indx

def molecule_image(mol, smile, fnam=None, label_with_num=True):
    """
    Create a molecule image of the smiles string
    :param mol: rdkit molecule object
    :param smile: smiles string
    :param fnam: file name to save the molecule image as
    """
    if fnam is None:
        mf = rdMolDescriptors.CalcMolFormula(mol)
        fnam = "molecule_{}".format(mf)
            
    if label_with_num is True:
         for atom in mol.GetAtoms():
             atom.SetProp('atomLabel',"{}_{}".format(atom.GetSymbol(),str(atom.GetIdx())))

    Draw.MolToFile(mol, "{}.png".format(fnam))

def smi2coordinates(smi, n_conformers=None, random_seed=17591):
    """
    smiles to 3D coordinates
    :param smi: smiles string
    :return: rdkit molecule object with Hydrogens added
    """
    mol = Chem.MolFromSmiles(smi)
    mol_with_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_with_h, randomSeed=random_seed)

    if n_conformers > 0:
        sample_conformations(mol_with_h, n_conformers, r_seed=random_seed)

    indx = energy_minimize_all_confs(mol_with_h)

    return mol_with_h, indx

def inchi2coordinates(inchi, n_conformers=None, random_seed=17591):
    """
    smiles to 3D coordinates
    :param smi: smiles string
    :return: rdkit molecule object with Hydrogens added
    """
    mol = Chem.MolFromInchi(inchi)
    mol_with_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_with_h, randomSeed=random_seed)

    if n_conformers > 0:
        sample_conformations(mol_with_h, n_conformers, r_seed=random_seed)

    indx = energy_minimize_all_confs(mol_with_h)

    return mol_with_h, indx


def xyz_representation(molecule, n_conf=-1, smiles=None):
    """
    make an xyz file representation
    :molecule: is a molecule object which has had Hydrogens added, has been embeded and has been energy minimized
    :n_conf: the molecule conformation number to make an xyz representation of
    """
    coords = Chem.MolToMolBlock(molecule, n_conf)
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
    coord_only = atomic_positions[["x", "y", "z"]]
    mf = rdMolDescriptors.CalcMolFormula(molecule)

    with open("{}-confromer{}.xyz".format(rdMolDescriptors.CalcMolFormula(molecule), n_conf), "w") as fout:
        fout.write("{}\n".format(molecule.GetNumAtoms()))
        if smiles is not None:
            fout.write("From RDKit: Molecule: {}: SMILES {}\n".format(mf, smiles))
        else:
            fout.write("From RDKit: Molecule: {}\n".format(mf))
        atomic_positions.to_csv(fout, header=None, index=None, sep=" ", mode="a")

def xyz(molecule, n_conf=-1):
    """
    Get elements and atomic coordinates
    :molecule: is a molecule object which has had Hydrogens added, has been embeded and has been energy minimized
    :n_conf: the molecule conformation number to make an xyz representation of
    return pandas dataframe of elements and coordinates
    """
    coords = Chem.MolToMolBlock(molecule, n_conf)
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

def gamess_input(mol, smi, method="SCFTYP=UHF MULT=1 RUNTYP=OPTIMIZE COORD=UNIQUE",
                 basis="GBASIS=N311 NGAUSS=6 NDFUNC=1 DIFFSP=.TRUE.", indx=-1, spin=None, charge=None):
    """
    Make Gamess input TODO: may need option for setting the charge and multiplicity if these are not defined in the
    TODO: method or basis lines
    :param mol: rdkit molecule class i.e. mol = Chem.MolFromSmiles(smile_string)
    :param method: The calculation method to use in Gamess
    :param basis: The calculation basis set set up
    :return: 
    """
    
    mf = rdMolDescriptors.CalcMolFormula(mol)
    coordinates = xyz(mol, n_conf=indx)
    coords_csv = coordinates.to_csv(header=False, index=False, sep=" ")

    with open("molecule.inp", "w") as fo:
        fo.write("$CONTRL {} $END\n".format(method))
        fo.write("$BASIS {} $END\n".format(basis))
        fo.write("$SYSTEM TIMLIM=1 $END\n")
        fo.write("$SCF DIRSCF=.TRUE. NCONV=9 $END\n")
        fo.write("$GUESS  GUESS=HUCKEL $END\n")
        fo.write("$DATA\n")
        fo.write("{} {}\n".format(mf, smi))
        fo.write("C1\n")
        for ent in coords_csv:
            fo.write(" {}\n".format(ent))
        fo.write("$END\n")

def gamess_input_from_template(mol, smi, name, template, indx=-1, spin=None, charge=None):
    """
    Function to make a GAMESS input from a template options file
    :param molecule: RDKit mol
    :param smi: smiles string
    :param name: molecule name
    :param template: template file to read and use the options from
    :param indx: molecule conformer index
    :return:
    """

    log = logging.getLogger(__name__)
    mf = rdMolDescriptors.CalcMolFormula(mol)
    coordinates = xyz(mol, n_conf=indx)
    atomic_masses = ["{:.1f}".format(round(atom.GetAtomicNum(),0)) for atom in mol.GetAtoms()]
    log.info("Atomic mass: {}".format(atomic_masses))
    coordinates.insert(loc=1, column="masses", value=atomic_masses)
    coords_csv = coordinates.to_csv(header=False, index=False, sep=" ")

    with open(template, "r") as fin:
        temp_data = [line.strip() for line in fin if line]

    out_filename = "molecule.inp"

    # Try to get data for charge and multiplicity (NOTE: currently not sure of the multiplicity
    # will be correct by this method needs some testing)
    elec_data = count_electons(mol)
    if spin is not None:
        elec_data[1] = spin.strip()
        log.info("Spin set to user defined value {}".format(elec_data[1]))

    if charge is not None:
        elec_data[0] = charge.strip()
        log.info("Charge set to user defined value {}".format(elec_data[0]))

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


def count_electons(mol, indx=-1):
    """
    This is a function to count the total number of electrons in a molecule and the number of valence electrons
    mol : RDKit molecule object
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


def get_cation(smile):
    """
    This is a function calculates the molecular charge and return +ve smiles and an index if the smile is multi
    molecule i.e. '.' separated
    smile : string - smiles string
    """

    log = logging.getLogger(__name__)

    for inx, smi in enumerate(smile.split(".")):
        log.info("checking index {} smiles {}".format(inx, smi))
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        total_charge = 0

        for atom in mol.GetAtoms():
            atom_chr = atom.GetFormalCharge()
            total_charge = total_charge + atom_chr

        if total_charge > 0:
            log.info("Smile {} is a cation. Molecular charge {}".format(inx, total_charge))
            # returns smile string or inx (assuming split on '.')
            return smi

    log.info("No cations found in {}.".format(smile))
    return None

def get_anion(smile):
    """
    his is a function calculates the molecular charge and return -ve smiles and an index if the smile is multi
    molecule i.e. '.' separated
    smile : string - smiles string
    """

    log = logging.getLogger(__name__)

    for inx, smi in enumerate(smile.split(".")):
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        total_charge = 0

        for atom in mol.GetAtoms():
            atom_chr = atom.GetFormalCharge()
            total_charge = total_charge + atom_chr

        if total_charge < 0:
            log.info("Smile {} is a anion. Molecular charge {}".format(inx, total_charge))
            # returns smile string or inx (assuming split on '.')
            return smi

    log.info("No anions found in {}.".format(smile))
    return None

def run():
    """
    Main run function to act based on user input
    """
    try:
        usage = "python {} Required Parameters {}\n".format(__title__, " options .....")

        parser = argparse.ArgumentParser(description="Command line binary points script",
                                         usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        #parser.add_argument("--smiles", metavar="smiles string(s)", action="store",
        #                    help="SMILES strings to parse and get 3D coordinates of", nargs="+", default=None)
        parser.add_argument("--input", metavar="smiles string(s)", action="store",
                            help="SMILES/InChI strings to parse and get 3D coordinates of or a csv file",
                            nargs="+", default=None)
        parser.add_argument("--spin", metavar="spin states for smiles", action="store",
                            help="SPin state (multiplicity) to use for the molecules", nargs="+", default=None)
        parser.add_argument("--charge", metavar="Charges for the smiles", action="store",
                            help="SPin state (multiplicity) to use for the molecules", nargs="+", default=None)

        group = parser.add_argument_group("Optional arguments")
        group.add_argument("--noimage", action="store_true", default=False, help="Save image of the molecule")
        group.add_argument("--noxyz", action="store_true", default=False, help="Save xyz file of the coordinates")
        group.add_argument("--method", action="store", default="SCFTYP=UHF MULT=1 RUNTYP=OPTIMIZE COORD=UNIQUE"
                           , help="The GAMESS method line")
        group.add_argument("--basis", action="store", default="GBASIS=N311 NGAUSS=6 NDFUNC=1 DIFFSP=.TRUE.",
                           help="The GAMESS basis line")
        group.add_argument("--template", action="store", default=None,
                           help="The GAMESS input file template specifying GAMESS options")
        group.add_argument("--names", action="store", default=None, nargs="+",
                           help="Molecule names")
        #group.add_argument("--inchi", metavar="InChI string(s)",action="store", default=None, nargs = "+",help = "InChI strings")
        group.add_argument("--random_seed", action="store", default=17591, type=int,
                           help="Random seed number for RDKit 3D structure set up")
        group.add_argument("--n_conformers", action="store", default=50, type=int,
                          help="Number of chemical conformtions to consider")
        group.add_argument("--row", action="store", default=None, type=int,
                          help="Row of csv file of smiles to make gamess input for")
        group.add_argument("--species", action="store", default="all",
                           help="If a multi-molecule smiles string is provided this species will be selected from the "
                                "string. Arguments can be cation, anion or int specifying the index (zero based). "
                                "default.")
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

    smiles_input = False
    inchi_input = False
   

    if len(op.input) == 1 and re.search(".[csv]+$", op.input[-1]):
        log.info("Input appears to be located in file {}\nDetermining file information\n".format(op.input[-1]))


        data = pd.read_csv(op.input[-1], header=0, dtype=str)
        data.columns = [h.strip().lower() for h in data.columns]
        
        if 'smiles' in data.columns:
            smiles_input = True
        elif 'inchi' in data.columns:
            inchi_input = True
        else:
            log.error("ERROR - file header smiles or inchi not found one must be given")
            raise RuntimeError
        
    else:
        if "inchi" in op.input.lower():
            inchi_input = True
            log.info("InChI given to command line:\n{}\nCycling through .....\n".format(op.input))
        else:
            smiles_input = True
            log.info("SMILES given to command line:\n{}\nCycling through .....\n".format(op.input))
        
    if smiles_input is True:
        try:
            if op.row is None:
                op.input = data["smiles"].values
            else:
                op.input = data["smiles"].values[op.row]

            if isinstance(op.input, str):
                log.debug("{} {}".format(type(op.input), op.input))
                op.input = [op.input]
            log.info("SMILES parsed from file:\n{}".format(op.input))

        except KeyError as err:
            log.error("ERROR - key 'SMILES' not found from csv file row 0 which are assumed to be headers")
            raise err

        except IndexError as err:
            log.error("ERROR - index error likely due to input option row. This value is assumed to be 0 based i.e. "
            "if there is one entry the row should be 0. {}".format(err))
            log.info("Data length {} index (row) {}".format(len(data), op.row))
            raise err

        try:
            if op.row is None:
                op.spin = data["spin"].values
            else:
                op.spin = data["spin"].values[op.row]

            log.info("spin parsed from file:\n{}".format(op.spin))

        except KeyError as err:
            log.info(
                "NOTE - No multiplicity (spin) states specified (no header spin) in file will set this automatically")
        
        except IndexError as err:
            log.error("ERROR - index error likely due to input option row. This value is assumed to be 0 based i.e. "
            "if there is one entry the row should be 0. {}".format(err))
            log.info("Data length {} index (row) {}".format(len(data), op.row))
            raise err

        try:
            if op.row is None:
                op.charge = data["charge"].values
            else:
                op.charge = data["charge"].values[op.row]

            log.info("charge parsed from file:\n{}".format(op.charge))

        except KeyError as err:
            log.info("NOTE - No charge states specified (no header charge) in file will set this automatically")
        
        except IndexError as err:
            log.error("ERROR - index error likely due to input option row. This value is assumed to be 0 based i.e. "
            "if there is one entry the row should be 0. {}".format(err))
            log.info("Data length {} index (row) {}".format(len(data), op.row))
            raise err

        if isinstance(op.spin, str):
            op.spin = [op.spin]

        if isinstance(op.charge, str):
            op.charge = [op.charge]

        if op.spin is None:
            op.spin = [None] * len(op.input)

        if op.charge is None:
            op.charge = [None] * len(op.input)


        for inx, smi in enumerate(op.input):
            if op.names is not None:
                name = op.names[inx]
            else:
                name = inx

            # Deal with user choices on species if smi is a multi-molecule smiles string(i.e. does it contain '.')
            if '.' in smi:
                log.info("Multi-molecule SMILES string detected.")
                if op.species.isdigit():
                    log.info("Selecting user defined index in the multi-molecule string {}.".format(op.species))
                    try:
                        smi = smi.split(".")[int(op.species)]
                        log.info("Smile {} is index {}.".format(smi, op.species))
                    except IndexError as err:
                        log.error("Index error for user specified molecule index ({}) "
                                  "in multi-molecule smiles".format(op.species))
                        raise err
                elif op.species.lower() == "cation":
                    log.info("Selecting the first cation found")
                    smi = get_cation(smi)
                    if smi is None:
                        log.error("ERROR - unable to locate cation as specified by the user through "
                                  "--species {}".format(op.species))
                        raise RuntimeError
                elif op.species.lower() == "anion":
                    log.info("Selecting the first anion found")
                    smi = get_anion(smi)
                    if smi is None:
                       log.error("ERROR - unable to locate anion as specified by the user through "
                                 "--species {}".format(op.species))
                       raise RuntimeError
                elif op.species.lower() == "all":
                    log.info("Running on multi-molecule system")


            molecule, indx = smi2coordinates(smi, n_conformers=op.n_conformers, random_seed=op.random_seed)
            mf = rdMolDescriptors.CalcMolFormula(molecule)
            with open("{}.sdf".format(mf), "w") as fo:
                fo.write(Chem.MolToMolBlock(molecule))

            if op.template is not None:
                gamess_input_from_template(molecule, smi, name, op.template, indx=indx, spin=op.spin[inx],
                                           charge=op.charge[inx])
            else:
                gamess_input(molecule, smi, method=op.method, basis=op.basis, indx=indx, spin=op.spin[inx],
                             charge=op.charge[inx])

            if not op.noxyz:
                xyz_representation(molecule, n_conf=indx, smiles=smi)
            if not op.noimage:
                molecule_image(molecule, smi)

    if inchi_input is True and smiles_input is False:

        if len(op.input) == 1 and re.search(".[csvtxt]+$", op.input[-1]):
            log.info("InChI appear to be located in file {}\nParsing .....\n".format(op.input[-1]))

            # Checks for the first row headers and assumes one is inchi case is set to lower. Also find the column
            # separators
            with open(op.input[-1]) as csvin:
                sniffer = csv.Sniffer().sniff(csvin.read(1024))
                csvin.seek(0)
                has_header = csv.Sniffer().has_header(csvin.read(1024))
                csvin.seek(0)
                log.info("file separator '{}' header present {}".format(str(sniffer.delimiter), has_header))
                if not has_header:
                    log.error("WARNING - Noheaders found in initial pass will try to continue but it is unlikely to "
                              "succeed")

            data = pd.read_csv(op.input[-1], header=0, dtype=str, sep=sniffer.delimiter)
            data.columns = [h.strip().lower() for h in data.columns]
            log.info("Data frame headers: {}".format(data.columns))
            log.info("Parsed csv file: {}".format(data))
            try:
                if op.row is None:
                    op.input = data["inchi"].values
                else:
                    op.input = data["inchi"].values[op.row]

                if isinstance(op.input, str):
                    log.debug("{} {}".format(type(op.input), op.input))
                    op.input = [op.input]
                log.info("InChI parsed from file:\n{}".format(op.input))

            except KeyError as err:
                log.error("ERROR - key 'inchi' not found from csv file row 0 which are assumed to be headers")
                raise err

            try:
                if op.row is None:
                    op.spin = data["spin"].values
                else:
                    op.spin = data["spin"].values[op.row]

                log.info("spin parsed from file:\n{}".format(op.spin))

                if isinstance(op.spin, str):
                    op.spin = [op.spin]

            except KeyError as err:
                log.info("NOTE - No multiplicity (spin) states specified (no header spin) in file will set this automatically")

            try:
                if op.row is None:
                    op.charge = data["charge"].values
                else:
                    op.charge = data["charge"].values[op.row]

                log.info("charge parsed from file:\n{}".format(op.charge))

            except KeyError as err:
                log.info("NOTE - No charge states specified (no header charge) in file will set this automatically")

            if isinstance(op.spin, str):
                op.spin = [op.spin]

            if isinstance(op.charge, str):
                op.charge = [op.charge]

            if op.spin is None:
                op.spin = [None] * len(op.input)

            if op.charge is None:
                op.charge = [None] * len(op.input)
        else:
            log.info("InChI given to command line:\n{}\nCycling through .....\n".format(op.input))

        if op.species != "all":
            # TODO: add functionality for Inchi molecule/somponent selection
            log.warning("WARNING - sub-structure/molecule selection is not available for InChi currently."
                        "The set up will proceed with all components/molecules")

        for inx, inch in enumerate(op.input):
            if op.names is not None:
                name = op.names[inx]
            else:
                name = inx

            molecule, indx = inchi2coordinates(inch, n_conformers=op.n_conformers, random_seed=op.random_seed)
            mf = rdMolDescriptors.CalcMolFormula(molecule)
            with open("{}.sdf".format(mf), "w") as fo:
                fo.write(Chem.MolToMolBlock(molecule, confId=indx))

            if op.template is not None:
                gamess_input_from_template(molecule, inch, name, op.template, indx=indx, spin=op.spin[inx],
                                           charge=op.charge[inx])
            else:
                gamess_input(molecule, inch, method=op.method, basis=op.basis, indx=indx, spin=op.spin[inx],
                                           charge=op.charge[inx])

            if not op.noxyz:
                xyz_representation(molecule, n_conf=indx, smiles=inch)
            if not op.noimage:
                molecule_image(molecule, inch)

if __name__ == "__main__":
    run()

# Example commandline data.csv -x -i
# data.csv
#label, SMILES
#0,c1ccccc1
#1,CCCCCC
#2,c1ccccc1C
#3,O
#4,HCl
# Example commandline CCCCCC c1ccccc1 c1ccccc1CCF
# Example commandline CCCCCC c1ccccc1 c1ccccc1CCF -x
# Example commandline CCCCCC c1ccccc1 c1ccccc1CCF -i
# Example commandline CCCCCC c1ccccc1 c1ccccc1CCF -x -i

# Direct command lines not reading a file
# --smiles O CO CCO --names water methanol ethanol --template gamess_template.txt
# Adds a test for spin state
# --smiles O CO CCO [CH2]Cl --names water methanol ethanol cl --template gamess_template.txt
# --inchi InChI=1S/H2O/h1H2 InChI=1S/CH4O/c1-2/h2H,1H3 InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3 --names water methanol ethanol --template gamess_template.txt
