#!/usr/bin/env python

# Copyright IBM Inc. 2022. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
# Author(s):
#   James McDonagh


"""
This module contains an XYZParser class to read and store XYZ file information
"""

# Python packages and utilities
import logging
import pandas as pd
import numpy as np
from rdkit import Chem

class XYZParser(object):
    """
    Class to parse and store xyz file information
    """

    def __init__(self, xyz_file: str = None):
        """
        Initialization function
        :param xyz_file: str - file name and path to read the xyz
        """
        self.xyz_df = None
        self.coordinates_np = None
        self.elements = None
        self.number_of_atoms = 0
        self.description = None
        self.input_file_name = xyz_file

        if self.input_file_name is not None:
            self.parse_xyz_file()


    def parse_xyz_file(self, xyz_file: str = None):
        """
        Function the read an xyz file and store the information
        :param xyz_file: str - path and xyz file name
        :return: None
        """
        log = logging.getLogger(__name__)

        if xyz_file is None:
            if self.input_file_name is not None:
                xyz_file = self.input_file_name
                log.info("Parsing {}".format(xyz_file))
            else:
                log.error("ERROR - no xyz file defined")

        else:
            log.info("Parsing {}".format(xyz_file))

        with open(xyz_file, "r") as xfin:
            self.number_of_atoms = xfin.readline().strip()
            try:
                self.number_of_atoms = int(self.number_of_atoms)
            except TypeError:
                log.warning("Could not convert line 0 to integer number of atoms please check. "
                            "Number of atoms = {}".format(self.number_of_atoms))

            self.description = xfin.readline().strip()

            log.debug("Number of atoms: {}\nDescription: {}".format(self.number_of_atoms, self.description))
            xyz_df_data = []
            xyz_df_with_masses_data = []
            xyz_df_with_charges_data = []
            per_tab = Chem.rdchem.GetPeriodicTable()
            self.masses = []
            self.nuclear_charge = []
            for ith, line in enumerate(xfin):
                if line:
                    try:
                        elm, x, y, z = line.strip().split()
                        log.debug("Parsed: '{}' '{}' '{}' '{}'".format(elm, x, y, z))
                        self.masses.append(per_tab.GetAtomicWeight(elm))
                        self.nuclear_charge.append(per_tab.GetAtomicNumber(elm))
                    except IndexError:
                        log.warning("Line {} should contain element x y z but index error encountered "
                                    "please check".format(ith + 2))

                    try:
                        x = float(x)
                    except TypeError:
                        log.warning("on line {} ({}) type error trying to convert x to float "
                                    "please check".format(ith + 2, line))

                    try:
                        y = float(y)
                    except TypeError:
                        log.warning("on line {} ({}) type error trying to convert y to float "
                                    "please check".format(ith + 2, line))

                    try:
                        z = float(z)
                    except TypeError:
                        log.warning("on line {} ({}) type error trying to convert z to float "
                                    "please check".format(ith + 2, line))

                    log.debug("Prepared: {} {} {} {}".format(elm, x, y, z))
                    xyz_df_data.append([elm, x, y, z])
                    xyz_df_with_masses_data.append([elm, self.masses[-1], x, y, z])
                    xyz_df_with_charges_data.append([elm, self.nuclear_charge[-1], x, y, z])


        self.xyz_df = pd.DataFrame(data=xyz_df_data, columns=["elements", "x", "y", "z"])
        log.debug(self.xyz_df)

        self.coordinates_np = self.xyz_df[["x", "y", "z"]].to_numpy()
        log.debug(self.coordinates_np)

        self.elements = self.xyz_df["elements"].to_list()
        log.debug(self.elements)

        self.xyz_with_masses_df = pd.DataFrame(data=xyz_df_with_masses_data, columns=["elements", "masses", "x", "y", "z"])
        self.xyz_with_charges_df = pd.DataFrame(data=xyz_df_with_charges_data, columns=["elements", "nuclear_charge", "x", "y", "z"])

        log.info("xyz {} parsed successfully".format(xyz_file))


    def get_xyz_df(self) -> pd.DataFrame:
        """
        Function to get the xyz elements and coordinates as a pandas dataframe
        :return: pd.DataFrame
        """
        return self.xyz_df


    def get_element_masses(self) -> list:
        """
        Function to get the masses of each atom in the molecule
        :return: list
        """
        return self.masses


    def get_coordinates(self) -> np.array:
        """
        Function to get the xyz coordinates as a numpy array
        :return: np.array
        """
        return self.coordinates_np


    def get_coordinates_df(self) -> pd.DataFrame:
        """
        Function to get the xyz coordinates as a pandas dataframe
        :return: pd.DataFrame
        """
        coordinates_df = pd.DataFrame(data=self.coordinates_np, columns=["x", "y", "z"])
        return coordinates_df

    def get_coordinates_and_symbols_with_atomic_weights(self) -> pd.DataFrame:
        """
        Function to get the dataframe of symbols atomic weights and coordinates
        :return:
        """
        return self.xyz_with_masses_df

    def get_coordinates_and_symbols_with_nuclear_charges(self) -> pd.DataFrame:
        """
        Function to get the dataframe of symbols nuclear charges (atomic number) and coordinates
        :return:
        """
        return self.xyz_with_charges_df


    def get_elements(self) -> list:
        """
        Function to get the xyz element symbols as a list
        :return: list
        """
        return self.elements


    def get_number_of_atoms(self) -> int:
        """
        Function to get the files number of atoms
        :return: list
        """
        return self.number_of_atoms


    def get_description(self) -> str:
        """
        Function to get the description
        :return: list
        """
        return self.description


#if __name__ == "__main__":
    # xyzparser = XYZParser()
    # xyzparser.parse_xyz_file("Test_data/glucose.xyz")
    # data = xyzparser.get_coordinates_and_symbols_with_atomic_weights()
    # print(data)
    # for ith, row in data.iterrows():
    #     print(" ".join([str(ent) for ent in row.tolist()]) + "\n")

