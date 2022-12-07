
# Copyright IBM Inc. 2015, 2019. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
# Author(s):
#   James McDonagh
#   Michael Johnston

from __future__ import annotations

import typing

import pandas

property_map = {'band-gap': 'gap', 'homo': 'homo', 'lumo': 'lumo', 'electric-moments': 'electric-moments',
                'total-energy': 'total-energy'}


def get_input_ids(input_id_file: str, variables: typing.Dict[str, str]) -> typing.List[str]:
    """Extracts input ids from input_id_file
       The input ids are expected to be in a column SMILES of filename

       Args:
          input_id_file: A path to the file containing the input ids.
            The value of the field interface.inputSpec.source.file
          variables: The global variables of the virtual experiment instance.

       Returns:
          Returns a list of the ids of the inputs in the naming specification defined by interface
    """

    df = pandas.read_csv(input_id_file)
    df.columns = df.columns.str.lower()
    smiles = df['smiles'].to_list()

    start_index = int(variables["startIndex"])
    number_molecules = int(variables["numberMolecules"])

    return smiles[start_index:(start_index + number_molecules)]


def get_properties(property_name: str, property_output_file: str, input_id_file: str) -> pandas.DataFrame:
    """This hook discovers the values of a property for all measured input ids.

       Args:
          property_output_file: A file path. The value of interface.propertiesSpec.name.source.output
            for the given property name (see next field)
          property_name: The value of interface.propertiesSpec.name.
          input_id_file: A file path. The path to the input id file.

       Returns:
          Returns a DataFrame with two columns (input-id, $PropertyName). Each row correspond to an input.
    """

    ids = pandas.read_csv(input_id_file)
    ids.columns = ids.columns.str.lower()
    ids = ids[['label', 'smiles']]
    ids.rename(columns={'smiles': 'input-id'}, inplace=True)

    df = pandas.read_csv(property_output_file)
    df = df[['label', property_map[property_name]]]
    df.rename(columns={property_map[property_name]: property_name}, inplace=True)
    df.columns = df.columns.str.lower()

    result = df.merge(right=ids, how='left', on='label')

    return result
