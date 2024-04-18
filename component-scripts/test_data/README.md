To test SMILES to GAMESS conversion

python ../rdkit_smiles2coordinates.py --input input_smiles.csv --row 0
python ../make_gamess_input_from_template_and_xyz.py -xf 0.xyz -g ../../dft/data-dft/input_molecule.txt -xp . -sf 0.sdf -sp . 
