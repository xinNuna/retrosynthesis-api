from cli import *

smiles = 'C1=CC(=C(C=C1/C=C/C(=O)O)O)O'
print(get_reactants_template_base(smiles, 50))