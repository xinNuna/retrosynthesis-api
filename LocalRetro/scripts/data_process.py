import csv
import rdkit
from rdkit import Chem

def canonicalize(smiles):
    try:        
        tmp = Chem.MolFromSmiles(smiles)
    except:
        print('no mol')
        return smiles
    if tmp is None:
        return smiles
    tmp = Chem.RemoveHs(tmp)
    [a.ClearProp('molAtomMapNumber') for a in tmp.GetAtoms()]
    return Chem.MolToSmiles(tmp)

input_file = '../outputs/decoded_prediction_train/LocalRetro_USPTO_50K.txt'
output_file = '../outputs/decoded_prediction_train/LocalRetro_USPTO_50K.csv'
rxn_file = '../data/USPTO_50K/raw_train.csv'
reactants = []

# 读取反应文件
with open(rxn_file, "r", encoding="utf-8") as infile:
    csv_reader = csv.reader(infile)
    next(csv_reader)
    for line_index, line in enumerate(csv_reader):
        react , prod = line[2].split(">>")
        react = canonicalize(react)
        reactants.append(react)

with open(input_file, "r", encoding="utf-8") as infile, open(output_file, "w", newline='', encoding="utf-8") as outfile:
        csv_writer = csv.writer(outfile)
        csv_writer.writerow(["sample_idx", "neg_reactants"])
        
        for line_index, line in enumerate(infile):
            line = line.strip()
            tuples = line.split("\t")[1:]

            for t in tuples:
                react = t.split(",")[0].strip("()").strip("'")
                react = canonicalize(react)
                if react == reactants[line_index]:
                    continue
                else:
                    csv_writer.writerow([line_index, react])