from argparse import ArgumentParser
import os
import torch
from Test import *
from Decode_predictions import *
import csv

def preprocess(smiles):
    smiles_full = '>>' + smiles
    data = [
        {
            "id": 0,
            "class": "UNK",
            "reactants>reagents>production": smiles_full
        }
    ]
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_file_path = os.path.join(script_dir, '..', 'data', 'USPTO_50K', 'raw_test.csv')
    with open(data_file_path, "w", newline="", encoding="utf-8") as csvfile:
        fieldnames = ["id", "class", "reactants>reagents>production"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

def get_edit():
    parser = ArgumentParser('LocalRetro input arguements')
    parser.add_argument('-g', '--gpu', default='cuda:0', help='GPU device to use')
    parser.add_argument('-d', '--dataset', default='USPTO_50K', help='Dataset to use')
    parser.add_argument('-c', '--config', default='default_config.json', help='Configuration of model')
    parser.add_argument('-b', '--batch-size', default=16, help='Batch size of dataloader')
    parser.add_argument('-k', '--top_num', default=100, help='Num. of predictions to write')
    parser.add_argument('-nw', '--num-workers', type=int, default=0, help='Number of processes for data loading')
    args = parser.parse_args().__dict__
    args['mode'] = 'test'
    args['device'] = torch.device(args['gpu']) if torch.cuda.is_available() else torch.device('cpu')
    print('Using device %s' % args['device'])
    main_test(args)

def decode_(top_k):
    parser = ArgumentParser('Decode Prediction')  
    parser.add_argument('-d', '--dataset', default='USPTO_50K', help='Dataset to use') 
    parser.add_argument('-m', '--model', default='default', help='Model to use')  
    parser.add_argument('-k', '--top-k', default=10, help='Number of top predictions')  
    args = parser.parse_args().__dict__ 
    args['top_k'] = top_k
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_file_path_1 = os.path.join(script_dir, '..', 'outputs', 'decoded_prediction')
    data_file_path_2 = os.path.join(script_dir, '..', 'outputs', 'decoded_prediction_class')
    mkdir_p(data_file_path_1)  
    mkdir_p(data_file_path_2)  

    reactants = main_decode(args) 
    return reactants

def get_reactants_template_base(smiles, top_k):
    preprocess(smiles)
    get_edit()
    reactants = decode_(top_k)
    print("\nThe output has been obtained")
    return reactants

if __name__ == '__main__':
    smiles = '[CH3:1][n:13]1[c:3]([CH3:2])[n:4][c:5]2[cH:6][c:7]([Cl:8])[c:9]([Cl:10])[cH:11][c:12]21'
    print(get_reactants_template_base(smiles, 10))

