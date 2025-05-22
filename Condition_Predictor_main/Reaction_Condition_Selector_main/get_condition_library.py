from ConditionClassifier import get_temp_condition, Classify_reaction_conditions
import argparse
from joblib import Parallel, delayed
import json,gzip
import os
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def main():
    parser = argparse.ArgumentParser(description='Classify reaction conditions')
    parser.add_argument('--data_path', type=str, default='data', help='path to data')
    parser.add_argument('--data_name', type=str, default='USPTO', help='name of data')
    parser.add_argument('--label_path', type=str, default='data/labels', help='path to data')
    parser.add_argument('--save_path', type=str, default='data/condition_library', help='path to save')
    parser.add_argument('--min_num_covered_rxns_by_rxn_centralized_template', type=int, default=5)
    parser.add_argument('--Inclusion', type=float, default=0.8, help='Inclusion of reaction conditions label')
    parser.add_argument('--data_set', type=str, default='train', help='data set to use,can be all, train, test or val.')
    parser.add_argument('--n_jobs', type=int, default=-1, help='number of jobs to run in parallel')
    parser.add_argument('--tpl_radius', type=str, default='0*', help='radius of template')
    args = parser.parse_args()
    # Gets the reaction conditions under each template
    condition_list = get_temp_condition(args)
    # Create the save path directory if it does not exist
    if not os.path.exists(args.save_path):
        os.makedirs(args.save_path)
    # *as a path will report an error, so _1 will be used instead
    if args.tpl_radius != '0*':
        with gzip.open('%s/condition_library_r%s.json.gz'%(args.save_path,args.tpl_radius), 'wt') as f:
            json.dump(condition_list, f)
    else:
        with gzip.open('%s/condition_library_r0_1.json.gz'%(args.save_path), 'wt') as f:
            json.dump(condition_list, f)
    print('start to classify conditions')
    # Classify the conditions in parallel
    classed_condition_list = Parallel(n_jobs=-1, verbose=4)(delayed(Classify_reaction_conditions)(condition_list[i]['conditions'],condition_list[i]['tpl'],condition_list[i]['tpl_smarts'],args) for i in list(condition_list.keys()))
    classed_condition_dic = {}
    for i in range(len(classed_condition_list)):
        classed_condition_dic[classed_condition_list[i][0]['tpl_smart']] = classed_condition_list[i]
    if args.tpl_radius != '0*':
        with gzip.open('%s/classed_conditions_library_r%s.json.gz'%(args.save_path,args.tpl_radius),'wt') as f:
            json.dump(classed_condition_dic,f)
    else:
        with gzip.open('%s/classed_conditions_library_r0_1.json.gz'%(args.save_path),'wt') as f:
            json.dump(classed_condition_dic,f)
    print('classify conditions done')

if __name__ == '__main__':
    main()
    
