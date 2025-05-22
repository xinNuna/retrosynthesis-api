import pandas as pd
import re
import numpy as np
from tqdm import tqdm
import template_extractor
from rxnmapper import BatchedMapper
from func_timeout import func_timeout, FunctionTimedOut
from joblib import Parallel, delayed
from rdkit.Chem import rdChemReactions
from rdkit.Chem.rdChemReactions import RemoveMappingNumbersFromReactions
from rdkit import Chem
from enum import Enum
from typing import List
import os
import argparse
import uuid
import warnings

class ErrorCode(Enum):
    InvalidSMILES = 2001 # 非法输入
    RDkitParseError = 2002 # rxnmapper处理映射失败
    ErrorRemoveChirality = 2003 # 去除分子手性失败
    SplitFailed = 2004 # 无法将反应物与产物分离
    EmptyComponent = 2005 # 反应物或产物为空
    TempalteExtractionFailed = 2006 # 模板生成失败
    MissingTemplateSmarts = 2007 # 模板中没有SMARTS
    MappingRemovalFailed = 2008  # 无法处理到模板SMARTS映射
    FunctionTimedOut = 2009 # 函数运行超时

class TemplateExtractionError(Exception):
    def __init__(self, error_code: ErrorCode, message: str, reaction_id: str = "N/A"):
        self.error_code = error_code
        self.message = message
        self.reaction_id = reaction_id
        super().__init__(f"[错误代码: {error_code.value}, 反应 ID: {reaction_id}] {message}")

def remove_chirality(reaction, reaction_id = 'N/A'):
    try:
        if not isinstance(reaction, str):
            raise TemplateExtractionError(ErrorCode.ErrorRemoveChirality, "输入不是有效字符串", reaction_id)
         
        reagent, product = reaction.split('>>')
        reagent = reagent.split('.')
        product = product.split('.')
        reagent = [Chem.MolFromSmiles(i) for i in reagent]
        product = [Chem.MolFromSmiles(i) for i in product]
        reagent = [Chem.MolToSmiles(i, isomericSmiles=False) for i in reagent]
        product = [Chem.MolToSmiles(i, isomericSmiles=False) for i in product]
        reagent = '.'.join(reagent)
        product = '.'.join(product)
        smiles = reagent + '>>' + product
        return smiles
    except TemplateExtractionError:
        raise
    except Exception as e:
        raise TemplateExtractionError(ErrorCode.ErrorRemoveChirality, f"无法处理分子手性: {e}", reaction_id)

def remove_mapping(reaction, reaction_id = 'N/A'):
    try:
        if reaction is None:
            raise TemplateExtractionError(ErrorCode.MAPPING_REMOVAL_FAILED, "反应式为None", reaction_id)
        reaction = rdChemReactions.ReactionFromSmarts(reaction)
        RemoveMappingNumbersFromReactions(reaction)    
        reaction = rdChemReactions.ReactionToSmiles(reaction)
        return reaction
    except TemplateExtractionError:
         raise
    except Exception as e:
        raise TemplateExtractionError(ErrorCode.MappingRemovalFailed, f"去除反应映射失败: {e}", reaction_id)

def _extract(reaction,radius=1):
    try:
        return template_extractor.extract_from_reaction(reaction,radius)
    except KeyboardInterrupt:
        # print('Interrupted')
        raise KeyboardInterrupt
    except Exception:
        return {'reaction_id': reaction['_id']}

def extract(reaction,radius):
    reaction_id = reaction.get('_id', 'N/A')
    try:
        return func_timeout(20, _extract, args=(reaction,radius))
    except FunctionTimedOut:
        # print('Timeout') 
        raise TemplateExtractionError(ErrorCode.FunctionTimedOut, f'模板抽取超时，超过{20}s',reaction['_id'])

def can_parse(reaction, reaction_id = 'N/A'):
    if not isinstance(reaction, str) or '>>' not in reaction:
        raise TemplateExtractionError(ErrorCode.InvalidSMILES, "输入不是有效的Smiles字符串 (缺失 '>>' 或 无效字符串)", reaction_id)
    try:
        reagenr, product = reaction.split('>>')
        if Chem.MolFromSmiles(reagenr) and Chem.MolFromSmiles(product):
            return True
        else:
            return False
    except TemplateExtractionError:
         raise
    except Exception as e:
        raise TemplateExtractionError(ErrorCode.RDkitParseError, f"RDkit解析反应时发生错误: {e}", reaction_id)

def get_data(args):
    data = pd.read_csv(args.data_path)
    can_parse_list = Parallel(n_jobs=-1)(delayed(can_parse)(reaction) for reaction in data[args.reaction_smiles_column])
    rxn_mapper = BatchedMapper(batch_size=32)
    # skip the reactions that cannot be parsed and don't map the reactions(but don't remove them,fill the mapped reactions with the original reactions)
    data = data[can_parse_list] 
    data[args.reaction_smiles_column] = [remove_chirality(reaction) for reaction in data[args.reaction_smiles_column]]
    if not args.skip_mapping:
        print('Mapping reactions...')
        data['Mapped_Reaction'] = list(rxn_mapper.map_reactions(list(data[args.reaction_smiles_column])))
        mapping_column = 'Mapped_Reaction'
    else:
        mapping_column = args.reaction_smiles_column
    print('Done')
    split_smiles = data[mapping_column].str.split('>>',expand=True)
    data['reactants'] = split_smiles[0]
    data['products'] = split_smiles[1]
    data['_id'] = data[args.id_column]
    reactions = data[['_id', 'reactants', 'products']].to_dict('records')
    print('Extracting templates...')
    templates_r0 = Parallel(n_jobs=-1, verbose=4)(delayed(extract)(reaction,0) for reaction in reactions)
    templates_r1 = Parallel(n_jobs=-1, verbose=4)(delayed(extract)(reaction,1) for reaction in reactions)
    print('Done')


    # Add the templates to the DataFrame
    templates_list0 = []
    for template in templates_r0:
        try:
            templates_list0.append(template['reaction_smarts'])
        except:
            templates_list0.append('')
        # Add the templates to the DataFrame
    templates_list1 = []
    for template in templates_r1:
        try:
            templates_list1.append(template['reaction_smarts'])
        except:
            templates_list1.append('')
    data['tpl_SMARTS_r0'] = templates_list0
    data['tpl_SMARTS_r1'] = templates_list1
    templates_r_1 = [remove_mapping(tpl) for tpl in data['tpl_SMARTS_r0']]
    data['tpl_SMARTS_r0*'] = templates_r_1

    # drop if the template is empty
    data = data[data['tpl_SMARTS_r0'] != '']
    data = data[data['tpl_SMARTS_r1'] != '']
    data = data[data['tpl_SMARTS_r0*'] != '']

    # number the templates
    data['template_r0'] = data['tpl_SMARTS_r0'].astype('category')
    data['template_r0'] = data['template_r0'].cat.codes
    data['template_r1'] = data['tpl_SMARTS_r1'].astype('category')
    data['template_r1'] = data['template_r1'].cat.codes
    data['template_r0*'] = data['tpl_SMARTS_r0*'].astype('category')
    data['template_r0*'] = data['template_r0*'].cat.codes

    data.to_csv(args.out_path, index=False)

def get_template(reaction_smiles: List[str], skip_mapping: bool = False) -> pd.DataFrame:
    if not reaction_smiles:
        raise TemplateExtractionError(ErrorCode.InvalidSMILES, f'输入反应式为空')
    # 为输入的每个smiles指定一个id
    reactions_initial_data = [{'_id': i + 1, 'reaction': smiles}
                            for i, smiles in enumerate(reaction_smiles)]

    # processed_reaction_data = []
    # errors = []
    parsable_inputs = []
    # ---解析反应方程式---
    for rxn_data in reactions_initial_data:
        reaction_id = rxn_data['_id']
        reaction_smiles = rxn_data['reaction']
        can_parse(reaction_smiles, reaction_id) # 如果反应不可解析会报错
        parsable_inputs.append(rxn_data)
    
    processed_df = pd.DataFrame(parsable_inputs)
    # ---去除分子手性---
    chirality_results = Parallel(n_jobs=-1)(
        delayed(remove_chirality)(row['reaction'], row['_id']) for index, row in processed_df.iterrows()
    )
    processed_df['reaction_no_chirality'] = chirality_results
    
    # ---使用rxnmapper进行映射
    mapping_column_source = 'reaction_no_chirality'
    mapped_df = processed_df.copy() # 使用处理过后的数据
    if not skip_mapping:
        rxn_mapper = BatchedMapper(batch_size=32)
        input_smiles_for_mapping = mapped_df['reaction_no_chirality'].tolist()
        mapped_reactions = list(rxn_mapper.map_reactions(input_smiles_for_mapping))

        mapped_results_data = []
        for i, mapped_smile in enumerate(mapped_reactions):
             original_id = mapped_df.iloc[i]['_id']
             mapped_results_data.append({'_id': original_id, 'Mapped_Reaction': mapped_smile})

        mapped_results_df = pd.DataFrame(mapped_results_data)
        if not mapped_results_df.empty:
            mapped_df = pd.merge(mapped_results_df, mapped_df[['_id', 'reaction', 'reaction_no_chirality']], on='_id', how='inner')
            mapping_column_source = 'Mapped_Reaction'
        else:
            mapping_column_source = 'reaction_no_chirality'

    processed_df = mapped_df

    # ---分离反应物与产物---
    split_results_data = []
    for index, row in processed_df.iterrows():
         reaction_id = row['_id']
         reaction_smiles = row[mapping_column_source]
         try:
             reactant_smiles, product_smiles = reaction_smiles.split('>>')
             if not reactant_smiles or not product_smiles:
                 raise TemplateExtractionError(ErrorCode.EmptyComponent, "反应物或产物为空", reaction_id)

             split_results_data.append({'_id': reaction_id, 'reactants': reactant_smiles, 'products': product_smiles})
         except TemplateExtractionError as e:
            raise TemplateExtractionError(ErrorCode.SplitFailed, "无效的反应字符串或是缺失 '>>'", reaction_id)
    
    split_df = pd.DataFrame(split_results_data)
    processed_df = pd.merge(split_df, processed_df[['_id', 'reaction']], on='_id', how='inner')


    # ---提取反应模板---
    reactions_for_extraction = processed_df[['_id', 'reactants', 'products']].to_dict('records')
    extraction_results_r0 = Parallel(n_jobs=-1, verbose=0)(delayed(extract)(reaction, 0) for reaction in reactions_for_extraction)
    extraction_results_r1 = Parallel(n_jobs=-1, verbose=0)(delayed(extract)(reaction, 1) for reaction in reactions_for_extraction)
    successful_extractions_data = []
    templates_list0 = []
    templates_list1 = []

    # 将模板添加到相应的反应中
    for template in extraction_results_r0:
        try:
            templates_list0.append(template['reaction_smarts'])
        except:
            templates_list0.append('')

    for template in extraction_results_r1:
        try:
            templates_list1.append(template['reaction_smarts'])
        except:
            templates_list1.append('')

    processed_df['tpl_SMARTS_r0'] = templates_list0
    processed_df['tpl_SMARTS_r1'] = templates_list1
    templates_r_1 = [remove_mapping(tpl) for tpl in processed_df['tpl_SMARTS_r0']]
    processed_df['tpl_SMARTS_r0*'] = templates_r_1
    # 去除没有模板的反应，并记录下相应的id
    # empty_ids = []
    # for i in range(len(reactions_for_extraction)):
    #     original_reaction_dict = reactions_for_extraction[i]
    #     reaction_id = original_reaction_dict['_id']
    #     if processed_df.iloc[i].loc['tpl_SMARTS_r0'] != '':
    #         successful_extractions_data = processed_df.iloc[i]
    #     else:
    #         empty_ids.append()



    processed_df = processed_df[processed_df['tpl_SMARTS_r0'] != '']
    processed_df = processed_df[processed_df['tpl_SMARTS_r1'] != '']
    processed_df = processed_df[processed_df['tpl_SMARTS_r0*'] != '']

    # number the templates
    processed_df['template_r0'] = processed_df['tpl_SMARTS_r0'].astype('category')
    processed_df['template_r0'] = processed_df['template_r0'].cat.codes
    processed_df['template_r1'] = processed_df['tpl_SMARTS_r1'].astype('category')
    processed_df['template_r1'] = processed_df['template_r1'].cat.codes
    processed_df['template_r0*'] = processed_df['tpl_SMARTS_r0*'].astype('category')
    processed_df['template_r0*'] = processed_df['template_r0*'].cat.codes

    return processed_df






    
            
    

if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--data_path', type=str, default='out.csv')
    # parser.add_argument('--reaction_smiles_column', type=str, default='Mapped_Reaction')
    # parser.add_argument('--id_column', type=str, default='_id')
    # parser.add_argument('--out_path', type=str, default='out.csv')
    # parser.add_argument('--skip_mapping', action='store_true')
    # args = parser.parse_args()
    # get_data(args)
    smiles = ['COC(=O)C1=CC(C)=C(Br)C=C1.COC(=O)C1=CC=C(B2OC(C)(C)C(C)(C)O2)C(C)=C1>>COC(=O)C1=CC(C)=C(C=C1)C1=CC=C(C=C1C)C(=O)OC', 
              'COC(=O)C1=CC(C)=C(Br)C=C1.COC(=O)C1=CC=C(B2OC(C)(C)C(C)(C)O2)C(C)=C1>>COC(=O)C1=CC(C)=C(C=C1)C1=CC=C(C=C1C)C(=O)OC', 
              'COC(=O)C1=CC(C)=C(Br)C=C1.COC(=O)C1=CC=C(B2OC(C)(C)C(C)(C)O2)C(C)=C1>>COC(=O)C1=CC(C)=C(C=C1)C1=CC=C(C=C1C)C(=O)OC'
            ]
    reaction = []
    get_template(reaction)
