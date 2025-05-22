from chemprop.args import PredictArgs
from typing import List
from chemprop.data import get_data,StandardScaler
#from chemprop.train.make_predictions import make_predictions
from chemprop.train.make_predictions import make_predictions
import os

# import inspect # 导入 inspect 模块
# print(f"DEBUG: Actual make_predictions loaded from: {inspect.getfile(make_predictions)}")


from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs

# from extract_template_from_data import get_template
from extract_template_from_data import get_template

import pandas as pd
import csv
import json,gzip
import time
import os
import argparse
from enum import Enum

# test_path = 'Reaction_Condition_Selector_main/data/data_test.csv'
# model_path = 'Reaction_Condition_Selector_main/models'
# label_path = 'Reaction_Condition_Selector_main/data/labels'
# library_path = 'Reaction_Condition_Selector_main/data/condition_library'
# save_path = '/data/prediction'

current_dir = os.path.dirname(os.path.abspath(__file__))
test_path = os.path.join(current_dir, 'data/data_test.csv')
model_path = os.path.join(current_dir, 'models')
label_path = os.path.join(current_dir, 'data/labels')
library_path = os.path.join(current_dir, 'data/condition_library')

class ErrorCode(Enum):
    InvalidSmiles = 1001 # 输入smiles不合法
    TemplateNotFound = 1002 # 找不到模板文件
    LoadLabelsError = 1003 # 获取反应条件标签时出错
    ConditionSelectorError = 1004 # 获取反应条件时出错

class ConditionPredictionError(Exception):
    def __init__(self, error_code: ErrorCode, message: str):
        self.error_code = error_code
        self.message = message
        super().__init__(f"[错误代码: {self.error_code.value}，{self.message}")

class condition_candidate():
    '''
    This class contains operations on condition candidates.
    Candidate here refers to the reaction cluster
    '''
    def __init__(self):
        self.rxn_smart = str() # String to store the reaction SMARTS pattern
        self.class_id = int() # Integer to store the cluster ID
        self.class_label = list() # List to store the cluster labels
        self.conditions = list() # List to store the conditions under the cluster
        self.condition_score = dict() # Dictionary to store the score for each condition
        # Template similarity is the degree of similarity between the templates recorded in the template-condition library and the templates of the predicted reactions, and in this work this value is constant at 1
        self.temp_similarity = float()
        self.max_score = float() # Highest condition score under this reaction cluster
        self.max_score_condition = list() # Corresponding conditions for obtaining the highest score

    def get_class_id(self,class_id):
        """
        Set the class ID for the condition candidate.
        """
        self.class_id = class_id

    def get_class_label(self,class_label):
        '''
        Set the class label for the condition candidate.
        '''
        self.class_label = class_label

    def get_rxn_smart(self,rxn_smart):
        '''
        Set the reaction SMARTS pattern for the condition candidate.
        '''
        self.rxn_smart = rxn_smart

    def get_conditions(self,conditions):
        '''
        Set the conditions for the condition candidate.
        '''
        self.conditions = conditions

    def get_temp_similarity(self,temp_similarity):
        self.temp_similarity = temp_similarity

    def cal_condition_score(self,condition,MPNN_out):
        '''
        计算条件组合的得分
        '''
        score = self.temp_similarity

        if len(condition) != len(MPNN_out):
             # print("Warning: Condition length does not match MPNN_out length.")
             return 0.0

        for i in range(len(condition)-1):
            try:
                score *= MPNN_out[i][condition[i]]
            except Exception as e:
                score *= 1e-10
        return score
    
    # 根据输入的MPNN_out计算每种条件组合的score
    def get_condition_score(self,MPNN_out,condition_key):
        cat_list,solv_list,reag_list = condition_key
        for condition in self.conditions:
            text_condition = decode_condition([condition],cat_list,solv_list,reag_list)
            self.condition_score[str(text_condition[0])] = self.cal_condition_score(condition,MPNN_out)
        self.condition_score = sorted(self.condition_score.items(), key=lambda x:x[1],reverse=True)
    
    # 得到最佳条件
    def get_max_score(self):
        self.max_score = self.condition_score[0][1]
        self.max_score_condition = self.condition_score[0][0]


def MPNN_prediction(model_dir,smiles):
    MPNN_args = ['--test_path', '%s'%test_path, '--checkpoint_dir', '%s'%model_dir, '--preds_path', './sample_preds.csv']
    MPNN_args = PredictArgs().parse_args(MPNN_args)

    preds = make_predictions(MPNN_args,smiles)

    return preds

def get_candidate_range(target_temp,classed_conditions_library):
        candidate_range = []
        for condition_class in classed_conditions_library[target_temp]:
            candidate = condition_candidate()
            candidate.get_class_id(condition_class['class_id'])
            candidate.get_temp_similarity(1)
            candidate.get_class_label(condition_class['class_label'])
            candidate.get_conditions(condition_class['encoded_conditions'])
            candidate_range.append(candidate)
        return candidate_range

def get_condition_score(candidate_range,pred,condition_key):
    if len(candidate_range) == 0:
        return []
    condition_pre = []
    for candidate in candidate_range:
        candidate.get_condition_score(pred,condition_key)
        candidate.get_max_score()
        condition_pre.append({
            'class_id':candidate.class_id,
            'class_label':candidate.class_label,
            'best condition':eval(candidate.max_score_condition),
            'score':candidate.max_score,
            'condition_score':candidate.condition_score
            })
    condition_pre = sorted(condition_pre, key=lambda x:x['score'],reverse=True)
    return condition_pre

def merge_same_lable_candidate(candidate_range):
    label_list = []
    out = []
    for candidate in candidate_range:
        class_label = candidate.class_label[0]+candidate.class_label[1]
        if set(class_label) not in label_list:
            label_list.append(set(class_label))
            out.append(candidate)
        else:
            for i in out:
                if set(i.class_label[0]+i.class_label[1]) == set(class_label):
                    i.conditions += candidate.conditions
    return out

def condition_selector(temp:str,pred:list,classed_conditions_library:dict):
    try:
        candidate_range = get_candidate_range(temp,classed_conditions_library)
        candidate_range = get_condition_score(candidate_range,pred,condition_key)
        return candidate_range
    except ConditionPredictionError:
        raise
    except Exception as e:
        raise ConditionPredictionError(ErrorCode.ConditionSelectorError, f'无法获取反应条件: {e}')
    
def encode_condition(condition_list:list,cat_list:list,solv_list:list,reag_list:list):
        '''
        encode condition list to index
        args:
            condition_list: condition list
            cat_list: catalyst list
            solv_list: solvent list
            reag_list: reagent list
        '''
        out = []
        for condition in condition_list:
            cat = cat_list.index(condition[0])
            solv1 = solv_list.index(condition[1])
            solv2 = solv_list.index(condition[2])
            reag1 = reag_list.index(condition[3])
            reag2 = reag_list.index(condition[4])
            reag3 = reag_list.index(condition[5])
            out.append([cat,solv1,solv2,reag1,reag2,reag3])
        return out
    
def decode_condition(condition_list:list,cat_list:list,solv_list:list,reag_list:list):
    '''
    decode condition list to index
    args:
        condition_list: condition list
        cat_list: catalyst list
        solv_list: solvent list
        reag_list: reagent list
    '''
    out = []
    for condition in condition_list:
        cat = cat_list[condition[0]]
        solv1 = solv_list[condition[1]]
        solv2 = solv_list[condition[2]]
        reag1 = reag_list[condition[3]]
        reag2 = reag_list[condition[4]]
        reag3 = reag_list[condition[5]]
        out.append([cat,solv1,solv2,reag1,reag2,reag3])
    return out

def cal_temp_similarity(tem1,tem2):
    '''
    calculate template similarity
    args:
        tem1: template1
        tem2: template2
    '''
    if tem1 == 'None' or tem2 == 'None':
        return 0
    tem1 = tem1.split('>>')
    tem2 = tem2.split('>>')
    mol1 = [Chem.MolFromSmarts(s) for s in tem1]
    mol2 = [Chem.MolFromSmarts(s) for s in tem2]
    fps1 = [FingerprintMols.FingerprintMol(m) for m in mol1]
    fps2 = [FingerprintMols.FingerprintMol(m) for m in mol2]
    score = 1
    for i in range(len(fps1)):
        score *= DataStructs.FingerprintSimilarity(fps1[i],fps2[i])
    return score

def get_condition_labels(path:str):
    try:
        with open('%s/cat_labels.csv'%path,'r') as f:
            reader = csv.DictReader(f)
            # 只是为了获取催化剂种类，与count无关
            cat_list_N = [row['cat'] for row in reader]
        with open('%s/solv_labels.csv'%path,'r') as f:
            reader = csv.DictReader(f)
            solv_list_N = [row['solv'] for row in reader]
        with open('%s/reag_labels.csv'%path,'r') as f:
            reader = csv.DictReader(f)
            reag_list_N = [row['reag'] for row in reader]
    except ConditionPredictionError:
        raise
    except Exception as e:
        raise ConditionPredictionError(ErrorCode.LoadLabelsError, f"无法获取反应条件标签: {e}")
    
    return cat_list_N,solv_list_N,reag_list_N
        
    
    # 预测函数
def get_rcs_prediction(smiles: List[str] = None, top_k = 5, label_path = label_path):

    global condition_key

    # t1 = time.time() # 计时
    condition_pred = {} # 用于存储预测结果

    # 从输入的SMILES中获得模板并做预测
    
    # 从外部获得的SMILES字符串

    test_data = get_template(smiles)

    # test_data = test_data[:1]
    ids = test_data['_id'].tolist() # Convert to list for easier iteration
    smiles = [[s] for s in test_data['reaction'].tolist()] # Format smiles for MPNN_prediction
    template_r0 = test_data['tpl_SMARTS_r0'].tolist()
    template_r1 = test_data['tpl_SMARTS_r1'].tolist()
    template_r_1 = test_data['tpl_SMARTS_r0*'].tolist()

    classed_conditions_library_r_1 = None # 模板库r0*
    classed_conditions_library_r0 = None # 模板库r0
    classed_conditions_library_r1 = None # 模板库r1


    condition_key = get_condition_labels(label_path)

    # 加载模板库
    try:
        with gzip.open(os.path.join(library_path, 'classed_conditions_library_r0_1.json.gz'), 'r') as f:
            classed_conditions_library_r_1 = json.load(f)
        with gzip.open(os.path.join(library_path, 'classed_conditions_library_r0.json.gz'), 'r') as f:
            classed_conditions_library_r0 = json.load(f)
        with gzip.open(os.path.join(library_path, 'classed_conditions_library_r1.json.gz'), 'r') as f:
            classed_conditions_library_r1 = json.load(f)
    except ConditionPredictionError:
        raise
    except Exception as e:
        raise ConditionPredictionError(ErrorCode.TemplateNotFound, f"无法找到模板文件：{e}")
    
    # 加载本地数据集
    # try:
    #     test_data = pd.read_csv(self.test_path)
    #     test_data = test_data[:1] # Keep or remove slicing as needed
    #     if test_data.empty:
    #          print("Warning: Test data file is empty.")
    #          return {}

    #     ids = test_data['_id'].tolist() # Convert to list for easier iteration
    #     smiles = [[s] for s in test_data['reaction'].tolist()] # Format smiles for MPNN_prediction
    #     template_r0 = test_data['tpl_SMARTS_r0'].tolist()
    #     template_r1 = test_data['tpl_SMARTS_r1'].tolist()
    #     template_r_1 = test_data['tpl_SMARTS_r0*'].tolist()

    #     if not (len(ids) == len(smiles) == len(template_r0) == len(template_r1) == len(template_r_1)):
    #          print("Error: Mismatch in the number of entries across columns in the test data.")
    #          return {}

    # except FileNotFoundError:
    #     print(f"Error loading test data: File not found at {self.test_path}")
    #     return {}
    # except KeyError as e:
    #     print(f"Error loading test data: Missing required column - {e}")
    #     return {}
    # except Exception as e:
    #     print(f"An error occurred while loading test data: {e}")
    #     return {}
    
    # MPNN预测条件
    MPNN_pred = {}

    for target in ['cat','solv0','solv1','reag0','reag1','reag2']:
        model_dir = "%s/MPNN_%s"%(model_path,target)
        MPNN_pred[target] = MPNN_prediction(model_dir,smiles)

    for i in range(test_data.shape[0]):
        if template_r1[i] in classed_conditions_library_r1:
            # MPNN_pred['cat'][i][0])这样可以得到催化剂关于第i个反应的待挑选标签概率列表，传入的参数分别为args、第i个反应的模板、第i个反应的催化剂、溶剂和试剂情况、反应模板库
            condition_pred[str(ids[i])] = ('r1:',condition_selector(template_r1[i],[list(MPNN_pred['cat'][i][0]),list(MPNN_pred['solv0'][i][0]),list(MPNN_pred['solv1'][i][0]),list(MPNN_pred['reag0'][i][0]),list(MPNN_pred['reag1'][i][0]),list(MPNN_pred['reag2'][i][0])],classed_conditions_library_r1))
        elif template_r0[i] in classed_conditions_library_r0:
            condition_pred[str(ids[i])] = ('r0:',condition_selector(template_r0[i],[list(MPNN_pred['cat'][i][0]),list(MPNN_pred['solv0'][i][0]),list(MPNN_pred['solv1'][i][0]),list(MPNN_pred['reag0'][i][0]),list(MPNN_pred['reag1'][i][0]),list(MPNN_pred['reag2'][i][0])],classed_conditions_library_r0))
        else:
            condition_pred[str(ids[i])] = ('r0*:',condition_selector(template_r_1[i],[list(MPNN_pred['cat'][i][0]),list(MPNN_pred['solv0'][i][0]),list(MPNN_pred['solv1'][i][0]),list(MPNN_pred['reag0'][i][0]),list(MPNN_pred['reag1'][i][0]),list(MPNN_pred['reag2'][i][0])],classed_conditions_library_r_1))
    


    final_condition_dict = {}
    for key, value in condition_pred.items():
        final_condition_list = []
        for i in range(len(value)):
            if not isinstance(value[i], str):
                for _, cond_list in enumerate(value[i]):
                    for _, cond_tup in enumerate(cond_list['condition_score']):
                        final_condition_list.append(cond_tup[0])
                temp = final_condition_list[:top_k]
        final_condition_dict[key] = temp

    return final_condition_dict

# if __name__ == '__main__':
#     smiles = ['OB(O)C1=C2C=CC=CC2=CC=C1.BrC1=CC=CC=C1C=O>>O=CC1=CC=CC=C1C1=C2C=CC=CC2=CC=C1', 'BrC1=CN=CS1.CC1=CC=CC(=C1)B(O)O>>CC1=CC=CC(=C1)C1=CN=CS1']
#     parser = argparse.ArgumentParser(description='reaction condition prediction')
#     parser.add_argument('--test_path', type=str, default='./data/data_test.csv', help='path to test data')
#     parser.add_argument('--model_path', type=str, default='./models', help='path to model')
#     parser.add_argument('--label_path', type=str, default='./data/labels', help='path to condition labels')
#     parser.add_argument('--library_path', type=str, default='./data/condition_library', help='path to classed conditions library')
#     parser.add_argument('--save_path', type=str, default='./data/prediction', help='path to save condition prediction')
#     args = parser.parse_args()
#     #predictor = cluster_predictor(args)
#     cond = get_rcs_prediction(smiles)
#     print(cond)
#     final_condition_dict = {}
#     for key, value in cond.items():
#         print(f"反应{key}的信息：\n")
#         final_condition_list = []
#         for i in range(len(value)):
#             if not isinstance(value[i], str):
#                 for _, cond_list in enumerate(value[i]):
#                     for _, cond_tup in enumerate(cond_list['condition_score']):
#                         final_condition_list.append(cond_tup[0])
#                 temp = final_condition_list[:5]
#         final_condition_dict[key] = temp
                





    

