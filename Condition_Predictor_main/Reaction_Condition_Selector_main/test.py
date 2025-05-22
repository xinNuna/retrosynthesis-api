import uuid
import pandas as pd
import numpy as np

class cluster_predictor:
    def __init__(self, args):
        self.args = args
        self.model_path = model_path # 模型加载路径
        self.label_path = label_path # 标签加载路径
        self.library_path = library_path # 模板库加载路径
        self.test_path = test_path # 测试集加载路径
        self.mpnn_model = []
        self.condition_key = self.get_condition_labels(label_path)
        # self.save_path = args.save_path # 保存路径，当前已关闭保存功能

        self.cat_list, self.solv_list, self.reag_list = self.condition_key
    
    def MPNN_prediction(self,model_dir,smiles):
        MPNN_args = ['--test_path', '%s'%self.test_path, '--checkpoint_dir', '%s'%model_dir, '--preds_path', './sample_preds.csv']
        MPNN_args = PredictArgs().parse_args(MPNN_args)
        preds = make_predictions(MPNN_args,smiles)
        return preds
    
    def get_candidate_range(self,target_temp,classed_conditions_library):
        candidate_range = []
        for condition_class in classed_conditions_library[target_temp]:
            candidate = condition_candidate()
            candidate.get_class_id(condition_class['class_id'])
            candidate.get_temp_similarity(1)
            candidate.get_class_label(condition_class['class_label'])
            candidate.get_conditions(condition_class['encoded_conditions'])
            candidate_range.append(candidate)
        return candidate_range
    
    
    def get_condition_score(self,candidate_range,pred,condition_key):
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
    
    def condition_selector(self, temp:str,pred:list,classed_conditions_library:dict):
        try:
            candidate_range = self.get_candidate_range(temp,classed_conditions_library)
            candidate_range = self.get_condition_score(candidate_range,pred,self.condition_key)
            return candidate_range
        except ConditionPredictionError:
            raise
        except Exception as e:
            raise ConditionPredictionError(ErrorCode.ConditionSelectorError, f'无法获取反应条件: {e}')
            # return []
    
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
    
    def get_condition_labels(self, path:str):
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
