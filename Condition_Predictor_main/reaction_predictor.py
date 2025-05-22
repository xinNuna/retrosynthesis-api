import warnings
import os

# 关闭进度条
os.environ['TQDM_DISABLE'] = '1'

# 忽略pytorch的参数警告
warnings.filterwarnings(
    "ignore",
    category=FutureWarning,
    message="You are using `torch.load` with `weights_only=False`"
)

from typing import List
from collections import OrderedDict
from decimal import Decimal
import ast
import sys
import os
import dataclasses

current_dir = os.path.dirname(os.path.abspath(__file__)) 
rxn_main_dir = os.path.join(current_dir, 'rxn_yield_context_main')
rcs_main_dir = os.path.join(current_dir, 'Reaction_Condition_Selector_main')
if rxn_main_dir not in sys.path:
    sys.path.insert(0, rcs_main_dir) 
    sys.path.append(rxn_main_dir)

from rxn_yield_context_main.api_rxn import get_rxn_predictions
from Reaction_Condition_Selector_main.api_rcs import get_rcs_prediction

# 自定义反应条件数据类
@dataclasses.dataclass
class Conditions:
     Catalyst: str
     Reagents: List[str]
     Solvents: List[str]
     Temperature: str

@dataclasses.dataclass
class Condition_Prediction:
     condition_prediction: List[Conditions]

def ReactionPredictor(reaction_smiles: List[str], top_k: int=5): # 通过top_k控制预测数量
    conditions = []
    temperature = []

    rcs_predictions = get_rcs_prediction(reaction_smiles, top_k) # 获取预测的溶剂，催化剂和试剂

    rxn_predictions = get_rxn_predictions(reaction_smiles, top_k) # 获取温度
    for j, rxn_conditions in enumerate(rxn_predictions.predictions):
        temp = []
        for i, condition in enumerate(rxn_conditions):
            temp.append(condition.Temperature)
        avg_temp = sum( float(t) for t in temp) / len(temp)
        temperature.append(Decimal(avg_temp).quantize(Decimal('0.0')))
    
    for key, value in rcs_predictions.items():
        cond = []
        for i in range(len(value)):
            temp = ast.literal_eval(value[i])
            cond.append(
                Conditions(
                    Catalyst=[temp[0]],
                    Reagents = [temp[1],temp[2]],
                    Solvents= [temp[3],temp[4],temp[5]],
                    Temperature=temperature[int(key)-1]
                )
            )
        conditions.append(cond)
    return Condition_Prediction(conditions)

if __name__ == "__main__":
    smiles = ['OB(O)C1=C2C=CC=CC2=CC=C1.BrC1=CC=CC=C1C=O>>O=CC1=CC=CC=C1C1=C2C=CC=CC2=CC=C1', 'BrC1=CN=CS1.CC1=CC=CC(=C1)B(O)O>>CC1=CC=CC(=C1)C1=CN=CS1']
    pred_cond = ReactionPredictor(smiles)
    for j, rec in enumerate(pred_cond.condition_prediction):
        for i, cond in enumerate(rec):
            print(f"反应 {j+1} 的top-5条件推荐如下所示 ") 
            print(f"Rank {i+1}:")
            print(f"催化剂: {cond.Catalyst}")
            print(f"  试剂: {cond.Reagents}")
            print(f"  溶剂: {cond.Solvents}")
        print(f"参考温度: {cond.Temperature} ℃\n")
    # print(pred_cond)
    # print("预测的条件：", cond)
    # print("所有的条件: ", cond_all)