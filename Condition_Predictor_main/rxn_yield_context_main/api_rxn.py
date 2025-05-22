# from rxn_yield_context.evaluate_model.eval_utils import ReactionContextPredictor
from rxn_yield_context_main.rxn_yield_context.evaluate_model.eval_utils import ReactionContextPredictor
from pathlib import Path
from typing import List
import dataclasses
from enum import Enum
import os

@dataclasses.dataclass
class InputSmiles:
    reactionSMILES: List[str]

@dataclasses.dataclass
class Conditions:
    Reagents: str
    Solvents: str
    Temperature: str


@dataclasses.dataclass
class Predictions:
    predictions: List[Conditions]

class ErrorCode(Enum):
    PredictorInitializedError = 3001 # 预测器初始化出现错误
    RecommandationError = 3002 # 获取推荐错误 

class RxnPredictionError(Exception):
    def __init__(self, error_code: ErrorCode, message: str):
        self.error_code = error_code
        self.message = message
        super().__init__(f"[错误代码: {self.error_code.value}，{self.message}")

def get_rxn_predictions(inputData: list, top_k: int=5) -> Predictions:
    current_dir = os.path.dirname(os.path.abspath(__file__))
    class_path = os.path.join(current_dir, 'data', 'reaxys_output')
    candidate_generation_model_path = os.path.join(current_dir, 'save_models', 'test_10R_first_local_10', 'multitask_model_epoch-80.checkpoint')
    ranking_model_path = os.path.join(current_dir, 'save_models', 'test_10R_second_7', 'rxn_model_relevance_listwise_morgan_epoch-80.checkpoint')

    chemicalFormula = [formula.strip() for formula in inputData]

    try:
        #  初始化 ReactionContextPredictor
        predictor = ReactionContextPredictor(
            class_path,
            candidate_generation_model_path,
            ranking_model_path,
            # 'rxn_yield_context_main/data/reaxys_output',  # class_data_path
            # 'rxn_yield_context_main/save_models/test_10R_first_local_10/multitask_model_epoch-80.checkpoint',  # candidate_generation_model_path
            # 'rxn_yield_context_main/save_models/test_10R_second_7/rxn_model_relevance_listwise_morgan_epoch-80.checkpoint',  # ranking_model_path
            # 'data/reaxys_output',  # class_data_path
            # 'save_models/test_10R_first_local_10/multitask_model_epoch-80.checkpoint',  # candidate_generation_model_path
            # 'save_models/test_10R_second_7/rxn_model_relevance_listwise_morgan_epoch-80.checkpoint',  # ranking_model_path
            cutoff_solv=0.3,
            cutoff_reag=0.3,
            verbose=False
        )
    except RxnPredictionError:
        raise
    except Exception as e:
        raise RxnPredictionError(ErrorCode.PredictorInitializedError, f'无法初始化预测器{e}')
    
    try:
        # 使用全局预测器进行预测
        rxn_smiles_list = chemicalFormula
        results = predictor.recommend_reaction_context(rxn_smiles_list, max_display=top_k) # 可以显示的推荐数量，最大为10
            
        # 转换结果格式
        recommendationsPerCondition = []
        recommendationsPerReaction = []
        for i in range(len(results)):
            for _, row in results[i].iterrows():
                recommendationsPerCondition.append(
                    Conditions(
                        Reagents=row.get('Reagent(s)'),
                        Solvents=row.get('Solvent(s)'),
                        Temperature=row.get('Temperature')
                    )
                )
            recommendationsPerReaction.append(recommendationsPerCondition[:top_k]) # 只返回前5条推荐
        return Predictions(predictions=recommendationsPerReaction)
    except RxnPredictionError:
        raise
    except Exception as e:
        raise RxnPredictionError(ErrorCode.RecommandationError, f'预测反应条件失败: {e}')



if __name__ == "__main__":
    input = ['OB(O)C1=CC=NC=C1.ClC1=NC=C(I)C(OC2=CC(Br)=CC=N2)=N1>>ClC1=NC=C(C(OC2=NC=CC(Br)=C2)=N1)C1=CC=NC=C1', 'CC(=O)OC1N(CC2=CC=CC=C2)C(=O)C2=C1C=CC=C2.C[Si](C)(C)CC=C>>C=CCC1N(CC2=CC=CC=C2)C(=O)C2=C1C=CC=C2']

    print("开始预测")

    results = get_rxn_predictions(input)
    if results.predictions:
        for j, rec in enumerate(results.predictions):
            print("\n\n")
            for i, cond in enumerate(rec):
                print(f"反应 {j+1} 的top-5条件推荐如下所示 ") 
                print(f"Rank {i+1}:")
                print(f"  试剂: {cond.Reagents}")
                print(f"  溶剂: {cond.Solvents}")
                print(f"  温度: {cond.Temperature} ℃\n")
    else:
        pass
    