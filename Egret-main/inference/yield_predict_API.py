import os
import pandas as pd
import sys
sys.path.append(os.path.join(os.path.abspath(os.path.dirname(__file__)), '..'))
from inference.models import SmilesClassificationModel
import numpy as np
import torch
from argparse import ArgumentParser

class SmilesClassificationModelSelf(SmilesClassificationModel):
    def load_and_cache_examples(self, examples, evaluate=False, no_cache=True, multi_label=False, verbose=True, silent=False):
        return super().load_and_cache_examples(examples, evaluate, no_cache, multi_label, verbose, silent)

# 定义 REACTION 类（作为输入结构）
class REACTION:
    def __init__(self, Product, Reactant, Condition):
        self.Product = Product                # 单个产物 SMILES
        self.Reactant = Reactant              # 列表：反应物 SMILES
        self.Condition = Condition            # 列表：条件 SMILES

# Yield 预测 API 类
class YieldPredictorAPI:
    def __init__(self, model_state_path, cuda_device=-1):
        self.model_state_path = model_state_path
        use_cuda = torch.cuda.is_available() if cuda_device != -1 else False
        self.model_class = SmilesClassificationModelSelf(
            "bert",
            self.model_state_path,
            num_labels=1,
            use_cuda=use_cuda,
            cuda_device=cuda_device
        )
        self.yield_mean = 79.29119663076209
        self.yield_std = 18.858441890553195

    def predict_yield(self, inputReaction: REACTION) -> float:
        # 拼接 Reaction SMILES：Reactants + Conditions >> Product
        reactants_and_conditions = ".".join(inputReaction.Reactant + inputReaction.Condition)
        rxn_smiles = f"{reactants_and_conditions}>>{inputReaction.Product}"
        
        # 预测
        y_pred_raw = self.model_class.predict([rxn_smiles])[0]
# 如果 y_pred_raw 是 0 维，直接取
        if isinstance(y_pred_raw, np.ndarray) and y_pred_raw.ndim == 0:
            y_pred = y_pred_raw.item()
# 如果是 1 维
        elif isinstance(y_pred_raw, (list, np.ndarray)):
            y_pred = y_pred_raw[0]
        else:
            raise ValueError(f"Unexpected prediction output: {y_pred_raw}")
        y_pred = y_pred * self.yield_std + self.yield_mean
        y_pred = np.clip(y_pred, 0, 100)
        
        return y_pred

# 示例用法
if __name__ == "__main__":
    # 定义反应
    reaction = REACTION(
        Product="CC(=O)c1ccccc1",
        Reactant=["CC(O)c1ccccc1", "O=Ic1ccccc1"],
        Condition=["CC1(C)CCCC(C)(C)N1[O]", "CCCCCCCCCCCCOS(=O)(=O)[O-]", "[Na+]", "O", "[Li+]", "[OH-]", "[Br-]", "[K+]"]
    )

    # 初始化预测器
    this_path = os.path.abspath(os.path.dirname(__file__))
    model_state_path = os.path.join(this_path, 'yield_prediction_model')
    predictor = YieldPredictorAPI(model_state_path=model_state_path)

    # 预测产率
    predicted_yield = predictor.predict_yield(reaction)
    print(f"Predicted Yield: {predicted_yield:.2f}%")
