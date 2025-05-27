import sys
import os
import pandas as pd
from typing import List
from tqdm import tqdm
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.dirname(CURRENT_DIR)
sys.path.append(CODE_DIR)

from ErrorCodes import ErrorCode
from ReactionData import Reactant, ReactionData
from API import SemiTemplateRetrosynthesis

import re

def normalize_smiles_set(smiles_str: str) -> set:
    """
    将反应物字符串拆分成集合，支持用 ';' 和 '.' 两种分隔符，
    并去除空字符串。
    """
    if not isinstance(smiles_str, str) or smiles_str.strip() == "":
        return set()
    # 用正则同时拆分 ; 和 .
    parts = re.split(r'[.;]', smiles_str)
    # 去除空字符串和多余空白，得到干净的smiles集合
    return set(p.strip() for p in parts if p.strip() != "")

def evaluate_top5_accuracy(csv_path: str) -> float:
    data = pd.read_excel(csv_path)
    total = 0
    correct = 0

    for idx, row in tqdm(data.iterrows(), total=len(data)):
        product = row["products"]
        true_reactants = normalize_smiles_set(row["reactants"])
        results: List[ReactionData] = []

        return_code = SemiTemplateRetrosynthesis(
            productSmiles=product,
            reactionResultsList=results,
            topNResults=5,
            rxnClass=None,
            tryAllClasses=True,
            maxSteps=9,
            beamSize=10
        )

        if return_code != ErrorCode.SUCCESS:
            continue

        total += 1
        for result in results:
            pred_reactants = set([r.smiles for r in result.reactants])
            if pred_reactants == true_reactants:
                correct += 1
                break

    accuracy = correct / total if total > 0 else 0.0
    print(f"\n共评估样本数: {total}")
    print(f"Top-5命中数: {correct}")
    print(f"Top-5准确率: {accuracy:.2%}")
    return accuracy


if __name__ == "__main__":
    file_path = "data_retro/reactions_flat.xlsx"  # 替换为你的CSV路径
    evaluate_top5_accuracy(file_path)
