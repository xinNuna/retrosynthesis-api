import sys
import os
import pandas as pd
from typing import List
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from tqdm import tqdm
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.dirname(CURRENT_DIR)
sys.path.append(CODE_DIR)

from ErrorCodes import ErrorCode
from ReactionData import Reactant, ReactionData
from API import TemplateRetrosynthesis
import re

def smiles_similarity(smiles1: str, smiles2: str) -> float:
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 is None or mol2 is None:
        return 0.0
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def is_similar_set(pred_set: set, true_set: set, threshold: float = 0.85) -> bool:
    """
    判断两个反应物集合是否结构相似。
    要求每个预测分子都在真实集合中至少有一个“足够相似”的匹配。
    """
    if len(pred_set) != len(true_set):
        return False

    matched = set()
    for pred in pred_set:
        for true in true_set:
            if true in matched:
                continue
            sim = smiles_similarity(pred, true)
            if sim >= threshold:
                matched.add(true)
                break
    return len(matched) == len(true_set)


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

        return_code = TemplateRetrosynthesis(
        productSmiles=product,
        reactionResultsList=results,
        topNResults=5
        )

        if return_code != ErrorCode.SUCCESS:
            continue

        total += 1
        for result in results:
            pred_reactants = set([r.smiles for r in result.reactants])
            if is_similar_set(pred_reactants, true_reactants, threshold=0.85):
                correct += 1
                break

    accuracy = correct / total if total > 0 else 0.0
    print(f"\n共评估样本数: {total}")
    print(f"Top-5命中数: {correct}")
    print(f"Top-5准确率: {accuracy:.2%}")
    return accuracy


if __name__ == "__main__":
    file_path = "data_retro/test.xlsx"  # 替换为你的CSV路径
    evaluate_top5_accuracy(file_path)
