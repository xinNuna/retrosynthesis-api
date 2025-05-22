"""
反应条件预测API测试程序
"""
import sys
import os
from typing import List

# 设置导入路径
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.dirname(CURRENT_DIR)  # 指向 code 文件夹
sys.path.append(CODE_DIR)

# 导入所需模块
from ErrorCodes import ErrorCode
from ReactionData import Reactant, ReactionData, ReactionCondition
from API import PredictReactionConditions

def test_condition_prediction(reaction_smiles_list: List[str]) -> None:
    """
    测试反应条件预测API
    
    参数:
        reaction_smiles_list: 待测试的反应SMILES列表
    """
    print("\n===== 测试反应条件预测 =====")
    
    # 创建ReactionData对象列表
    reaction_data_list = []
    
    for i, reaction_smiles in enumerate(reaction_smiles_list):
        # 分割反应物和产物
        if ">>" not in reaction_smiles:
            print(f"错误: 反应SMILES {i+1} 格式不正确，缺少'>>'")
            continue
            
        parts = reaction_smiles.split(">>")
        if len(parts) != 2:
            print(f"错误: 反应SMILES {i+1} 格式不正确")
            continue
            
        reactants_smiles, product_smiles = parts
        
        # 将反应物SMILES分割为单独的分子
        reactant_molecules = reactants_smiles.split(".")
        
        # 创建反应物对象
        reactants = []
        for j, r_smiles in enumerate(reactant_molecules):
            reactant = Reactant(
                smiles=r_smiles,
                rank=j+1
            )
            reactants.append(reactant)
            
        # 创建反应对象
        reaction_data = ReactionData(
            product=product_smiles,
            reactants=reactants
        )
        
        reaction_data_list.append(reaction_data)
        
        print(f"反应 {i+1}:")
        print(f"  产物: {product_smiles}")
        print(f"  反应物: {reactants_smiles}")
    
    # 调用条件预测API
    print("\n调用条件预测API...")
    return_code = PredictReactionConditions(reaction_data_list, topNResults=5)
    
    # 检查返回结果
    if return_code == ErrorCode.SUCCESS:
        print(f"成功: {return_code.name} (代码: {return_code.value})")
        print(f"获得 {len(reaction_data_list)} 个反应的条件预测")
        
        # 打印预测结果
        for i, reaction in enumerate(reaction_data_list):
            print(f"\n反应 {i+1} 的条件预测结果:")
            
            if not reaction.conditions:
                print("  没有预测到条件")
                continue
                
            for j, condition in enumerate(reaction.conditions):
                print(f"  条件 {j+1}:")
                print(f"    催化剂: {', '.join(condition.catalyst) if condition.catalyst else '无'}")
                print(f"    溶剂: {', '.join(condition.solvents) if condition.solvents else '无'}")
                print(f"    试剂: {', '.join(condition.reagents) if condition.reagents else '无'}")
                
                if condition.temperature:
                    print(f"    温度: {condition.temperature}")
                else:
                    print(f"    温度: 未知")
    else:
        print(f"错误: {return_code.name} (代码: {return_code.value})")
    
    print("="*50)

if __name__ == "__main__":
    # 测试反应SMILES
    reaction_smiles = [
        'OB(O)C1=C2C=CC=CC2=CC=C1.BrC1=CC=CC=C1C=O>>O=CC1=CC=CC=C1C1=C2C=CC=CC2=CC=C1',
        'BrC1=CN=CS1.CC1=CC=CC(=C1)B(O)O>>CC1=CC=CC(=C1)C1=CN=CS1'
    ]
    
    test_condition_prediction(reaction_smiles)
