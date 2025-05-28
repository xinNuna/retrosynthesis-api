#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
快速测试：验证 Flask 运行时的环境问题
"""

import os
import sys

# 设置路径，模拟 Flask 应用的环境
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from ReactionData import ReactionData, Reactant
from API.ConditionAPI import PredictReactionConditions
from ErrorCodes import ErrorCode

def test_in_flask_context():
    """在模拟 Flask 环境中测试条件预测"""
    print("=== 在模拟 Flask 环境中测试条件预测 ===")
    print(f"当前工作目录: {os.getcwd()}")
    print(f"脚本目录: {CURRENT_DIR}")
    
    # 创建测试数据（与您的 Flask 请求相同）
    reaction_data_list = []
    
    # 创建反应物
    reactants = [
        Reactant(smiles="OB(O)C1=C2C=CC=CC2=CC=C1"),
        Reactant(smiles="BrC1=CC=CC=C1C=O")
    ]
    
    # 创建反应数据
    reaction_data = ReactionData(
        product="O=CC1=CC=CC=C1C1=C2C=CC=CC2=CC=C1",
        reactants=reactants
    )
    
    reaction_data_list.append(reaction_data)
    
    print(f"测试数据准备完成，包含 {len(reaction_data_list)} 个反应")
    
    # 调用 API
    try:
        print("调用 PredictReactionConditions...")
        result_code = PredictReactionConditions(reaction_data_list, topNResults=1)
        
        if result_code == ErrorCode.SUCCESS:
            print("✓ 调用成功！")
            print(f"预测结果数量: {len(reaction_data_list[0].conditions)}")
            for i, condition in enumerate(reaction_data_list[0].conditions):
                print(f"条件 {i+1}:")
                print(f"  催化剂: {condition.catalyst}")
                print(f"  溶剂: {condition.solvents}")
                print(f"  试剂: {condition.reagents}")
                print(f"  温度: {condition.temperature}")
        else:
            print(f"✗ 调用失败，错误码: {result_code.name} ({result_code.value})")
            
    except Exception as e:
        print(f"✗ 调用过程中出现异常: {e}")
        import traceback
        print(f"详细错误: {traceback.format_exc()}")

if __name__ == "__main__":
    test_in_flask_context()