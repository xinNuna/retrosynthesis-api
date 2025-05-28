#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
详细调试测试：检查预测结果的详细结构
"""

import os
import sys

# 设置路径
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

# 添加反应条件预测模块目录到路径
CONDITION_PREDICTOR_DIR = os.path.join(CURRENT_DIR, 'Condition_Predictor_main')
if CONDITION_PREDICTOR_DIR not in sys.path:
    sys.path.insert(0, CONDITION_PREDICTOR_DIR)

def test_reaction_predictor_directly():
    """直接测试 ReactionPredictor 查看原始返回结果"""
    print("=== 直接测试 ReactionPredictor ===")
    
    try:
        # 保存当前工作目录
        original_cwd = os.getcwd()
        
        # 切换到反应条件预测器目录
        os.chdir(CONDITION_PREDICTOR_DIR)
        print(f"切换到工作目录: {os.getcwd()}")
        
        # 导入并调用预测器
        from reaction_predictor import ReactionPredictor
        
        # 构建测试反应
        reaction_smiles_list = [
            'OB(O)C1=C2C=CC=CC2=CC=C1.BrC1=CC=CC=C1C=O>>O=CC1=CC=CC=C1C1=C2C=CC=CC2=CC=C1'
        ]
        
        print(f"测试反应: {reaction_smiles_list[0]}")
        
        # 调用预测器
        predictions = ReactionPredictor(reaction_smiles_list, top_k=5)
        
        print(f"预测结果类型: {type(predictions)}")
        print(f"预测结果属性: {dir(predictions)}")
        
        if hasattr(predictions, 'condition_prediction'):
            print(f"condition_prediction类型: {type(predictions.condition_prediction)}")
            print(f"condition_prediction长度: {len(predictions.condition_prediction)}")
            
            if predictions.condition_prediction:
                first_reaction_conditions = predictions.condition_prediction[0]
                print(f"第一个反应的条件数量: {len(first_reaction_conditions)}")
                
                for i, condition in enumerate(first_reaction_conditions):
                    print(f"\n--- 条件 {i+1} ---")
                    print(f"条件对象类型: {type(condition)}")
                    print(f"条件对象属性: {dir(condition)}")
                    
                    # 检查每个属性
                    for attr_name in ['Catalyst', 'Solvents', 'Reagents', 'Temperature']:
                        if hasattr(condition, attr_name):
                            attr_value = getattr(condition, attr_name)
                            print(f"  {attr_name}: {attr_value} (类型: {type(attr_value)})")
                            
                            # 如果是列表，显示详细内容
                            if isinstance(attr_value, list):
                                print(f"    列表长度: {len(attr_value)}")
                                for j, item in enumerate(attr_value):
                                    print(f"    [{j}]: {item} (类型: {type(item)})")
                        else:
                            print(f"  {attr_name}: 属性不存在")
                    
                    # 检查是否有其他可能的属性名
                    print("  所有属性:")
                    for attr in dir(condition):
                        if not attr.startswith('_'):
                            try:
                                value = getattr(condition, attr)
                                if not callable(value):
                                    print(f"    {attr}: {value}")
                            except:
                                pass
                                
                    # 只显示前2个条件的详细信息
                    if i >= 1:
                        break
        else:
            print("预测结果没有 condition_prediction 属性")
            
    except Exception as e:
        print(f"测试失败: {e}")
        import traceback
        print(f"详细错误: {traceback.format_exc()}")
    finally:
        # 恢复原工作目录
        try:
            os.chdir(original_cwd)
        except:
            pass

def test_api_processing():
    """测试 API 的处理逻辑"""
    print("\n=== 测试 API 处理逻辑 ===")
    
    # 导入模拟条件对象
    import dataclasses
    from typing import List
    
    # 创建一个模拟的条件对象
    @dataclasses.dataclass
    class MockCondition:
        Catalyst: List[str]
        Reagents: List[str] 
        Solvents: List[str]
        Temperature: str
    
    # 创建测试条件
    test_condition = MockCondition(
        Catalyst=['Pd(PPh3)4'],
        Reagents=['K2CO3', 'H2O'],
        Solvents=['DMF', 'THF'], 
        Temperature='80'
    )
    
    print(f"测试条件对象: {test_condition}")
    
    # 模拟 API 中的处理逻辑
    catalysts = []
    if hasattr(test_condition, 'Catalyst') and test_condition.Catalyst:
        if isinstance(test_condition.Catalyst, list):
            catalysts = [cat for cat in test_condition.Catalyst if cat not in ["nan", "None", None, ""]]
        elif test_condition.Catalyst not in ["nan", "None", None, ""]:
            catalysts = [test_condition.Catalyst]
    
    solvents = []
    if hasattr(test_condition, 'Solvents') and test_condition.Solvents:
        solvents = [s for s in test_condition.Solvents if s not in ["nan", "None", None, ""]]
    
    reagents = []
    if hasattr(test_condition, 'Reagents') and test_condition.Reagents:
        reagents = [r for r in test_condition.Reagents if r not in ["nan", "None", None, ""]]
    
    temperature = None
    if hasattr(test_condition, 'Temperature') and test_condition.Temperature:
        temp_str = str(test_condition.Temperature).strip()
        if "°C" in temp_str or "K" in temp_str or "℃" in temp_str or "F" in temp_str:
            temperature = temp_str
        else:
            temperature = f"{temp_str} °C"
    
    print(f"处理后结果:")
    print(f"  催化剂: {catalysts}")
    print(f"  溶剂: {solvents}")
    print(f"  试剂: {reagents}")
    print(f"  温度: {temperature}")

if __name__ == "__main__":
    test_reaction_predictor_directly()
    test_api_processing()