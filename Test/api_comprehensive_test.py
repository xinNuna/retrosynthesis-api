"""
API综合测试程序：逆合成预测、条件预测、产率预测的完整流程
支持不带原子映射的SMILES输入
"""
import os
import sys
import time
from typing import List

# 设置导入路径
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.dirname(CURRENT_DIR)  # 指向 code 文件夹
sys.path.append(CODE_DIR)

# 导入需要的模块
from API import TemplateRetrosynthesis
from API import SemiTemplateRetrosynthesis
from API import PredictReactionConditions
from API import PredictReactionYield
from ReactionData import ReactionData, Reactant, ReactionCondition
from ErrorCodes import ErrorCode

def print_header(title, char="="):
    """打印标题"""
    print(f"\n{title}")
    print(char * len(title))

def print_reaction_data(reaction: ReactionData, full_detail: bool = False):
    """打印单条反应数据"""
    print(f"产物: {reaction.product}")
    
    print("\n反应物:")
    for i, reactant in enumerate(reaction.reactants):
        print(f"  {i+1}. {reactant.smiles}")
    
    if reaction.conditions and (full_detail or len(reaction.conditions) <= 3):
        print("\n反应条件:")
        for i, condition in enumerate(reaction.conditions[:3 if not full_detail else None]):
            print(f"  条件 {i+1}:")
            # 添加字段存在性检查
            print(f"    催化剂: {condition.catalyst if condition.catalyst else '无'}")
            print(f"    溶剂: {', '.join(condition.solvents) if condition.solvents else '无'}")
            print(f"    试剂: {', '.join(condition.reagents) if condition.reagents else '无'}")
            print(f"    温度: {condition.temperature if condition.temperature else '无'}")
    elif reaction.conditions:
        print(f"\n反应条件: {len(reaction.conditions)}个条件可用 (使用 -d 查看详情)")
    
    if reaction.yields:
        print("\n产率预测:")
        for i, yield_data in enumerate(reaction.yields):
            print(f"  预测产率: {yield_data.value:.2f}%")


def predict_reactants(product_smiles: str) -> List[ReactionData]:
    """预测反应物"""
    print_header("反应物预测", "-")
    
    # 选择逆合成方法
    while True:
        choice = input("选择逆合成预测方法 [1: 模板法 (默认), 2: 半模板法]: ").strip()
        if not choice:
            choice = "1"
        
        if choice in ["1", "2"]:
            method = "template" if choice == "1" else "semi-template"
            break
        else:
            print(" 无效选择，请输入1或2")
    
    # 获取返回结果数量
    topn = get_integer_input("返回结果数量", 5, min_value=1, max_value=20)
    
    # 保存预测结果的列表
    results = []
    
    # 记录开始时间
    start_time = time.time()
    
    print("\n 正在预测反应物，请稍候...")
    
    # 根据选择调用不同的API
    if method == "template":
        # 模板法 - 不需要原子映射，API内部会处理
        error_code = TemplateRetrosynthesis(product_smiles, results, topn)
    else:  # semi-template
        # 半模板法 - 可选择指定反应类别
        use_class = input("是否指定反应类别? [y/n]: ").strip().lower() == 'y'
        
        if not use_class:
            # 不指定类别，尝试所有类别
            error_code = SemiTemplateRetrosynthesis(product_smiles, results, topn)
        else:
            # 指定反应类别
            rxn_class = get_integer_input("反应类别 (1-10)", 1, min_value=1, max_value=10)
            error_code = SemiTemplateRetrosynthesis(
                product_smiles, results, topn, 
                rxn_class=rxn_class, 
                tryAllClasses=False
            )
    
    # 计算用时
    elapsed_time = time.time() - start_time
    
    # 输出结果
    if error_code == ErrorCode.SUCCESS:
        print(f"\n 反应物预测成功！")
        print(f"  预测用时: {elapsed_time:.2f}秒")
        print(f"  找到 {len(results)} 个可能的反应路径")
        
        # 打印预测结果
        for i, reaction in enumerate(results[:min(3, len(results))]):
            print(f"\n结果 {i+1}/{len(results)}:")
            print_reaction_data(reaction)
        
        if len(results) > 3:
            print(f"\n... 还有 {len(results) - 3} 个结果未显示 ...")
        
        return results
    else:
        print(f"\n 反应物预测失败!")
        print(f"  错误代码: {error_code.name}")
        return []


def predict_conditions(reactions: List[ReactionData]) -> bool:
    """预测反应条件"""
    if not reactions:
        print("\n 无法预测反应条件: 没有有效的反应数据")
        return False
    
    print_header("反应条件预测", "-")
    
    # 获取返回条件数量
    topn = get_integer_input("返回条件数量", 3, min_value=1, max_value=10)
    
    # 记录开始时间
    start_time = time.time()
    
    print("\n 正在预测反应条件，请稍候...")
    
    # 调用条件预测API
    error_code = PredictReactionConditions(reactions, topn)
    
    # 计算用时
    elapsed_time = time.time() - start_time
    
    # 输出结果
    if error_code == ErrorCode.SUCCESS:
        print(f"\n 反应条件预测成功！")
        print(f"  预测用时: {elapsed_time:.2f}秒")
        
        # 打印条件预测结果
        for i, reaction in enumerate(reactions[:min(3, len(reactions))]):
            if i > 0:
                print("\n" + "-" * 30)
            print(f"\n反应 {i+1} 的条件预测结果:")
            if reaction.conditions:
                for j, condition in enumerate(reaction.conditions[:min(3, len(reaction.conditions))]):
                    print(f"  条件 {j+1}:")
                    if condition.catalyst:
                        print(f"    催化剂: {condition.catalyst}")
                    if condition.solvents:
                        print(f"    溶剂: {', '.join(condition.solvents)}")
                    if condition.reagents:
                        print(f"    试剂: {', '.join(condition.reagents)}")
                    if condition.temperature:
                        print(f"    温度: {condition.temperature}")
                
                if len(reaction.conditions) > 3:
                    print(f"  ... 还有 {len(reaction.conditions) - 3} 个条件未显示 ...")
            else:
                print("  未预测到条件")
        
        return True
    else:
        print(f"\n 反应条件预测失败!")
        print(f"  错误代码: {error_code.name}")
        return False


def predict_yield(reaction: ReactionData) -> bool:
    """预测反应产率"""
    if not reaction or not reaction.reactants:
        print("\n 无法预测产率: 没有有效的反应数据")
        return False
    
    if not reaction.conditions:
        print("\n 无法预测产率: 没有有效的反应条件")
        return False
    
    print_header("产率预测", "-")
    print("\n 正在预测产率，请稍候...")
    
    # 记录开始时间
    start_time = time.time()
    
    # 调用产率预测API
    error_code = PredictReactionYield(reaction)
    
    # 计算用时
    elapsed_time = time.time() - start_time
    
    # 输出结果
    if error_code == ErrorCode.SUCCESS:
        print(f"\n 产率预测成功！")
        print(f"  预测用时: {elapsed_time:.2f}秒")
        
        if reaction.yields:
            predicted_yield = reaction.yields[0].value
            print(f"  预测产率: {predicted_yield:.2f}%")
        else:
            print("   错误: 没有返回产率结果")
        
        return True
    else:
        print(f"\n 产率预测失败!")
        print(f"  错误代码: {error_code.name}")
        return False


def get_integer_input(prompt_text, default_value=None, min_value=None, max_value=None):
    """获取整数输入，支持范围验证"""
    while True:
        prompt = f"{prompt_text} "
        if default_value is not None:
            prompt += f"[默认: {default_value}]: "
        else:
            prompt += ": "
            
        user_input = input(prompt).strip()
        
        if not user_input and default_value is not None:
            return default_value
            
        try:
            value = int(user_input)
            
            # 范围验证
            if min_value is not None and value < min_value:
                print(f" 输入值必须大于或等于 {min_value}")
                continue
                
            if max_value is not None and value > max_value:
                print(f" 输入值必须小于或等于 {max_value}")
                continue
                
            return value
        except ValueError:
            print(" 请输入有效的整数")


def get_reaction_choice(reactions, prompt_text="选择要使用的反应"):
    """从多个反应中选择一个"""
    if len(reactions) == 1:
        return reactions[0]
        
    print(f"\n{prompt_text}:")
    for i, reaction in enumerate(reactions):
        # 仅显示概述
        print(f"  [{i+1}] 产物: {reaction.product}")
        reactants_text = ", ".join([r.smiles[:30] + ('...' if len(r.smiles) > 30 else '') 
                                   for r in reaction.reactants[:2]])
        if len(reaction.reactants) > 2:
            reactants_text += f" ... (共{len(reaction.reactants)}个)"
        print(f"      反应物: {reactants_text}")
    
    # 获取用户选择
    choice = get_integer_input("\n请输入选择", 1, min_value=1, max_value=len(reactions))
    return reactions[choice - 1]


def get_condition_choice(reaction, prompt_text="选择要使用的条件"):
    """从多个条件中选择一个"""
    if len(reaction.conditions) == 1:
        return reaction.conditions[0]
        
    print(f"\n{prompt_text}:")
    for i, condition in enumerate(reaction.conditions):
        # 改进条件摘要显示，确保显示所有可用信息
        summary_parts = []
        
        # 检查催化剂
        if condition.catalyst:
            catalyst = condition.catalyst[:30] + ('...' if len(condition.catalyst) > 30 else '')
            summary_parts.append(f"催化剂: {catalyst}")
        
        # 检查溶剂
        if condition.solvents:
            solvents_str = ", ".join(condition.solvents[:2])
            if len(condition.solvents) > 2:
                solvents_str += f" ... (共{len(condition.solvents)}个)"
            summary_parts.append(f"溶剂: {solvents_str}")
        
        # 检查试剂
        if condition.reagents:
            reagents_str = ", ".join(condition.reagents[:2])
            if len(condition.reagents) > 2:
                reagents_str += f" ... (共{len(condition.reagents)}个)"
            summary_parts.append(f"试剂: {reagents_str}")
        
        # 检查温度
        if condition.temperature:
            summary_parts.append(f"温度: {condition.temperature}")
        
        # 合并所有摘要部分
        condition_text = " | ".join(summary_parts) if summary_parts else "无详细信息"
        print(f"  [{i+1}] {condition_text}")
    
    # 获取用户选择
    choice = get_integer_input("\n请输入选择", 1, min_value=1, max_value=len(reaction.conditions))
    return reaction.conditions[choice - 1]


def show_example_smiles():
    """显示示例SMILES"""
    examples = [
        ("苯乙酮", "CC(=O)c1ccccc1"),
        ("阿司匹林", "CC(=O)Oc1ccccc1C(=O)O"),
        ("联苯", "c1ccc(-c2ccccc2)cc1"),
        ("苯酚", "Oc1ccccc1"),
        ("对氯苯酚", "Oc1ccc(Cl)cc1")
    ]
    
    print("\n示例SMILES:")
    for name, smiles in examples:
        print(f"  {name}: {smiles}")


def main():
    """主函数"""
    print_header("API综合测试程序: 逆合成 + 条件预测 + 产率预测", "=")
    print("提示: 在任何输入提示处输入 'q' 退出程序，输入 'examples' 显示示例SMILES")
    
    while True:
        # 获取产物SMILES
        product_smiles = input("\n请输入产物SMILES: ").strip()
        
        # 处理特殊命令
        if product_smiles.lower() in ('q', 'quit', 'exit'):
            break
        elif product_smiles.lower() in ('examples', 'ex', 'sample'):
            show_example_smiles()
            continue
            
        if not product_smiles:
            print(" 错误: 产物SMILES不能为空")
            continue
        
        # 步骤1: 预测反应物
        reactions = predict_reactants(product_smiles)
        if not reactions:
            print(" 反应物预测失败，请尝试另一个产物结构或检查SMILES格式")
            continue
        
        # 选择要继续的反应
        print_header("选择反应路径", "-")
        selected_reaction = get_reaction_choice(reactions)
        print("\n选择的反应路径:")
        print_reaction_data(selected_reaction, full_detail=True)
        
        # 步骤2: 为选择的反应预测条件
        if not predict_conditions([selected_reaction]):
            continue
        
        # 步骤3: 为选择的反应和条件预测产率
        if not selected_reaction.conditions:
            print(" 无法预测产率: 未成功预测反应条件")
            continue
            
        # 选择要用于产率预测的条件
        selected_condition = get_condition_choice(selected_reaction)
        
        # 创建只包含选定条件的反应对象
        yield_reaction = ReactionData(
            product=selected_reaction.product,
            reactants=selected_reaction.reactants,
            conditions=[selected_condition]
        )
        
        # 预测产率
        predict_yield(yield_reaction)
        
        # 询问是否继续
        if input("\n是否继续测试其他产物? [y/n]: ").strip().lower() == 'n':
            break
    
    print("\n测试程序已结束。感谢使用！")


if __name__ == "__main__":
    main()
