"""
模板法逆合成预测测试程序（无原子映射版本）
"""
import sys
import os

# 设置导入路径
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.dirname(CURRENT_DIR)  # 指向 code 文件夹
sys.path.append(CODE_DIR)

from API import TemplateRetrosynthesis
from ReactionData import ReactionData
from ErrorCodes import ErrorCode


def main():
    """测试模板法逆合成预测（不带原子映射的SMILES）"""
    
    # 测试SMILES - 不带原子映射
    smiles = 'Cn1c(C)nc2cc(Cl)c(Cl)cc21'  # 与之前相同的分子，但没有原子映射
    print(f"测试产物（无原子映射）: {smiles}")
    print("-" * 50)
    
    # 创建结果列表
    results = []
    
    # 调用API
    status = TemplateRetrosynthesis(
        productSmiles=smiles,
        reactionResultsList=results,
        topNResults=5
    )
    
    # 输出结果
    if status == ErrorCode.SUCCESS:
        print(f"\n预测成功! 找到 {len(results)} 个反应路径:")
        
        for i, reaction in enumerate(results, 1):
            print(f"\n路径 {i}:")
            print(f"  产物: {reaction.product}")
            print(f"  反应物:")
            for j, reactant in enumerate(reaction.reactants, 1):
                print(f"    {j}. {reactant.smiles}")
    else:
        print(f"\n预测失败! 错误代码: {status.name} ({status.value})")


def test_multiple_molecules():
    """测试多个分子"""
    print("\n\n=== 测试多个分子（无原子映射） ===")
    
    test_molecules = [
        'CC(=O)Oc1ccccc1C(=O)O',  # 阿司匹林
        'CC(C)Cc1ccc(C(C)C(=O)O)cc1',  # 布洛芬类似物
        'O=C(O)c1ccc2ccccc2c1',  # 萘甲酸
    ]
    
    for idx, smiles in enumerate(test_molecules, 1):
        print(f"\n分子 {idx}: {smiles}")
        results = []
        
        status = TemplateRetrosynthesis(
            productSmiles=smiles,
            reactionResultsList=results,
            topNResults=3
        )
        
        if status == ErrorCode.SUCCESS:
            print(f"成功! 找到 {len(results)} 个路径")
            for i, reaction in enumerate(results, 1):
                print(f"  路径{i}: {' + '.join([r.smiles for r in reaction.reactants])}")
        else:
            print(f"失败: {status.name}")


if __name__ == "__main__":
    main()
    test_multiple_molecules()
