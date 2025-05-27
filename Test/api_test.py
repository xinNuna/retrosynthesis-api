"""
逆合成API测试程序 - 测试半模板法
"""
import sys
import os
# 设置导入路径
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.dirname(CURRENT_DIR)  # 指向 code 文件夹
sys.path.append(CODE_DIR)
# 导入所需模块
from ErrorCodes import ErrorCode
from ReactionData import Reactant, ReactionData
from API import SemiTemplateRetrosynthesis

def test_semi_template(product_smiles: str, rxn_class: int = None, try_all_classes: bool = True) -> ErrorCode:
    """
    测试半模板逆合成预测
    参数:
        product_smiles: 产物SMILES
        rxn_class: 反应类别
        try_all_classes: 是否尝试所有反应类别
    返回:
        状态码
    """
    results = []
    print(f"\n===== 测试半模板逆合成 =====")
    print(f"产物: {product_smiles}")
    print(f"参数: rxnClass={rxn_class}, tryAllClasses={try_all_classes}")
    # 调用半模板逆合成预测函数
    returnCode = SemiTemplateRetrosynthesis(
        productSmiles=product_smiles,
        reactionResultsList=results,
        topNResults=5,
        rxnClass=rxn_class,
        tryAllClasses=try_all_classes,
        maxSteps=9,
        beamSize=10
    )
    
    # 检查返回结果
    print_results(returnCode, results)
    return returnCode


def print_results(returnCode: ErrorCode, results):
    """打印预测结果"""
    if returnCode == ErrorCode.SUCCESS:
        print(f"成功: {returnCode.name} (代码: {returnCode.value})")
        print(f"获得 {len(results)} 个预测结果")
        
        for i, result in enumerate(results):
            print(f"\n结果 {i+1}:")
            print(f"产物: {result.product}")
            print(f"反应物数量: {len(result.reactants)}")
            
            # 打印每个反应物的信息
            for j, reactant in enumerate(result.reactants):
                print(f"  反应物 {j+1}:")
                print(f"    SMILES: {reactant.smiles}")
                if reactant.rxnClass is not None:
                    print(f"    反应类别: {reactant.rxnClass}")
                print(f"    排名: {reactant.rank}")
    else:
        print(f"错误: {returnCode.name} (代码: {returnCode.value})")
    
    print("="*50)


def main():
    """主函数，运行所有测试"""
    #测试用例 - 咖啡因
    caffeine = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    
    #测试用例 - 可能生成多个反应物的化合物
    benzyl_phenyl_ether = "c1ccccc1OCc2ccccc2"  # 苄基苯基醚
    biphenyl = "c1ccccc1-c2ccccc2"  # 联苯
    ethyl_benzoate = "CCOC(=O)c1ccccc1"  # 苯甲酸乙酯
    
    print("\n====================================================")
    print("          开始测试半模板法逆合成预测")
    print("====================================================")
    
    # 测试1: 半模板法 - 指定反应类别
    #test_semi_template(caffeine, rxn_class=7, try_all_classes=False)
    
    # 测试2: 半模板法 - 尝试所有反应类别
    test_semi_template(biphenyl, rxn_class=None, try_all_classes=True)
    
    # 测试3: 半模板法 - 错误情况(未指定反应类别且不尝试所有类别)
    #test_semi_template(benzyl_phenyl_ether, rxn_class=None, try_all_classes=False)


if __name__ == "__main__":
    main()
