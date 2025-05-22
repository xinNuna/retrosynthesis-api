"""
产率预测API简单测试程序
"""
import os
import sys
import time

# 设置导入路径
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.dirname(CURRENT_DIR)  # 指向 code 文件夹
sys.path.append(CODE_DIR)

# 导入需要的模块
from API.YieldPredictAPI import PredictReactionYield
from ReactionData import ReactionData, Reactant, ReactionCondition
from ErrorCodes import ErrorCode

def run_test(name, reaction):
    """运行单个测试"""
    print(f"\n测试: {name}")
    print("-" * 30)
    
    # 记录开始时间
    start_time = time.time()
    
    # 调用产率预测API
    error_code = PredictReactionYield(reaction)
    
    # 记录结束时间
    elapsed_time = time.time() - start_time
    
    # 输出结果
    if error_code == ErrorCode.SUCCESS:
        print(f"✓ 预测成功！")
        print(f"  执行时间: {elapsed_time:.3f}秒")
        if reaction.yields:
            predicted_yield = reaction.yields[0].value
            print(f"  预测产率: {predicted_yield:.2f}%")
        else:
            print("  错误: 没有返回产率结果")
    else:
        print(f"✗ 预测失败!")
        print(f"  错误代码: {error_code.name}")

def main():
    """主测试函数"""
    print("产率预测API测试")
    print("=" * 30)
    
    # 测试1: 苯乙酮氧化反应
    test1 = ReactionData(
        product="CC(=O)c1ccccc1",  # 产物: 苯乙酮
        reactants=[
            Reactant(smiles="CC(O)c1ccccc1"),  # 反应物1: 苯乙醇
            Reactant(smiles="O=Ic1ccccc1")     # 反应物2: 碘苯甲醛
        ],
        conditions=[
            ReactionCondition(
                catalyst=["CC1(C)CCCC(C)(C)N1[O]"],  # TEMPO催化剂
                solvents=["O"],  # 水
                reagents=["CCCCCCCCCCCCOS(=O)(=O)[O-]", "[Na+]", "[Li+]", "[OH-]", "[Br-]", "[K+]"]  
            )
        ]
    )
    run_test("苯乙醇氧化反应", test1)
    
    # 测试2: Suzuki偶联反应
    test2 = ReactionData(
        product="c1ccc(-c2ccccc2)cc1",  # 产物: 联苯
        reactants=[
            Reactant(smiles="c1ccccc1B(O)O"),    # 反应物1: 苯硼酸
            Reactant(smiles="Brc1ccccc1")        # 反应物2: 溴苯
        ],
        conditions=[
            ReactionCondition(
                catalyst=["[Pd]"],  # 钯催化剂
                solvents=["O", "CC(=O)OC"],  # 水和乙酸乙酯
                reagents=["[K+][OH-]"]  # 氢氧化钾
            )
        ]
    )
    run_test("Suzuki偶联反应", test2)
    
    # 测试3: Wittig反应
    test3 = ReactionData(
        product="CC(/C=C/c1ccccc1)=O",  # 产物: 4-苯基-3-丁烯-2-酮
        reactants=[
            Reactant(smiles="CC(=O)C[P+](c1ccccc1)(c1ccccc1)c1ccccc1"),  # 反应物1: Wittig试剂
            Reactant(smiles="O=Cc1ccccc1")         # 反应物2: 苯甲醛
        ],
        conditions=[
            ReactionCondition(
                catalyst=[],
                solvents=["CN1CCCC1=O"],  # N-甲基吡咯烷酮
                reagents=["[Na+][O-]CC"]   # 乙醇钠
            )
        ]
    )
    run_test("Wittig反应", test3)
    
    # 测试4: 酰胺化反应
    test4 = ReactionData(
        product="O=C(NCc1ccccc1)c1ccccc1",  # 产物: N-苄基苯甲酰胺
        reactants=[
            Reactant(smiles="O=C(O)c1ccccc1"),    # 反应物1: 苯甲酸
            Reactant(smiles="NCc1ccccc1")         # 反应物2: 苄胺
        ],
        conditions=[
            ReactionCondition(
                catalyst=[],
                solvents=["Cc1ccccc1"],  # 甲苯
                reagents=["O=C1OCCO1", "N#CC(C)(C)C"]  # 1,3-二氧戊环、叔丁基异腈
            )
        ]
    )
    run_test("酰胺化反应", test4)
    
    # 测试5: 烷基化反应
    test5 = ReactionData(
        product="CCN(CC)Cc1ccccc1",  # 产物: N,N-二乙基苄胺
        reactants=[
            Reactant(smiles="CCN(CC)C"),  # 反应物1: 三乙胺
            Reactant(smiles="BrCc1ccccc1")  # 反应物2: 溴化苄
        ],
        conditions=[
            ReactionCondition(
                catalyst=[],
                solvents=["CN(C)C=O"],  # DMF
                reagents=["[K+][OH-]"]  # 氢氧化钾
            )
        ]
    )
    run_test("烷基化反应", test5)

    print("\n测试完成。")


if __name__ == "__main__":
    main()
