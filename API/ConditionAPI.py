"""
反应条件预测 API
"""
import os
import sys
from typing import List

# 获取当前脚本的绝对路径
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)  # 指向 code 文件夹

# 添加根目录到路径以便导入其他模块
sys.path.append(ROOT_DIR)

# 添加反应条件预测模块目录到路径
CONDITION_PREDICTOR_DIR = os.path.join(ROOT_DIR, 'Condition_Predictor_main')
if CONDITION_PREDICTOR_DIR not in sys.path:
    sys.path.insert(0, CONDITION_PREDICTOR_DIR)

# 导入所需的模块
from ErrorCodes import ErrorCode, RetroApiError
from ReactionData import ReactionData, ReactionCondition
from .UtilityTools import SuppressAllOutput

# 尝试导入反应条件预测模块，抑制所有输出
with SuppressAllOutput():
    try:
        from reaction_predictor import ReactionPredictor
        CONDITION_PREDICTION_AVAILABLE = True
    except ImportError:
        CONDITION_PREDICTION_AVAILABLE = False

# 反应条件预测
def PredictReactionConditions(
    reactionDataList: List[ReactionData],
    topNResults: int = 5
) -> ErrorCode:
    """
    反应条件预测函数
    
    参数:
        reactionDataList: 包含反应信息的ReactionData对象列表，每个对象必须包含产物和反应物信息
        topNResults: 指定返回前几个预测结果
        
    返回:
        状态码枚举，ErrorCode.SUCCESS表示成功
    """
    try:
        # 检查模型是否可用
        if not CONDITION_PREDICTION_AVAILABLE:
            raise RetroApiError("反应条件预测模型不可用", ErrorCode.MODEL_LOAD_FAILED)
        
        # 验证输入参数
        if not reactionDataList or not isinstance(reactionDataList, list):
            raise RetroApiError("无效的反应数据列表", ErrorCode.INPUT_INVALID_RESULTS_LIST)
        
        # 验证每个反应数据对象
        for reaction in reactionDataList:
            if not isinstance(reaction, ReactionData):
                raise RetroApiError("反应数据列表中含有非ReactionData类型对象", ErrorCode.INPUT_INVALID_RESULTS_LIST)
            
            # 检查产物和反应物是否存在
            if not reaction.product or not reaction.reactants:
                raise RetroApiError("反应数据中缺少产物或反应物信息", ErrorCode.INPUT_INVALID_SMILES)
        
        # 构建反应SMILES列表
        reaction_smiles_list = []
        for reaction in reactionDataList:
            # 收集所有反应物的SMILES
            reactants_smiles = ".".join([r.smiles for r in reaction.reactants if r.smiles])
            if not reactants_smiles:
                raise RetroApiError("反应物SMILES无效", ErrorCode.INPUT_INVALID_SMILES)
            
            # 构造反应SMILES: 反应物>>产物
            reaction_smiles = reactants_smiles + ">>" + reaction.product
            reaction_smiles_list.append(reaction_smiles)
        
        if not reaction_smiles_list:
            raise RetroApiError("没有有效的反应SMILES", ErrorCode.INPUT_INVALID_SMILES)
        
        # 使用上下文管理器抑制所有输出
        with SuppressAllOutput():
            try:
                predictions = ReactionPredictor(reaction_smiles_list, top_k=topNResults)
            except Exception as e:
                raise RetroApiError(f"调用反应条件预测器失败: {str(e)}", ErrorCode.PREDICTION_EXECUTION_ERROR)
        
        # 如果没有得到有效的预测结果
        if not predictions or not predictions.condition_prediction:
            raise RetroApiError("没有获取到有效的反应条件预测", ErrorCode.PREDICTION_NO_VALID_RESULT)
        
        # 将预测结果填充到输入的ReactionData对象中
        for i, reaction_conditions in enumerate(predictions.condition_prediction):
            if i >= len(reactionDataList):
                break  # 防止索引越界
            
            # 清空现有条件列表
            reactionDataList[i].conditions.clear()
            
            # 添加预测的条件
            for j, condition in enumerate(reaction_conditions):
                if j >= topNResults:
                    break  # 只保留前topNResults个结果
                
                # 创建反应条件对象
                catalysts = []
                if condition.Catalyst:
                    # 处理催化剂数据，可能是字符串或列表
                    if isinstance(condition.Catalyst, list):
                        # 收集所有有效的催化剂
                        catalysts = [cat for cat in condition.Catalyst if cat not in ["nan", "None", None, ""]]
                    elif condition.Catalyst not in ["nan", "None", None, ""]:
                        catalysts = [condition.Catalyst]
                
                # 过滤无效的溶剂和试剂
                solvents = []
                if condition.Solvents:
                    solvents = [s for s in condition.Solvents if s not in ["nan", "None", None, ""]]
                
                reagents = []
                if condition.Reagents:
                    reagents = [r for r in condition.Reagents if r not in ["nan", "None", None, ""]]
                
                # 创建条件对象
                temperature = None
                if condition.Temperature:
                    # 检查温度是否已经包含单位
                    temp_str = str(condition.Temperature).strip()
                    if "°C" in temp_str or "K" in temp_str or "℃" in temp_str or "F" in temp_str:
                        temperature = temp_str
                    else:
                        # 没有单位，添加摄氏度单位
                        temperature = f"{temp_str} °C"
                
                condition_obj = ReactionCondition(
                    catalyst=catalysts,  # 使用催化剂列表
                    solvents=solvents,
                    reagents=reagents,
                    temperature=temperature
                )
                
                # 添加到反应数据中
                reactionDataList[i].conditions.append(condition_obj)
        
        return ErrorCode.SUCCESS
        
    except RetroApiError as e:
        # 如果是已经包装好的RetroApiError，直接返回其状态码
        return e.status
    except Exception as e:
        # 对于未预期的错误，返回未知错误状态码
        return ErrorCode.UNKNOWN_ERROR
