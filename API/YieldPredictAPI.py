"""
产率预测 API
"""
import os
import sys
from typing import List, Optional, Any

import torch

# 获取当前脚本的绝对路径
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)  # 指向 code 文件夹

# 添加根目录到路径以便导入其他模块
sys.path.append(ROOT_DIR)

# 添加Egret-main目录到路径
EGRET_DIR = os.path.join(ROOT_DIR, 'Egret-main')
sys.path.append(EGRET_DIR)

# 导入所需的模块
from ErrorCodes import ErrorCode, RetroApiError
from ReactionData import ReactionData, ReactionYield
from .UtilityTools import SuppressAllOutput

# 导入原始的产率预测API
sys.path.append(os.path.join(EGRET_DIR, 'inference'))
from yield_predict_API import YieldPredictorAPI, REACTION

# 全局变量 - 模型路径
MODEL_STATE_PATH = os.path.join(EGRET_DIR, 'inference', 'yield_prediction_model')

# 全局预测器实例（单例模式）
_YIELD_PREDICTOR = None


def GetYieldPredictor() -> YieldPredictorAPI:
    """获取产率预测器实例（单例模式）"""
    global _YIELD_PREDICTOR
    
    if _YIELD_PREDICTOR is not None:
        return _YIELD_PREDICTOR
    
    try:
        # 检查模型目录是否存在
        if not os.path.exists(MODEL_STATE_PATH):
            raise RetroApiError(
                f"模型目录不存在: {MODEL_STATE_PATH}", 
                ErrorCode.MODEL_FILE_NOT_FOUND
            )
        
        # 检查必要的模型文件
        required_files = ["config.json", "pytorch_model.bin", "vocab.txt"]
        for file in required_files:
            file_path = os.path.join(MODEL_STATE_PATH, file)
            if not os.path.exists(file_path):
                raise RetroApiError(
                    f"模型文件缺失: {file_path}", 
                    ErrorCode.MODEL_FILE_NOT_FOUND
                )
        
        print(f"模型路径: {MODEL_STATE_PATH}")
        
        # 初始化预测器（自动检测CUDA）
        cuda_device = 0 if torch.cuda.is_available() else -1
        print(f"CUDA可用: {torch.cuda.is_available()}, 设备ID: {cuda_device}")
        
        # 尝试加载模型
        _YIELD_PREDICTOR = YieldPredictorAPI(
            model_state_path=MODEL_STATE_PATH,
            cuda_device=cuda_device
        )
        print("产率预测模型加载成功")
        return _YIELD_PREDICTOR
    except RetroApiError as e:
        raise e
    except Exception as e:
        import traceback
        print(f"初始化产率预测器失败: {str(e)}")
        print(f"详细错误信息: {traceback.format_exc()}")
        raise RetroApiError(f"初始化产率预测器失败: {str(e)}", ErrorCode.MODEL_INIT_ERROR)


def ConvertReactionDataToREACTION(reactionData: ReactionData) -> REACTION:
    """将ReactionData转换为原始API需要的REACTION格式"""
    try:
        # 验证输入
        if not reactionData.product:
            raise RetroApiError("产物SMILES不能为空", ErrorCode.INPUT_INVALID_SMILES)
        
        if not reactionData.reactants:
            raise RetroApiError("反应物列表不能为空", ErrorCode.INPUT_INVALID_SMILES)
        
        # 提取反应物SMILES
        reactant_smiles = []
        for reactant in reactionData.reactants:
            if reactant.smiles:
                # 确保 SMILES 是有效的字符串
                smiles_str = reactant.smiles.strip()
                if smiles_str:
                    reactant_smiles.append(smiles_str)
        
        if not reactant_smiles:
            raise RetroApiError("没有有效的反应物SMILES", ErrorCode.INPUT_INVALID_SMILES)
        
        # 提取条件SMILES（催化剂、溶剂、试剂）
        condition_smiles = []
        for condition in reactionData.conditions:
            # 处理催化剂
            if condition.catalyst:
                for cat in condition.catalyst:
                    if cat and isinstance(cat, str):
                        condition_smiles.append(cat.strip())
            
            # 处理溶剂
            if condition.solvents:
                for solv in condition.solvents:
                    if solv and isinstance(solv, str):
                        condition_smiles.append(solv.strip())
            
            # 处理试剂
            if condition.reagents:
                for reag in condition.reagents:
                    if reag and isinstance(reag, str):
                        condition_smiles.append(reag.strip())
        
        # 过滤空字符串
        condition_smiles = [c for c in condition_smiles if c]
        
        # 如果没有条件，添加一个默认条件以确保模型正常工作
        if not condition_smiles:
            print("警告: 没有有效的反应条件，添加默认占位符")
            condition_smiles = ["O"]  # 水作为默认条件
        
        # 创建REACTION对象
        print(f"创建REACTION对象: 产物={reactionData.product}, 反应物数量={len(reactant_smiles)}, 条件数量={len(condition_smiles)}")
        reaction = REACTION(
            Product=reactionData.product.strip(),
            Reactant=reactant_smiles,
            Condition=condition_smiles
        )
        
        return reaction
        
    except Exception as e:
        if isinstance(e, RetroApiError):
            raise
        import traceback
        print(f"转换反应数据时出错: {str(e)}")
        print(f"详细错误信息: {traceback.format_exc()}")
        raise RetroApiError(f"转换反应数据时出错: {str(e)}", ErrorCode.INPUT_INVALID_SMILES)


def PredictReactionYield(reactionData: ReactionData) -> ErrorCode:
    """
    预测反应产率的主函数
    
    参数:
        reactionData: 包含反应信息的ReactionData对象
        
    返回:
        状态码枚举，ErrorCode.SUCCESS表示成功
    """
    try:
        # 参数验证
        if reactionData is None:
            raise RetroApiError("反应数据不能为None", ErrorCode.INPUT_INVALID_RESULTS_LIST)
        
        if not isinstance(reactionData, ReactionData):
            raise RetroApiError("反应数据必须是ReactionData类型", ErrorCode.INPUT_INVALID_RESULTS_LIST)
        
        # 准备清空原有的产率数据
        reactionData.yields = []
        
        try:
            # 获取预测器实例
            print("尝试获取产率预测器实例...")
            predictor = GetYieldPredictor()
            
            # 转换数据格式
            print("转换数据格式...")
            reaction = ConvertReactionDataToREACTION(reactionData)
            
            # 记录输入数据用于调试
            print(f"产物: {reaction.Product}")
            print(f"反应物: {reaction.Reactant}")
            print(f"条件: {reaction.Condition}")
            
            # 调用原始API进行预测
            print("调用产率预测函数...")
            predicted_yield = predictor.predict_yield(reaction)
            print(f"预测结果: {predicted_yield}")
            
            # 创建产率对象
            yield_result = ReactionYield()
            yield_result.value = float(predicted_yield)
            
            # 将产率结果添加到反应数据中
            reactionData.yields.append(yield_result)
            
            print("产率预测完成")
            return ErrorCode.SUCCESS
            
        except Exception as e:
            import traceback
            print(f"产率预测过程中出现异常: {str(e)}")
            print(f"详细错误信息: {traceback.format_exc()}")
            raise RetroApiError(f"产率预测失败: {str(e)}", ErrorCode.PREDICTION_EXECUTION_ERROR)
        
    except RetroApiError as e:
        print(f"产率预测API错误: [{e.status}] {e.message}")
        return e.status
    except Exception as e:
        import traceback
        print(f"产率预测中出现未预期的异常: {str(e)}")
        print(f"详细错误信息: {traceback.format_exc()}")
        return ErrorCode.UNKNOWN_ERROR