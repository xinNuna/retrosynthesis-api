"""
半模板法逆合成预测 API
"""
import os
import sys
from typing import List, Optional

# 获取当前脚本的绝对路径
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)  # 指向 code 文件夹

# 添加根目录到路径以便导入其他模块
sys.path.append(ROOT_DIR)

# 添加Graph2Edits-master目录到路径
GRAPH2EDITS_DIR = os.path.join(ROOT_DIR, 'Graph2Edits-master')
sys.path.append(GRAPH2EDITS_DIR)

# 导入所需的模块
from ErrorCodes import ErrorCode, RetroApiError
from ReactionData import ReactionData, Reactant
from .UtilityTools import SuppressAllOutput

# 导入半模板法的预测函数
from api.RetroApi import PredictRetrosynthesis as _PredictRetrosynthesis
from api.RetroApi import Canonicalize

# 半模板法
def SemiTemplateRetrosynthesis(
    productSmiles: str, 
    reactionResultsList: List[ReactionData],
    topNResults: int = 5,
    rxnClass: Optional[int] = None,
    tryAllClasses: bool = True,
    maxSteps: int = 9,
    beamSize: int = 10
) -> ErrorCode:
    """
    半模板逆合成预测函数
    
    参数:
        productSmiles: 产物的SMILES字符串
        reactionResultsList: 结果列表，用于存储预测结果
        topNResults: 指定返回前几个预测结果
        rxnClass: 反应类别（可选，None表示不使用反应类别）
        tryAllClasses: 是否尝试所有反应类别
        maxSteps: 最大编辑步骤数
        beamSize: 束搜索宽度
        
    返回:
        状态码枚举，ErrorCode.SUCCESS表示成功
    """
    # 使用上下文管理器抑制所有输出
    with SuppressAllOutput():
        return _PredictRetrosynthesis(
            productSmiles=productSmiles,
            reactionResultsList=reactionResultsList,
            topNResults=topNResults,
            rxnClass=rxnClass,
            tryAllClasses=tryAllClasses,
            maxSteps=maxSteps,
            beamSize=beamSize
        )
