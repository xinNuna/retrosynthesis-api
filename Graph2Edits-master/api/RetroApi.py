import os
import sys
import torch
import glob
import re
from datetime import datetime
from rdkit import Chem, RDLogger
from typing import List, Dict, Union, Optional, Tuple, Any, Callable


# 添加code目录到Python路径，以便能够找到ErrorCodes和ReactionData模块
CODE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(CODE_DIR)


from ErrorCodes import ErrorCode, RetroApiError
from ReactionData import Reactant, ReactionCondition, ReactionYield, ReactionData

# 导入模型和束搜索
sys.path.append('.')  # 将当前目录添加到路径
try:
    from models import Graph2Edits, BeamSearch
except ImportError:
    print("错误: 无法导入 Graph2Edits。请确保您正在项目根目录下运行。")
    sys.exit(1)

# 抑制RDKit警告
RDLogger.logger().setLevel(RDLogger.CRITICAL)

# 全局变量
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
DEFAULT_MAX_STEPS = 9
DEFAULT_BEAM_SIZE = 10
DEFAULT_TOP_K = 5

# 固定的模型路径
MODEL_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 
                        "experiments/uspto_50k/with_rxn_class/02-04-2025--17-34-39/epoch_121.pt")


def Canonicalize(smi: str) -> str:
    """规范化SMILES字符串"""
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return smi
        mol = Chem.RemoveHs(mol)
        [a.ClearProp('molAtomMapNumber') for a in mol.GetAtoms()]
        return Chem.MolToSmiles(mol)
    except Exception as e:
        raise RetroApiError(f"规范化SMILES时出错: {str(e)}", ErrorCode.INPUT_INVALID_SMILES)


def LoadModel(modelPath: str, useRxnClass: bool = False) -> Any:
    """加载训练好的模型"""
    try:
        if not os.path.exists(modelPath):
            raise RetroApiError(f"模型文件不存在: {modelPath}", ErrorCode.MODEL_FILE_NOT_FOUND)
            
        checkpoint = torch.load(modelPath, map_location=torch.device(DEVICE))
        config = checkpoint['saveables']
        
        model = Graph2Edits(**config, device=DEVICE)
        model.load_state_dict(checkpoint['state'])
        model.to(DEVICE)
        model.eval()
        
        beamModel = BeamSearch(model=model, step_beam_size=10, beam_size=10, use_rxn_class=useRxnClass)
        return beamModel
    except Exception as e:
        # 如果异常已经是RetroApiError，直接重新抛出
        if isinstance(e, RetroApiError):
            raise
        # 否则，包装为RetroApiError
        raise RetroApiError(f"加载模型时出错: {str(e)}", ErrorCode.MODEL_LOAD_FAILED)


def PrepareProduct(smiles: str) -> str:
    """准备用于预测的产物SMILES"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise RetroApiError(f"无法解析SMILES: {smiles}", ErrorCode.INPUT_INVALID_SMILES)
            
        # 芳香化并在没有原子映射时添加映射
        try:
            Chem.Kekulize(mol)
        except:
            pass
            
        # 检查分子是否有原子映射，如果没有则添加
        hasMapping = False
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() > 0:
                hasMapping = True
                break
                
        if not hasMapping:
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx() + 1)
                
        return Chem.MolToSmiles(mol)
    except Exception as e:
        # 如果异常已经是RetroApiError，直接重新抛出
        if isinstance(e, RetroApiError):
            raise
        # 否则，包装为RetroApiError
        raise RetroApiError(f"准备产物时出错: {str(e)}", ErrorCode.INPUT_INVALID_SMILES)


def SplitReactants(reactants_smiles: str) -> List[str]:
    """将点分隔的SMILES拆分为单个反应物SMILES"""
    reactants = reactants_smiles.split('.')
    return [r for r in reactants if r.strip()]


def PredictReactantsWithModel(beamModel: Any, productSmiles: str, 
                            rxnClass: Optional[int] = None, 
                            maxSteps: int = DEFAULT_MAX_STEPS) -> List[Any]:
    """使用已加载的模型从产物SMILES预测反应物"""
    try:
        with torch.no_grad():
            topKResults = beamModel.run_search(
                prod_smi=productSmiles, 
                max_steps=maxSteps, 
                rxn_class=rxnClass
            )
            
        
        class PredictionResult:
            def __init__(self, rank, probability, reactants, reactantsCanonical, rxnClass=None):
                self.rank = rank
                self.probability = probability
                self.reactants = reactants
                self.reactantsCanonical = reactantsCanonical
                self.rxnClass = rxnClass
                
        results = []
        for i, path in enumerate(topKResults):
            if 'final_smi' in path and path['final_smi'] != 'final_smi_unmapped':
                prob = path['prob']
                result = PredictionResult(
                    rank=i + 1,
                    probability=round(prob, 4),
                    reactants=path['final_smi'],
                    reactantsCanonical=Canonicalize(path['final_smi']),
                    rxnClass=rxnClass
                )
                results.append(result)
            
        if not results:
            raise RetroApiError("没有预测到有效的反应物", ErrorCode.PREDICTION_NO_VALID_REACTANTS)
            
        return results
    except Exception as e:
        # 如果异常已经是RetroApiError，直接重新抛出
        if isinstance(e, RetroApiError):
            raise
        # 否则，包装为RetroApiError
        raise RetroApiError(f"预测过程中出错: {str(e)}", ErrorCode.PREDICTION_EXECUTION_ERROR)


def PredictRetrosynthesis(
    productSmiles: str, 
    reactionResultsList: List[ReactionData],
    topNResults: int = DEFAULT_TOP_K,
    rxnClass: Optional[int] = None,
    tryAllClasses: bool = False,
    maxSteps: int = DEFAULT_MAX_STEPS,
    beamSize: int = DEFAULT_BEAM_SIZE
) -> ErrorCode:
    """
    执行逆合成预测的主函数，返回状态码枚举
    
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
    try:
        # ===== 参数验证 =====
        # 验证产物SMILES
        if not productSmiles or not isinstance(productSmiles, str):
            raise RetroApiError("产物SMILES不能为空且必须是字符串", ErrorCode.INPUT_INVALID_SMILES)
            
        # 验证结果列表
        if reactionResultsList is None:
            raise RetroApiError("结果列表不能为None", ErrorCode.INPUT_INVALID_RESULTS_LIST)
            
        if not isinstance(reactionResultsList, list):
            raise RetroApiError("结果列表必须是列表类型", ErrorCode.INPUT_INVALID_RESULTS_LIST)
        
        # 验证topNResults
        if topNResults <= 0:
            topNResults = DEFAULT_TOP_K
            
        MAX_RESULTS = 100
        if topNResults > MAX_RESULTS:
            raise RetroApiError(f"返回结果数量不能超过{MAX_RESULTS}", ErrorCode.INPUT_INVALID_TOP_N)
            
        # 验证maxSteps
        if maxSteps <= 0:
            raise RetroApiError("最大编辑步骤数必须大于0", ErrorCode.INPUT_INVALID_MAX_STEPS)
            
        if maxSteps > 20:
            raise RetroApiError("最大编辑步骤数不能超过20", ErrorCode.INPUT_INVALID_MAX_STEPS)
            
        # 验证beamSize
        if beamSize < 2:
            raise RetroApiError("束搜索宽度必须至少为2", ErrorCode.INPUT_INVALID_BEAM_SIZE)
            
        if beamSize > 50:
            raise RetroApiError("束搜索宽度不能超过50", ErrorCode.INPUT_INVALID_BEAM_SIZE)
        
        # 验证反应类别
        if rxnClass is not None and (rxnClass < 0 or rxnClass > 9):
            raise RetroApiError(f"无效的反应类别: {rxnClass}，有效范围为0-9", ErrorCode.INPUT_INVALID_RXN_CLASS)
            
        # 验证反应类别设置组合
        if rxnClass is None and not tryAllClasses:
            raise RetroApiError("未指定反应类别(rxnClass=None)时，必须设置tryAllClasses=True", 
                                ErrorCode.INPUT_MISSING_RXN_CLASS)
            
        # 清空结果列表
        reactionResultsList.clear()
            
        # 确定是否使用反应类别
        useRxnClass = (rxnClass is not None) or tryAllClasses
        
        # 加载模型
        beamModel = LoadModel(MODEL_PATH, useRxnClass)
        
        # 准备产物SMILES
        preparedSmiles = PrepareProduct(productSmiles)
        productCanonical = Canonicalize(productSmiles)  # 规范化产物SMILES
        
        # 处理反应类别
        results = []
        if tryAllClasses:
            allResults = []
            
            for rxnClassIdx in range(10):  # USPTO-50k 有10个反应类别
                try:
                    classResults = PredictReactantsWithModel(
                        beamModel, preparedSmiles, 
                        rxnClass=rxnClassIdx, 
                        maxSteps=maxSteps
                    )
                    if classResults:
                        for result in classResults:
                            # 拆分反应物SMILES
                            reactant_smiles_list = SplitReactants(result.reactantsCanonical)
                            
                            # 创建反应物对象列表
                            reactants = []
                            for i, r_smiles in enumerate(reactant_smiles_list):
                                reactant = Reactant(
                                    smiles=r_smiles,
                                    rxnClass=rxnClassIdx,  # 使用当前成功的反应类别
                                    rank=result.rank,
                                    probability=result.probability / len(reactant_smiles_list)  # 简单地将概率平分
                                )
                                reactants.append(reactant)
                            
                            # 创建反应条件和产率对象列表
                            conditions = []
                            yields = []
                            
                            # 创建反应结果对象 - 使用新的多反应物模式
                            reactionData = ReactionData(
                                product=productCanonical,
                                reactants=reactants,  # 现在传递反应物列表
                                conditions=conditions,
                                yields=yields
                            )
                            
                            allResults.append(reactionData)
                except RetroApiError:
                    # 某个反应类别没有预测结果，继续尝试其他类别
                    continue
            
            if not allResults:
                raise RetroApiError("所有反应类别均未预测到有效的反应物", ErrorCode.PREDICTION_ALL_CLASSES_FAILED)
                
            # 按照反应物的概率排序
            allResults.sort(key=lambda x: x.reactants[0].probability if x.reactants else 0, reverse=True)
            results = allResults[:topNResults]  # 使用topNResults参数
                
        else:
            # 指定了反应类别的情况
            try:
                oldResults = PredictReactantsWithModel(
                    beamModel, preparedSmiles, 
                    rxnClass=rxnClass, 
                    maxSteps=maxSteps
                )
                
                # 只保留指定数量的结果
                topResults = oldResults[:topNResults]  # 使用topNResults参数
                
                # 转换为新的ReactionData对象 - 使用多反应物模式
                for result in topResults:
                    # 拆分反应物SMILES
                    reactant_smiles_list = SplitReactants(result.reactantsCanonical)
                    
                    # 创建反应物对象列表
                    reactants = []
                    for i, r_smiles in enumerate(reactant_smiles_list):
                        reactant = Reactant(
                            smiles=r_smiles,
                            rxnClass=rxnClass,  # 用户指定的反应类别
                            rank=result.rank,
                            probability=result.probability / len(reactant_smiles_list)  # 简单地将概率平分
                        )
                        reactants.append(reactant)
                    
                    # 创建反应条件和产率对象列表
                    conditions = []
                    yields = []
                    
                    # 创建反应结果对象
                    reactionData = ReactionData(
                        product=productCanonical,
                        reactants=reactants,  # 现在传递反应物列表
                        conditions=conditions,
                        yields=yields
                    )
                    
                    results.append(reactionData)
            except RetroApiError as e:
                # 如果指定反应类别失败，直接报错，不再尝试其他类别
                raise RetroApiError(f"使用指定反应类别 {rxnClass} 预测失败: {e.message}", ErrorCode.PREDICTION_NO_VALID_REACTANTS)
        
        # 将结果保存到结果列表中
        reactionResultsList.extend(results)
        
        return ErrorCode.SUCCESS
        
    except RetroApiError as e:
        # 只返回错误状态，不打印错误信息
        return e.status
    except Exception as e:
        # 处理未预期的错误，同样不打印错误信息
        return ErrorCode.UNKNOWN_ERROR