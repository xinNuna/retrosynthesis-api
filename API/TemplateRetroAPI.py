"""
模板法逆合成预测 API
"""
import os
import sys
from typing import List, Optional

# 获取当前脚本的绝对路径
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)  # 指向 code 文件夹

# 添加根目录到路径以便导入其他模块
sys.path.append(ROOT_DIR)

# LocalRetro目录
LOCALRETRO_DIR = os.path.join(ROOT_DIR, 'LocalRetro')
SCRIPTS_DIR = os.path.join(LOCALRETRO_DIR, 'scripts')

# 导入所需的模块
from ErrorCodes import ErrorCode, RetroApiError
from ReactionData import ReactionData, Reactant
from .UtilityTools import SuppressAllOutput


def AddAtomMapping(smiles: str) -> str:
    """为SMILES添加原子映射"""
    from rdkit import Chem
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise RetroApiError(f"无效的SMILES: {smiles}", ErrorCode.INPUT_INVALID_SMILES)
        
        # 检查是否已有原子映射
        has_mapping = any(atom.GetAtomMapNum() > 0 for atom in mol.GetAtoms())
        
        if has_mapping:
            # 如果已有映射，先清除
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
        
        # 添加新的原子映射
        for i, atom in enumerate(mol.GetAtoms()):
            atom.SetAtomMapNum(i + 1)
        
        return Chem.MolToSmiles(mol)
    except Exception as e:
        raise RetroApiError(f"处理SMILES时出错: {str(e)}", ErrorCode.INPUT_INVALID_SMILES)


def RunLocalRetro(smiles: str, top_k: int) -> List[str]:
    """运行LocalRetro"""
    # 保存原始环境
    original_cwd = os.getcwd()
    original_path = sys.path.copy()
    
    try:
        # 1. 切换到scripts目录，这样相对路径才能正确工作
        os.chdir(SCRIPTS_DIR)
        
        # 2. 将scripts目录加入Python路径的最前面
        if SCRIPTS_DIR not in sys.path:
            sys.path.insert(0, SCRIPTS_DIR)
        if LOCALRETRO_DIR not in sys.path:
            sys.path.insert(0, LOCALRETRO_DIR)
        
        # 3. 清理可能冲突的模块
        modules_to_clean = ['utils', 'Test', 'get_edit', 'Decode_predictions', 
                          'models', 'model_utils', 'dataset', 'cli']
        for module in modules_to_clean:
            if module in sys.modules:
                del sys.modules[module]
        
        # 4. 导入并运行
        from cli import get_reactants_template_base
        results = get_reactants_template_base(smiles, top_k)
        
        return results if results else []
        
    except Exception as e:
        #print(f"RunLocalRetro错误: {e}")
        import traceback
        traceback.print_exc()
        return []
    finally:
        # 恢复原始环境
        os.chdir(original_cwd)
        sys.path = original_path


def ConvertToReactionData(productSmiles: str, reactants_list: List[str]) -> List[ReactionData]:
    """将预测结果转换为ReactionData对象"""
    results = []
    
    for i, reactants_str in enumerate(reactants_list):
        if not reactants_str:
            continue
            
        try:
            # 解析反应物SMILES
            reactant_smiles_list = reactants_str.split('.')
            
            # 创建反应物对象列表
            reactants = []
            for r_smiles in reactant_smiles_list:
                r_smiles = r_smiles.strip()
                if r_smiles:
                    reactant = Reactant(
                        smiles=r_smiles,
                        rxnClass=None,  # 模板法不使用反应类别
                        rank=i + 1,
                        probability=1.0 / (i + 1)  # 简单的概率估计
                    )
                    reactants.append(reactant)
            
            if reactants:
                # 创建反应结果对象
                reactionData = ReactionData(
                    product=productSmiles,
                    reactants=reactants,
                    conditions=[],
                    yields=[]
                )
                results.append(reactionData)
                
        except Exception:
            continue
    
    return results


def TemplateRetrosynthesis(
    productSmiles: str, 
    reactionResultsList: List[ReactionData],
    topNResults: int = 10,
    **kwargs  # 接受其他参数但不使用
) -> ErrorCode:
    """
    模板法逆合成预测函数
    
    参数:
        productSmiles: 产物的SMILES字符串（不需要原子映射）
        reactionResultsList: 结果列表，用于存储预测结果
        topNResults: 指定返回前几个预测结果（默认10个）
        
    返回:
        状态码枚举，ErrorCode.SUCCESS表示成功
    """
    try:
        # ===== 参数验证 =====
        if not productSmiles or not isinstance(productSmiles, str):
            raise RetroApiError("产物SMILES不能为空且必须是字符串", ErrorCode.INPUT_INVALID_SMILES)
            
        if reactionResultsList is None:
            raise RetroApiError("结果列表不能为None", ErrorCode.INPUT_INVALID_RESULTS_LIST)
            
        if not isinstance(reactionResultsList, list):
            raise RetroApiError("结果列表必须是列表类型", ErrorCode.INPUT_INVALID_RESULTS_LIST)
        
        if topNResults <= 0:
            topNResults = 10
            
        # 清空结果列表
        reactionResultsList.clear()
        
        # 使用上下文管理器抑制所有输出
        with SuppressAllOutput():
            # 为SMILES添加原子映射
            mapped_smiles = AddAtomMapping(productSmiles)
            
            # 运行LocalRetro
            reactants_list = RunLocalRetro(mapped_smiles, topNResults)
            
            # 转换为ReactionData对象
            results = ConvertToReactionData(productSmiles, reactants_list)
            
            if not results:
                raise RetroApiError("没有预测到有效的反应物", ErrorCode.PREDICTION_NO_VALID_REACTANTS)
            
            # 将结果保存到结果列表中
            reactionResultsList.extend(results)
            
            return ErrorCode.SUCCESS
        
    except RetroApiError as e:
        return e.status
    except Exception as e:
        #print(f"TemplateRetrosynthesis错误: {e}")
        import traceback
        traceback.print_exc()
        return ErrorCode.UNKNOWN_ERROR