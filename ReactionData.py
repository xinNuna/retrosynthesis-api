from typing import List, Optional, Union

class Reactant:
    """反应物类"""
    
    def __init__(self, 
                 smiles: str, 
                 rxnClass: Optional[int] = None,
                 rank: int = 0,
                 probability: float = 0.0):
        self.smiles = smiles                      # 反应物SMILES字符串
        self.rxnClass = rxnClass                  # 反应类别
        self.rank = rank                          # 结果排名
        self.probability = probability            # 预测概率


class ReactionCondition:
    """反应条件类"""
    def __init__(self, 
                 catalyst: List[str] = None,
                 solvents: List[str] = None, 
                 reagents: List[str] = None,
                 temperature: str = None):
        self.catalyst = catalyst or []      # 催化剂列表
        self.solvents = solvents or []      # 溶剂列表
        self.reagents = reagents or []      # 试剂列表
        self.temperature = temperature      # 温度（包含单位）


class ReactionYield:
    """反应产率类"""
    
    def __init__(self, value: float = 0.0):
        self.value = value  # 产率值（百分比，0-100）


class ReactionData:
    """反应数据类，存储反应相关的所有数据，支持多个反应物"""
    
    def __init__(self, 
                 product: str = None, 
                 reactants: List[Reactant] = None,  
                 conditions: List[ReactionCondition] = None, 
                 yields: List[ReactionYield] = None):
        self.product = product                     # 产物SMILES
        self.reactants = reactants or []           # 反应物对象列表
        self.conditions = conditions or []         # 反应条件对象列表
        self.yields = yields or []                 # 产率对象列表