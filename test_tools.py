import sys

sys.stdout.reconfigure(encoding='utf-8')
sys.stdin.reconfigure(encoding='utf-8')

import asyncio
import json
import logging
import os
from typing import Any, List, Optional, Dict

from mcp.server import Server
import mcp.types as types
from mcp.server.stdio import stdio_server

# 确保可以导入API模块
try:
    # 直接尝试导入模块
    from API.SemiTemplateRetroAPI import SemiTemplateRetrosynthesis
    from API.TemplateRetroAPI import TemplateRetrosynthesis
    from API.ConditionAPI import PredictReactionConditions
    from API.YieldPredictAPI import PredictReactionYield
    from ReactionData import ReactionData, Reactant, ReactionCondition, ReactionYield
    from ErrorCodes import ErrorCode
except ImportError:
    # 如果导入失败，尝试添加当前目录的父目录到路径
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(current_dir)
    if parent_dir not in sys.path:
        sys.path.append(parent_dir)
    
    # 再次尝试导入
    from API.SemiTemplateRetroAPI import SemiTemplateRetrosynthesis
    from API.TemplateRetroAPI import TemplateRetrosynthesis
    from API.ConditionAPI import PredictReactionConditions
    from API.YieldPredictAPI import PredictReactionYield
    from ReactionData import ReactionData, Reactant, ReactionCondition, ReactionYield
    from ErrorCodes import ErrorCode

logger = logging.getLogger(__name__)

app = Server("retrosynthesis-server")


@app.list_tools()
async def ListTools() -> list[types.Tool]:
    """列出可用的药物逆合成预测工具"""
    return [
        types.Tool(
            name="deep_research",
            description="深度研究工具，提供深入的研究和分析功能，可以对任何主题进行迭代式深度研究",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "研究查询内容，您想要研究的主题或问题"
                    },
                    "breadth": {
                        "type": "integer",
                        "description": "研究广度，控制生成多少个搜索查询（推荐值：3-10）",
                        "default": 4,
                        "minimum": 1,
                        "maximum": 20
                    },
                    "depth": {
                        "type": "integer",
                        "description": "研究深度，控制递归探索的层级数量（推荐值：1-5）",
                        "default": 2,
                        "minimum": 1,
                        "maximum": 5
                    },
                    "concurrency": {
                        "type": "integer",
                        "description": "并发限制，控制同时处理的搜索查询数量（建议设为1）",
                        "default": 3,
                        "minimum": 1,
                        "maximum": 10
                    }
                },
                "required": ["query"]
            },
        ),
        # 半模板法工具
        types.Tool(
            name="semi_template_retrosynthesis",
            description="半模板法逆合成预测工具，预测产物的合成路径和反应物",
            inputSchema={
                "type": "object",
                "properties": {
                    "product_smiles": {
                        "type": "string",
                        "description": "产物的SMILES字符串"
                    },
                    "top_n_results": {
                        "type": "integer",
                        "description": "返回前N个预测结果",
                        "default": 5,
                        "minimum": 1,
                        "maximum": 20
                    },
                    "rxn_class": {
                        "type": ["integer", "null"],
                        "description": "反应类别（可选，不指定则自动判断）",
                        "default": None
                    },
                    "try_all_classes": {
                        "type": "boolean",
                        "description": "是否尝试所有反应类别",
                        "default": True
                    },
                    "max_steps": {
                        "type": "integer",
                        "description": "最大编辑步骤数",
                        "default": 9,
                        "minimum": 1,
                        "maximum": 20
                    },
                    "beam_size": {
                        "type": "integer",
                        "description": "束搜索宽度",
                        "default": 10,
                        "minimum": 1,
                        "maximum": 50
                    }
                },
                "required": ["product_smiles"]
            },
        ),
        # 模板法工具
        types.Tool(
            name="template_retrosynthesis",
            description="模板法逆合成预测工具，使用反应模板预测产物的合成路径和反应物",
            inputSchema={
                "type": "object",
                "properties": {
                    "product_smiles": {
                        "type": "string",
                        "description": "产物的SMILES字符串"
                    },
                    "top_n_results": {
                        "type": "integer",
                        "description": "返回前N个预测结果",
                        "default": 10,
                        "minimum": 1,
                        "maximum": 50
                    }
                },
                "required": ["product_smiles"]
            },
        ),
        # 反应条件预测工具
        types.Tool(
            name="predict_reaction_conditions",
            description="反应条件预测工具，预测反应所需的条件（催化剂、溶剂、试剂等）",
            inputSchema={
                "type": "object",
                "properties": {
                    "product_smiles": {
                        "type": "string",
                        "description": "产物的SMILES字符串"
                    },
                    "reactants_smiles": {
                        "type": "array",
                        "description": "反应物的SMILES字符串列表",
                        "items": {
                            "type": "string"
                        }
                    },
                    "top_n_results": {
                        "type": "integer",
                        "description": "返回前N个预测结果",
                        "default": 5,
                        "minimum": 1,
                        "maximum": 20
                    }
                },
                "required": ["product_smiles", "reactants_smiles"]
            },
        ),
        # 产率预测工具
        types.Tool(
            name="predict_reaction_yield",
            description="反应产率预测工具，预测化学反应的产率",
            inputSchema={
                "type": "object",
                "properties": {
                    "product_smiles": {
                        "type": "string",
                        "description": "产物的SMILES字符串"
                    },
                    "reactants_smiles": {
                        "type": "array",
                        "description": "反应物的SMILES字符串列表",
                        "items": {
                            "type": "string"
                        }
                    },
                    "conditions": {
                        "type": "object",
                        "description": "反应条件",
                        "properties": {
                            "catalyst": {
                                "type": ["string", "null"],
                                "description": "催化剂SMILES"
                            },
                            "solvents": {
                                "type": "array",
                                "description": "溶剂SMILES列表",
                                "items": {
                                    "type": "string"
                                }
                            },
                            "reagents": {
                                "type": "array",
                                "description": "试剂SMILES列表",
                                "items": {
                                    "type": "string"
                                }
                            },
                            "temperature": {
                                "type": ["string", "null"],
                                "description": "反应温度（包含单位）"
                            }
                        }
                    }
                },
                "required": ["product_smiles", "reactants_smiles"]
            },
        )
    ]


# 辅助函数：将ReactionData对象转换为字典
def ReactionDataToDict(reaction: ReactionData) -> Dict:
    """将ReactionData对象转换为可JSON序列化的字典"""
    result = {
        "product": reaction.product,
        "reactants": []
    }
    
    # 转换反应物
    for reactant in reaction.reactants:
        reactant_dict = {
            "smiles": reactant.smiles,
            "rank": reactant.rank,
            "probability": reactant.probability
        }
        if reactant.rxnClass is not None:
            reactant_dict["rxnClass"] = reactant.rxnClass
        result["reactants"].append(reactant_dict)
      # 转换反应条件
    if reaction.conditions:
        result["conditions"] = []
        for condition in reaction.conditions:
            condition_dict = {
                "solvents": condition.solvents,
                "reagents": condition.reagents
            }
            if condition.catalyst:
                condition_dict["catalyst"] = condition.catalyst
            if condition.temperature:
                condition_dict["temperature"] = condition.temperature
            result["conditions"].append(condition_dict)
    
    # 转换产率
    if reaction.yields:
        result["yields"] = [{"value": y.value} for y in reaction.yields]
    
    return result


@app.call_tool()
async def CallTool(name: str, arguments: dict) -> list[types.TextContent]:
    """处理工具调用请求"""
    try:
        if name == "deep_research":
            # 获取参数
            query = arguments["query"]
            breadth = arguments.get("breadth", 4)  # 默认广度为4
            depth = arguments.get("depth", 2)  # 默认深度为2
            concurrency = arguments.get("concurrency", 3)  # 默认并发为3
            
            # 构建Docker命令
            dockerCmd = [
                "docker", "exec", "-i", "deep-research", 
                "node", "src/cli.js", 
                "--query", query,
                "--breadth", str(breadth),
                "--depth", str(depth),
                "--concurrency", str(concurrency)
            ]
            
            # 执行Docker命令
            try:
                result = await RunDockerCommand(dockerCmd)
                content_type = "markdown"
                return [types.TextContent(type=content_type, text=result)]
            except Exception as e:
                logger.exception(f"Deep Research服务执行错误: {str(e)}")
                return [types.TextContent(
                    type="text",
                    text=f"Deep Research服务执行错误: {str(e)}",
                    isError=True
                )]
                
        # ==== 半模板法逆合成预测 ====
        elif name == "semi_template_retrosynthesis":
            # 获取参数
            productSmiles = arguments["product_smiles"]
            topNResults = arguments.get("top_n_results", 5)
            rxnClass = arguments.get("rxn_class")
            tryAllClasses = arguments.get("try_all_classes", True)
            maxSteps = arguments.get("max_steps", 9)
            beamSize = arguments.get("beam_size", 10)
            
            # 准备结果列表
            reactionResults = []
            
            # 调用API
            status = SemiTemplateRetrosynthesis(
                productSmiles=productSmiles,
                reactionResultsList=reactionResults,
                topNResults=topNResults,
                rxnClass=rxnClass,
                tryAllClasses=tryAllClasses,
                maxSteps=maxSteps,
                beamSize=beamSize
            )
            
            # 检查状态
            if status != ErrorCode.SUCCESS:
                return [types.TextContent(
                    type="text",
                    text=f"半模板法逆合成预测失败，错误码: {status.name} ({status.value})",
                    isError=True
                )]
            
            # 将结果转换为字典列表
            resultsDict = [ReactionDataToDict(r) for r in reactionResults]
            
            # 返回结果
            return [types.TextContent(
                type="json",
                text=json.dumps(resultsDict, ensure_ascii=False)
            )]
            
        # ==== 模板法逆合成预测 ====
        elif name == "template_retrosynthesis":
            # 获取参数
            productSmiles = arguments["product_smiles"]
            topNResults = arguments.get("top_n_results", 10)
            
            # 准备结果列表
            reactionResults = []
            
            # 调用API
            status = TemplateRetrosynthesis(
                productSmiles=productSmiles,
                reactionResultsList=reactionResults,
                topNResults=topNResults
            )
            
            # 检查状态
            if status != ErrorCode.SUCCESS:
                return [types.TextContent(
                    type="text",
                    text=f"模板法逆合成预测失败，错误码: {status.name} ({status.value})",
                    isError=True
                )]
            
            # 将结果转换为字典列表
            resultsDict = [ReactionDataToDict(r) for r in reactionResults]
            
            # 返回结果
            return [types.TextContent(
                type="json",
                text=json.dumps(resultsDict, ensure_ascii=False)
            )]
            
        # ==== 反应条件预测 ====
        elif name == "predict_reaction_conditions":
            # 获取参数
            productSmiles = arguments["product_smiles"]
            reactantsSmiles = arguments["reactants_smiles"]
            topNResults = arguments.get("top_n_results", 5)
            
            # 创建反应数据对象
            reactionData = ReactionData(
                product=productSmiles,
                reactants=[Reactant(smiles=s) for s in reactantsSmiles]
            )
            
            # 准备反应数据列表
            reactionDataList = [reactionData]
            
            # 调用API
            status = PredictReactionConditions(
                reactionDataList=reactionDataList,
                topNResults=topNResults
            )
            
            # 检查状态
            if status != ErrorCode.SUCCESS:
                return [types.TextContent(
                    type="text",
                    text=f"反应条件预测失败，错误码: {status.name} ({status.value})",
                    isError=True
                )]
            
            # 将结果转换为字典列表
            resultsDict = [ReactionDataToDict(r) for r in reactionDataList]
            
            # 返回结果
            return [types.TextContent(
                type="json",
                text=json.dumps(resultsDict, ensure_ascii=False)
            )]
            
        # ==== 产率预测 ====
        elif name == "predict_reaction_yield":
            # 获取参数
            productSmiles = arguments["product_smiles"]
            reactantsSmiles = arguments["reactants_smiles"]
            conditions = arguments.get("conditions", {})
            
            # 创建反应物对象
            reactants = [Reactant(smiles=s) for s in reactantsSmiles]
            
            # 创建反应条件对象
            reactionConditions = []
            if conditions:
                reactionCondition = ReactionCondition(
                    catalyst=conditions.get("catalyst"),
                    solvents=conditions.get("solvents", []),
                    reagents=conditions.get("reagents", []),
                    temperature=conditions.get("temperature")
                )
                reactionConditions = [reactionCondition]
            
            # 创建反应数据对象
            reactionData = ReactionData(
                product=productSmiles,
                reactants=reactants,
                conditions=reactionConditions
            )
            
            # 调用API
            status = PredictReactionYield(reactionData)
            
            # 检查状态
            if status != ErrorCode.SUCCESS:
                return [types.TextContent(
                    type="text",
                    text=f"产率预测失败，错误码: {status.name} ({status.value})",
                    isError=True
                )]
            
            # 将结果转换为字典
            resultDict = ReactionDataToDict(reactionData)
            
            # 返回结果
            return [types.TextContent(
                type="json",
                text=json.dumps(resultDict, ensure_ascii=False)
            )]
            
        else:
            return [types.TextContent(
                type="text",
                text=f"未知工具: {name}",
                isError=True
            )]

    except Exception as e:
        logger.exception("调用工具时发生错误")
        return [types.TextContent(
            type="text",
            text=f"错误: {str(e)}",
            isError=True
        )]


# 辅助函数：运行Docker命令
async def RunDockerCommand(cmd):
    # 此处应实现Docker命令执行逻辑
    # 由于未提供具体实现，暂时返回模拟结果
    return "执行deep research查询成功。\n\n查询结果将在此显示。"


async def Main():
    """运行MCP服务器"""
    async with stdio_server() as (read_stream, write_stream):
        await app.run(
            read_stream,
            write_stream,
            app.create_initialization_options()
        )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    asyncio.run(Main())