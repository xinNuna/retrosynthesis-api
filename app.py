#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Flask API服务器
提供对四个主要API函数的网络访问:
1. PredictReactionConditions - 反应条件预测
2. SemiTemplateRetrosynthesis - 半模板法逆合成预测  
3. TemplateRetrosynthesis - 模板法逆合成预测
4. PredictReactionYield - 产率预测
"""

import os
import sys
import json
import traceback
from typing import List, Dict, Any, Optional
from flask import Flask, request, jsonify, abort

# 获取当前脚本的绝对路径
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

# 确保当前目录在PATH中
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

# 导入所需的模块
from ReactionData import ReactionData, Reactant, ReactionCondition, ReactionYield
from ErrorCodes import ErrorCode, RetroApiError
from API.ConditionAPI import PredictReactionConditions
from API.SemiTemplateRetroAPI import SemiTemplateRetrosynthesis
from API.TemplateRetroAPI import TemplateRetrosynthesis
from API.YieldPredictAPI import PredictReactionYield

# 创建Flask应用
app = Flask(__name__)

# JSON序列化和反序列化函数
def reaction_data_to_dict(reaction_data: ReactionData) -> Dict[str, Any]:
    """将ReactionData对象转换为字典，用于JSON序列化"""
    result = {
        "product": reaction_data.product,
        "reactants": [],
        "conditions": [],
        "yields": []
    }
    
    # 转换反应物
    for reactant in reaction_data.reactants:
        result["reactants"].append({
            "smiles": reactant.smiles,
            "rxnClass": reactant.rxnClass,
            "rank": reactant.rank,
            "probability": reactant.probability
        })
    
    # 转换反应条件
    for condition in reaction_data.conditions:
        result["conditions"].append({
            "catalyst": condition.catalyst,
            "solvents": condition.solvents,
            "reagents": condition.reagents,
            "temperature": condition.temperature
        })
    
    # 转换产率
    for yield_obj in reaction_data.yields:
        result["yields"].append({
            "value": yield_obj.value
        })
    
    return result

def dict_to_reaction_data(data_dict: Dict[str, Any]) -> ReactionData:
    """将字典转换为ReactionData对象，用于JSON反序列化"""
    try:
        print(f"正在转换数据: {data_dict}")  # 调试信息
        
        # 创建反应物对象列表
        reactants = []
        if "reactants" in data_dict and isinstance(data_dict["reactants"], list):
            for i, reactant_dict in enumerate(data_dict["reactants"]):
                print(f"处理反应物 {i}: {reactant_dict}")  # 调试信息
                
                if isinstance(reactant_dict, dict) and "smiles" in reactant_dict:
                    # 确保所有字段都有默认值
                    smiles = reactant_dict.get("smiles")
                    if not smiles or not isinstance(smiles, str):
                        raise ValueError(f"反应物 {i} 的 smiles 字段无效: {smiles}")
                    
                    reactant = Reactant(
                        smiles=smiles.strip(),
                        rxnClass=reactant_dict.get("rxnClass"),  # 可以是None
                        rank=int(reactant_dict.get("rank", 0)),  # 默认值为0
                        probability=float(reactant_dict.get("probability", 0.0))  # 默认概率为0
                    )
                    reactants.append(reactant)
                    print(f"成功创建反应物: {reactant.smiles}")  # 调试信息
                else:
                    raise ValueError(f"反应物 {i} 格式错误或缺少 smiles 字段")
        
        # 创建反应条件对象列表
        conditions = []
        if "conditions" in data_dict and isinstance(data_dict["conditions"], list):
            for i, condition_dict in enumerate(data_dict["conditions"]):
                print(f"处理反应条件 {i}: {condition_dict}")  # 调试信息
                
                if isinstance(condition_dict, dict):
                    # 确保catalyst是列表
                    catalyst = condition_dict.get("catalyst")
                    if catalyst is not None and not isinstance(catalyst, list):
                        if isinstance(catalyst, str):
                            catalyst = [catalyst]
                        else:
                            catalyst = []
                    
                    # 确保solvents和reagents是列表
                    solvents = condition_dict.get("solvents", [])
                    if not isinstance(solvents, list):
                        solvents = [solvents] if solvents else []
                    
                    reagents = condition_dict.get("reagents", [])
                    if not isinstance(reagents, list):
                        reagents = [reagents] if reagents else []
                    
                    condition = ReactionCondition(
                        catalyst=catalyst,
                        solvents=solvents,
                        reagents=reagents,
                        temperature=condition_dict.get("temperature")
                    )
                    conditions.append(condition)
        
        # 创建产率对象列表
        yields = []
        if "yields" in data_dict and isinstance(data_dict["yields"], list):
            for i, yield_dict in enumerate(data_dict["yields"]):
                print(f"处理产率 {i}: {yield_dict}")  # 调试信息
                
                if isinstance(yield_dict, dict):
                    yield_obj = ReactionYield(
                        value=float(yield_dict.get("value", 0.0))
                    )
                    yields.append(yield_obj)
        
        # 创建ReactionData对象
        product = data_dict.get("product")
        if not product or not isinstance(product, str):
            raise ValueError(f"产物 SMILES 无效: {product}")
        
        reaction_data = ReactionData(
            product=product.strip(),
            reactants=reactants,
            conditions=conditions,
            yields=yields
        )
        
        print(f"成功创建ReactionData对象，产物: {reaction_data.product}, 反应物数量: {len(reaction_data.reactants)}")  # 调试信息
        return reaction_data
        
    except Exception as e:
        print(f"转换ReactionData时出错: {e}")  # 调试信息
        print(f"错误详情: {traceback.format_exc()}")  # 详细错误信息
        raise ValueError(f"转换反应数据时出错: {str(e)}")

def error_code_to_message(error_code: ErrorCode) -> str:
    """将ErrorCode转换为错误消息"""
    error_messages = {
        # 成功状态
        ErrorCode.SUCCESS: "成功",
        
        # 输入类错误 (100-199)
        ErrorCode.INPUT_INVALID_SMILES: "无效的SMILES输入",
        ErrorCode.INPUT_INVALID_RESULTS_LIST: "无效的结果列表",
        ErrorCode.INPUT_INVALID_TOP_N: "无效的topN参数值",
        ErrorCode.INPUT_INVALID_MAX_STEPS: "无效的最大步骤参数",
        ErrorCode.INPUT_INVALID_BEAM_SIZE: "无效的束宽参数",
        ErrorCode.INPUT_INVALID_RXN_CLASS: "无效的反应类别",
        ErrorCode.INPUT_MISSING_RXN_CLASS: "缺少反应类别设置",
        
        # 模型加载和执行类错误 (200-299)
        ErrorCode.MODEL_FILE_NOT_FOUND: "模型文件未找到",
        ErrorCode.MODEL_LOAD_FAILED: "模型加载失败",
        ErrorCode.MODEL_INIT_ERROR: "模型初始化错误",
        
        # 预测执行类错误 (300-399)
        ErrorCode.PREDICTION_NO_VALID_REACTANTS: "未预测到有效反应物",
        ErrorCode.PREDICTION_ALL_CLASSES_FAILED: "所有反应类别预测均失败",
        ErrorCode.PREDICTION_EXECUTION_ERROR: "预测执行错误",
        ErrorCode.PREDICTION_NO_VALID_RESULT: "没有有效的预测结果",
        
        # 未知错误
        ErrorCode.UNKNOWN_ERROR: "未知错误",
    }
    return error_messages.get(error_code, f"未定义的错误码: {error_code.value}")

# API错误处理
@app.errorhandler(400)
def bad_request(error):
    return jsonify({"status": "error", "message": "请求参数错误"}), 400

@app.errorhandler(404)
def not_found(error):
    return jsonify({"status": "error", "message": "请求的资源不存在"}), 404

@app.errorhandler(500)
def server_error(error):
    return jsonify({"status": "error", "message": "服务器内部错误"}), 500

# 健康检查API
@app.route('/api/health', methods=['GET'])
def health_check():
    return jsonify({"status": "ok", "message": "服务正常运行"})

# 反应条件预测API
@app.route('/api/condition', methods=['POST'])
def predict_reaction_conditions():
    try:
        print("收到反应条件预测请求")  # 调试信息
        
        # 获取请求数据
        request_data = request.get_json()
        if not request_data:
            print("请求体为空")  # 调试信息
            return jsonify({"status": "error", "message": "请求体不能为空"}), 400
        
        print(f"请求数据: {json.dumps(request_data, indent=2, ensure_ascii=False)}")  # 调试信息
        
        # 获取topNResults参数，默认为5
        top_n_results = request_data.get("topNResults", 5)
        if not isinstance(top_n_results, int) or top_n_results <= 0:
            print(f"topNResults参数无效: {top_n_results}")  # 调试信息
            return jsonify({"status": "error", "message": "topNResults必须是正整数"}), 400
        
        # 获取反应数据列表
        reaction_data_list_json = request_data.get("reactionDataList", [])
        if not isinstance(reaction_data_list_json, list):
            print(f"reactionDataList不是列表: {type(reaction_data_list_json)}")  # 调试信息
            return jsonify({"status": "error", "message": "reactionDataList必须是列表"}), 400
        
        if not reaction_data_list_json:
            print("reactionDataList为空")  # 调试信息
            return jsonify({"status": "error", "message": "reactionDataList不能为空"}), 400
        
        # 将JSON反应数据列表转换为ReactionData对象列表
        reaction_data_list = []
        for i, reaction_data_json in enumerate(reaction_data_list_json):
            try:
                print(f"转换第{i}个反应数据")  # 调试信息
                reaction_data = dict_to_reaction_data(reaction_data_json)
                reaction_data_list.append(reaction_data)
            except Exception as e:
                print(f"转换第{i}个反应数据时出错: {e}")  # 调试信息
                return jsonify({
                    "status": "error", 
                    "message": f"转换第{i}个反应数据时出错: {str(e)}"
                }), 400
        
        print(f"成功转换{len(reaction_data_list)}个反应数据对象")  # 调试信息
        
        # 调用API
        print("开始调用PredictReactionConditions")  # 调试信息
        result_code = PredictReactionConditions(
            reactionDataList=reaction_data_list,
            topNResults=top_n_results
        )
        print(f"API调用完成，结果码: {result_code}")  # 调试信息
        
        # 检查结果
        if result_code != ErrorCode.SUCCESS:
            error_msg = error_code_to_message(result_code)
            print(f"API调用失败: {error_msg}")  # 调试信息
            return jsonify({
                "status": "error", 
                "code": result_code.value, 
                "message": error_msg
            }), 400
        
        # 将ReactionData对象列表转换回JSON
        response_data = []
        for reaction_data in reaction_data_list:
            response_data.append(reaction_data_to_dict(reaction_data))
        
        print(f"成功返回{len(response_data)}个预测结果")  # 调试信息
        return jsonify({
            "status": "success",
            "code": ErrorCode.SUCCESS.value,
            "message": "反应条件预测成功",
            "data": response_data
        })
    
    except Exception as e:
        print(f"反应条件预测API发生异常: {e}")  # 调试信息
        print(f"异常详情: {traceback.format_exc()}")  # 详细错误信息
        return jsonify({"status": "error", "message": f"服务器错误: {str(e)}"}), 500

# 半模板法逆合成预测API
@app.route('/api/semi-template', methods=['POST'])
def semi_template_retrosynthesis():
    try:
        # 获取请求数据
        request_data = request.get_json()
        if not request_data:
            return jsonify({"status": "error", "message": "请求体不能为空"}), 400
        
        # 获取必需参数
        product_smiles = request_data.get("productSmiles")
        if not product_smiles or not isinstance(product_smiles, str):
            return jsonify({"status": "error", "message": "productSmiles是必需的，且必须是字符串"}), 400
        
        # 获取可选参数
        top_n_results = request_data.get("topNResults", 5)
        rxn_class = request_data.get("rxnClass")
        try_all_classes = request_data.get("tryAllClasses", True)
        max_steps = request_data.get("maxSteps", 9)
        beam_size = request_data.get("beamSize", 10)
        
        # 参数类型检查
        if not isinstance(top_n_results, int) or top_n_results <= 0:
            return jsonify({"status": "error", "message": "topNResults必须是正整数"}), 400
        
        if rxn_class is not None and not isinstance(rxn_class, int):
            return jsonify({"status": "error", "message": "rxnClass必须是整数"}), 400
        
        if not isinstance(try_all_classes, bool):
            return jsonify({"status": "error", "message": "tryAllClasses必须是布尔值"}), 400
        
        if not isinstance(max_steps, int) or max_steps <= 0:
            return jsonify({"status": "error", "message": "maxSteps必须是正整数"}), 400
        
        if not isinstance(beam_size, int) or beam_size <= 0:
            return jsonify({"status": "error", "message": "beamSize必须是正整数"}), 400
        
        # 创建结果列表
        reaction_results_list = []
        
        # 调用API
        result_code = SemiTemplateRetrosynthesis(
            productSmiles=product_smiles,
            reactionResultsList=reaction_results_list,
            topNResults=top_n_results,
            rxnClass=rxn_class,
            tryAllClasses=try_all_classes,
            maxSteps=max_steps,
            beamSize=beam_size
        )
        
        # 检查结果
        if result_code != ErrorCode.SUCCESS:
            return jsonify({
                "status": "error", 
                "code": result_code.value, 
                "message": error_code_to_message(result_code)
            }), 400
        
        # 将ReactionData对象列表转换为JSON
        response_data = []
        for reaction_data in reaction_results_list:
            response_data.append(reaction_data_to_dict(reaction_data))
        
        return jsonify({
            "status": "success",
            "code": ErrorCode.SUCCESS.value,
            "message": "半模板法逆合成预测成功",
            "data": response_data
        })
    
    except Exception as e:
        return jsonify({"status": "error", "message": f"服务器错误: {str(e)}"}), 500

# 模板法逆合成预测API
@app.route('/api/template', methods=['POST'])
def template_retrosynthesis():
    try:
        # 获取请求数据
        request_data = request.get_json()
        if not request_data:
            return jsonify({"status": "error", "message": "请求体不能为空"}), 400
        
        # 获取必需参数
        product_smiles = request_data.get("productSmiles")
        if not product_smiles or not isinstance(product_smiles, str):
            return jsonify({"status": "error", "message": "productSmiles是必需的，且必须是字符串"}), 400
        
        # 获取可选参数
        top_n_results = request_data.get("topNResults", 10)
        
        # 参数类型检查
        if not isinstance(top_n_results, int) or top_n_results <= 0:
            return jsonify({"status": "error", "message": "topNResults必须是正整数"}), 400
        
        # 创建结果列表
        reaction_results_list = []
        
        # 调用API
        result_code = TemplateRetrosynthesis(
            productSmiles=product_smiles,
            reactionResultsList=reaction_results_list,
            topNResults=top_n_results
        )
        
        # 检查结果
        if result_code != ErrorCode.SUCCESS:
            return jsonify({
                "status": "error", 
                "code": result_code.value, 
                "message": error_code_to_message(result_code)
            }), 400
        
        # 将ReactionData对象列表转换为JSON
        response_data = []
        for reaction_data in reaction_results_list:
            response_data.append(reaction_data_to_dict(reaction_data))
        
        return jsonify({
            "status": "success",
            "code": ErrorCode.SUCCESS.value,
            "message": "模板法逆合成预测成功",
            "data": response_data
        })
    
    except Exception as e:
        return jsonify({"status": "error", "message": f"服务器错误: {str(e)}"}), 500

# 产率预测API
@app.route('/api/yield', methods=['POST'])
def predict_reaction_yield():
    try:
        # 获取请求数据
        request_data = request.get_json()
        if not request_data:
            return jsonify({"status": "error", "message": "请求体不能为空"}), 400
        
        print(f"产率预测 - 接收到请求数据: {json.dumps(request_data, ensure_ascii=False)}")
        
        # 将JSON转换为ReactionData对象
        try:
            reaction_data = dict_to_reaction_data(request_data)
        except Exception as e:
            print(f"JSON转换为ReactionData失败: {str(e)}")
            return jsonify({
                "status": "error", 
                "code": ErrorCode.INPUT_INVALID_RESULTS_LIST.value, 
                "message": f"请求数据格式错误: {str(e)}"
            }), 400
        
        # 验证必要数据
        if not reaction_data.product:
            return jsonify({
                "status": "error", 
                "code": ErrorCode.INPUT_INVALID_SMILES.value,
                "message": "产物SMILES不能为空"
            }), 400
        
        if not reaction_data.reactants:
            return jsonify({
                "status": "error", 
                "code": ErrorCode.INPUT_INVALID_SMILES.value,
                "message": "反应物列表不能为空"
            }), 400
        
        if not reaction_data.conditions:
            return jsonify({
                "status": "error", 
                "code": ErrorCode.INPUT_INVALID_SMILES.value,
                "message": "反应条件不能为空"
            }), 400
        
        print("开始调用产率预测API")
        # 调用API
        result_code = PredictReactionYield(reaction_data)
        
        # 检查结果
        if result_code != ErrorCode.SUCCESS:
            print(f"产率预测失败，错误代码: {result_code.name} ({result_code.value})")
            return jsonify({
                "status": "error", 
                "code": result_code.value, 
                "message": error_code_to_message(result_code)
            }), 400
        
        # 将ReactionData对象转换回JSON
        response_data = reaction_data_to_dict(reaction_data)
        
        print("产率预测成功")
        return jsonify({
            "status": "success",
            "code": ErrorCode.SUCCESS.value,
            "message": "产率预测成功",
            "data": response_data
        })
    
    except Exception as e:
        import traceback
        error_details = traceback.format_exc()
        print(f"产率预测API出现未捕获的异常: {str(e)}")
        print(f"详细错误信息: {error_details}")
        return jsonify({
            "status": "error", 
            "code": ErrorCode.UNKNOWN_ERROR.value,
            "message": f"服务器错误: {str(e)}"
        }), 500

# 主函数


if __name__ == '__main__':
    # 主程序入口
    print("启动Flask API服务器...")
    
    # 设置环境变量避免多进程问题
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'
    
    # 获取服务器运行配置
    host = os.environ.get('API_HOST', '0.0.0.0')  # 默认绑定到所有网络接口
    port = int(os.environ.get('API_PORT', 5000))  # 默认端口5000
    debug = False  # 关闭调试模式以避免多进程问题
    
    print(f"服务器将监听在 {host}:{port}")
    # 关键：使用单线程模式
    app.run(host=host, port=port, debug=debug, threaded=False, processes=1)