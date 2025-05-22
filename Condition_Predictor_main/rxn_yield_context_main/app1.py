import sys
from pathlib import Path
project_root = Path(__file__).parent.resolve()  # 获取 app.py 所在目录
sys.path.append(str(project_root))

from flask import Flask, render_template, request, jsonify, url_for, redirect
import os
import uuid
from werkzeug.utils import secure_filename
from utils.utils import extract_formula_from_image
from rxn_yield_context_main.rxn_yield_context.evaluate_model.eval_utils import ReactionContextPredictor

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['ALLOWED_EXTENSIONS'] = {'txt', 'png', 'jpg', 'jpeg'}

# 确保上传目录存在
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

# 初始化 ReactionContextPredictor
predictor = ReactionContextPredictor(
    'rxn_yield_context_main/data/reaxys_output',  # class_data_path
    'rxn_yield_context_main/save_models/test_10R_first_local_10/multitask_model_epoch-80.checkpoint',  # candidate_generation_model_path
    'rxn_yield_context_main/save_models/test_10R_second_7/rxn_model_relevance_listwise_morgan_epoch-80.checkpoint',  # ranking_model_path
    cutoff_solv=0.3,
    cutoff_reag=0.3,
    verbose=False
)

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in app.config['ALLOWED_EXTENSIONS']

def process_chemical_formula(formula):
    """处理化学式并获取推荐"""
    try:
        # 使用全局预测器进行预测
        rxn_smiles_list = [formula]
        results = predictor.recommend_reaction_context(rxn_smiles_list, max_display=5)
        if not results or len(results) == 0:
            return get_default_recommendations()
            
        # 转换结果格式
        recommendations = []
        for _, row in results[0].iterrows():
            recommendations.append({
                "试剂": row['Reagent(s)'],
                "溶剂": row['Solvent(s)'],
                "温度": f"{row['Temperature']}°C"
            })
        return recommendations[:5]  # 只返回前5条推荐
    except Exception as e:
        print(f"处理化学式失败: {str(e)}")
        return get_default_recommendations()

def get_default_recommendations():
    """默认推荐"""
    return [
        {"试剂": "无", "溶剂": "水", "温度": "25°C"},
        {"试剂": "H2SO4", "溶剂": "甲醇", "温度": "80°C"}
    ]

@app.route('/')
def index():
    return render_template('index1.html')

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return jsonify({'error': '没有文件部分'}), 400
        
    file = request.files['file']
        
    if file.filename == '':
        return jsonify({'error': '没有选择文件'}), 400
        
    if file and allowed_file(file.filename):
        # 保存文件
        filename = secure_filename(file.filename)
        unique_filename = f"{uuid.uuid4()}_{filename}"
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], unique_filename)
        file.save(filepath)
                
        # 根据文件类型处理
        file_extension = filename.rsplit('.', 1)[1].lower()
                
        if file_extension in ['png', 'jpg', 'jpeg']:
            # 处理图像文件
            chemical_formula = extract_formula_from_image(filepath)
        else:
            # 处理文本文件
            with open(filepath, 'r') as f:
                chemical_formula = f.read().strip()
                
        # 通过工具函数处理化学式获取推荐
        recommendations = process_chemical_formula(chemical_formula)
                
        return jsonify({
            'formula': chemical_formula,
            'recommendations': recommendations
        })
        
    return jsonify({'error': '不允许的文件类型'}), 400

@app.route('/api/recommend', methods=['POST'])
def recommend():
    data = request.json
    formula = data.get('formula', '').strip()
    
    if not formula:
        return jsonify({'error': '请提供化学式'}), 400
    
    try:
        # 使用预处理函数处理化学式
        recommendations = process_chemical_formula(formula)
        
        return jsonify({
            'formula': formula,
            'recommendations': recommendations
        })
    except Exception as e:
        return jsonify({
            'error': f'预测失败: {str(e)}',
            'recommendations': get_default_recommendations()
        }), 500

# 添加历史记录页面路由
@app.route('/history')
def history():
    return render_template('history.html')

if __name__ == '__main__':
    app.run(debug=True)