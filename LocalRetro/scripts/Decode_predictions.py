import os, sys, re  # 导入操作系统、系统和正则表达式模块
import pandas as pd  # 导入Pandas库，用于数据处理和CSV文件操作
import multiprocessing  # 导入多进程模块，用于并行处理
from tqdm import tqdm  # 导入tqdm库，用于显示进度条
from functools import partial  # 导入partial函数，用于固定函数参数
from collections import defaultdict  # 导入defaultdict，用于创建默认字典
from argparse import ArgumentParser  # 导入命令行参数解析模块

sys.path.append('../')  # 将父目录添加到系统路径以导入模块

import rdkit  # 导入RDKit化学信息学库
from rdkit import Chem, RDLogger  # 从RDKit导入化学分子处理和日志模块
from rdkit.Chem import rdChemReactions  # 导入RDKit化学反应模块

from utils import mkdir_p  # 从utils模块导入目录创建函数
from LocalTemplate.template_decoder import *  # 从LocalTemplate模块导入模板解码相关函数
import csv

# 函数功能：为指定测试ID获取前k个解码预测结果
def get_k_predictions(test_id, args):
    raw_prediction = args['raw_predictions'][test_id]  # 获取指定测试ID的原始预测数据
    all_prediction = []  # 初始化所有预测结果列表
    class_prediction = []  # 初始化基于反应类别的预测结果列表
    product = raw_prediction[0]  # 获取产物SMILES（原始预测的第一个元素）
    predictions = raw_prediction[1:]  # 获取预测编辑列表
    for prediction in predictions:  # 遍历所有预测编辑
        mol, pred_site, template, template_info, score = read_prediction(product, prediction, args['atom_templates'],
                                                                         args['bond_templates'],
                                                                         args['template_infos'])  # 解析预测编辑
        local_template = '>>'.join(
            ['(%s)' % smarts for smarts in template.split('_')[0].split('>>')])  # 提取局部模板并格式化为SMARTS
        decoded_smiles = decode_localtemplate(mol, pred_site, local_template, template_info)  # 解码生成SMILES
        try:
            decoded_smiles = decode_localtemplate(mol, pred_site, local_template,
                                                  template_info)  # 再次尝试解码SMILES（重复代码，可能为冗余）
            if decoded_smiles == None or str((decoded_smiles, score)) in all_prediction:  # 如果解码结果为空或已存在
                continue  # 跳过当前预测
        except Exception as e:  # 捕获解码过程中的异常
            continue  # 跳过当前预测
        all_prediction.append(str((decoded_smiles)))  # 将解码结果和得分添加到所有预测列表

        if args['rxn_class_given']:  # 如果提供了反应类别信息
            rxn_class = args['test_rxn_class'][test_id]  # 获取测试ID对应的反应类别
            if template in args['templates_class'][str(rxn_class)].values:  # 如果模板属于该反应类别
                class_prediction.append(str((decoded_smiles)))  # 添加到类别预测列表
            if len(class_prediction) >= args['top_k']:  # 如果类别预测数量达到top_k
                break  # 退出循环

        elif len(all_prediction) >= args['top_k']:  # 如果未提供类别信息且所有预测数量达到top_k
            break  # 退出循环
    return (test_id, (all_prediction, class_prediction))  # 返回测试ID和预测结果元组


# 函数功能：主函数，加载模板和预测数据，执行多进程解码并保存结果
def main_decode(args):
    atom_templates = pd.read_csv('../data/%s/atom_templates.csv' % args['dataset'])  # 读取原子模板CSV文件
    bond_templates = pd.read_csv('../data/%s/bond_templates.csv' % args['dataset'])  # 读取键模板CSV文件
    template_infos = pd.read_csv('../data/%s/template_infos.csv' % args['dataset'])  # 读取模板信息CSV文件
    class_test = '../data/%s/class_test.csv' % args['dataset']  # 设置测试反应类别文件路径
    if os.path.exists(class_test):  # 如果反应类别文件存在
        args['rxn_class_given'] = True  # 设置反应类别可用标志
        args['templates_class'] = pd.read_csv('../data/%s/template_rxnclass.csv' % args['dataset'])  # 读取模板与反应类别的映射
        args['test_rxn_class'] = pd.read_csv(class_test)['class']  # 读取测试数据的反应类别
    else:
        args['rxn_class_given'] = False  # 设置反应类别不可用标志
    args['atom_templates'] = {atom_templates['Class'][i]: atom_templates['Template'][i] for i in
                              atom_templates.index}  # 创建原子模板的类别到模板映射
    args['bond_templates'] = {bond_templates['Class'][i]: bond_templates['Template'][i] for i in
                              bond_templates.index}  # 创建键模板的类别到模板映射
    args['template_infos'] = {template_infos['Template'][i]: {'edit_site': eval(template_infos['edit_site'][i]),
                                                              'change_H': eval(template_infos['change_H'][i]),
                                                              'change_C': eval(template_infos['change_C'][i]),
                                                              'change_S': eval(template_infos['change_S'][i])} for i in
                              template_infos.index}  # 创建模板信息的字典

    if args['model'] == 'default':  # 如果模型为默认值
        result_name = 'LocalRetro_%s.txt' % args['dataset']  # 设置结果文件名为数据集名称
    else:
        result_name = 'LocalRetro_%s.txt' % args['model']  # 设置结果文件名为模型名称

    prediction_file = '../outputs/raw_prediction/' + result_name  # 设置原始预测文件路径
    raw_predictions = {}  # 初始化原始预测字典
    with open(prediction_file, 'r') as f:  # 打开原始预测文件
        for line in f.readlines():  # 逐行读取文件
            seps = line.split('\t')  # 以制表符分割行
            if seps[0] == 'Test_id':  # 如果是标题行
                continue  # 跳过标题行
            raw_predictions[int(seps[0])] = seps[1:]  # 将测试ID和预测数据存入字典
    args['raw_predictions'] = raw_predictions  # 将原始预测字典添加到参数
    # 多进程处理
    result_dict = {}  # 初始化结果字典
    partial_func = partial(get_k_predictions, args=args)  # 创建部分应用函数，固定args参数
    with multiprocessing.Pool(processes=8) as pool:  # 创建8个进程的进程池
        tasks = range(len(raw_predictions))  # 创建任务索引列表
        for result in tqdm(pool.imap_unordered(partial_func, tasks), total=len(tasks),
                           desc='Decoding LocalRetro predictions'):  # 并行执行任务并显示进度
            result_dict[result[0]] = result[1]  # 将结果存入字典
    reactants = []  # 初始化所有预测数据列表
    for i in sorted(result_dict.keys()):  # 按测试ID排序遍历结果
        all_prediction, class_prediction = result_dict[i]  # 获取所有预测和类别预测
        reactants = all_prediction  # 添加测试ID和所有预测到数据列表

        print('\rDecoding LocalRetro predictions %d/%d' % (i, len(raw_predictions)), end='', flush=True)  # 打印解码进度
    return reactants


# 主程序入口
if __name__ == '__main__':
    parser = ArgumentParser('Decode Prediction')  # 创建命令行参数解析器
    parser.add_argument('-d', '--dataset', default='USPTO_50K', help='Dataset to use')  # 添加数据集参数
    parser.add_argument('-m', '--model', default='default', help='Model to use')  # 添加模型参数
    parser.add_argument('-k', '--top-k', default=50, help='Number of top predictions')  # 添加前k个预测数量参数
    args = parser.parse_args().__dict__  # 解析参数为字典
    mkdir_p('../outputs/decoded_prediction')  # 创建解码预测输出目录
    mkdir_p('../outputs/decoded_prediction_class')  # 创建基于类别的解码预测输出目录
    main_decode(args)  # 调用主函数