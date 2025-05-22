import os, time  # 导入操作系统和时间模块，用于文件操作和时间记录
import numpy as np  # 导入NumPy库，用于数值计算

from rdkit import Chem  # 从RDKit导入Chem模块，用于处理化学分子

import torch  # 导入PyTorch库，用于深度学习模型操作
import torch.nn as nn  # 导入PyTorch的神经网络模块
import dgl  # 导入DGL库，用于图神经网络处理

from utils import predict  # 从utils模块导入predict函数，用于模型预测


# 函数功能：将预测输出转换为编辑位点索引和模板编号
def get_id_template(a, class_n):
    class_n = class_n  # 获取模板类别数量（注释中无实际作用，可能为占位符）
    edit_idx = a // class_n  # 计算编辑位点索引（整数除法）
    template = a % class_n  # 计算模板编号（取余）
    return (edit_idx, template)  # 返回编辑位点索引和模板编号的元组


# 函数功能：从模型输出中提取排名靠前的编辑位点和概率
def output2edit(out, top_num):
    class_n = out.size(-1)  # 获取输出张量的最后一个维度大小（模板类别数）
    readout = out.cpu().detach().numpy()  # 将张量移到CPU并转换为NumPy数组
    readout = readout.reshape(-1)  # 将数组展平为一维
    output_rank = np.flip(np.argsort(readout))  # 按值降序排序并获取索引
    output_rank = [r for r in output_rank if get_id_template(r, class_n)[1] != 0][:top_num]  # 过滤非零模板并取前top_num个

    selected_edit = [get_id_template(a, class_n) for a in output_rank]  # 获取排名靠前的编辑位点和模板
    selected_proba = [readout[a] for a in output_rank]  # 获取对应的概率值

    return selected_edit, selected_proba  # 返回编辑位点列表和概率列表


# 函数功能：结合原子和键的编辑预测，选择排名靠前的综合编辑
def combined_edit(graph, atom_out, bond_out, top_num):
    edit_id_a, edit_proba_a = output2edit(atom_out, top_num)  # 获取原子编辑位点和概率
    edit_id_b, edit_proba_b = output2edit(bond_out, top_num)  # 获取键编辑位点和概率
    edit_id_c = edit_id_a + edit_id_b  # 合并原子和键的编辑位点
    edit_type_c = ['a'] * top_num + ['b'] * top_num  # 创建编辑类型列表（原子为'a'，键为'b'）
    edit_proba_c = edit_proba_a + edit_proba_b  # 合并原子和键的概率
    edit_rank_c = np.flip(np.argsort(edit_proba_c))[:top_num]  # 按概率降序排序并取前top_num个
    edit_type_c = [edit_type_c[r] for r in edit_rank_c]  # 获取排序后的编辑类型
    edit_id_c = [edit_id_c[r] for r in edit_rank_c]  # 获取排序后的编辑位点
    edit_proba_c = [edit_proba_c[r] for r in edit_rank_c]  # 获取排序后的概率

    return edit_type_c, edit_id_c, edit_proba_c  # 返回综合编辑类型、位点和概率


# 函数功能：将批处理图分解为单个图，并获取节点和边的分隔点
def get_bg_partition(bg):
    sg = bg.remove_self_loop()  # 移除批处理图中的自环
    gs = dgl.unbatch(sg)  # 将批处理图分解为单个图列表
    nodes_sep = [0]  # 初始化节点分隔点列表
    edges_sep = [0]  # 初始化边分隔点列表
    for g in gs:  # 遍历每个单图
        nodes_sep.append(nodes_sep[-1] + g.num_nodes())  # 添加节点分隔点
        edges_sep.append(edges_sep[-1] + g.num_edges())  # 添加边分隔点
    return gs, nodes_sep[1:], edges_sep[1:]  # 返回单图列表、节点分隔点和边分隔点


# 函数功能：从SMILES字符串提取原子和键的SMARTS模式，用于逆合成分析
def get_edit_smarts_retro(smiles):
    mol = Chem.MolFromSmiles(smiles)  # 将SMILES字符串转换为RDKit分子对象
    A = [a.GetSymbol() for a in mol.GetAtoms()]  # 获取所有原子的化学符号列表
    B = []  # 初始化键的SMARTS模式列表
    for bond in mol.GetBonds():  # 遍历分子中的所有化学键
        bond_smart = bond.GetSmarts()  # 获取键的SMARTS表示
        if bond_smart == '':  # 如果键的SMARTS为空
            if bond.GetIsAromatic():  # 检查是否为芳香键
                bond_smart = '@'  # 将芳香键表示为'@'
            else:
                bond_smart = '-'  # 将普通单键表示为'-'
        u, v = bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()  # 获取键连接的两个原子的索引
        B += ['%s%s%s' % (A[u], bond_smart, A[v]), '%s%s%s' % (A[v], bond_smart, A[u])]  # 添加双向的键SMARTS表示
    return A, B  # 返回原子符号列表和键SMARTS列表


# 函数功能：根据位点模板掩码模型预测的原子和键预测结果
def mask_prediction(smiles, atom_logits, bond_logits, site_templates):
    atom_smarts, bond_smarts = get_edit_smarts_retro(smiles)  # 获取分子中的原子和键SMARTS
    atom_mask, bond_mask = torch.zeros_like(atom_logits), torch.zeros_like(bond_logits)  # 初始化原子和键的掩码张量
    for i, smarts in enumerate(atom_smarts):  # 遍历原子SMARTS
        if smarts in site_templates:  # 如果原子SMARTS在位点模板中
            atom_mask[i][site_templates[smarts]] = 1  # 设置对应模板类别的掩码为1
    for i, smarts in enumerate(bond_smarts):  # 遍历键SMARTS
        if smarts in site_templates:  # 如果键SMARTS在位点模板中
            bond_mask[i][site_templates[smarts]] = 1  # 设置对应模板类别的掩码为1
    return atom_logits * atom_mask, bond_logits * bond_mask  # 返回掩码后的原子和键预测结果


# 函数功能：使用模型对测试数据进行预测，并将编辑预测结果写入文件
def write_edits(args, model, test_loader):
    model.eval()  # 将模型设置为评估模式
    with open(args['result_path'], 'w') as f:  # 打开结果文件进行写入
        f.write(
            'Test_id\tProduct\t%s\n' % '\t'.join(['Prediction %s' % (i + 1) for i in range(args['top_num'])]))  # 写入文件头
        with torch.no_grad():  # 禁用梯度计算以节省内存
            for batch_id, data in enumerate(test_loader):  # 遍历测试数据加载器
                # smiles_list, bg, rxns = data  # 获取SMILES列表、批处理图和反应字符串
                smiles_list, bg = data  # 获取SMILES列表、批处理图和反应字符串
                batch_atom_logits, batch_bond_logits, _ = predict(args, model, bg)  # 使用模型进行预测
                sg = bg.remove_self_loop()  # 移除批处理图中的自环
                graphs = dgl.unbatch(sg)  # 将批处理图分解为单图
                batch_atom_logits = nn.Softmax(dim=1)(batch_atom_logits)  # 对原子预测结果应用Softmax
                batch_bond_logits = nn.Softmax(dim=1)(batch_bond_logits)  # 对键预测结果应用Softmax
                graphs, nodes_sep, edges_sep = get_bg_partition(bg)  # 获取单图和节点、边分隔点
                start_node = 0  # 初始化节点起始索引
                start_edge = 0  # 初始化边起始索引
                print('\rWriting test molecule batch %s/%s' % (batch_id, len(test_loader)), end='',
                      flush=True)  # 打印处理进度
                for single_id, (graph, end_node, end_edge) in enumerate(zip(graphs, nodes_sep, edges_sep)):  # 遍历单图
                    smiles = smiles_list[single_id]  # 获取当前分子的SMILES
                    test_id = (batch_id * args['batch_size']) + single_id  # 计算测试ID
                    atom_logits = batch_atom_logits[start_node:end_node]  # 提取当前分子的原子预测
                    bond_logits = batch_bond_logits[start_edge:end_edge]  # 提取当前分子的键预测
                    atom_logits, bond_logits = mask_prediction(smiles, atom_logits, bond_logits,
                                                               args['site_templates'])  # 应用位点模板掩码
                    pred_types, pred_sites, pred_scores = combined_edit(graph, atom_logits, bond_logits,
                                                                        args['top_num'])  # 合并原子和键编辑预测
                    start_node = end_node  # 更新节点起始索引
                    start_edge = end_edge  # 更新边起始索引
                    f.write('%s\t%s\t%s\n' % (test_id, smiles, '\t'.join(
                        ['(%s, %s, %s, %.3f)' % (pred_types[i], pred_sites[i][0], pred_sites[i][1], pred_scores[i]) for
                         i in range(args['top_num'])])))  # 写入预测结果

    print()  # 打印空行
    return  # 返回（无返回值）
