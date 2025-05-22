import os  # 导入操作系统模块，用于文件路径操作
import pandas as pd  # 导入Pandas库，用于数据处理和CSV文件操作

from rdkit import Chem  # 从RDKit导入Chem模块，用于处理化学分子

import torch  # 导入PyTorch库，用于深度学习相关操作
import numpy as np  # 导入NumPy库，用于数值计算
import sklearn  # 导入Scikit-learn库（未直接使用，可能为依赖）
import dgl.backend as F  # 导入DGL后端模块（用于图神经网络操作）
from dgl.data.utils import save_graphs, load_graphs  # 从DGL导入图保存和加载工具函数


# 函数功能：规范化反应SMILES字符串，去除原子映射并返回标准化的反应表达式
def canonicalize_rxn(rxn):
    canonicalized_smiles = []  # 初始化规范化SMILES列表
    r, p = rxn.split('>>')  # 将反应字符串按'>>'分割为反应物和产物
    for s in [r, p]:  # 遍历反应物和产物
        mol = Chem.MolFromSmiles(s)  # 将SMILES字符串转换为RDKit分子对象
        [atom.SetAtomMapNum(0) for atom in mol.GetAtoms()]  # 将所有原子的映射编号设置为0
        canonicalized_smiles.append(Chem.MolToSmiles(mol))  # 将分子转换回规范化SMILES并添加到列表
    return '>>'.join(canonicalized_smiles)  # 将规范化SMILES以'>>'连接并返回


# 类功能：USPTODataset类，用于加载和处理USPTO数据集，生成分子图并支持训练、验证和测试数据
class USPTODataset(object):
    def __init__(self, args, smiles_to_graph, node_featurizer, edge_featurizer, load=True, log_every=1000):
        # 初始化函数，加载数据并设置参数
        df = pd.read_csv('%s/labeled_data.csv' % args['data_dir'])  # 读取标注数据的CSV文件
        self.train_ids = df.index[df['Split'] == 'train'].values  # 获取训练集的索引
        self.val_ids = df.index[df['Split'] == 'val'].values  # 获取验证集的索引
        self.test_ids = df.index[df['Split'] == 'test'].values  # 获取测试集的索引
        self.smiles = df['Products'].tolist()  # 获取产物SMILES列表
        self.masks = df['Mask'].tolist()  # 获取掩码列表
        self.labels = [eval(t) for t in df['Labels']]  # 将标签字符串解析为Python对象
        self.cache_file_path = '../data/saved_graphs/%s_dglgraph.bin' % args['dataset']  # 设置分子图缓存文件路径
        self._pre_process(smiles_to_graph, node_featurizer, edge_featurizer, load, log_every)  # 调用预处理函数生成或加载分子图

    # 函数功能：预处理SMILES数据，生成或加载DGL分子图
    def _pre_process(self, smiles_to_graph, node_featurizer, edge_featurizer, load, log_every):
        if os.path.exists(self.cache_file_path) and load:  # 如果缓存文件存在且允许加载
            print('Loading previously saved dgl graphs...')  # 打印加载缓存提示
            self.graphs, label_dict = load_graphs(self.cache_file_path)  # 加载缓存的DGL图
        else:
            print('Processing dgl graphs from scratch...')  # 打印从头处理提示
            self.graphs = []  # 初始化分子图列表
            for i, s in enumerate(self.smiles):  # 遍历所有SMILES
                if (i + 1) % log_every == 0:  # 每处理log_every条记录打印进度
                    print('\rProcessing molecule %d/%d' % (i + 1, len(self.smiles)), end='', flush=True)  # 打印处理进度
                self.graphs.append(smiles_to_graph(s, node_featurizer=node_featurizer,
                                                   edge_featurizer=edge_featurizer,
                                                   canonical_atom_order=False))  # 生成分子图并添加到列表
            print()  # 打印空行
            save_graphs(self.cache_file_path, self.graphs)  # 保存分子图到缓存文件

    # 函数功能：实现数据集的索引访问，返回SMILES、分子图、标签和掩码
    def __getitem__(self, item):
        return self.smiles[item], self.graphs[item], self.labels[item], self.masks[item]  # 返回指定索引的SMILES、图、标签和掩码

    # 函数功能：返回数据集的长度
    def __len__(self):
        return len(self.smiles)  # 返回SMILES列表的长度（数据集大小）


# 类功能：USPTOTestDataset类，用于加载和处理USPTO测试数据集，生成分子图并支持反应预测
class USPTOTestDataset(object):
    def __init__(self, args, smiles_to_graph, node_featurizer, edge_featurizer, load=True, log_every=1000):
        # 初始化函数，加载测试数据并设置参数
        df = pd.read_csv('../data/%s/raw_test.csv' % (args['dataset']))  # 读取原始测试数据的CSV文件
        self.rxns = df['reactants>reagents>production'].tolist()  # 获取反应字符串列表
        # self.rxns = [canonicalize_rxn(rxn) for rxn in self.rxns]  # 规范化所有反应字符串
        self.smiles = [rxn.split('>>')[-1] for rxn in self.rxns]  # 提取反应中的产物SMILES
        self.cache_file_path = '../data/saved_graphs/%s_dglgraph.bin' % args['dataset']  # 设置分子图缓存文件路径（使用训练数据集的缓存）
        self._pre_process(smiles_to_graph, node_featurizer, edge_featurizer, load, log_every)  # 调用预处理函数生成或加载分子图

    # 函数功能：预处理测试数据的SMILES，生成或加载DGL分子图
    def _pre_process(self, smiles_to_graph, node_featurizer, edge_featurizer, load, log_every):
        # if os.path.exists(self.cache_file_path) and load:  # 如果缓存文件存在且允许加载
        #     print('Loading previously saved test dgl graphs...')  # 打印加载缓存提示
        #     self.graphs, label_dict = load_graphs(self.cache_file_path)  # 加载缓存的DGL图
        # else:
        print('Processing test dgl graphs from scratch...')  # 打印从头处理提示
        self.graphs = []  # 初始化分子图列表
        for i, s in enumerate(self.smiles):  # 遍历所有SMILES
            if (i + 1) % log_every == 0:  # 每处理log_every条记录打印进度
                print('Processing molecule %d/%d' % (i + 1, len(self.smiles)))  # 打印处理进度
            self.graphs.append(smiles_to_graph(s, node_featurizer=node_featurizer,
                                               edge_featurizer=edge_featurizer,
                                               canonical_atom_order=False))  # 生成分子图并添加到列表
        save_graphs(self.cache_file_path, self.graphs)  # 保存分子图到缓存文件

    # 函数功能：实现测试数据集的索引访问，返回SMILES、分子图和反应字符串
    def __getitem__(self, item):
        # return self.smiles[item], self.graphs[item], self.rxns[item]  # 返回指定索引的SMILES、图和反应字符串
        return self.smiles[item], self.graphs[item]  # 返回指定索引的图

    # 函数功能：返回测试数据集的长度
    def __len__(self):
        return len(self.smiles)  # 返回SMILES列表的长度（数据集大小）