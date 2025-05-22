from argparse import ArgumentParser
import os
import torch
import json 

from utils import * 
from get_edit import write_edits 


# ��ȡ Test.py �ű����ڵ�Ŀ¼ (�� .../LocalRetro/scripts/)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# ��ȡ��Ŀ�ĸ�Ŀ¼ (�� .../LocalRetro/)
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)


def main_test(args):
    model_name = 'LocalRetro_%s.pth' % args['dataset']

    args['model_path'] = os.path.join(PROJECT_ROOT, 'models', model_name)
    args['config_path'] = os.path.join(PROJECT_ROOT, 'data', 'configs', args['config'])
    args['data_dir'] = os.path.join(PROJECT_ROOT, 'data', args['dataset']) 
    args['result_path'] = os.path.join(PROJECT_ROOT, 'outputs', 'raw_prediction', model_name.replace('.pth', '.txt'))

    outputs_main_dir = os.path.join(PROJECT_ROOT, 'outputs')
    outputs_raw_pred_dir = os.path.join(outputs_main_dir, 'raw_prediction')
    mkdir_p(outputs_main_dir)
    mkdir_p(outputs_raw_pred_dir)

    args = init_featurizer(args) # ���������������Ҫ�޸�·���������ڲ�Ҳʹ������PROJECT_ROOT���߼�
    args = get_site_templates(args) 
    model = load_model(args)
    test_loader = load_dataloader(args)
    write_edits(args, model, test_loader)
    return
