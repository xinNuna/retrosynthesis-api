#!/usr/bin/env python
import os
import sys
import argparse
import torch
import glob
import re
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

# Import the model and beam search
sys.path.append('.')  # Add current directory to path
try:
    from models import Graph2Edits, BeamSearch
except ImportError:
    print("Error: Cannot import Graph2Edits. Make sure you're running from the project root directory.")
    sys.exit(1)

# Suppress RDKit warnings
RDLogger.logger().setLevel(RDLogger.CRITICAL)

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

def find_latest_model(use_rxn_class=False):
    """Find the latest and best model checkpoint"""
    base_dir = os.path.join('experiments', 'uspto_50k')
    rxn_dir = 'with_rxn_class' if use_rxn_class else 'without_rxn_class'
    model_dir = os.path.join(base_dir, rxn_dir)
    
    if not os.path.exists(model_dir):
        print(f"Error: Model directory {model_dir} not found")
        return None
        
    # Find the latest date directory
    date_dirs = sorted([d for d in os.listdir(model_dir) if os.path.isdir(os.path.join(model_dir, d))])
    if not date_dirs:
        print(f"Error: No date directories found in {model_dir}")
        return None
        
    latest_date_dir = date_dirs[-1]
    model_path = os.path.join(model_dir, latest_date_dir)
    
    # Find the highest epoch model
    model_files = glob.glob(os.path.join(model_path, "epoch_*.pt"))
    if not model_files:
        print(f"Error: No model files found in {model_path}")
        return None
        
    # Extract epoch numbers and find the highest
    epoch_numbers = []
    for model_file in model_files:
        match = re.search(r'epoch_(\d+)\.pt', model_file)
        if match:
            epoch_numbers.append((int(match.group(1)), model_file))
    
    if not epoch_numbers:
        print(f"Error: Could not parse epoch numbers from model files in {model_path}")
        return None
        
    # Sort by epoch number and get the highest
    epoch_numbers.sort(key=lambda x: x[0])
    highest_epoch_model = epoch_numbers[-1][1]
    
    print(f"Found latest model: {highest_epoch_model}")
    return highest_epoch_model

def canonicalize(smi):
    """Canonicalize a SMILES string"""
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return smi
        mol = Chem.RemoveHs(mol)
        [a.ClearProp('molAtomMapNumber') for a in mol.GetAtoms()]
        return Chem.MolToSmiles(mol)
    except:
        print('Error canonicalizing SMILES', flush=True)
        return smi

def load_model(model_path, use_rxn_class=False):
    """Load the trained model"""
    try:
        print(f"正在加载模型: {model_path}")
        checkpoint = torch.load(model_path, map_location=torch.device(DEVICE))
        config = checkpoint['saveables']
        
        model = Graph2Edits(**config, device=DEVICE)
        model.load_state_dict(checkpoint['state'])
        model.to(DEVICE)
        model.eval()
        
        print(f"模型加载成功，运行设备: {DEVICE}")
        beam_model = BeamSearch(model=model, step_beam_size=10, beam_size=10, use_rxn_class=use_rxn_class)
        return beam_model
    except Exception as e:
        print(f"加载模型时出错: {str(e)}")
        sys.exit(1)

def prepare_product(smiles):
    """Prepare product SMILES for prediction"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"错误: 无法解析SMILES: {smiles}")
            return None
            
        # Kekulize and add atom mapping if not present
        try:
            Chem.Kekulize(mol)
        except:
            print("警告: 无法将分子转化为Kekulé结构，继续处理")
            
        # Check if molecule has atom mapping, if not add it
        has_mapping = False
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() > 0:
                has_mapping = True
                break
                
        if not has_mapping:
            print("添加原子映射编号...")
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx() + 1)
        else:
            print("分子已有原子映射编号")
                
        return Chem.MolToSmiles(mol)
    except Exception as e:
        print(f"准备产物时出错: {str(e)}")
        return None

def predict_reactants(beam_model, product_smiles, rxn_class=None, max_steps=9):
    """Predict reactants from product SMILES"""
    try:
        with torch.no_grad():
            top_k_results = beam_model.run_search(
                prod_smi=product_smiles, 
                max_steps=max_steps, 
                rxn_class=rxn_class
            )
            
        results = []
        for i, path in enumerate(top_k_results):
            if 'final_smi' in path and path['final_smi'] != 'final_smi_unmapped':
                prob = path['prob']
                results.append({
                    'rank': i + 1,
                    'probability': round(prob, 4),
                    'reactants': path['final_smi'],
                    'reactants_canonical': canonicalize(path['final_smi'])
                })
        
        return results
    except Exception as e:
        print(f"Error during prediction: {str(e)}")
        return []

def main():
    parser = argparse.ArgumentParser(description='Retrosynthesis prediction using Graph2Edits')
    parser.add_argument('--model', type=str, help='Path to the model checkpoint')
    parser.add_argument('--product', type=str, help='Product SMILES string')
    parser.add_argument('--rxn_class', type=int, default=None, help='Reaction class (optional)')
    parser.add_argument('--beam_size', type=int, default=10, help='Beam search width')
    parser.add_argument('--max_steps', type=int, default=9, help='Maximum number of edit steps')
    parser.add_argument('--use_rxn_class', action='store_true', help='Whether to use reaction class information')
    parser.add_argument('--top_k', type=int, default=5, help='Number of top predictions to show')
    
    args = parser.parse_args()
    use_rxn_class = args.use_rxn_class
    
    # Check if arguments provided, otherwise switch to interactive mode
    if args.product is None:
        # Find the latest model if not specified
        model_path = args.model if args.model else find_latest_model(use_rxn_class=use_rxn_class)
        
        # Check if model exists
        if not model_path or not os.path.exists(model_path):
            print("No model found.")
            print("Please specify a model path using --model")
            return
            
        print(f"使用模型: {model_path}")
        beam_model = load_model(model_path, use_rxn_class)
        
        # Interactive mode
        print("Graph2Edits 反合成预测系统")
        print("输入产物的SMILES结构式来预测反应物 (输入 'q' 退出)")
        
        while True:
            try:
                product_smiles = input("\n产物SMILES> ")
                if product_smiles.lower() in ['q', 'quit', 'exit']:
                    break
                    
                if not product_smiles.strip():
                    continue
                    
                prepared_smiles = prepare_product(product_smiles)
                if prepared_smiles is None:
                    continue
                    
                print(f"预处理后的SMILES: {prepared_smiles}")
                print("正在预测反应物...")
                
                results = predict_reactants(beam_model, prepared_smiles, 
                                           rxn_class=args.rxn_class, 
                                           max_steps=args.max_steps)
                
                if not results:
                    print("没有预测到有效的反应物。")
                    continue
                    
                print("\n预测的反应物:")
                top_k = min(args.top_k, len(results))
                for i in range(top_k):
                    result = results[i]
                    print(f"排名 {result['rank']} (概率: {result['probability']:.4f}):")
                    print(f"  SMILES: {result['reactants']}")
                    print(f"  标准化SMILES: {result['reactants_canonical']}")
                    print()
                    
            except KeyboardInterrupt:
                print("\n正在退出...")
                break
            except Exception as e:
                print(f"错误: {str(e)}")
    else:
        # Command line mode
        if args.model is None:
            model_path = find_latest_model(use_rxn_class=args.use_rxn_class)
            if not model_path:
                print("错误: 未找到模型。请使用 --model 指定模型路径。")
                return
            args.model = model_path
            
        beam_model = load_model(args.model, args.use_rxn_class)
        prepared_smiles = prepare_product(args.product)
        
        if prepared_smiles is None:
            return
            
        print(f"为以下产物预测反应物: {prepared_smiles}")
        results = predict_reactants(beam_model, prepared_smiles, args.rxn_class, args.max_steps)
        
        if not results:
            print("没有预测到有效的反应物。")
            return
            
        top_k = min(args.top_k, len(results))
        for i in range(top_k):
            result = results[i]
            print(f"排名 {result['rank']} (概率: {result['probability']:.4f}): {result['reactants_canonical']}")

if __name__ == "__main__":
    main()