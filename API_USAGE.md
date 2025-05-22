# Flask API服务使用说明

本文档说明如何启动API服务器以及如何调用各个API端点。

## 启动服务器

```bash
# 进入项目目录
cd /home/zhangxiaohong/publicenv/code

# 可选：设置环境变量来配置服务器
# export API_HOST=0.0.0.0  # 监听所有网络接口（默认值）
# export API_PORT=5000     # 监听端口（默认值）
# export API_DEBUG=False   # 调试模式（默认关闭）

# 启动服务器
python app.py
```

启动后，服务器将默认监听在 `0.0.0.0:5000`，其他电脑可通过 `http://192.168.3.138:5000/` 来访问API。

## API端点

### 1. 健康检查

验证API服务器是否正常运行。

- **URL**: `/api/health`
- **方法**: `GET`
- **响应示例**:
  ```json
  {
    "status": "ok",
    "message": "服务正常运行"
  }
  ```

### 2. 反应条件预测

预测反应的最佳条件。

- **URL**: `/api/condition`
- **方法**: `POST`
- **请求体**:
  ```json
  {
  "topNResults": 1,
  "reactionDataList": [
    {
      "product": "O=CC1=CC=CC=C1C1=C2C=CC=CC2=CC=C1",
      "reactants": [
        {
          "smiles": "OB(O)C1=C2C=CC=CC2=CC=C1"
        },
        {
          "smiles": "BrC1=CC=CC=C1C=O"
        }
      ]
    }
    ]
  }
  

  ```
- **响应示例**:
  ```json
  {
    "status": "success",
    "code": 0,
    "message": "反应条件预测成功",
    "data": [
      {
        "product": "O=CC1=CC=CC=C1C1=C2C=CC=CC2=CC=C1",
        "reactants": [
          {
            "smiles": "OB(O)C1=C2C=CC=CC2=CC=C1",
            "rxnClass": null,
            "rank": 0,
            "probability": 0.0
          },
          {
            "smiles": "BrC1=CC=CC=C1C=O",
            "rxnClass": null,
            "rank": 0,
            "probability": 0.0
          }
        ],
        "conditions": [
          {
            "catalyst": ["Pd(PPh3)4"],
            "solvents": ["DMF", "H2O"],
            "reagents": ["K2CO3"],
            "temperature": "80 °C"
          },
          {
            "catalyst": ["Pd(PPh3)2Cl2"],
            "solvents": ["1,4-dioxane", "H2O"],
            "reagents": ["K3PO4"],
            "temperature": "100 °C"
          }
        ],
        "yields": []
      },
      {
        "product": "CC1=CC=CC(=C1)C1=CN=CS1",
        "reactants": [
          {
            "smiles": "BrC1=CN=CS1",
            "rxnClass": null,
            "rank": 0,
            "probability": 0.0
          },
          {
            "smiles": "CC1=CC=CC(=C1)B(O)O",
            "rxnClass": null,
            "rank": 0,
            "probability": 0.0
          }
        ],
        "conditions": [
          {
            "catalyst": ["Pd(PPh3)4"],
            "solvents": ["THF", "H2O"],
            "reagents": ["K2CO3"],
            "temperature": "80 °C"
          },
          {
            "catalyst": ["Pd(dppf)Cl2"],
            "solvents": ["DMF"],
            "reagents": ["K3PO4"],
            "temperature": "100 °C"
          }
        ],
        "yields": []
      }
    ]
  }
  ```

### 3. 半模板法逆合成预测

使用半模板法预测合成路径。

- **URL**: `/api/semi-template`
- **方法**: `POST`
- **请求体**:
  ```json
  {
    "productSmiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "topNResults": 5,
    "tryAllClasses": true,
    "maxSteps": 9,
    "beamSize": 10
  }
  ```
- **响应示例**:
  ```json
  {
    "status": "success",
    "code": 0,
    "message": "半模板法逆合成预测成功",
    "data": [
      {
        "product": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "reactants": [
          {
            "smiles": "CN1C(=O)N(C)C2=C1N=CN2C",
            "rxnClass": 7,
            "rank": 0,
            "probability": 0.0
          },
          {
            "smiles": "O",
            "rxnClass": 7,
            "rank": 0,
            "probability": 0.0
          }
        ],
        "conditions": [],
        "yields": []
      },
      {
        "product": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "reactants": [
          {
            "smiles": "CC(=O)N(C)c1nc(C)n(C)c1C(=O)N(C)C=O",
            "rxnClass": 9,
            "rank": 0,
            "probability": 0.0
          }
        ],
        "conditions": [],
        "yields": []
      }
    ]
  }
  ```

### 4. 模板法逆合成预测

使用模板法预测合成路径。

- **URL**: `/api/template`
- **方法**: `POST`
- **请求体**:
  ```json
  {
    "productSmiles": "CC(=O)c1ccc(Br)cc1",
    "topNResults": 5
  }
  ```
- **响应示例**:
  ```json
  {
    "status": "success",
    "code": 0,
    "message": "模板法逆合成预测成功",
    "data": [
      {
        "product": "CC(=O)c1ccc(Br)cc1",
        "reactants": [
          {
            "smiles": "CC(=O)c1ccccc1",
            "rxnClass": null,
            "rank": 0,
            "probability": 0.0
          },
          {
            "smiles": "Br2",
            "rxnClass": null,
            "rank": 0,
            "probability": 0.0
          }
        ],
        "conditions": [],
        "yields": []
      },
      {
        "product": "CC(=O)c1ccc(Br)cc1",
        "reactants": [
          {
            "smiles": "Brc1ccc(Br)cc1",
            "rxnClass": null,
            "rank": 0,
            "probability": 0.0
          },
          {
            "smiles": "CC(=O)Cl",
            "rxnClass": null,
            "rank": 0,
            "probability": 0.0
          }
        ],
        "conditions": [],
        "yields": []
      }
    ]
  }
  ```

### 5. 产率预测

预测反应产率。

- **URL**: `/api/yield`
- **方法**: `POST`
- **请求体**:
  ```json
  {
    "product": "CC(=O)c1ccccc1",
    "reactants": [
      {
        "smiles": "CC(O)c1ccccc1"
      },
      {
        "smiles": "O=Ic1ccccc1"
      }
    ],
    "conditions": [
      {
        "catalyst": ["CC1(C)CCCC(C)(C)N1[O]"],
        "solvents": ["O"],
        "reagents": ["CCCCCCCCCCCCOS(=O)(=O)[O-]", "[Na+]", "[Li+]", "[OH-]", "[Br-]", "[K+]"],
        "temperature": "25 °C"
      }
    ]
  }
  ```
- **响应示例**:
  ```json
  {
    "status": "success",
    "code": 0,
    "message": "产率预测成功",
    "data": {
      "product": "CC(=O)c1ccccc1",
      "reactants": [
        {
          "smiles": "CC(O)c1ccccc1",
          "rxnClass": null,
          "rank": 0,
          "probability": 0.0
        },
        {
          "smiles": "O=Ic1ccccc1",
          "rxnClass": null,
          "rank": 0,
          "probability": 0.0
        }
      ],
      "conditions": [
        {
          "catalyst": ["CC1(C)CCCC(C)(C)N1[O]"],
          "solvents": ["O"],
          "reagents": ["CCCCCCCCCCCCOS(=O)(=O)[O-]", "[Na+]", "[Li+]", "[OH-]", "[Br-]", "[K+]"],
          "temperature": "25 °C"
        }
      ],
      "yields": [
        {
          "value": 82.5
        }
      ]
    }
  }
  ```

## 错误代码

接口可能返回的错误代码及其含义：

### 成功状态
- **0**: 成功 (SUCCESS)

### 输入类错误 (100-199)
- **101**: 无效的SMILES字符串 (INPUT_INVALID_SMILES)
- **102**: 无效的结果列表 (INPUT_INVALID_RESULTS_LIST)
- **103**: 无效的topN参数值 (INPUT_INVALID_TOP_N)
- **104**: 无效的最大步骤参数 (INPUT_INVALID_MAX_STEPS)
- **105**: 无效的束宽参数 (INPUT_INVALID_BEAM_SIZE)
- **106**: 无效的反应类别 (INPUT_INVALID_RXN_CLASS)
- **107**: 缺少反应类别设置 (INPUT_MISSING_RXN_CLASS)

### 模型加载和执行类错误 (200-299)
- **201**: 模型文件未找到 (MODEL_FILE_NOT_FOUND)
- **202**: 模型加载失败 (MODEL_LOAD_FAILED)
- **203**: 模型初始化失败 (MODEL_INIT_ERROR)

### 预测执行类错误 (300-399)
- **301**: 未预测到有效反应物 (PREDICTION_NO_VALID_REACTANTS)
- **302**: 所有反应类别预测均失败 (PREDICTION_ALL_CLASSES_FAILED)
- **303**: 预测执行过程错误 (PREDICTION_EXECUTION_ERROR)
- **304**: 没有有效的预测结果 (PREDICTION_NO_VALID_RESULT)

### 未知错误
- **999**: 未知错误 (UNKNOWN_ERROR)

## 客户端调用示例（Python）

```python
import requests
import json

# 服务器地址
API_URL = "http://服务器IP:5000"

# 反应条件预测示例
def test_condition_prediction():
    url = f"{API_URL}/api/condition"
    data = {
        "topNResults": 3,
        "reactionDataList": [
            {
                "product": "O=CC1=CC=CC=C1C1=C2C=CC=CC2=CC=C1",
                "reactants": [
                    {
                        "smiles": "OB(O)C1=C2C=CC=CC2=CC=C1"
                    },
                    {
                        "smiles": "BrC1=CC=CC=C1C=O"
                    }
                ]
            }
        ]
    }
    response = requests.post(url, json=data)
    print(json.dumps(response.json(), indent=2))

# 半模板法逆合成示例
def test_semi_template_retrosynthesis():
    url = f"{API_URL}/api/semi-template"
    data = {
        "productSmiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "topNResults": 3,
        "tryAllClasses": True,
        "maxSteps": 9,
        "beamSize": 10
    }
    response = requests.post(url, json=data)
    print(json.dumps(response.json(), indent=2))

# 模板法逆合成示例
def test_template_retrosynthesis():
    url = f"{API_URL}/api/template"
    data = {
        "productSmiles": "CC(=O)c1ccc(Br)cc1",
        "topNResults": 3
    }
    response = requests.post(url, json=data)
    print(json.dumps(response.json(), indent=2))

# 产率预测示例
def test_yield_prediction():
    url = f"{API_URL}/api/yield"
    data = {
        "product": "CC(=O)c1ccccc1",
        "reactants": [
            {
                "smiles": "CC(O)c1ccccc1"
            },
            {
                "smiles": "O=Ic1ccccc1"
            }
        ],
        "conditions": [
            {
                "catalyst": ["CC1(C)CCCC(C)(C)N1[O]"],
                "solvents": ["O"],
                "reagents": ["CCCCCCCCCCCCOS(=O)(=O)[O-]", "[Na+]", "[Li+]", "[OH-]", "[Br-]", "[K+]"],
                "temperature": "25 °C"
            }
        ]
    }
    response = requests.post(url, json=data)
    print(json.dumps(response.json(), indent=2))

if __name__ == "__main__":
    test_condition_prediction()
    test_semi_template_retrosynthesis()
    test_template_retrosynthesis()
    test_yield_prediction()
```
