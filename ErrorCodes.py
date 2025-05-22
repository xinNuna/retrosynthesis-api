from enum import Enum

class ErrorCode(Enum):
    # 成功
    SUCCESS = 0
    
    # SMILES / 输入类错误 (100-199)
    INPUT_INVALID_SMILES = 101        # 无效的SMILES字符串
    INPUT_INVALID_RESULTS_LIST = 102  # 无效的结果列表
    INPUT_INVALID_TOP_N = 103         # 无效的topN参数值
    INPUT_INVALID_MAX_STEPS = 104     # 无效的最大步骤参数
    INPUT_INVALID_BEAM_SIZE = 105     # 无效的束宽参数
    INPUT_INVALID_RXN_CLASS = 106     # 无效的反应类别
    INPUT_MISSING_RXN_CLASS = 107     # 缺少反应类别设置
    
    # 模型加载和执行类错误 (200-299)
    MODEL_FILE_NOT_FOUND = 201        # 模型文件未找到
    MODEL_LOAD_FAILED = 202           # 模型加载失败
    MODEL_INIT_ERROR = 203            # 模型初始化失败
    
    # 预测执行类错误 (300-399)
    PREDICTION_NO_VALID_REACTANTS = 301  # 未预测到有效反应物
    PREDICTION_ALL_CLASSES_FAILED = 302  # 所有反应类别预测均失败
    PREDICTION_EXECUTION_ERROR = 303     # 预测执行过程错误
    PREDICTION_NO_VALID_RESULT = 304     # 没有有效的预测结果
    
    # 未知错误
    UNKNOWN_ERROR = 999


class RetroApiError(Exception):
    """逆合成预测API的异常类"""
    def __init__(self, message: str, status: ErrorCode):
        self.message = message
        self.status = status
        super().__init__(f"[错误代码 {status.value}] {message}")
