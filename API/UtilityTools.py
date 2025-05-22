"""
API模块的公共工具函数
"""
import os
import sys
import warnings

class SuppressAllOutput:
    """上下文管理器：抑制所有的标准输出、错误和警告"""
    
    def __enter__(self):
        # 抑制标准输出和错误
        self.stdout_fd = sys.stdout
        self.stderr_fd = sys.stderr
        self.null = open(os.devnull, "w")
        sys.stdout = self.null
        sys.stderr = self.null
        
        # 抑制警告
        self.catch_warnings = warnings.catch_warnings()
        self.catch_warnings.__enter__()
        warnings.simplefilter("ignore")
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        # 恢复标准输出和错误
        sys.stdout = self.stdout_fd
        sys.stderr = self.stderr_fd
        self.null.close()
        
        # 恢复警告
        self.catch_warnings.__exit__(exc_type, exc_val, exc_tb)
