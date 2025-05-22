#!/bin/bash

# 启动API服务器的shell脚本

# 默认参数
HOST="0.0.0.0"
PORT=5000
DEBUG=false

# 帮助信息
show_help() {
    echo "用法: $0 [选项]"
    echo
    echo "选项:"
    echo "  -h, --host HOST      设置监听主机 (默认: 0.0.0.0)"
    echo "  -p, --port PORT      设置监听端口 (默认: 5000)"
    echo "  -d, --debug          启用调试模式 (默认: 禁用)"
    echo "  --help               显示此帮助信息"
    echo
    echo "示例:"
    echo "  $0"
    echo "  $0 --host 127.0.0.1 --port 8080"
    echo "  $0 -h 0.0.0.0 -p 5000 -d"
    exit 1
}

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--host)
            HOST="$2"
            shift 2
            ;;
        -p|--port)
            PORT="$2"
            shift 2
            ;;
        -d|--debug)
            DEBUG=true
            shift
            ;;
        --help)
            show_help
            ;;
        *)
            echo "未知选项: $1"
            show_help
            ;;
    esac
done

echo "启动API服务器..."
echo "监听地址: $HOST:$PORT"
echo "调试模式: $DEBUG"

# 设置环境变量
export API_HOST="$HOST"
export API_PORT="$PORT"
export API_DEBUG="$DEBUG"

# 启动服务器
cd "$(dirname "$0")"
python app.py
