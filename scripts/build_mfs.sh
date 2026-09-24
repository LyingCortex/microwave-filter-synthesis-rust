#!/usr/bin/env bash
# =============================================================================
# mfs — 多平台 Python 扩展构建脚本
#
# 构建 pyo3 扩展并产出 wheel，或安装到当前 venv。
#
#   Linux (WSL)   -> mfs*.so  (绑定本 venv 的 CPython 版本)
#   Windows 宿主  -> mfs*.pyd (绑定本 venv 的 CPython 版本)
#
# 为什么"各平台 native 编译"而非在 WSL 内交叉编译 Windows 产物：
#   - pyo3 产物绑定具体 CPython ABI，两端 venv Python 版本不同，
#     native 编译各自匹配最干净。
#   - WSL 缺 mingw-w64 且 sudo 需密码，交叉产 Windows .pyd 不可行；
#     WSL 也不可直接执行 Windows .exe。
#   - Windows 宿主已有 MSVC 工具链 (cargo.exe + stable-x86_64-pc-windows-msvc)
#     和 .venv (CPython 3.11 + maturin.exe)，在该端 native 编译即可。
#
# 用法：
#   scripts/build_mfs.sh                  # 当前平台：wheel 构建 + 装进当前 venv
#   scripts/build_mfs.sh --wheel-only      # 只构建 wheel 到 dist/，不安装
#   scripts/build_mfs.sh --linux           # 当前平台构建（Linux/WSL）
#   scripts/build_mfs.sh --windows         # 在当前 Windows 端构建（git-bash）
#   scripts/build_mfs.sh --print-commands  # 只打印两端的构建命令，不执行
#
# 可选环境变量：
#   MFS_VENV   安装目标 venv（默认 $VIRTUAL_ENV，再默认当前目录 .venv）
#   MFS_NO_INSTALL  设为 1 时跳过安装，仅构建
# =============================================================================
set -euo pipefail

MFS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$MFS_DIR"

log()  { printf '\033[1;32m[build-mfs]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[build-mfs]\033[0m %s\n' "$*" >&2; }
die()  { printf '\033[1;31m[build-mfs]\033[0m %s\n' "$*" >&2; exit 1; }

platform() {
  case "$(uname -s)" in
    Linux)  echo linux ;;
    MINGW*|MSYS*|CYGWIN*) echo windows ;;
    Darwin) echo macos ;;
    *) die "unsupported platform: $(uname -s)" ;;
  esac
}

resolve_venv() {
  local v="${MFS_VENV:-${VIRTUAL_ENV:-$MFS_DIR/.venv}}"
  [ -d "$v" ] || die "venv 不存在: $v (设 MFS_VENV 指向目标 venv)"
  echo "$v"
}

build_wheel() {
  local venv="$1"
  command -v uv >/dev/null 2>&1 || die "需要 uv (curl -LsSf https://astral.sh/uv/install.sh | sh)"
  log "uv build --wheel -> dist/"
  uv build --wheel 2>&1 | tail -6
  ls -la dist/*.whl 2>/dev/null || die "未找到产物 wheel"
}

install_into_venv() {
  local venv="$1"
  if [ "${MFS_NO_INSTALL:-0}" = "1" ]; then
    warn "MFS_NO_INSTALL=1，跳过安装"
    return 0
  fi
  log "安装到 $venv"
  VIRTUAL_ENV="$venv" uv pip install --reinstall . 2>&1 | tail -6
  "$venv/bin/python" -c "import mfs; print('  verify OK:', mfs.__file__)" \
    || "$venv/Scripts/python.exe" -c "import mfs; print('  verify OK:', mfs.__file__)"
}

do_linux() {
  [ "$(platform)" = "linux" ] || warn "当前不是 Linux/WSL；若要在本端编 Linux 产物请换环境"
  local venv; venv="$(resolve_venv)"
  build_wheel "$venv"
  install_into_venv "$venv"
}

do_windows() {
  if [ "$(platform)" = "windows" ]; then
    # git-bash 原生 Windows 端
    local venv="${MFS_VENV:-$MFS_DIR/.venv}"
    [ -f "$venv/Scripts/maturin.exe" ] || die "maturin.exe 不在 $venv/Scripts/ (先 pip install maturin)"
    log "maturin develop --release -> $venv"
    "$venv/Scripts/maturin.exe" develop --release 2>&1 | tail -15
    install_into_venv "$venv"
  else
    # WSL 无法执行 Windows .exe：只打印命令
    warn "当前是 WSL，无法直接构建 Windows 产物。在 Windows 端执行："
    cat <<'EOF'
  cd %USERPROFILE%\dev\mfs
  .venv\Scripts\maturin.exe develop --release
  # 只产 wheel：.venv\Scripts\python.exe -m uv build --wheel
EOF
  fi
}

print_commands() {
  cat <<EOF
== Linux (WSL) ==
  cd $MFS_DIR
  scripts/build_mfs.sh                 # wheel + 装进当前 venv
  MFS_VENV=<venv路径> scripts/build_mfs.sh

== Windows (PowerShell / git-bash) ==
  cd $MFS_DIR (Windows 路径)
  .venv\\Scripts\\maturin.exe develop --release
  # 或 wheel: uv build --wheel
EOF
}

mode="${1:-auto}"
case "$mode" in
  --linux)          do_linux ;;
  --windows)        do_windows ;;
  --wheel-only)     MFS_NO_INSTALL=1; do_linux ;;
  --print-commands) print_commands ;;
  auto)             [ "$(platform)" = "linux" ] && do_linux || do_windows ;;
  *) die "未知参数: $mode (支持: --linux | --windows | --wheel-only | --print-commands)" ;;
esac
log "完成。"
