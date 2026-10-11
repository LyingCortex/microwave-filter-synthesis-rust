#!/usr/bin/env bash
# =============================================================================
# mfs — 构建 vendored 分发布 wheel（跨 glibc 兼容）
#
# 背景（fuqtgui CI 2026-10-08 全红事故，详见 fuqtgui 的
# docs/architecture/circuit-implementation-analysis.md §13.4）：
#   maturin 按**构建机 glibc** 给 wheel 打 manylinux 标签，并按构建机符号表
#   绑定版本节点。在 glibc 较新的机器（如 Ubuntu 24.04 = 2.39）构建后，
#   分发给 glibc 较旧的消费方（如 CI 的 python:3.12-slim = bookworm 2.36）
#   两层都挂：
#     1. 安装层：manylinux_2_39 标签 -> uv/pip 直接拒装
#        ("not a supported wheel on this platform")
#     2. 加载层：即使手改标签，.so 里 weak 符号（pidfd_spawnp/pidfd_getpid）
#        在 .gnu.version 绑定 GLIBC_2.39 版本节点 + .gnu.version_r 版本需求
#        表含该条目 -> glibc 加载器对 DT_VERNEED 表**整体校验**，import 报
#        "version `GLIBC_2.39' not found"。weak 也照样卡。
#   本脚本在构建后对 wheel 内 .so 做 ELF 修复（纯字节级，二进制逻辑零改动）：
#     a. WHEEL 标签 manylinux_X -> MFS_MANYLINUX（默认 2_35）
#     b. .gnu.version：绑定版本高于目标的 weak undefined 符号 / STT_GNU_IFUNC
#        分发桩 versym -> 1 (VER_NDX_GLOBAL 无版本)；目标 glibc 上解析不到时
#        loader 给 NULL，Rust 弱依赖回退到 posix_spawn
#     c. .gnu.version_r：摘掉不再被任何 dynsym 引用的版本节点（前驱
#        vna_next 续接 + vn_cnt 递减）
#   修复器有断言保护：只解除 weak 符号的版本绑定；若存在**强**符号绑定高于
#   目标的 glibc 版本，直接 REFUSED 退出（那说明真不能降到该标签）。
#   之后做本地验证 + 目标 glibc docker 镜像验证（安装 + import + 数值对照）。
#
#   目标下限为什么是 2_35：硬符号（非 weak）真实需求最高 GLIBC_2.35
#   (hypot)；2.39 只出现在两个 weak 的 pidfd_* 上。
#
# 用法：
#   scripts/build_wheel.sh                  # 构建 + 修复 + 本地/docker 验证
#   scripts/build_wheel.sh --no-verify      # 构建 + 修复，跳过验证
#   scripts/build_wheel.sh --verify-only W  # 不构建，只验证已有 wheel
#
# 环境变量：
#   MFS_MANYLINUX   目标 manylinux 标签，默认 2_35（消费方 glibc 2.36 可装）
#   MFS_VERIFY_IMG  验证用 docker 镜像，默认 python:3.12-slim-bookworm
#   MFS_PY          本地验证用的 python（默认本仓库 .venv，无则系统 python3）
#
# 产出（dist/）：
#   mfs-<ver>-cp3XX-cp3XX-manylinux_<MFS_MANYLINUX>_x86_64.whl   # 分发用
#   mfs-<ver>-cp3XX-cp3XX-manylinux_<原标签>_x86_64.whl          # 构建原样，留档
#   分发用 wheel 拷进消费方（如 ../fuqtgui/src/vendor/，删掉旧标签文件）后，
#   在消费方仓库跑 uv lock 刷新 mfs 的 wheel path（wheel 文件本身无 hash）。
#   Windows wheel 无 manylinux 问题，继续用 build_mfs.sh --windows 构建。
#
# 依赖：uv + cargo（构建）、python3（修复/本地验证，纯标准库）、docker（验证）
# =============================================================================
set -euo pipefail

MFS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$MFS_DIR"

log()  { printf '\033[1;32m[build-wheel]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[build-wheel]\033[0m %s\n' "$*" >&2; }
die()  { printf '\033[1;31m[build-wheel]\033[0m %s\n' "$*" >&2; exit 1; }

TARGET="${MFS_MANYLINUX:-2_35}"
IMG="${MFS_VERIFY_IMG:-python:3.12-slim-bookworm}"
NEW_WHEEL=""
NEW_REPAIRED=""

# ── ELF 修复器（纯 python3 标准库）────────────────────────────────────────────
# 符号版本映射：VERSYM 的 15 位值 == 对应 VERNAMAUX 的 vna_other（低 15 位），
# 按值匹配，不是链内位置。
PATCHER='
import base64, hashlib, re, struct, sys, zipfile

WH, TAG = sys.argv[1], sys.argv[2]
z = zipfile.ZipFile(WH)
names = z.namelist()
data = {n: z.read(n) for n in names}
z.close()
SO = next(n for n in names if n.endswith(".so"))
WP = next(n for n in names if n.endswith("dist-info/WHEEL"))
RP = next(n for n in names if n.endswith("dist-info/RECORD"))
# 输入标签：有 patchelf 时 manylinux_X_Y（maturin 自审计）；无 patchelf 时
# 裸 linux（构建机 glibc 的 weak 符号绑定照样存在，ELF 修复照样能剥）。
src_tag = re.search(r"(?:manylinux_\d+_\d+|linux)_x86_64", WH).group(0)
dst_tag = "manylinux_" + TAG + "_x86_64"
assert src_tag != dst_tag, "wheel 已是 manylinux_" + TAG

def b64(d):
    return base64.urlsafe_b64encode(hashlib.sha256(d).digest()).rstrip(b"=").decode()

def vkey(s):
    # glibc 版本逐分量数值比较（2.9 < 2.30！元组字典序会错），补零到 4 位统一长度
    n = [int(x) for x in re.findall(r"\d+", s)] + [0, 0, 0, 0]
    return tuple(n[:4])

floor = vkey("GLIBC_" + TAG.replace("_", "."))

so = bytearray(data[SO])
e_shoff, = struct.unpack_from("<Q", so, 0x28)
e_shentsize, e_shnum, e_shstrndx = struct.unpack_from("<HHH", so, 0x3A)
shs = []
for i in range(e_shnum):
    o = e_shoff + i * e_shentsize
    n, t, f, a, off, sz, lk, info, al, es = struct.unpack_from("<IIQQQQIIQQ", so, o)
    shs.append(dict(name=n, offset=off, size=sz, entsize=es))
shstr = shs[e_shstrndx]
def secname(i):
    p = shstr["offset"] + shs[i]["name"]
    return so[p:so.index(b"\0", p)].decode()
sec = {secname(i): shs[i] for i in range(e_shnum)}
dstr, dsym, ver, vr = sec[".dynstr"], sec[".dynsym"], sec[".gnu.version"], sec[".gnu.version_r"]
dstrb = so[dstr["offset"]:dstr["offset"] + dstr["size"]]
N = dsym["size"] // dsym["entsize"]

def sname(k):
    st_name = struct.unpack_from("<I", so, dsym["offset"] + k * dsym["entsize"])[0]
    e = dstrb.index(b"\0", st_name)
    return dstrb[st_name:e].decode()

def vname(vi15):
    # 按 vna_other 值匹配（不是链内位置）
    i = 0
    while i < vr["size"]:
        v, cnt, vf, vaux, vnext = struct.unpack_from("<HHIII", so, vr["offset"] + i)
        a = vr["offset"] + i + vaux
        for _ in range(cnt):
            _, _, oth, vn, nxt = struct.unpack_from("<IHHII", so, a)
            if (oth & 0x7FFF) == vi15:
                e = dstrb.index(b"\0", vn)
                return dstrb[vn:e].decode()
            a += nxt
        if vnext == 0:
            break
        i += vnext
    raise AssertionError("verneed 节点不存在: vi15=" + str(vi15))

# 1) 收集：undefined 且绑定版本高于目标的 weak / STT_GNU_IFUNC 符号
targets = []
for k in range(N):
    o = dsym["offset"] + k * dsym["entsize"]
    st_name, st_info, st_other, st_shndx = struct.unpack_from("<IBBH", so, o)
    if st_shndx != 0:
        continue  # 只看 undefined 符号
    vi = struct.unpack_from("<H", so, ver["offset"] + k * 2)[0]
    if vi < 2:
        continue
    vn = vname(vi & 0x7FFF)
    # 只对 glibc 提供的版本节点（GLIBC_*，含 libm 的节点）套 glibc 下限；
    # GCC_3.x 属 libgcc_s.so.1，不是 glibc，不参与 manylinux 判定
    if vn.startswith("GLIBC_") and vkey(vn) > floor:
        # 可解除绑定：STB_WEAK(bind=2, st_info&0xF0==0x20) 或 STT_GNU_IFUNC
        # (type=10, glibc 的 ifunc 分发桩)。两者在缺版本时都安全回退。
        # STB_GLOBAL(bind=1, 0x10) 真强符号则 REFUSED
        if (st_info & 0xF0) != 0x20 and (st_info & 0x0F) != 0x0A:
            print("REFUSED: 强符号 " + sname(k) + " 绑定 " + vn + " 高于目标，"
                  "不能降到 manylinux_" + TAG, file=sys.stderr)
            sys.exit(1)
        targets.append((k, sname(k), vi & 0x7FFF))
print("weak 符号解除版本绑定: " + (str(targets) if targets else "无"))

# 2) 解除后不再被任何 dynsym 引用的版本节点
ref = set()
tk = {k for k, _, _ in targets}
for j in range(N):
    if j in tk:
        continue
    vi = struct.unpack_from("<H", so, ver["offset"] + j * 2)[0]
    if vi >= 2:
        ref.add(vi & 0x7FFF)
unref = {vi for _, _, vi in targets if vi not in ref}

# 3) .gnu.version -> VER_NDX_GLOBAL(1)
for k, nm_, vi in targets:
    o = ver["offset"] + k * 2
    so[o:o + 2] = struct.pack("<H", 1)
    print("  .gnu.version[%d] %s: %d -> 1 (VER_NDX_GLOBAL)" % (k, nm_, vi))

# 4) .gnu.version_r：摘掉无引用节点（游标先推进，删除不回跳）
removed = 0
i = 0
while i < vr["size"]:
    e = vr["offset"] + i
    v, cnt, vf, vaux, vnext = struct.unpack_from("<HHIII", so, e)
    a = e + vaux
    prev = None
    while a < vr["offset"] + vr["size"]:
        h, fl, oth, vn, nxt = struct.unpack_from("<IHHII", so, a)
        cur = a
        a += nxt
        if (oth & 0x7FFF) in unref:
            assert prev is not None, "VERNAUX 链首（基版本）不能无引用"
            so[prev + 12:prev + 16] = struct.pack("<I", nxt)  # 前驱 vna_next 续接
            so[e + 2:e + 4] = struct.pack("<H", cnt - 1)      # vn_cnt 递减
            removed += 1
            fn = dstrb[vf:dstrb.index(b"\0", vf)].decode()
            print("  .gnu.version_r: 摘 %s 节点 other=%d (vn_cnt %d->%d)"
                  % (fn, oth & 0x7FFF, cnt, cnt - 1))
        else:
            prev = cur
        if nxt == 0:
            break
    if vnext == 0:
        break
    i += vnext
print("  .gnu.version_r 摘除节点 %d 个" % removed)

# 5) 结构验证：无残留引用 + 链长 == vn_cnt + 无重复节点
for k in range(N):
    vi = struct.unpack_from("<H", so, ver["offset"] + k * 2)[0]
    assert vi < 2 or (vi & 0x7FFF) not in unref, "dynsym[%d] 仍引用已摘节点" % k
i = 0
while i < vr["size"]:
    e = vr["offset"] + i
    v, cnt, vf, vaux, vnext = struct.unpack_from("<HHIII", so, e)
    a = e + vaux
    seen = set()
    c = 0
    while a < vr["offset"] + vr["size"]:
        h, fl, oth, vn, nxt = struct.unpack_from("<IHHII", so, a)
        assert (oth & 0x7FFF) not in unref, "链中仍有已摘节点"
        assert (oth & 0x7FFF) not in seen, "链中节点重复"
        seen.add(oth & 0x7FFF)
        c += 1
        if nxt == 0:
            break
        a += nxt
    assert c == cnt, "verneed 链长 %d != vn_cnt %d" % (c, cnt)
    if vnext == 0:
        break
    i += vnext

# 6) 重打标签 + 重建 RECORD
wh = data[WP].replace(src_tag.encode(), dst_tag.encode())
assert src_tag.encode() not in wh and dst_tag.encode() in wh, "WHEEL 标签替换失败"
data[WP] = wh
data[SO] = bytes(so)
rec = []
for ln in data[RP].decode().splitlines():
    p = ln.split(",", 1)[0]
    rec.append("%s,sha256=%s,%d" % (p, b64(data[p]), len(data[p])) if p != RP else RP + ",,")
data[RP] = "\n".join(rec).encode()

out = WH.rsplit("/", 1)[0] + "/" + WH.rsplit("/", 1)[1].replace(src_tag, dst_tag)
with zipfile.ZipFile(out, "w") as zo:
    for n in names:
        zo.writestr(n, data[n], compress_type=zipfile.ZIP_DEFLATED)
print("RETAG: " + WH + " -> " + out)
'

k12_probe() {
  # 固定数值探针：4 阶 Chebyshev RL=20dB 的 k12。纯 Rust 计算 + glibc 数学库，
  # 与构建机/消费机 glibc 无关（libm 2.35->2.36 无算法变更），必须一致
  cat <<'EOF'
import mfs
d = mfs.bandpass(order=4, rl=20.0, center=900e6, bw=30e6)
print("MFSK12=" + str(round(float(d.matrix()[0][1]), 9)))
EOF
}

verify_local() {
  local wheel="$1"
  local py="${MFS_PY:-$MFS_DIR/.venv/bin/python}"
  [ -x "$py" ] || py="$(command -v python3 || true)"
  [ -n "$py" ] || die "找不到 python3 做本地验证"
  local tmp out
  tmp="$(mktemp -d)"
  k12_probe > "$tmp/p.py"
  "$py" -m venv "$tmp/v"
  "$tmp/v/bin/pip" install --no-cache-dir -q "$wheel" --no-deps
  out="$("$tmp/v/bin/python" "$tmp/p.py" 2>&1)" || { rm -rf "$tmp"; die "本地验证失败: $wheel"$'\n'"$out"; }
  rm -rf "$tmp"
  local num
  num="$(printf '%s\n' "$out" | grep -oE 'MFSK12=[0-9.]+' | head -1 | cut -d= -f2)"
  [ -n "$num" ] || die "本地验证: 未取到数值"$'\n'"$out"
  printf '%s' "$num"
}

verify_docker() {
  local wheel="$1" ref="$2"
  command -v docker >/dev/null 2>&1 || die "docker 验证需要 docker（--no-verify 可跳过）"
  log "docker 验证: $IMG"
  local tmp
  tmp="$(mktemp -d)"
  k12_probe > "$tmp/p.py"
  cp "$wheel" "$tmp/$(basename "$wheel")"
  if ! docker run --rm -v "$tmp":/w "$IMG" bash -c '
    set -e
    pip install --no-cache-dir -q /w/*.whl --no-deps
    python3 /w/p.py
  ' > "$tmp/out.txt" 2>&1; then
    cat "$tmp/out.txt" >&2
    rm -rf "$tmp"
    die "docker 验证失败（目标 glibc 装不上或 import 失败）"
  fi
  local got
  got="$(grep -oE 'MFSK12=[0-9.]+' "$tmp/out.txt" | head -1 | cut -d= -f2)"
  rm -rf "$tmp"
  [ -n "$got" ] || die "docker 验证: 未取到数值"
  awk -v a="$ref" -v b="$got" 'BEGIN { d = a - b; if (d < 0) d = -d; exit (d > 1e-6) ? 1 : 0 }' \
    || die "数值不一致: docker=$got 构建机=$ref"
  log "docker 验证通过 ($IMG, k12=$got)"
}

build() {
  command -v uv >/dev/null 2>&1 || die "需要 uv (curl -LsSf https://astral.sh/uv/install.sh | sh)"
  command -v cargo >/dev/null 2>&1 || die "需要 cargo (maturin 构建后端)"
  log "uv build --wheel -> dist/"
  local before f b
  before="$(ls dist/ 2>/dev/null || true)"
  uv build --wheel 2>&1 | tail -6
  # 构建机有 patchelf 时 uv/maturin 打 manylinux_*_x86_64；没有 patchelf 时
  # maturin 只打 linux_x86_64（无 manylinux 标签）。两种都要能辨认，否则
  # 构建成功却在这里 die（review 2026-10-11：WSL2 无 patchelf 的实测路径）。
  for f in dist/mfs-*-cp3*-cp3*-manylinux_*_x86_64.whl \
           dist/mfs-*-cp3*-cp3*-linux_x86_64.whl; do
    [ -e "$f" ] || continue
    b="${f##*/}"
    case "$before" in
      *"$b"*) : ;;   # 旧产物，跳过
      *) NEW_WHEEL="$f" ;;
    esac
  done
  [ -n "$NEW_WHEEL" ] || die "未能从 dist/ 辨认新构建的 wheel"
  log "构建产物: ${NEW_WHEEL##*/}"
}

repair() {
  # 只修 build() 辨认出的**新** wheel（NEW_WHEEL）。
  # 不要遍历 dist/ 里所有非目标标签的 wheel：那里还留档着历次构建的
  # manylinux_<原标签> 原件，全量重打标签会改写这些留档文件，且 NEW_REPAIRED
  # 会落在 glob 排序的最后一个文件上（未必是新产物）。
  local w="$NEW_WHEEL" out rc
  case "$w" in
    *manylinux_${TARGET}_x86_64.whl)
      NEW_REPAIRED="$w"   # 构建机本来就是目标 glibc，无需 ELF 修复
      ;;
    *)
      out="$(python3 -c "$PATCHER" "$w" "$TARGET")"
      rc=$?
      printf '%s\n' "$out" | sed 's/^/  /'
      [ "$rc" -eq 0 ] || die "ELF 修复失败: ${w##*/}"
      NEW_REPAIRED="$(printf '%s\n' "$out" | grep '^RETAG:' | awk '{print $NF}')"
      [ -e "$NEW_REPAIRED" ] || die "修复未产出 manylinux_${TARGET} wheel"
      ;;
  esac
  log "分发用 wheel: ${NEW_REPAIRED##*/}"
}

post_notes() {
  log "完成。"
  log "分发: 把 ${NEW_REPAIRED##*/} 拷入消费方 ../fuqtgui/src/vendor/（删掉旧标签文件），"
  log "      在消费方仓库跑 uv lock 刷新 mfs 的 wheel path，uv lock --check 确认。"
}

mode="${1:-build}"
case "$mode" in
  --no-verify)
    build
    repair
    warn "已跳过验证（--no-verify）。分发前建议跑一次: scripts/build_wheel.sh --verify-only <wheel>"
    post_notes
    ;;
  --verify-only)
    [ -f "${2:-}" ] || die "用法: scripts/build_wheel.sh --verify-only <wheel>"
    k12="$(verify_local "$2")"
    log "本地验证通过 (k12=$k12)"
    verify_docker "$2" "$k12"
    post_notes
    ;;
  build)
    build
    repair
    k12="$(verify_local "$NEW_REPAIRED")"
    log "本地验证通过 (k12=$k12)"
    verify_docker "$NEW_REPAIRED" "$k12"
    post_notes
    ;;
  *) die "未知参数: $mode (支持: build | --no-verify | --verify-only <wheel>)" ;;
esac
