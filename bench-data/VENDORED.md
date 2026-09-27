被 vendored 的数据 — 请勿手改。

这里的文件是从权威源 filter-bench 逐字节复制的副本：

    python filter-bench/tools/vendor.py --into /mnt/c/Users/eynulai/dev/mfs/bench-data

用途：在没有兄弟目录 filter-bench 时仍能独立构建与测试。
校验分两级：
  1. 兄弟目录可用 -> 与权威源 data/ 逐字节比对；
  2. 不可用       -> 与本目录 MANIFEST.sha256 比对（抓本地手改）。

要改数据，改权威源后重新运行 vendor.py。手改这里会被校验测试抓到。
