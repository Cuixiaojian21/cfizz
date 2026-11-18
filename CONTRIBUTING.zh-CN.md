# 为 cfizz 贡献代码

> Chromosome Flow Integration with Zone-Terminal Visualization
> 端到端 Hi-C 差异分析与区段末端可视化。

感谢你关注 cfizz。本指南覆盖开发环境、代码与文档规范、示例编写和 Pull Request 流程。

## 快速开始

```bash
git clone https://github.com/Cuixiaojian21/cfizz.git
cd cfizz

micromamba env create -f environment.yml
micromamba activate cfizz

pip install -r requirements.txt
pip install -e ".[all]"
```

`examples/` 下的脚本会把 `src/` 加入 `sys.path`，所以不安装也可以直接运行；日常开发建议执行 `pip install -e ".[all]"`。

项目没有独立的 `pytest` 测试套件。示例就是冒烟测试：相关示例能跑通并生成预期图表，就说明功能正常。

## 设计原则

cfizz 围绕四个承诺构建。

### 1. 端到端流程

每个分析步骤都应该可以不用手动拼装中间文件就能跑通。目标链路是：

```text
FASTQ（可选） → .mcool → Compartment / TAD / Loop → 变化特征 → 出版级图表
```

FASTQ 预处理脚本是可选的参考实现。已有 `.mcool` 的用户可以直接从这里开始。

### 2. 多组学整合

cfizz 把 Hi-C、BigWig、GTF 和 BED 都视为一等输入。轨道通过字典组合：

```python
tracks = [
    {"file": "genes.gtf", "gtf_style": "flybase", "labels": True},
    {"file": "CTCF.bw", "color": "#264653", "name": "CTCF"},
    {"file": "enhancers.bed", "color": "purple", "name": "Enhancers"},
]
```

Hi-C 面板使用同样的方式：

```python
hics = [
    {
        "file": "sample.mcool",
        "name": "sample",
        "balance": True,
        "resolution": 10_000,
        "triangle_ratio": 1.0,
        "insulation_path": "tad/1_0.sample.10000.insulation.tsv",
        "window_size": 100_000,
        "loops_path": "loops/sample.10k.loops.txt",
    }
]
```

### 3. Zone-Terminal 可视化

Zone-Terminal 模块是 cfizz 的可视化层：

- **Zonal（分区式适配）**：叠加层适配 Compartment、TAD 和 Loop 区域。
- **Terminal（末端模块）**：位于流水线末端，把分析结果转化为出版级图表。

图表暴露物理尺寸（`width_cm` / `height_cm` 或 `plot_size`）和高 dpi 输出。SVG 中内嵌的热图通过 `svg_dpi` 控制（默认 1000），放大查看时仍然清晰。

### 4. 示例 = 测试 = 文档

每个公开功能都应该在 `examples/` 下至少有一个可运行示例。一个示例同时是：

- **冒烟测试**：能跑通就说明功能正常；
- **文档**：它是标准用法；
- **演示**：生成的图就是视觉证明。

新增或修改函数时，要新增或更新对应示例。

## 项目结构

```text
cfizz/
├── src/cfizz/              # 可安装包（src layout）
│   ├── analyze/            # Compartment / TAD / Loop / O/E / distance
│   ├── api/integrated/     # quick_plot / saddle / APA / tracks
│   ├── io/                 # cooler / BED / GTF / BigWig / 路径工具
│   ├── preprocessing/      # 可选的 FASTQ → mcool 包装
│   ├── utils/              # 坐标 / 区间 / BEDPE
│   ├── viz/                # heatmap / compartment / TAD / loop / pileup
│   └── config.py
├── examples/               # 可运行示例（冒烟测试 + 文档）
├── demo/                   # 15 MB demo 数据集与产物
├── assets/figures/         # README 引用的图片
├── README.md               # 英文主文档
├── README.zh-CN.md         # 简体中文版
├── CONTRIBUTING.md         # 本文件英文版
├── CONTRIBUTING.zh-CN.md   # 本文件
└── LICENSE
```

不要把包直接放在 `cfizz/`、`io/`、`analyze/` 等位置。包位于 `src/cfizz/`；仓库根目录只放项目文件。

## 代码规范

### 命名

- 函数：`verb_noun`，例如 `compute_insulation`、`plot_multi_heatmap`。
- 类：`PascalCase`。
- 常量：`UPPER_SNAKE_CASE`。
- 配置字典键：`snake_case`。
- 内部辅助函数：以下划线开头，例如 `_helper_function`。

### 公开 API

- 公开函数放在模块层级。
- 每个公开函数都要有 docstring：一行摘要、参数、返回值，必要时给出可调用示例。
- 内部辅助函数以下划线开头。
- 优先追加可选关键字参数，避免改变已有位置参数顺序，保持向后兼容。

### 导入与类型

- 导入顺序：标准库 → 第三方 → cfizz。
- 避免 `import *`。
- 公开 API 使用类型标注。
- 兼容 Python 3.9，使用 `Optional[X]`，不要用 `X | None`。

### 图表与输出

- 用 `width_cm` / `plot_size` 控制物理尺寸，不要为出版图硬编码像素尺寸。
- 面向用户的图同时保存 PNG 和 SVG。
- 区分矢量文字/线条和内嵌栅格热图。
- SVG 中的热图优先使用 `svg_dpi=1000` 或更高；如有例外需在文档中说明。

## 示例规范

示例应：

- 尽量使用相对路径和仓库自带 demo 数据；
- 能从仓库根目录直接运行；
- 打印产物位置；
- 写入 `demo/output/`；
- 必需输入缺失时明确报错；
- 提交的示例不要包含机器相关的绝对路径（用 `CFIZZ_FASTA` 等环境变量说明覆盖方式）。

新增示例时：

1. 把脚本放到 `examples/`（或 `examples/diff/`、`examples/integrated/`）。
2. 默认使用 `demo/data/` 和 `demo/output/`，除非示例明确说明使用外部数据。
3. 运行脚本。
4. 检查生成的图或表。
5. 在 `README.md` 和 `README.zh-CN.md` 中加一小节。
6. 如果 README 要展示该图，把图片放进 `assets/figures/` 或引用已提交的 demo 产物。

## 文档规范

- `README.md` 是英文主文档；`README.zh-CN.md` 是简体中文版。
- `README.md` 使用 GitHub raw 绝对图片链接，保证 GitHub 和 PyPI 都能显示。
- `README.zh-CN.md` 使用相对路径，方便 GitHub 和本地预览。
- 两个语言版本保持结构对齐。
- 标题使用句首大写；每个文档只有一个 H1。
- 代码块标注语言（`bash`、`python`、`text`、`mermaid`）。
- 项目名统一写小写 `cfizz`，包括命令和代码。
- 正文中用 `altered feature` / “变化特征”描述 6_x；代码和文件名保留 `differential` / `diff`。
- 新增或重命名标题时同步更新目录。

## 测试与验证

项目没有 `pytest` 测试套件。通过运行相关示例验证改动：

```bash
python examples/5_1_primary_analysis_template.py
python examples/5_2_primary_analysis_visualization.py
python examples/diff/6_1_differential_compute.py
```

快速语法检查：

```bash
python -m py_compile path/to/changed_file.py
```

文档改动请检查：

- 图片路径存在；
- 链接和锚点可解析；
- 代码块成对；
- 数字和文件名与当前产物一致。

## 提交与 PR 流程

1. Fork 并 clone 仓库。
2. 创建分支：`git checkout -b feat/your-feature`。
3. 修改代码。
4. 运行并检查相关示例。
5. 用清晰的提交信息提交。
6. 发起 Pull Request，说明改了什么、为什么改、哪个示例可以验证。

提交信息：第一行简短（≤72 字符），正文说明 **为什么**。

示例：

```text
Add SVG rasterization dpi control for quick_plot_integrated

- Pass svg_dpi through to HeatmapTracks.plot
- Keep PNG dpi independent
- Verified with examples/5_2_primary_analysis_visualization.py
```

## 缺陷报告

请包含：

- cfizz 版本（`python -c "import cfizz; print(cfizz.__version__)"`）；
- Python 版本和操作系统；
- 最小复现代码或所运行的示例；
- 预期输出与实际输出；
- traceback 或相关日志片段。

Issue 地址：<https://github.com/Cuixiaojian21/cfizz/issues>。

## 社区

- Issues：<https://github.com/Cuixiaojian21/cfizz/issues>
- 维护者：cfizz developers

讨论请保持尊重和建设性。

## 许可证

提交贡献即表示你同意贡献按 MIT License 授权。详见 [LICENSE](LICENSE)。
