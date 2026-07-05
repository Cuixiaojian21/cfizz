"""

借鉴NeoLoopFinder的简洁设计 + pyGenomeTracks的配置系统

快速、简单、高效的基因组轨道可视化
"""

from dataclasses import dataclass
from typing import Dict, Any, Optional, List, Union
import logging
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# 配置matplotlib
import matplotlib
matplotlib.use('Agg')


@dataclass
class GenomeRange:
    """
    基因组区域 - 轻量级数据结构
    """
    chrom: str
    start: int
    end: int

    @property
    def size(self) -> int:
        """区域大小"""
        return self.end - self.start

    def __str__(self) -> str:
        """字符串表示"""
        return f"{self.chrom}:{self.start:,}-{self.end:,}"

    def overlaps(self, other: 'GenomeRange') -> bool:
        """检查是否与其他区域重叠"""
        return (self.chrom == other.chrom and
                self.start < other.end and
                self.end > other.start)


@dataclass
class TrackConfig:
    """简化的Track配置"""
    file: str
    type: str  # 'bigwig' | 'gtf' | 'bed'
    color: str = 'blue'
    alpha: float = 0.8
    number_of_bins: int = 700
    summary_method: str = 'mean'
    plot_type: str = 'fill'  # 'fill' | 'line' | 'points'
    line_width: float = 0  # Line width for 'line' plot type
    min_value: Optional[float] = None
    max_value: Optional[float] = None
    labels: bool = True
    fontsize: int = 5
    name: Optional[str] = None  # 自定义track名称

    # GTF-specific parameters (following pyGenomeTracks design)
    gtf_style: str = 'flybase'  # 'flybase' | 'UCSC' | 'tssarrow' | 'exonarrows'
    border_color: str = 'black'
    color_utr: str = '#666666'  # UTR color (default blue, same as exon)
    color_backbone: str = 'black'  # Backbone color
    color_arrow: str = 'black'  # Arrow color
    line_width: float = 0.5
    height_utr: float = 1.0  # UTR height relative to coding region
    height_intron: float = 0.5  # Intron height relative to exon
    arrowhead_fraction: float = 0.004  # Fraction of region width for arrowhead
    arrowhead_included: bool = False  # Whether arrowhead is included in interval
    arrow_interval: int = 2  # Distance between arrows
    prefered_name: str = 'transcript_name'  # Use transcript_name or gene_name


class SimpleTrack:
    """简化的Track - 借鉴NeoLoopFinder的设计"""

    def __init__(self, file_path: str, track_type: str, **kwargs):
        """
        初始化Track

        Args:
            file_path: 数据文件路径
            track_type: track类型 ('bigwig', 'gtf', 'bed')
            **kwargs: 配置参数
        """
        # 确保file_path是字符串
        self.file_path = str(file_path)
        self.track_type = track_type
        self.config = TrackConfig(
            file=self.file_path,
            type=track_type,
            **kwargs
        )

        # 日志
        self.log = logging.getLogger(__name__)

        # 打开文件
        if track_type == 'bigwig':
            self._open_bigwig()
            # 处理plot_type参数，支持 'line:0.5' 格式
            self._process_plot_type()
        elif track_type in ['gtf', 'bed']:
            self._read_intervals()

    def _open_bigwig(self):
        """打开BigWig文件"""
        try:
            import pyBigWig
            self.bw = pyBigWig.open(self.file_path)
        except Exception as e:
            self.log.error(f"Failed to open BigWig file: {e}")
            raise

    def _process_plot_type(self):
        """
        处理plot_type参数，支持 'line:0.5' 格式
        参考pyGenomeTracks的设计
        """
        if hasattr(self.config, 'plot_type') and isinstance(self.config.plot_type, str):
            # 检查是否包含 ':' 分隔符
            if self.config.plot_type.find(":") > 0:
                try:
                    plot_type, size_str = self.config.plot_type.split(":")
                    self.config.plot_type = plot_type
                    # 尝试将size解析为float
                    size = float(size_str)
                    # 根据plot_type设置对应的size
                    if plot_type == 'line':
                        self.config.line_width = size
                    elif plot_type == 'points':
                        # 对于points，可以添加point_size参数（如果需要）
                        pass
                except (ValueError, AttributeError) as e:
                    self.log.warning(f"Invalid plot_type format: {self.config.plot_type}. "
                                   f"Using default. Error: {e}")
                    # 保持默认值

    def _read_intervals(self):
        """读取区间数据 (GTF/BED)"""
        try:
            if self.track_type == 'gtf':
                self.interval_data = self._read_gtf()
            elif self.track_type == 'bed':
                self.interval_data = self._read_bed()
        except Exception as e:
            self.log.error(f"Failed to read {self.track_type} file: {e}")
            raise

    def _read_gtf(self) -> List[Dict[str, Any]]:
        """
        读取GTF文件
        参考 pyGenomeTracks 的实现，默认合并转录本和合并重叠外显子
        """
        import gzip

        # 第一步：读取所有特征
        raw_data = []
        open_func = gzip.open if self.file_path.endswith('.gz') else open
        mode = 'rt'

        with open_func(self.file_path, mode) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 8:
                    continue

                chrom = parts[0]
                feature = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]

                # 解析attributes
                attrs = {}
                if len(parts) > 8:
                    attr_str = parts[8]
                    pairs = attr_str.split(';')
                    for pair in pairs:
                        pair = pair.strip()
                        if not pair:
                            continue
                        if ' ' in pair:
                            key, value = pair.split(' ', 1)
                            value = value.strip('"')
                            attrs[key] = value

                raw_data.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'feature': feature,
                    'attrs': attrs
                })

        # 第二步：按基因合并转录本（默认行为）
        genes = self._merge_transcripts(raw_data)

        # 第三步：合并重叠外显子（默认行为）
        for gene in genes:
            if 'exons' in gene:
                gene['exons'] = self._merge_overlapping_exons(gene['exons'])

        return genes

    def _merge_transcripts(self, data: List[Dict]) -> List[Dict[str, Any]]:
        """
        合并同一基因的所有转录本为一个基因
        参考 pyGenomeTracks readGtf.py 的实现
        """
        # 按基因分组
        genes = {}
        for feature in data:
            # 获取基因ID
            gene_id = feature['attrs'].get('gene_id', '')
            if not gene_id:
                continue

            if gene_id not in genes:
                genes[gene_id] = {
                    'gene_id': gene_id,
                    'chrom': feature['chrom'],
                    'strand': feature['strand'],
                    'gene_name': feature['attrs'].get('gene_name', ''),
                    'gene_biotype': feature['attrs'].get('gene_biotype', ''),
                    'exons': [],
                    'cds': [],
                    'utr_5': [],
                    'utr_3': [],
                    'start_codon': [],
                    'stop_codon': []
                }

            # 收集所有特征
            if feature['feature'] == 'exon':
                genes[gene_id]['exons'].append({
                    'start': feature['start'],
                    'end': feature['end'],
                    'strand': feature['strand'],
                    'attrs': feature['attrs']
                })
            elif feature['feature'] == 'CDS':
                genes[gene_id]['cds'].append({
                    'start': feature['start'],
                    'end': feature['end'],
                    'strand': feature['strand']
                })
            elif feature['feature'] == 'five_prime_utr':
                genes[gene_id]['utr_5'].append({
                    'start': feature['start'],
                    'end': feature['end'],
                    'strand': feature['strand']
                })
            elif feature['feature'] == 'three_prime_utr':
                genes[gene_id]['utr_3'].append({
                    'start': feature['start'],
                    'end': feature['end'],
                    'strand': feature['strand']
                })
            elif feature['feature'] == 'start_codon':
                genes[gene_id]['start_codon'].append({
                    'start': feature['start'],
                    'end': feature['end'],
                    'strand': feature['strand']
                })
            elif feature['feature'] == 'stop_codon':
                genes[gene_id]['stop_codon'].append({
                    'start': feature['start'],
                    'end': feature['end'],
                    'strand': feature['strand']
                })

        # 计算基因范围
        for gene_id, gene in genes.items():
            all_features = gene['exons'] + gene['cds'] + gene['utr_5'] + gene['utr_3']
            if all_features:
                gene['start'] = min(f['start'] for f in all_features)
                gene['end'] = max(f['end'] for f in all_features)
            else:
                # 如果没有特征，跳过这个基因
                continue

        # 转换为列表
        result = [gene for gene in genes.values() if 'start' in gene]
        result.sort(key=lambda g: g['start'])

        return result

    def _merge_overlapping_exons(self, exons: List[Dict]) -> List[Dict]:
        """
        合并重叠的外显子
        参考 pyGenomeTracks readGtf.py 的实现
        """
        if not exons:
            return []

        # 按位置排序
        exons_sorted = sorted(exons, key=lambda e: e['start'])

        merged = []
        current = exons_sorted[0].copy()

        for exon in exons_sorted[1:]:
            if exon['start'] <= current['end']:
                # 重叠，合并
                current['end'] = max(current['end'], exon['end'])
            else:
                # 不重叠，保存当前并开始新的
                merged.append(current)
                current = exon.copy()

        merged.append(current)
        return merged

    def _read_bed(self) -> List[Dict[str, Any]]:
        """读取BED文件"""
        import gzip

        data = []
        open_func = gzip.open if self.file_path.endswith('.gz') else open
        mode = 'rt'

        with open_func(self.file_path, mode) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                name = parts[3] if len(parts) > 3 else ''
                score = float(parts[4]) if len(parts) > 4 and parts[4] != '.' else 0
                strand = parts[5] if len(parts) > 5 else '+'

                data.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': name,
                    'score': score,
                    'strand': strand
                })

        return data

    def plot(self, ax, region: GenomeRange):
        """
        绘制track - 借鉴NeoLoopFinder的直接matplotlib调用

        Args:
            ax: matplotlib axes
            region: 基因组区域
        """
        if self.track_type == 'bigwig':
            self._plot_bigwig(ax, region)
        elif self.track_type == 'gtf':
            self._plot_gtf(ax, region)
        elif self.track_type == 'bed':
            self._plot_bed(ax, region)

    def _plot_bigwig(self, ax, region: GenomeRange):
        """绘制BigWig数据 - 借鉴NeoLoopFinder"""
        # 读取数据
        try:
            scores_raw = self.bw.stats(
                region.chrom,
                region.start,
                region.end,
                nBins=self.config.number_of_bins,
                type=self.config.summary_method
            )
            scores = np.array(scores_raw, dtype=float)
        except Exception as e:
            self.log.warning(f"Failed to read BigWig data: {e}")
            return

        # 处理NaN和None值
        scores = np.nan_to_num(scores, nan=0.0)

        # X轴坐标 - 使用相对坐标（0-1范围）避免极大数值问题
        x = np.linspace(0, 1, len(scores))

        # 绘制 (借鉴NeoLoopFinder)
        if self.config.plot_type == 'fill':
            ax.fill_between(
                x, scores,
                color=self.config.color,
                alpha=self.config.alpha,
                edgecolor='none'  # fill模式默认无边框
            )
        elif self.config.plot_type == 'line':
            ax.plot(
                x, scores,
                color=self.config.color,
                linewidth=self.config.line_width,
                alpha=self.config.alpha
            )
        elif self.config.plot_type == 'points':
            ax.scatter(
                x, scores,
                c=self.config.color,
                alpha=self.config.alpha,
                s=1
            )

        # 设置坐标轴 - 使用相对坐标
        ax.set_xlim(0, 1)

        # Y轴范围 - 需要同时设置min和max，或者其中一个 + 自动计算另一个
        # 先获取当前ylim（如果有数据的话）
        current_ylim = ax.get_ylim()

        # 如果用户指定了min_value和max_value，都设置
        if self.config.min_value is not None and self.config.max_value is not None:
            ax.set_ylim(self.config.min_value, self.config.max_value)
        # 如果只指定了min_value，设置min，让max自动计算
        elif self.config.min_value is not None:
            ax.set_ylim(self.config.min_value, current_ylim[1])
        # 如果只指定了max_value，设置max，让min自动计算
        elif self.config.max_value is not None:
            ax.set_ylim(current_ylim[0], self.config.max_value)
        # 如果都没指定，让matplotlib自动计算

    def _plot_gtf(self, ax, region: GenomeRange):
        """绘制GTF基因 - 基于pyGenomeTracks flybase样式设计
        遵循pyGenomeTracks的设计规范
        """
        from matplotlib.patches import Rectangle, Polygon

        # 过滤在region内的基因
        genes_in_region = [
            g for g in self.interval_data
            if g['chrom'] == region.chrom and
               g['end'] > region.start and
               g['start'] < region.end
        ]

        if not genes_in_region:
            return

        # 计算Y轴位置（堆叠显示）
        y_pos = 0
        row_scale = 2.3  # pyGenomeTracks标准行间距
        max_genes_per_row = 7  # 每行最多显示7个基因，超出后换行

        for idx, gene in enumerate(genes_in_region):
            # 每7个基因后重置y_pos到0，实现多行循环显示
            if idx > 0 and idx % max_genes_per_row == 0:
                y_pos = 0
            
            # DEBUG: 打印基因和y_pos信息
            gene_name = self._get_preferred_gene_name(gene) if self.config.labels else f"gene_{idx}"
            print(f"  [DEBUG GTF] idx={idx}, gene={gene_name}, y_pos={y_pos:.2f}, row={idx // max_genes_per_row}")
            
            # 将基因组坐标转换为相对坐标
            start_rel = (gene['start'] - region.start) / (region.end - region.start)
            end_rel = (gene['end'] - region.start) / (region.end - region.start)

            # 基因高度 - 遵循pyGenomeTracks标准
            gene_height = 1.0
            half_height = gene_height / 2

            # 1. 绘制基因骨架线（backbone）
            # 遵循pyGenomeTracks: color_backbone, linewidth, zorder=-1
            ax.plot([start_rel, end_rel], [y_pos + half_height, y_pos + half_height],
                    color=self.config.color_backbone,
                    linewidth=self.config.line_width,
                    zorder=-1)

            # 2. 处理外显子和CDS
            # pyGenomeTracks使用_split_bed_to_blocks方法
            # 我们模拟这个过程，将exons分成coding和UTR区块
            blocks = self._split_gene_to_blocks(gene, region)

            if blocks:
                # 根据strand决定箭头方向
                if gene['strand'] == '-':
                    blocks = blocks[::-1]  # 反转顺序

                # 绘制最后一个块（可能带箭头）- 箭头在转录方向末端
                last_block = blocks.pop()
                self._draw_block_with_arrow(ax, last_block, y_pos, gene_height,
                                          gene['strand'], region_length=1.0)

                # 绘制剩余块（全部矩形）
                for block in blocks:
                    self._draw_block_rectangle(ax, block, y_pos, gene_height)

            # 3. 绘制基因标签 - 参考pyGenomeTracks，放在基因右侧
            if self.config.labels:
                gene_name = self._get_preferred_gene_name(gene)
                if gene_name:
                    # 标签放在基因右侧，垂直居中
                    # 参考pyGenomeTracks BedTrack.py:613-618的实现
                    ax.text(
                        end_rel + 0.01,  # 基因结束位置右侧偏移（减小间距）
                        y_pos + half_height,  # 垂直居中
                        gene_name,
                        fontsize=self.config.fontsize,
                        va='center',
                        ha='left',  # 左对齐
                        fontstyle=self.config.fontstyle if hasattr(self.config, 'fontstyle') else 'normal'
                    )

            y_pos += row_scale  # pyGenomeTracks标准行间距

        print(f"  [DEBUG GTF] Final y_pos={y_pos:.2f}, max_genes_per_row={max_genes_per_row}, ylim={-0.5} to {y_pos:.2f}")

        # 设置坐标轴 - 使用相对坐标
        # 动态设置ylim，适应实际基因数量（恢复旧版本逻辑）
        ax.set_xlim(0, 1)
        # T-7.18: ylim 给最后一行 gene 留 row_scale + gene_height 空间
        # 原代码 ax.set_ylim(-0.5, y_pos) 算错,导致 idx 5/6/12/13 这 4 个 gene body 被裁
        # y_pos 循环结束 = 11.5;idx 6/13 画 body 在 y=[13.8, 14.8] 需要 ylim >= 14.8
        # 公式: ylim 上限 = y_pos + row_scale + gene_height = 11.5 + 2.3 + 1.0 = 14.8
        ax.set_ylim(-0.5, y_pos + row_scale + gene_height)

    def _split_gene_to_blocks(self, gene, region):
        """将基因分割成coding和UTR区块 - 参考pyGenomeTracks _split_bed_to_blocks"""
        blocks = []

        # 获取所有exons
        if 'exons' in gene and gene['exons']:
            for exon in gene['exons']:
                exon_start = (exon['start'] - region.start) / (region.end - region.start)
                exon_end = (exon['end'] - region.start) / (region.end - region.start)

                # 检查exon是否与区域重叠
                if exon_end <= 0 or exon_start >= 1:
                    continue

                exon_start = max(0, exon_start)
                exon_end = min(1, exon_end)

                # 检查是否有CDS信息
                cds_in_exon = False
                cds_start_in_exon = exon_start
                cds_end_in_exon = exon_end

                if 'cds' in gene:
                    for cds in gene['cds']:
                        cds_start = (cds['start'] - region.start) / (region.end - region.start)
                        cds_end = (cds['end'] - region.start) / (region.end - region.start)

                        # 检查CDS是否在exon范围内
                        if cds_start < exon_end and cds_end > exon_start:
                            cds_in_exon = True
                            # 分割exon为UTR和CDS部分
                            if cds_start > exon_start:
                                # 前UTR
                                blocks.append({
                                    'start': exon_start,
                                    'end': cds_start,
                                    'type': 'UTR',
                                    'color': self.config.color_utr
                                })
                            # CDS部分
                            blocks.append({
                                'start': max(exon_start, cds_start),
                                'end': min(exon_end, cds_end),
                                'type': 'coding',
                                'color': self.config.color
                            })
                            if cds_end < exon_end:
                                # 后UTR
                                blocks.append({
                                    'start': cds_end,
                                    'end': exon_end,
                                    'type': 'UTR',
                                    'color': self.config.color_utr
                                })
                            break

                if not cds_in_exon:
                    # 整个exon都是UTR或未知类型
                    block_type = 'UTR'
                    color = self.config.color_utr

                    # 检查是否是起始/终止密码子
                    if 'start_codon' in gene:
                        for codon in gene['start_codon']:
                            if codon['start'] >= exon['start'] and codon['end'] <= exon['end']:
                                block_type = 'start_codon'
                                color = '#00FF00'
                                break

                    if 'stop_codon' in gene and block_type == 'UTR':
                        for codon in gene['stop_codon']:
                            if codon['start'] >= exon['start'] and codon['end'] <= exon['end']:
                                block_type = 'stop_codon'
                                color = '#FF0000'
                                break

                    blocks.append({
                        'start': exon_start,
                        'end': exon_end,
                        'type': block_type,
                        'color': color
                    })

        return blocks

    def _draw_block_with_arrow(self, ax, block, y_pos, gene_height, strand, region_length=1.0):
        """绘制带箭头的区块 - 参考pyGenomeTracks draw_gene_with_introns_flybase_style"""
        from matplotlib.patches import Polygon

        # 计算箭头大小 - 遵循pyGenomeTracks: arrowhead_fraction是相对于整个region的
        # 参考pyGenomeTracks: self.current_small_relative = self.properties['arrowhead_fraction'] * (end_region - start_region)
        arrowhead_length = self.config.arrowhead_fraction * region_length
        if self.config.arrowhead_included:
            arrowhead_length = min(arrowhead_length, (block['end'] - block['start']) / 2)

        # 计算区块位置和大小
        if block['type'] == 'UTR':
            height = gene_height * self.config.height_utr
            y0 = y_pos + (gene_height - height) / 2
            color = block['color']
        else:
            height = gene_height
            y0 = y_pos
            color = block['color']

        # 绘制箭头 - 遵循BED约定
        if strand != '.':
            if strand == '+':
                # 正向箭头 (+)
                # 箭头从block['start']开始，到block['end'] + arrowhead_length结束
                x0 = block['start']
                x1 = block['end']
                x2 = block['end'] + arrowhead_length

                # 5个顶点：矩形主体 + 箭头
                # 遵循pyGenomeTracks的顶点顺序
                vertices = [
                    (x0, y0),                # 左下
                    (x0, y0 + height),      # 左上
                    (x1, y0 + height),      # 右上
                    (x2, y0 + height / 2),  # 箭头尖
                    (x1, y0)                # 右下
                ]
            else:
                # 反向箭头 (-)
                # 箭头从block['start'] - arrowhead_length开始，到block['end']结束
                x0 = block['start'] - arrowhead_length
                x1 = block['start']
                x2 = block['end']

                # 5个顶点：箭头 + 矩形主体
                vertices = [
                    (x0, y0 + height / 2),  # 箭头尖
                    (x1, y0 + height),      # 左上
                    (x2, y0 + height),      # 右上
                    (x2, y0),               # 右下
                    (x1, y0)                # 左下
                ]

            polygon = Polygon(vertices, closed=True, fill=True,
                            facecolor=color,
                            edgecolor=self.config.border_color,
                            linewidth=self.config.line_width,
                            alpha=self.config.alpha)
            ax.add_patch(polygon)

    def _draw_block_rectangle(self, ax, block, y_pos, gene_height):
        """绘制矩形区块 - 参考pyGenomeTracks"""
        from matplotlib.patches import Rectangle

        # 计算区块位置和大小
        if block['type'] == 'UTR':
            height = gene_height * self.config.height_utr
            y0 = y_pos + (gene_height - height) / 2
        else:
            height = gene_height
            y0 = y_pos

        # 绘制矩形
        rect = Rectangle(
            (block['start'], y0),
            block['end'] - block['start'],
            height=height,
            facecolor=block['color'],
            edgecolor=self.config.border_color,
            linewidth=self.config.line_width,
            alpha=self.config.alpha
        )
        ax.add_patch(rect)

    def _get_preferred_gene_name(self, gene):
        """获取首选的基因名 - 参考pyGenomeTracks prefered_name"""
        if self.config.prefered_name == 'transcript_name':
            # 尝试获取转录本名
            if 'transcripts' in gene and gene['transcripts']:
                transcript = gene['transcripts'][0]
                if 'transcript_name' in transcript.get('attrs', {}):
                    return transcript['attrs']['transcript_name']
            # 回退到基因名
            return gene.get('gene_name', gene.get('gene_id', ''))
        else:
            # 使用基因名
            return gene.get('gene_name', gene.get('gene_id', ''))

    def _adjust_color(self, color: str, factor: float) -> str:
        """
        调整颜色亮度
        参考 pyGenomeTracks 的颜色处理
        """
        import matplotlib.colors as mcolors
        try:
            rgb = mcolors.to_rgb(color)
            # 将颜色调暗
            rgb_adjusted = tuple(factor * c for c in rgb)
            return mcolors.to_hex(rgb_adjusted)
        except:
            return color

    def _plot_bed(self, ax, region: GenomeRange):
        """绘制BED区间"""
        from matplotlib.patches import Rectangle

        # 过滤在region内的区间
        intervals_in_region = [
            i for i in self.interval_data
            if i['chrom'] == region.chrom and
               i['end'] > region.start and
               i['start'] < region.end
        ]

        if not intervals_in_region:
            return

        # 计算Y轴位置
        y_pos = 0
        for interval in intervals_in_region:
            # 将基因组坐标转换为相对坐标
            start_rel = (interval['start'] - region.start) / (region.end - region.start)
            end_rel = (interval['end'] - region.start) / (region.end - region.start)

            # 绘制区间
            rect = Rectangle(
                (start_rel, y_pos),
                end_rel - start_rel,
                height=1.0,
                facecolor=self.config.color,
                edgecolor='black',
                linewidth=0.5,
                alpha=self.config.alpha
            )
            ax.add_patch(rect)

            # 区间标签 - 使用相对坐标
            if self.config.labels and interval['name']:
                ax.text(
                    start_rel,
                    y_pos + 1.5,
                    interval['name'],
                    fontsize=self.config.fontsize,
                    va='bottom'
                )

            y_pos += 2.5

        # 设置坐标轴 - 使用相对坐标
        ax.set_xlim(0, 1)
        ax.set_ylim(-1, y_pos)


def calculate_track_layout(
    n_tracks: int,
    total_width_cm: float,
    total_height_cm: float,
    width: float,
    left_margin: float,
    right_margin: float,
    track_heights: List[float],
    margin_top: float = 0.8,
    margin_bottom: float = 1.2  # 增加默认底部边距以容纳x轴标签
) -> List[tuple]:
    """
    计算每个track的布局位置 - 独立布局函数

    参数:
        n_tracks: track数量
        total_width_cm: 总宽度（cm）
        total_height_cm: 总高度（cm）
        width: 主绘图区宽度（cm）
        left_margin: 左侧留白（cm）
        right_margin: 右侧留白（cm）
        track_heights: 每个track的高度列表（cm）
        margin_top: 顶部边距（cm）
        margin_bottom: 底部边距（cm，用于容纳x轴标签）

    返回:
        List[tuple]: 每个track的布局信息 (x_left, y_bottom, x_width, y_height)
    """
    # 计算X位置和宽度（相对于figure）
    x_left = left_margin / total_width_cm
    x_width = width / total_width_cm

    # 计算每个track的Y位置和高度（从上往下排列）
    axes_layout = []
    for i in range(n_tracks):
        # 从顶部开始计算Y位置
        # track 0在顶部，track 1在中间，track 2在底部
        y_bottom = (margin_top + sum(track_heights[:i])) / total_height_cm
        y_height = track_heights[i] / total_height_cm

        axes_layout.append((x_left, y_bottom, x_width, y_height))

    return axes_layout


def plot_bw_tracks(
    tracks: List[SimpleTrack],
    region: GenomeRange,
    output: Optional[str] = None,
    width: Optional[float] = None,
    left_margin: float = 0.8,
    right_margin: float = 1.6,
    track_heights: Optional[List[float]] = None,
    dpi: int = 300,
    **kwargs
):
    """
    BigWig轨道绘图API - 专门处理BigWig文件

    此函数专门用于绘制BigWig轨道（信号强度数据），会自动：
    - 显示Y轴刻度和数值
    - 处理min_value和max_value
    - 添加坐标标签

    参数:
        tracks: SimpleTrack对象列表（BigWig类型）
        region: 基因组区域
        output: 输出文件路径 (可选)
        width: 用户指定的主绘图区域宽度（cm）
        left_margin: 左侧留白宽度（cm），用于y轴、y轴刻度、tracks的title
        right_margin: 右侧留白宽度（cm）
        track_heights: 每个track的高度列表（cm），如果不指定则默认为每个1cm
        dpi: 图像分辨率（默认300 DPI）
    """
    # 1. 用户友好的BigWig tracks顺序处理
    # 用户输入的tracks顺序就是从上到下的视觉顺序
    # 例如：用户输入[track1, track2] → track1在顶部，track2在底部
    # 内部处理：将列表倒序，从下往上绘制
    tracks = tracks[::-1]  # 倒序处理

    # 2. 计算图片尺寸
    # 高度 = 每个BigWig track高度 + 上下边距
    if track_heights is None:
        track_heights = [1.0] * len(tracks)  # 默认每个BigWig track 1cm
    total_track_height = sum(track_heights)
    margin_top = 0.8
    margin_bottom = 1.2  # 增加底部边距以容纳x轴标签
    total_height_cm = total_track_height + margin_top + margin_bottom

    # 宽度 = 用户指定宽度 + 左侧留白 + 右侧留白
    if width is None:
        width = 5.0  # 默认5 cm
    total_width_cm = width + left_margin + right_margin

    # 3. 设置matplotlib参数
    # 使用Arial字体，5pt
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 5

    # 3. 创建subplot
    # 使用add_axes精确控制axes位置和大小，实现100%精确的尺寸控制
    # 注意：matplotlib的figsize参数是英寸，需要转换
    cm_to_inch = 1 / 2.54  # 1英寸 = 2.54厘米
    fig = plt.figure(figsize=(total_width_cm * cm_to_inch, total_height_cm * cm_to_inch))

    # 计算每个axes的位置和大小（使用独立布局函数）
    axes_layout = calculate_track_layout(
        n_tracks=len(tracks),
        total_width_cm=total_width_cm,
        total_height_cm=total_height_cm,
        width=width,
        left_margin=left_margin,
        right_margin=right_margin,
        track_heights=track_heights,
        margin_top=margin_top,
        margin_bottom=margin_bottom
    )

    axes = []
    for layout in axes_layout:
        # 使用add_axes精确控制axes位置和大小
        ax = fig.add_axes(layout)
        axes.append(ax)

    # 4. 绘制每个track
    for track_idx, (track, ax) in enumerate(zip(tracks, axes)):
        track.plot(ax, region)

        # 5. 为每个track添加标题（右侧，垂直居中）
        # 获取标题文本（优先使用自定义name，否则使用文件名）
        if track.config.name:
            track_name = track.config.name
        else:
            track_name = track.config.file.split('/')[-1] if hasattr(track.config, 'file') else f'Track {track_idx + 1}'

        # 使用ax.text手动设置标题位置在右侧，垂直居中
        # 获取当前y轴范围用于垂直居中
        ylim = ax.get_ylim()
        y_center = (ylim[0] + ylim[1]) / 2

        # 在tracks右侧添加标题（使用相对坐标）
        ax.text(
            1.02,  # 右侧偏移
            0.5,   # 垂直居中（相对于axes）
            track_name,
            fontsize=5,
            ha='left',  # 左对齐
            va='center',  # 垂直居中
            transform=ax.transAxes  # 使用axes相对坐标
        )

        # 6. 隐藏所有默认坐标轴
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)  # 隐藏底部边框

        # 7. Y轴处理：根据track类型区别对待
        track_config = track.config

        # GTF tracks不应该有Y轴数值（它们显示基因结构，不是数值信号）
        if track.track_type == 'gtf':
            # GTF tracks：完全隐藏Y轴
            ax.spines['left'].set_visible(False)  # 隐藏左侧边框
            ax.set_yticks([])  # 隐藏Y轴刻度

            # 设置合适的Y轴范围（不显示数值，但要有合适的显示空间）
            # GTF tracks的Y轴范围通常比较小，比如-0.3到1.5
            ax.set_ylim(-0.3, 1.5)

            # 不添加手动Y轴标签（因为不需要数值）
        else:
            # BigWig/Bed tracks：显示Y轴数值
            ax.spines['left'].set_visible(True)  # 显示左侧边框（y轴刻度线）

            # 设置y轴范围（只显示最小值和最大值）
            # 获取当前y轴范围
            ylim = ax.get_ylim()
            y_min, y_max = ylim

            # 检查用户是否指定了自定义范围
            has_custom_range = track_config.min_value is not None or track_config.max_value is not None

            # 使用用户指定的min/max，或保持当前范围
            new_min = track_config.min_value if track_config.min_value is not None else y_min
            new_max = track_config.max_value if track_config.max_value is not None else y_max

            # 只有在用户没有指定范围的情况下，才扩展y轴为坐标标签留出空间
            if not has_custom_range:
                # 智能扩展y轴范围，避免出现不必要的负值
                y_range = y_max - y_min
                y_padding = y_range * 0.2

                # 只向上扩展最大值，除非最小值本身小于0
                if y_min >= 0:
                    # 最小值 >= 0，只向上扩展，避免负值
                    new_min = y_min
                    new_max = y_max + y_padding
                else:
                    # 最小值 < 0，上下都扩展
                    new_min = y_min - y_padding
                    new_max = y_max + y_padding

            # 设置y轴范围
            ax.set_ylim(new_min, new_max)

            # 设置y轴刻度 - 只显示最小值和最大值位置的刻度线
            ax.set_yticks([new_min, new_max])  # 只在最小值和最大值位置显示刻度线
            ax.tick_params(axis='y', which='major', length=3, labelsize=5, left=True, labelleft=False)
            # 隐藏默认刻度标签，使用手动黑色标签（带偏移量）

            # 手动添加y轴标签（黑色）- 最小值向上移，最大值向下移
            # 计算y轴偏移量（相对于y轴数据范围的百分比）
            y_range = new_max - new_min
            y_offset = y_range * 0.2  # y轴偏移量：20%的数据范围

            # 最小值标签 - 向上偏移（数值增大），但不超过最大值
            label_y_min = new_min + y_offset
            ax.text(
                -0.05,  # 左侧偏移（相对于坐标轴的百分比）
                label_y_min,
                format_y_axis_value(new_min),  # 使用智能格式化函数
                ha='right', va='center',
                fontsize=5,
                color='black',  # 黑色显示手动标签
                transform=ax.get_yaxis_transform()
            )

            # 最大值标签 - 向下偏移（数值减小），但不低于最小值
            label_y_max = new_max - y_offset
            ax.text(
                -0.05,  # 左侧偏移（相对于坐标轴的百分比）
                label_y_max,
                format_y_axis_value(new_max),  # 使用智能格式化函数
                ha='right', va='center',
                fontsize=5,
                color='black',  # 黑色显示手动标签
                transform=ax.get_yaxis_transform()
            )

        # 10. 先隐藏x轴刻度和标签（所有BigWig track都隐藏）
        ax.set_xticks([])
        ax.set_xlabel('')

        # 11. 然后在最下方的BigWig track上添加坐标标签
        # 注意：由于tracks列表被倒序处理，track_idx=0对应最下方的track
        if track_idx == 0:
            # 最下方的BigWig track：添加坐标标签
            # 直接使用当前y轴范围（已经在上面设置好了）
            ylim = ax.get_ylim()
            ymin, ymax = ylim

            # 添加坐标标签 - 借鉴add_coordinate_labels的设计
            # x轴独立偏移量：5%（与y轴偏移量分开配置）
            add_coordinate_labels_for_tracks(
                ax=ax,
                xmin=0, xmax=1,
                ymin=ymin, ymax=ymax,
                start_pos=region.start,
                end_pos=region.end,
                chrom=region.chrom,
                fontsize=5,
                offset_ratio=0.2  # x轴独立偏移量：2%
            )

    plt.subplots_adjust(hspace=0.2)

    # 6. 保存
    if output:
        # 获取输出路径的基名
        from pathlib import Path
        output_path = Path(output)

        # 保存为PNG格式（高分辨率）
        png_path = output_path.with_suffix('.png')
        # 使用bbox_inches=None保持精确的尺寸控制
        plt.savefig(png_path, dpi=dpi, bbox_inches=None)
        print(f"✓ Saved PNG: {png_path}")

        # 保存为SVG格式（矢量格式）
        svg_path = output_path.with_suffix('.svg')
        plt.savefig(svg_path, bbox_inches=None)
        print(f"✓ Saved SVG: {svg_path}")

        plt.close()
    else:
        plt.show()


def add_gtf_legend(ax, fontsize=7, framealpha=0.8):
    """
    为GTF track添加图例

    参数:
        ax: matplotlib的axes对象
        fontsize: 图例字体大小
        framealpha: 图例背景透明度
    """
    import matplotlib.pyplot as plt

    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor='#2E86AB', edgecolor='black', label='Exon'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#1B4F72', edgecolor='black', label='CDS'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#CCCCCC', edgecolor='black', label='UTR')
    ]

    ax.legend(
        handles=legend_elements,
        loc='upper left',
        fontsize=fontsize,
        framealpha=framealpha
    )


def plot_gtf_tracks(
    tracks: List[SimpleTrack],
    region: GenomeRange,
    output: Optional[str] = None,
    width: Optional[float] = None,
    left_margin: float = 0.8,
    right_margin: float = 1.6,
    track_heights: Optional[List[float]] = None,
    dpi: int = 300,
    # GTF专用参数（可选）
    gene_name_offset: Optional[float] = None,  # gene名字距离基因的距离 (相对于axes的偏移)
    gene_name_vertical: Optional[float] = None,  # gene名字垂直位置 (0-1, 0.5是中间)
    **kwargs
):
    """
    GTF轨道绘图API - 专门处理GTF文件

    此函数专门用于绘制GTF轨道（基因结构数据），会自动：
    - 隐藏Y轴（基因结构是分类特征，不是数值）
    - 添加网格
    - 设置合适的坐标轴范围
    - 可选显示gene name（通过gene_name_offset和gene_name_vertical参数控制）

    参数:
        tracks: SimpleTrack对象列表（GTF类型）
        region: 基因组区域
        output: 输出文件路径 (可选)
        width: 用户指定的主绘图区域宽度（cm）
        left_margin: 左侧留白宽度（cm）
        right_margin: 右侧留白宽度（cm）
        track_heights: 每个track的高度列表（cm），如果不指定则默认为每个1cm
        dpi: 图像分辨率（默认300 DPI）
        gene_name_offset: gene名字距离基因的距离（相对于axes的x轴偏移，可选）
        gene_name_vertical: gene名字的垂直位置（0-1，可选）
    """
    # 1. 用户友好的GTF tracks顺序处理
    tracks = tracks[::-1]  # 倒序处理

    # 2. 计算图片尺寸
    if track_heights is None:
        track_heights = [1.0] * len(tracks)  # 默认每个GTF track 1cm
    total_track_height = sum(track_heights)
    margin_top = 0.8
    margin_bottom = 1.2  # 底部边距以容纳x轴标签
    total_height_cm = total_track_height + margin_top + margin_bottom

    if width is None:
        width = 5.0  # 默认5 cm
    total_width_cm = width + left_margin + right_margin

    # 3. 设置matplotlib参数
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 5

    # 4. 创建subplot
    cm_to_inch = 1 / 2.54
    fig = plt.figure(figsize=(total_width_cm * cm_to_inch, total_height_cm * cm_to_inch))

    # 5. 计算每个axes的位置
    axes_layout = calculate_track_layout(
        n_tracks=len(tracks),
        total_width_cm=total_width_cm,
        total_height_cm=total_height_cm,
        width=width,
        left_margin=left_margin,
        right_margin=right_margin,
        track_heights=track_heights,
        margin_top=margin_top,
        margin_bottom=margin_bottom
    )

    axes = []
    for layout in axes_layout:
        ax = fig.add_axes(layout)
        axes.append(ax)

    # 6. 绘制每个GTF track
    for track_idx, (track, ax) in enumerate(zip(tracks, axes)):
        # 绘制GTF
        track.plot(ax, region)

        # 7. GTF特殊处理：

        # 7.1 添加网格（基于demo_gtf_simple.py）
        ax.grid(True, axis='x', alpha=0.2)

        # 7.2 隐藏Y轴（基因结构是分类特征）
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_yticks([])

        # 7.3 设置合适的坐标轴范围
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.3, 1.5)

        # 7.4 可选显示gene name（如果用户提供了偏移量参数）
        if gene_name_offset is not None or gene_name_vertical is not None:
            # 用户想要显示gene name，使用track的name属性
            track_name = track.config.name if track.config.name else None
            if track_name:
                # 设置默认值
                offset = gene_name_offset if gene_name_offset is not None else 1.02
                vertical = gene_name_vertical if gene_name_vertical is not None else 0.5

                ax.text(
                    offset,  # 横向偏移
                    vertical,  # 垂直位置
                    track_name,
                    fontsize=5,
                    ha='left',
                    va='center',
                    transform=ax.transAxes
                )

        # 8. 隐藏x轴刻度和标签
        ax.set_xticks([])
        ax.set_xlabel('')

    # 9. 在最下方的track上添加坐标标签（基于demo_gtf_simple.py）
    # 注意：由于tracks列表被倒序处理，track_idx=0对应最下方的track
    bottom_ax = axes[0]
    ylim = bottom_ax.get_ylim()
    ymin, ymax = ylim

    add_coordinate_labels_for_tracks(
        ax=bottom_ax,
        xmin=0, xmax=1,
        ymin=ymin, ymax=ymax,
        start_pos=region.start,
        end_pos=region.end,
        chrom=region.chrom,
        fontsize=5,
        offset_ratio=0.2
    )

    # 10. 保存图像
    if output:
        output_path = Path(output)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # 保存为PNG格式
        png_path = output_path.with_suffix('.png')
        plt.savefig(png_path, dpi=dpi, bbox_inches=None)
        print(f"✓ Saved PNG: {png_path}")

        # 保存为SVG格式（矢量格式）
        svg_path = output_path.with_suffix('.svg')
        plt.savefig(svg_path, bbox_inches=None)
        print(f"✓ Saved SVG: {svg_path}")

        plt.close()
    else:
        plt.show()


def plot_bed_tracks(
    tracks: List[SimpleTrack],
    region: GenomeRange,
    output: Optional[str] = None,
    width: Optional[float] = None,
    left_margin: float = 0.8,
    right_margin: float = 1.6,
    track_heights: Optional[List[float]] = None,
    dpi: int = 300,
    **kwargs
):
    """
    BED轨道绘图API - 专门处理BED文件

    此函数专门用于绘制BED轨道（区间数据），类似于GTF但更简单，会自动：
    - 隐藏Y轴（区间是分类特征，不是数值）
    - 添加网格
    - 设置合适的坐标轴范围

    参数:
        tracks: SimpleTrack对象列表（BED类型）
        region: 基因组区域
        output: 输出文件路径 (可选)
        width: 用户指定的主绘图区域宽度（cm）
        left_margin: 左侧留白宽度（cm）
        right_margin: 右侧留白宽度（cm）
        track_heights: 每个track的高度列表（cm），如果不指定则默认为每个1cm
        dpi: 图像分辨率（默认300 DPI）
    """
    # 1. 用户友好的BED tracks顺序处理
    tracks = tracks[::-1]  # 倒序处理

    # 2. 计算图片尺寸
    if track_heights is None:
        track_heights = [1.0] * len(tracks)  # 默认每个BED track 1cm
    total_track_height = sum(track_heights)
    margin_top = 0.8
    margin_bottom = 1.2  # 底部边距以容纳x轴标签
    total_height_cm = total_track_height + margin_top + margin_bottom

    if width is None:
        width = 5.0  # 默认5 cm
    total_width_cm = width + left_margin + right_margin

    # 3. 设置matplotlib参数
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 5

    # 4. 创建subplot
    cm_to_inch = 1 / 2.54
    fig = plt.figure(figsize=(total_width_cm * cm_to_inch, total_height_cm * cm_to_inch))

    # 5. 计算每个axes的位置
    axes_layout = calculate_track_layout(
        n_tracks=len(tracks),
        total_width_cm=total_width_cm,
        total_height_cm=total_height_cm,
        width=width,
        left_margin=left_margin,
        right_margin=right_margin,
        track_heights=track_heights,
        margin_top=margin_top,
        margin_bottom=margin_bottom
    )

    axes = []
    for layout in axes_layout:
        ax = fig.add_axes(layout)
        axes.append(ax)

    # 6. 绘制每个BED track
    for track_idx, (track, ax) in enumerate(zip(tracks, axes)):
        # 绘制BED
        track.plot(ax, region)

        # 7. BED特殊处理：

        # 7.1 添加网格
        ax.grid(True, axis='x', alpha=0.2)

        # 7.2 隐藏Y轴（区间是分类特征）
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_yticks([])

        # 7.3 设置合适的坐标轴范围（与GTF类似）
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.3, 1.5)

        # 7.4 添加标题（右侧，垂直居中）
        if track.config.name:
            track_name = track.config.name
        else:
            track_name = track.config.file.split('/')[-1] if hasattr(track.config, 'file') else f'Track {track_idx + 1}'

        ylim = ax.get_ylim()
        y_center = (ylim[0] + ylim[1]) / 2

        ax.text(
            1.02,  # 右侧偏移
            0.5,   # 垂直居中（相对于axes）
            track_name,
            fontsize=5,
            ha='left',
            va='center',
            transform=ax.transAxes
        )

        # 7.5 隐藏x轴刻度和标签
        ax.set_xticks([])
        ax.set_xlabel('')

    # 8. 在最下方的track上添加坐标标签
    # 注意：由于tracks列表被倒序处理，track_idx=0对应最下方的track
    bottom_ax = axes[0]
    ylim = bottom_ax.get_ylim()
    ymin, ymax = ylim

    add_coordinate_labels_for_tracks(
        ax=bottom_ax,
        xmin=0, xmax=1,
        ymin=ymin, ymax=ymax,
        start_pos=region.start,
        end_pos=region.end,
        chrom=region.chrom,
        fontsize=5,
        offset_ratio=0.2
    )

    # 9. 保存图像
    if output:
        output_path = Path(output)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # 保存为PNG格式
        png_path = output_path.with_suffix('.png')
        plt.savefig(png_path, dpi=dpi, bbox_inches=None)
        print(f"✓ Saved PNG: {png_path}")

        # 保存为SVG格式（矢量格式）
        svg_path = output_path.with_suffix('.svg')
        plt.savefig(svg_path, bbox_inches=None)
        print(f"✓ Saved SVG: {svg_path}")

        plt.close()
    else:
        plt.show()


def plot_mixed_tracks(
    tracks: List[SimpleTrack],
    region: GenomeRange,
    output: Optional[str] = None,
    width: Optional[float] = None,
    left_margin: float = 0.8,
    right_margin: float = 1.6,
    track_heights: Optional[List[float]] = None,
    dpi: int = 300,
    **kwargs
):
    """
    混合类型轨道绘图API - 专门处理多种类型混合的tracks

    此函数专门用于绘制混合类型的tracks（BigWig + GTF + BED），会根据每个track的类型
    自动选择正确的绘制方法，并确保所有tracks在同一个figure中正确显示。

    注意：此函数使用与plot_bw_tracks和plot_gtf_tracks相同的精确布局方法，
    确保绘图风格的一致性。

    参数:
        tracks: SimpleTrack对象列表
        region: GenomeRange对象，指定要绘制的基因组区域
        output: 输出文件路径（可选）
        width: 图像宽度（可选）
        left_margin: 左边距，默认0.8
        right_margin: 右边距，默认1.6
        track_heights: 每个track的高度列表（可选）
        dpi: 图像DPI，默认300
        **kwargs: 其他绘图参数

    返回:
        无，直接保存图像或显示
    """
    # 1. 用户友好的混合tracks顺序处理
    # 用户输入的tracks顺序就是从上到下的视觉顺序
    # 例如：用户输入[track1, track2, track3] → track1在顶部，track3在底部
    # 内部处理：将列表倒序，从下往上绘制
    tracks = tracks[::-1]  # 倒序处理
    # 同时倒序track_heights，确保高度与tracks同步
    if track_heights is not None:
        track_heights = track_heights[::-1]  # 同步倒序

    # 注意：所有参数（min_value, max_value等）已经在_build_tracks_params_matrix中
    # 设置到track.config中，会随着tracks的倒序而自动同步，无需额外处理

    # 2. 设置默认参数 - 与plot_bw_tracks保持一致
    if width is None:
        width = 5.0

    if track_heights is None:
        track_heights = [1.0] * len(tracks)
    elif len(track_heights) != len(tracks):
        raise ValueError(f"track_heights length ({len(track_heights)}) must match tracks length ({len(tracks)})")

    total_track_height = sum(track_heights)
    margin_top = 0.8
    margin_bottom = 1.2  # 增加底部边距以容纳x轴标签
    total_height_cm = total_track_height + margin_top + margin_bottom

    # 宽度 = 用户指定宽度 + 左侧留白 + 右侧留白
    total_width_cm = width + left_margin + right_margin

    # 2. 设置matplotlib参数 - 与plot_bw_tracks保持一致
    # 使用Arial字体，5pt
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 5

    # 3. 创建subplot - 使用精确布局方法
    # 使用add_axes精确控制axes位置和大小，实现100%精确的尺寸控制
    cm_to_inch = 1 / 2.54  # 1英寸 = 2.54厘米
    fig = plt.figure(figsize=(total_width_cm * cm_to_inch, total_height_cm * cm_to_inch))

    # 计算每个axes的位置和大小（使用独立布局函数）
    axes_layout = calculate_track_layout(
        n_tracks=len(tracks),
        total_width_cm=total_width_cm,
        total_height_cm=total_height_cm,
        width=width,
        left_margin=left_margin,
        right_margin=right_margin,
        track_heights=track_heights,
        margin_top=margin_top,
        margin_bottom=margin_bottom
    )

    axes = []
    for layout in axes_layout:
        # 使用add_axes精确控制axes位置和大小
        ax = fig.add_axes(layout)
        axes.append(ax)

    # 4. 绘制每个track - 与专用绘制函数保持一致
    for track_idx, (track, ax) in enumerate(zip(tracks, axes)):
        # 4.1 绘制track
        track.plot(ax, region)

        # 4.2 添加track标题（右侧，垂直居中）
        # 只有在用户明确提供了name（非None）时才显示标题
        if track.config.name is not None:
            track_name = track.config.name
        else:
            track_name = None  # 不显示标题

        ylim = ax.get_ylim()
        y_center = (ylim[0] + ylim[1]) / 2

        # 只有在有track_name时才添加标题
        if track_name is not None:
            ax.text(
                1.02,  # 右侧偏移
                0.5,   # 垂直居中（相对于axes）
                track_name,
                fontsize=5,
                ha='left',
                va='center',
                transform=ax.transAxes
            )

        # 4.3 根据track类型处理坐标轴
        track_type = track.config.type

        # 隐藏所有默认坐标轴
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        # Y轴处理：根据track类型区别对待
        if track_type == 'bigwig':
            # BigWig tracks：显示Y轴数值
            ax.spines['left'].set_visible(True)  # 显示左侧边框（y轴刻度线）

            # 添加Y轴网格
            ax.grid(True, axis='y', alpha=0.2)

            # 设置y轴范围（只显示最小值和最大值）
            # 获取当前y轴范围
            ylim = ax.get_ylim()
            y_min, y_max = ylim

            # 检查用户是否指定了自定义范围
            has_custom_range = track.config.min_value is not None or track.config.max_value is not None

            # 使用用户指定的min/max，或保持当前范围
            new_min = track.config.min_value if track.config.min_value is not None else y_min
            new_max = track.config.max_value if track.config.max_value is not None else y_max

            # 只有在用户没有指定范围的情况下，才扩展y轴为坐标标签留出空间
            if not has_custom_range:
                # 智能扩展y轴范围，避免出现不必要的负值
                y_range = y_max - y_min
                y_padding = y_range * 0.2

                # 只向上扩展最大值，除非最小值本身小于0
                if y_min >= 0:
                    # 最小值 >= 0，只向上扩展，避免负值
                    new_min = y_min
                    new_max = y_max + y_padding
                else:
                    # 最小值 < 0，上下都扩展
                    new_min = y_min - y_padding
                    new_max = y_max + y_padding

            # 设置y轴范围
            ax.set_ylim(new_min, new_max)

            # 设置y轴刻度 - 只显示最小值和最大值位置的刻度线
            ax.set_yticks([new_min, new_max])  # 只在最小值和最大值位置显示刻度线
            ax.tick_params(axis='y', which='major', length=3, labelsize=5, left=True, labelleft=False)
            # 隐藏默认刻度标签，使用手动黑色标签（带偏移量）

            # 手动添加y轴标签（黑色）- 最小值向上移，最大值向下移
            # 计算y轴偏移量（相对于y轴数据范围的百分比）
            y_range = new_max - new_min
            y_offset = y_range * 0.2  # y轴偏移量：20%的数据范围

            # 最小值标签 - 向上偏移（数值增大），但不超过最大值
            label_y_min = new_min + y_offset
            ax.text(
                -0.05,  # 左侧偏移（相对于坐标轴的百分比）
                label_y_min,
                format_y_axis_value(new_min),  # 使用智能格式化函数
                ha='right', va='center',
                fontsize=5,
                color='black',  # 黑色显示手动标签
                transform=ax.get_yaxis_transform()
            )

            # 最大值标签 - 向下偏移（数值减小），但不低于最小值
            label_y_max = new_max - y_offset
            ax.text(
                -0.05,  # 左侧偏移（相对于坐标轴的百分比）
                label_y_max,
                format_y_axis_value(new_max),  # 使用智能格式化函数
                ha='right', va='center',
                fontsize=5,
                color='black',  # 黑色显示手动标签
                transform=ax.get_yaxis_transform()
            )

        elif track_type == 'gtf':
            # GTF tracks：不显示Y轴数值（它们显示基因结构，不是数值信号）
            ax.spines['left'].set_visible(False)  # 隐藏左侧边框
            ax.set_yticks([])  # 隐藏Y轴刻度

            # 获取当前ylim（由_plot_gtf设置）- 不要覆盖！
            # GTF tracks使用_plot_gtf中设置的ylim，支持多行循环显示
            current_ylim = ax.get_ylim()
            # 如果ylim是默认的小范围（-0.3到1.5），则扩展以适应多行显示
            if current_ylim[1] <= 2.0:
                # _plot_gtf已经设置了ylim，保持不变
                pass
            
            # 添加X轴网格
            ax.grid(True, axis='x', alpha=0.2)

        elif track_type == 'bed':
            # BED tracks：不显示Y轴数值（它们显示区间，不是数值信号）
            ax.spines['left'].set_visible(False)  # 隐藏左侧边框
            ax.set_yticks([])  # 隐藏Y轴刻度

            # 设置合适的Y轴范围（不显示数值，但要有合适的显示空间）
            ax.set_ylim(-0.3, 1.5)

            # 添加X轴网格
            ax.grid(True, axis='x', alpha=0.2)

        else:
            raise ValueError(f"Unsupported track type: {track_type}")

        # 4.4 设置X轴范围和隐藏刻度
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.set_xlabel('')

    # 5. 在最下方的track上添加坐标标签 - 与专用绘制函数保持一致
    # 注意：由于tracks列表被倒序处理，track_idx=0对应最下方的track
    bottom_ax = axes[0]
    ylim = bottom_ax.get_ylim()
    ymin, ymax = ylim

    add_coordinate_labels_for_tracks(
        ax=bottom_ax,
        xmin=0, xmax=1,
        ymin=ymin, ymax=ymax,
        start_pos=region.start,
        end_pos=region.end,
        chrom=region.chrom,
        fontsize=5,
        offset_ratio=0.2
    )

    # 6. 保存图像 - 与专用绘制函数保持一致
    if output:
        # 获取输出路径的基名
        from pathlib import Path
        output_path = Path(output)

        # 保存为PNG格式（高分辨率）
        png_path = output_path.with_suffix('.png')
        # 使用bbox_inches=None保持精确的尺寸控制
        plt.savefig(png_path, dpi=dpi, bbox_inches=None)
        print(f"✓ Saved PNG: {png_path}")

        # 保存为SVG格式（矢量格式）
        svg_path = output_path.with_suffix('.svg')
        plt.savefig(svg_path, bbox_inches=None)
        print(f"✓ Saved SVG: {svg_path}")

        plt.close()
    else:
        plt.show()


def format_y_axis_value(value: float) -> str:
    """
    格式化y轴标签显示 - 智能选择精度

    参数:
        value: y轴数值

    返回:
        str: 格式化后的字符串
    """
    # 获取数值的绝对值
    abs_value = abs(value)

    # 判断是否应该显示小数
    if abs_value >= 1000:
        # 大于等于1000，显示整数
        return f'{value:.0f}'
    elif abs_value >= 100:
        # 100-1000之间，保留一位小数
        return f'{value:.1f}'
    elif abs_value >= 10:
        # 10-100之间，保留两位小数
        return f'{value:.2f}'
    else:
        # 小于10，保留三位小数
        return f'{value:.3f}'


def format_genomic_coord(pos: int) -> str:
    """
    格式化基因组坐标显示

    参数:
        pos: 基因组位置

    返回:
        格式化后的字符串（如：1.0 Mb, 500 kb, 1000）
    """
    if pos >= 1_000_000:
        return f"{pos / 1_000_000:.1f} Mb"
    elif pos >= 1_000:
        return f"{pos / 1_000:.0f} kb"
    else:
        return str(pos)


def create_track(file_path: str, track_type: Optional[str] = None, **kwargs) -> SimpleTrack:
    """
    工厂函数 - 自动检测文件类型

    参数:
        file_path: 文件路径
        track_type: track类型 (可选，自动检测)
        **kwargs: 配置参数

    返回:
        SimpleTrack对象
    """
    # 1. 自动检测类型（如果未提供）
    if track_type is None:
        file_path_str = str(file_path)  # 确保是字符串
        ext = Path(file_path_str).suffix.lower()

        # 支持多种文件扩展名
        if ext in ['.bw', '.bigwig']:
            track_type = 'bigwig'
        elif ext in ['.gtf', '.gtf.gz']:
            track_type = 'gtf'
        elif ext in ['.bed', '.bed.gz']:
            track_type = 'bed'
        else:
            # 如果无法自动检测，尝试根据文件内容判断
            # 或者使用默认值
            raise ValueError(
                f"无法自动检测文件类型: {file_path}\n"
                f"请在字典中明确指定 'type': 'bigwig'/'gtf'/'bed'\n"
                f"支持的扩展名: .bw, .bigwig, .gtf, .gtf.gz, .bed, .bed.gz"
            )

    # 2. 创建track
    return SimpleTrack(file_path, track_type, **kwargs)


def _process_tracks_dict_list(
    tracks_dict_list: List[dict],
    **kwargs
) -> tuple[List[SimpleTrack], dict]:
    """
    处理tracks字典列表 - 新推荐的API格式

    将字典列表格式的tracks转换为SimpleTrack对象列表

    参数:
        tracks_dict_list: tracks字典列表，每个字典包含一个track的配置
            例如: [{'file': 'file1.bw', 'color': 'blue', 'min_value': 0}, ...]
        **kwargs: 全局参数（不按track索引分配）

    返回:
        tuple: (tracks列表, 全局参数字典)
    """
    # 1. TrackConfig支持的参数
    supported_params = {
        'file', 'type',  # 特殊处理：file和type是必需的
        'color', 'alpha', 'number_of_bins', 'summary_method',
        'plot_type', 'min_value', 'max_value', 'labels', 'fontsize', 'name',
        # GTF-specific parameters
        'gtf_style', 'border_color', 'color_utr', 'color_backbone',
        'color_arrow', 'line_width', 'height_utr', 'height_intron',
        'arrowhead_fraction', 'arrowhead_included', 'arrow_interval',
        'prefered_name'
    }

    # 2. 创建tracks列表
    tracks = []
    for track_dict in tracks_dict_list:
        # 2.1 提取必需参数
        if 'file' not in track_dict:
            raise ValueError("每个track字典必须包含'file'参数")

        file_path = track_dict['file']

        # 2.2 提取可选参数
        track_type = track_dict.get('type', None)
        track_kwargs = {}

        # 处理字典中的所有参数
        for param_name, param_value in track_dict.items():
            # 跳过已处理的特殊参数
            if param_name in ['file', 'type']:
                continue

            # 验证参数是否支持
            if param_name in supported_params:
                track_kwargs[param_name] = param_value
            else:
                # 未知参数，忽略或警告
                # 为了向后兼容，我们忽略未知参数
                pass

        # 2.3 创建track
        track = create_track(file_path, track_type, **track_kwargs)
        tracks.append(track)

    # 3. 所有kwargs都是全局参数（字典格式不需要全局参数）
    global_params = kwargs

    return tracks, global_params


def _build_tracks_params_matrix(
    file_paths: Optional[List[str]] = None,
    tracks: Optional[List[dict]] = None,
    track_type: Optional[List[str]] = None,
    colors: Optional[List[str]] = None,
    names: Optional[List[str]] = None,
    track_heights: Optional[List[float]] = None,
    **kwargs
) -> tuple[List[SimpleTrack], dict]:
    """
    构建tracks参数矩阵 - 系统性的参数处理函数

    支持两种输入格式：
    1. 新推荐：tracks字典列表格式
    2. 传统：file_paths + 平行列表参数格式

    参数:
        file_paths: 文件路径列表（传统格式）
        tracks: tracks字典列表（新推荐格式）
        track_type: track类型列表
        colors: 颜色列表
        names: 名称列表
        track_heights: 高度列表
        **kwargs: 其他track参数 (如min_value, max_value, fontsize等)

    返回:
        tuple: (tracks列表, 全局参数字典)
    """
    # 检查输入格式
    if tracks is not None:
        # 新推荐格式：tracks字典列表
        return _process_tracks_dict_list(tracks, **kwargs)
    elif file_paths is not None:
        # 传统格式：file_paths + 平行列表参数
        return _process_traditional_params(
            file_paths=file_paths,
            track_type=track_type,
            colors=colors,
            names=names,
            track_heights=track_heights,
            **kwargs
        )
    else:
        raise ValueError("必须提供file_paths参数（传统格式）或tracks参数（新推荐格式）")


def _process_traditional_params(
    file_paths: List[str],
    track_type: Optional[List[str]],
    colors: Optional[List[str]],
    names: Optional[List[str]],
    track_heights: Optional[List[float]],
    **kwargs
) -> tuple[List[SimpleTrack], dict]:
    """
    处理传统格式的参数 - file_paths + 平行列表

    参数:
        file_paths: 文件路径列表
        track_type: track类型列表
        colors: 颜色列表
        names: 名称列表
        track_heights: 高度列表
        **kwargs: 其他track参数

    返回:
        tuple: (tracks列表, 全局参数字典)
    """
    # 1. TrackConfig支持的参数
    supported_params = {
        'color', 'alpha', 'number_of_bins', 'summary_method',
        'plot_type', 'min_value', 'max_value', 'labels', 'fontsize', 'name',
        # GTF-specific parameters
        'gtf_style', 'border_color', 'color_utr', 'color_backbone',
        'color_arrow', 'line_width', 'height_utr', 'height_intron',
        'arrowhead_fraction', 'arrowhead_included', 'arrow_interval',
        'prefered_name'
    }

    # 2. 提取所有需要按track索引分配的参数
    track_params = {}
    for param_name in supported_params:
        if param_name in kwargs:
            param_value = kwargs[param_name]
            # 处理列表形式的参数（如min_value=[0, 0]）
            if isinstance(param_value, list):
                track_params[param_name] = param_value
            else:
                # 单值参数（如所有track共用）
                track_params[param_name] = param_value

    # 3. 创建tracks列表
    tracks = []
    for i, file_path in enumerate(file_paths):
        print(file_path)

        # 3.1 确定当前track的类型
        current_track_type = None
        if track_type is not None:
            if isinstance(track_type, list):
                if i < len(track_type):
                    current_track_type = track_type[i]
            else:
                current_track_type = track_type

        # 3.2 提取当前track的参数
        track_kwargs = {}

        # 处理track_type
        if current_track_type:
            # 类型已通过参数传递，不需要存储到config中
            pass

        # 处理颜色
        if colors and i < len(colors):
            track_kwargs['color'] = colors[i]

        # 处理名称
        if names and i < len(names):
            name_value = names[i]
            if name_value:  # 非空字符串才设置name
                track_kwargs['name'] = name_value
            # 空字符串或不提供name -> 不设置name属性（保持为None）

        # 处理其他track参数
        for param_name, param_value in track_params.items():
            if param_name == 'name':
                continue  # name已经处理过了

            if isinstance(param_value, list):
                # 列表形式的参数，按索引分配
                if i < len(param_value):
                    track_kwargs[param_name] = param_value[i]
            else:
                # 单值参数，所有track共用
                track_kwargs[param_name] = param_value

        # 3.3 创建track
        track = create_track(file_path, current_track_type, **track_kwargs)
        tracks.append(track)

    # 4. 分离track参数和全局参数
    global_params = {}
    for param_name, param_value in kwargs.items():
        if param_name not in supported_params:
            global_params[param_name] = param_value

    return tracks, global_params


def add_coordinate_labels_for_tracks(ax, xmin, xmax, ymin, ymax, start_pos, end_pos, chrom, fontsize=5, offset_ratio=0.02):
    """
    为tracks添加坐标标签 - 借鉴add_coordinate_labels的设计

    参数:
        ax: matplotlib的axes对象
        xmin, xmax, ymin, ymax: 坐标轴范围
        start_pos, end_pos: 起始和结束位置（bp）
        chrom: 染色体号
        fontsize: 字体大小，默认5
        offset_ratio: 偏移比例，默认0.2
    """
    # 使用axes相对坐标，确保位置稳定
    # 将标签放在axes内部，靠近底部边缘的位置
    offset_axes = 0.05  # axes高度的5%，让标签靠近底部边缘但仍在内部

    # 将坐标标签添加到axes内部底部
    # 使用transform=ax.get_xaxis_transform()将标签放在axes底部
    ax.text(
        xmin, -offset_axes,
        format_genomic_coord(start_pos),
        va='top', ha='left', fontsize=fontsize,
        transform=ax.get_xaxis_transform()
    )
    ax.text(
        xmax, -offset_axes,
        format_genomic_coord(end_pos),
        va='top', ha='right', fontsize=fontsize,
        transform=ax.get_xaxis_transform()
    )

    # 添加染色体标签（在坐标标签下方）
    ax.text(
        (xmin + xmax) / 2, -offset_axes,
        chrom,
        va='top', ha='center', fontsize=fontsize,
        transform=ax.get_xaxis_transform()
    )


def quick_plot(
    file_paths: Optional[List[str]] = None,
    tracks: Optional[List[dict]] = None,
    region: Optional[GenomeRange] = None,
    output: Optional[str] = None,
    track_type: Optional[List[str]] = None,
    colors: Optional[List[str]] = None,
    names: Optional[List[str]] = None,
    width: Optional[float] = None,
    left_margin: float = 0.8,
    right_margin: float = 1.6,
    track_heights: Optional[List[float]] = None,
    dpi: int = 300,
    **kwargs
):
    """
    快速绘图 - 最简单的API

    系统性重构：所有tracks（包括单个）都使用plot_mixed_tracks处理，
    实现真正的统一架构和扩展性。

    支持两种输入格式：
    1. 新推荐：tracks字典列表格式（每个track的配置独立组织）
       ```python
       quick_plot(
           tracks=[
               {'file': 'genes.gtf', 'color': 'blue', 'name': 'Gene'},
               {'file': 'ctcf.bw', 'color': 'orange', 'min_value': 0, 'max_value': 100}
           ],
           region=region,
           output='output.png'
       )
       ```

    2. 传统：file_paths + 平行列表参数格式（保持向后兼容）
       ```python
       quick_plot(
           file_paths=['genes.gtf', 'ctcf.bw'],
           colors=['blue', 'orange'],
           min_value=[0],
           max_value=[100],
           region=region,
           output='output.png'
       )
       ```

    参数:
        file_paths: 文件路径列表（传统格式，向后兼容）
        tracks: tracks字典列表（新推荐格式，每个字典包含一个track的完整配置）
            例如: [{'file': 'file1.bw', 'color': 'blue', 'min_value': 0}, ...]
        region: 基因组区域（必需）
        output: 输出文件路径
        track_type: track类型列表 (可选，自动检测)
        colors: 颜色列表（传统格式）
        names: 自定义track名称列表（可选，显示在tracks右侧）
        width: 主绘图区域宽度（cm）
        left_margin: 左侧留白宽度（cm）
        right_margin: 右侧留白宽度（cm）
        track_heights: 每个track的高度列表（cm）
        dpi: 图像分辨率（默认300 DPI）
        **kwargs: 其他配置参数（包括GTF专用参数等）

    注意:
        - 必须提供region参数
        - 必须提供file_paths（传统格式）或tracks（新推荐格式）之一
        - tracks字典格式中，每个字典支持的所有参数：
          * 'file': 文件路径（必需）
          * 'type': track类型 'bigwig'/'gtf'/'bed'（可选，自动检测）
          * 'color': 颜色（可选，默认'blue'）
          * 'alpha': 透明度（可选，默认0.8）
          * 'min_value': 最小值（可选）
          * 'max_value': 最大值（可选）
          * 'name': 自定义track名称（可选）
          * 'fontsize': 字体大小（可选，默认5）
          * 以及其他TrackConfig支持的参数
    """
    # 1. 验证必需参数
    if region is None:
        raise ValueError("必须提供region参数")

    # 2. 构建tracks参数矩阵 - 系统性参数处理
    tracks, global_params = _build_tracks_params_matrix(
        file_paths=file_paths,
        tracks=tracks,
        track_type=track_type,
        colors=colors,
        names=names,
        track_heights=track_heights,
        **kwargs
    )

    # 3. 检查tracks数量和类型
    track_types = [track.track_type for track in tracks]
    print(f"Tracks detected: {track_types}")

    # 4. 所有tracks（包括单个）都使用plot_mixed_tracks
    # 这实现了真正的统一架构和扩展性
    print("Using unified mixed tracks plotting function...")

    # 5. 传递参数给plot_mixed_tracks
    # 分离全局参数和track专用参数
    plot_mixed_tracks(
        tracks=tracks,
        region=region,
        output=output,
        width=width,
        left_margin=left_margin,
        right_margin=right_margin,
        track_heights=track_heights,
        dpi=dpi,
        **global_params  # 全局参数（如GTF专用参数等）
    )


# 向后兼容别名
plot_tracks = plot_bw_tracks  # 保留旧的plot_tracks名称作为plot_bw_tracks的别名

