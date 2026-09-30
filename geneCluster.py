#!/usr/bin/env python3
import os
import sys
import re
import json
import gzip
import subprocess
import threading
import traceback
from datetime import datetime
from collections import defaultdict, Counter

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

import pandas as pd
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    HAS_MPL = True
except Exception:
    HAS_MPL = False


CLUSTER_WINDOW_KB = 500
MAX_INTERVENING_GENES = 10
BGC_WINDOW_KB = 200
BGC_MIN_GENES = 3
BGC_MIN_CORE = 2
BGC_CORE_PLUS_ACCESSORY = 4
TANDEM_WINDOW_KB = 50

RESISTANCE_CLASSES = {
    'NLR', 'TNL', 'TNJ', 'TNJ_like', 'CC_NLR', 'NLR_Jacalin',
    'TIR_NLR', 'RPW8', 'Jacalin', 'LRR_repeat', 'NB_ARC', 'TIR', 'CC'
}

RESISTANCE_PFAMS = {
    'PF00931': 'NB_ARC',
    'PF01582': 'TIR',
    'PF13676': 'TIR',
    'PF00560': 'LRR',
    'PF07723': 'LRR',
    'PF13855': 'LRR',
    'PF01419': 'Jacalin',
    'PF05659': 'RPW8',
    'PF05660': 'RPW8',
}

BIOSYNTHESIS_CORE_PFAMS = {
    'PF01397': 'Terpene_synthase',
    'PF03936': 'Terpene_synthase_N',
    'PF13243': 'Terpene_synthase_C',
    'PF13249': 'Terpene_synthase_C2',
    'PF08491': 'Squalene_synthase',
    'PF00494': 'Prenyltransferase',
    'PF00195': 'Chalcone_synthase',
    'PF02797': 'Flavonoid_3_hydroxylase',
    'PF01596': 'SAM_methyltransferase',
    'PF00891': 'O_methyltransferase',
    'PF00682': 'HMGL_like',
}

BIOSYNTHESIS_ACCESSORY_PFAMS = {
    'PF00067': 'Cytochrome_P450',
    'PF00201': 'UDP_glucosyltransferase',
    'PF00155': 'Aminotransferase',
    'PF00282': 'Pyridoxal_enzyme',
    'PF00501': 'AMP_binding',
    'PF00141': 'Peroxidase',
    'PF00933': 'Glycosyl_hydrolase',
    'PF01565': 'FAD_binding_4',
    'PF13460': 'NAD_binding',
    'PF13561': 'Enoyl_reductase',
}

BIOSYNTHESIS_PFAMS = {}
BIOSYNTHESIS_PFAMS.update(BIOSYNTHESIS_CORE_PFAMS)
BIOSYNTHESIS_PFAMS.update(BIOSYNTHESIS_ACCESSORY_PFAMS)


BG_MAIN = "#dbeafe"
BG_CARD = "#eff6ff"
BG_INPUT = "#f8fbff"
FG_MAIN = "#1e3a5f"
FG_MUTED = "#5b7ea8"
ACCENT = "#2563eb"
ACCENT_ACTIVE = "#1d4ed8"
ACCENT_DARK = "#1e40af"
BORDER = "#93c5fd"
LOG_BG = "#f0f7ff"
LOG_FG = "#1e3a5f"


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def timestamp():
    return datetime.now().strftime('%Y%m%d_%H%M%S')


def open_maybe_gzip(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def normalize_gene_id(gid):
    gid = str(gid).strip()
    gid = re.sub(r'\.t\d+$', '', gid)
    gid = re.sub(r'\.mrna\d+$', '', gid)
    gid = re.sub(r'\.\d+$', '', gid)
    gid = re.sub(r'^gene[-_]', '', gid)
    gid = re.sub(r'^transcript[-_]', '', gid)
    gid = re.sub(r'^mrna[-_]', '', gid)
    return gid


def parse_gff3(gff_file):
    genes = {}
    gene_order = defaultdict(list)
    gene_re = re.compile(r'(?:^|;)ID=([^;]+)')
    with open_maybe_gzip(gff_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = fields
            if ftype.lower() != 'gene':
                continue
            m = gene_re.search(attrs)
            if not m:
                continue
            gene_id = m.group(1).strip()
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            genes[gene_id] = {
                'seqid': seqid,
                'start': start_i,
                'end': end_i,
                'strand': strand,
                'length': end_i - start_i + 1,
                'attributes': attrs,
            }
            gene_order[seqid].append(gene_id)
    for seqid in gene_order:
        gene_order[seqid].sort(key=lambda g: genes[g]['start'])
    return genes, gene_order


def parse_transcript_to_gene(gff_file):
    tx_to_gene = {}
    parent_re = re.compile(r'(?:^|;)Parent=([^;]+)')
    id_re = re.compile(r'(?:^|;)ID=([^;]+)')
    with open_maybe_gzip(gff_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            ftype = fields[2].lower()
            if ftype not in ('mrna', 'transcript'):
                continue
            attrs = fields[8]
            m_id = id_re.search(attrs)
            m_par = parent_re.search(attrs)
            if m_id and m_par:
                tx_id = m_id.group(1).strip()
                gene_id = m_par.group(1).strip()
                tx_to_gene[tx_id] = gene_id
                base = normalize_gene_id(tx_id)
                if base not in tx_to_gene:
                    tx_to_gene[base] = gene_id
    return tx_to_gene


def load_annotation_table(path):
    if not path or not os.path.exists(path):
        return pd.DataFrame()
    if path.endswith('.xlsx'):
        return pd.read_excel(path)
    if path.endswith('.tsv') or path.endswith('.txt'):
        return pd.read_csv(path, sep='\t')
    return pd.read_csv(path, sep=None, engine='python')


def build_id_map(df, tx_to_gene, id_col='Gene_ID'):
    exact = {}
    normalized = {}
    if df is None or df.empty or id_col not in df.columns:
        return exact, normalized
    for _, row in df.iterrows():
        gid = str(row[id_col]).strip()
        exact[gid] = row
        mapped = tx_to_gene.get(gid)
        if mapped is None:
            mapped = normalize_gene_id(gid)
        if mapped not in normalized:
            normalized[mapped] = row
        base = normalize_gene_id(gid)
        if base not in normalized:
            normalized[base] = row
    return exact, normalized


def classify_gene(gene_id, res_exact, res_norm, pharma_exact, pharma_norm):
    gene_classes = set()
    gene_pfams = set()
    source = []
    norm = normalize_gene_id(gene_id)

    r_row = res_exact.get(gene_id)
    if r_row is None:
        r_row = res_norm.get(norm)
    if r_row is not None:
        if 'Resistance_Classes' in r_row.index and pd.notna(r_row['Resistance_Classes']):
            for item in str(r_row['Resistance_Classes']).split(';'):
                item = item.strip()
                if item:
                    gene_classes.add(item)
        if 'Resistance_PFAMs' in r_row.index and pd.notna(r_row['Resistance_PFAMs']):
            for item in str(r_row['Resistance_PFAMs']).split(';'):
                item = item.strip()
                if item:
                    gene_pfams.add(item)
        if 'Resistance_Score' in r_row.index and pd.notna(r_row['Resistance_Score']):
            try:
                if float(r_row['Resistance_Score']) > 0 and not gene_classes:
                    gene_classes.add('NLR')
            except Exception:
                pass
        source.append('resistance')

    p_row = pharma_exact.get(gene_id)
    if p_row is None:
        p_row = pharma_norm.get(norm)
    if p_row is not None:
        if 'PFAM_Domains' in p_row.index and pd.notna(p_row['PFAM_Domains']):
            for item in str(p_row['PFAM_Domains']).split(';'):
                item = item.strip()
                if item:
                    gene_pfams.add(item)
        source.append('pharma')

    for pfam in gene_pfams:
        if pfam in RESISTANCE_PFAMS:
            gene_classes.add(RESISTANCE_PFAMS[pfam])

    return {
        'gene_id': gene_id,
        'classes': gene_classes,
        'pfams': gene_pfams,
        'source': source,
    }


def build_gene_table(genes, gene_order, resistance_df, pharma_df, tx_to_gene):
    res_exact, res_norm = build_id_map(resistance_df, tx_to_gene)
    ph_exact, ph_norm = build_id_map(pharma_df, tx_to_gene)
    print(f"[INFO] Resistance exact keys: {len(res_exact)}, normalized keys: {len(res_norm)}")
    print(f"[INFO] Pharma exact keys: {len(ph_exact)}, normalized keys: {len(ph_norm)}")
    print(f"[INFO] Sample RES normalized: {list(res_norm.keys())[:5]}")
    print(f"[INFO] Sample GFF IDs: {list(genes.keys())[:5]}")
    intersect = set(res_norm.keys()) & set(genes.keys())
    print(f"[INFO] Normalized intersection: {len(intersect)}")

    rows = []
    for seqid, gids in gene_order.items():
        for idx, gid in enumerate(gids):
            info = genes[gid]
            cls = classify_gene(gid, res_exact, res_norm, ph_exact, ph_norm)
            has_core = any(p in BIOSYNTHESIS_CORE_PFAMS for p in cls['pfams'])
            has_acc = any(p in BIOSYNTHESIS_ACCESSORY_PFAMS for p in cls['pfams'])
            rows.append({
                'gene_id': gid,
                'seqid': seqid,
                'start': info['start'],
                'end': info['end'],
                'strand': info['strand'],
                'length': info['length'],
                'index': idx,
                'classes': ';'.join(sorted(cls['classes'])),
                'pfams': ';'.join(sorted(cls['pfams'])),
                'sources': ';'.join(cls['source']),
                'is_resistance': int(bool(cls['classes'] & RESISTANCE_CLASSES)),
                'is_biosynthesis_core': int(has_core),
                'is_biosynthesis_accessory': int(has_acc),
                'is_biosynthesis': int(has_core or has_acc),
            })
    df = pd.DataFrame(rows)
    print(f"[INFO] Total genes: {len(df)}")
    print(f"[INFO] Resistance flagged: {int(df['is_resistance'].sum())}")
    print(f"[INFO] Biosynthesis core flagged: {int(df['is_biosynthesis_core'].sum())}")
    print(f"[INFO] Biosynthesis accessory flagged: {int(df['is_biosynthesis_accessory'].sum())}")
    print(f"[INFO] Biosynthesis total flagged: {int(df['is_biosynthesis'].sum())}")
    return df


def find_clusters(df, class_flag, window_kb, max_intervening, min_genes=2):
    clusters = []
    for seqid, sub in df.groupby('seqid'):
        sub = sub.sort_values('start').reset_index(drop=True)
        flagged = sub[sub[class_flag] == 1].reset_index(drop=True)
        if len(flagged) < min_genes:
            continue
        current = [flagged.iloc[0]]
        for i in range(1, len(flagged)):
            prev = current[-1]
            cur = flagged.iloc[i]
            dist_kb = (cur['start'] - prev['end']) / 1000.0
            intervening = sub[(sub['start'] > prev['end']) & (sub['end'] < cur['start'])]
            n_intervening = len(intervening)
            if dist_kb <= window_kb and n_intervening <= max_intervening:
                current.append(cur)
            else:
                if len(current) >= min_genes:
                    clusters.append(current)
                current = [cur]
        if len(current) >= min_genes:
            clusters.append(current)
    return clusters


def summarize_clusters(clusters, label):
    rows = []
    for cid, cluster in enumerate(clusters, 1):
        seqid = cluster[0]['seqid']
        start = min(g['start'] for g in cluster)
        end = max(g['end'] for g in cluster)
        span_kb = (end - start) / 1000.0
        genes_in = [g['gene_id'] for g in cluster]
        classes = set()
        for g in cluster:
            for c in str(g['classes']).split(';'):
                c = c.strip()
                if c:
                    classes.add(c)
        rows.append({
            'cluster_id': f"{label}_{cid}",
            'seqid': seqid,
            'start': start,
            'end': end,
            'span_kb': round(span_kb, 2),
            'n_genes': len(cluster),
            'gene_ids': ';'.join(genes_in),
            'classes': ';'.join(sorted(classes)),
        })
    return pd.DataFrame(rows)


def find_tandem_duplicates(df, class_flag, window_kb=50):
    rows = []
    for seqid, sub in df.groupby('seqid'):
        sub = sub.sort_values('start').reset_index(drop=True)
        flagged = sub[sub[class_flag] == 1].reset_index(drop=True)
        for i in range(len(flagged) - 1):
            a = flagged.iloc[i]
            b = flagged.iloc[i + 1]
            dist_kb = (b['start'] - a['end']) / 1000.0
            if 0 <= dist_kb <= window_kb:
                rows.append({
                    'seqid': seqid,
                    'gene1': a['gene_id'],
                    'gene2': b['gene_id'],
                    'distance_kb': round(dist_kb, 3),
                })
    return pd.DataFrame(rows)


def find_bgc_candidates(df, window_kb=200, min_core=2, core_plus_acc=4):
    clusters = []
    for seqid, sub in df.groupby('seqid'):
        sub = sub.sort_values('start').reset_index(drop=True)
        bio = sub[sub['is_biosynthesis'] == 1].reset_index(drop=True)
        if len(bio) < min_core:
            continue
        current = [bio.iloc[0]]
        for i in range(1, len(bio)):
            prev = current[-1]
            cur = bio.iloc[i]
            dist_kb = (cur['start'] - prev['end']) / 1000.0
            if dist_kb <= window_kb:
                current.append(cur)
            else:
                if len(current) >= min_core:
                    clusters.append(current)
                current = [cur]
        if len(current) >= min_core:
            clusters.append(current)

    validated = []
    for cluster in clusters:
        n_core = sum(1 for g in cluster if g['is_biosynthesis_core'] == 1)
        n_accessory = sum(1 for g in cluster if g['is_biosynthesis_accessory'] == 1)
        if n_core >= min_core:
            validated.append(cluster)
        elif n_core >= 1 and (n_core + n_accessory) >= core_plus_acc:
            validated.append(cluster)
    return validated


def summarize_bgc(clusters, label='BGC'):
    rows = []
    for cid, cluster in enumerate(clusters, 1):
        seqid = cluster[0]['seqid']
        start = min(g['start'] for g in cluster)
        end = max(g['end'] for g in cluster)
        span_kb = (end - start) / 1000.0
        pfams = set()
        core_domains = set()
        accessory_domains = set()
        for g in cluster:
            for p in str(g['pfams']).split(';'):
                p = p.strip()
                if p:
                    pfams.add(p)
        for p in pfams:
            if p in BIOSYNTHESIS_CORE_PFAMS:
                core_domains.add(BIOSYNTHESIS_CORE_PFAMS[p])
            elif p in BIOSYNTHESIS_ACCESSORY_PFAMS:
                accessory_domains.add(BIOSYNTHESIS_ACCESSORY_PFAMS[p])
        n_core = sum(1 for g in cluster if g['is_biosynthesis_core'] == 1)
        n_acc = sum(1 for g in cluster if g['is_biosynthesis_accessory'] == 1)
        rows.append({
            'cluster_id': f"{label}_{cid}",
            'seqid': seqid,
            'start': start,
            'end': end,
            'span_kb': round(span_kb, 2),
            'n_genes': len(cluster),
            'n_core': n_core,
            'n_accessory': n_acc,
            'gene_ids': ';'.join(g['gene_id'] for g in cluster),
            'pfams': ';'.join(sorted(pfams)),
            'core_domains': ';'.join(sorted(core_domains)),
            'accessory_domains': ';'.join(sorted(accessory_domains)),
        })
    return pd.DataFrame(rows)


def plot_cluster_map(df, resistance_clusters, bgc_clusters, output_prefix):
    if not HAS_MPL:
        return
    seqids = sorted(df['seqid'].unique())
    n = len(seqids)
    fig, ax = plt.subplots(figsize=(14, max(4, n * 0.6)))
    y_positions = {s: i for i, s in enumerate(seqids)}
    for s in seqids:
        sub = df[df['seqid'] == s]
        if sub.empty:
            continue
        ax.plot([sub['start'].min(), sub['end'].max()], [y_positions[s], y_positions[s]],
                color='#94a3b8', linewidth=2, zorder=1)
    for _, c in resistance_clusters.iterrows():
        y = y_positions.get(c['seqid'])
        if y is None:
            continue
        ax.plot([c['start'], c['end']], [y, y], color='#dc2626', linewidth=6, alpha=0.85, zorder=3)
    for _, c in bgc_clusters.iterrows():
        y = y_positions.get(c['seqid'])
        if y is None:
            continue
        ax.plot([c['start'], c['end']], [y, y], color='#2563eb', linewidth=6, alpha=0.85, zorder=2)
    ax.set_yticks(range(n))
    ax.set_yticklabels(seqids, fontsize=9)
    ax.set_xlabel('Position (bp)')
    ax.set_title('Genomic Distribution of Resistance and Biosynthesis Clusters')
    handles = [
        Patch(facecolor='#dc2626', label='Resistance cluster'),
        Patch(facecolor='#2563eb', label='Biosynthesis cluster'),
        Patch(facecolor='#94a3b8', label='Scaffold'),
    ]
    ax.legend(handles=handles, loc='upper right')
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    for ext in ('pdf', 'png', 'svg'):
        fig.savefig(f"{output_prefix}_cluster_map.{ext}", bbox_inches='tight')
    plt.close(fig)


def plot_cluster_sizes(resistance_clusters, bgc_clusters, output_prefix):
    if not HAS_MPL:
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    if not resistance_clusters.empty:
        axes[0].hist(resistance_clusters['n_genes'],
                     bins=range(2, resistance_clusters['n_genes'].max() + 3),
                     color='#dc2626', alpha=0.8, edgecolor='black')
        axes[0].set_title('Resistance Cluster Sizes')
        axes[0].set_xlabel('Number of genes')
        axes[0].set_ylabel('Count')
    else:
        axes[0].text(0.5, 0.5, 'No resistance clusters', ha='center', va='center')
        axes[0].set_title('Resistance Cluster Sizes')
    if not bgc_clusters.empty:
        axes[1].hist(bgc_clusters['n_genes'],
                     bins=range(3, bgc_clusters['n_genes'].max() + 3),
                     color='#2563eb', alpha=0.8, edgecolor='black')
        axes[1].set_title('Biosynthesis Cluster Sizes')
        axes[1].set_xlabel('Number of genes')
        axes[1].set_ylabel('Count')
    else:
        axes[1].text(0.5, 0.5, 'No biosynthesis clusters', ha='center', va='center')
        axes[1].set_title('Biosynthesis Cluster Sizes')
    plt.tight_layout()
    for ext in ('pdf', 'png', 'svg'):
        fig.savefig(f"{output_prefix}_cluster_sizes.{ext}", bbox_inches='tight')
    plt.close(fig)


def write_summary(output_dir, stats, resistance_clusters, bgc_clusters,
                  tandem_df, df):
    path = os.path.join(output_dir, 'cluster_summary.txt')
    with open(path, 'w') as f:
        f.write("=" * 78 + "\n")
        f.write("GFF3 CLUSTER ANALYSIS SUMMARY\n")
        f.write("=" * 78 + "\n")
        f.write(f"Generated: {datetime.now().isoformat()}\n\n")
        f.write(f"Total genes parsed: {len(df)}\n")
        f.write(f"Scaffolds: {df['seqid'].nunique()}\n")
        f.write(f"Resistance genes: {int(df['is_resistance'].sum())}\n")
        f.write(f"Biosynthesis core genes: {int(df['is_biosynthesis_core'].sum())}\n")
        f.write(f"Biosynthesis accessory genes: {int(df['is_biosynthesis_accessory'].sum())}\n")
        f.write(f"Biosynthesis total genes: {int(df['is_biosynthesis'].sum())}\n\n")
        f.write("-" * 78 + "\n")
        f.write("CLUSTER PARAMETERS\n")
        f.write("-" * 78 + "\n")
        for k, v in stats.items():
            f.write(f"{k}: {v}\n")
        f.write("\n")
        f.write("-" * 78 + "\n")
        f.write("RESISTANCE CLUSTERS\n")
        f.write("-" * 78 + "\n")
        if resistance_clusters.empty:
            f.write("No resistance clusters detected\n")
        else:
            f.write(f"Total clusters: {len(resistance_clusters)}\n")
            f.write(f"Total genes in clusters: {resistance_clusters['n_genes'].sum()}\n")
            f.write(f"Mean cluster size: {resistance_clusters['n_genes'].mean():.2f}\n")
            f.write(f"Largest cluster: {resistance_clusters['n_genes'].max()}\n\n")
            for _, row in resistance_clusters.sort_values('n_genes', ascending=False).iterrows():
                f.write(f"  {row['cluster_id']}\t{row['seqid']}:{row['start']}-{row['end']}\t"
                        f"n={row['n_genes']}\tspan={row['span_kb']}kb\n")
                f.write(f"    Genes: {row['gene_ids']}\n")
                f.write(f"    Classes: {row['classes']}\n")
        f.write("\n")
        f.write("-" * 78 + "\n")
        f.write("BIOSYNTHESIS CLUSTERS (CORE-VALIDATED)\n")
        f.write("-" * 78 + "\n")
        if bgc_clusters.empty:
            f.write("No biosynthesis clusters detected\n")
        else:
            f.write(f"Total clusters: {len(bgc_clusters)}\n")
            f.write(f"Total genes in clusters: {bgc_clusters['n_genes'].sum()}\n")
            f.write(f"Mean cluster size: {bgc_clusters['n_genes'].mean():.2f}\n")
            f.write(f"Total core genes: {bgc_clusters['n_core'].sum()}\n")
            f.write(f"Total accessory genes: {bgc_clusters['n_accessory'].sum()}\n\n")
            for _, row in bgc_clusters.sort_values('n_core', ascending=False).iterrows():
                f.write(f"  {row['cluster_id']}\t{row['seqid']}:{row['start']}-{row['end']}\t"
                        f"n={row['n_genes']}\tcore={row['n_core']}\tacc={row['n_accessory']}\t"
                        f"span={row['span_kb']}kb\n")
                f.write(f"    Genes: {row['gene_ids']}\n")
                f.write(f"    Core domains: {row['core_domains']}\n")
                f.write(f"    Accessory domains: {row['accessory_domains']}\n")
        f.write("\n")
        f.write("-" * 78 + "\n")
        f.write("TANDEM DUPLICATES (RESISTANCE)\n")
        f.write("-" * 78 + "\n")
        if tandem_df.empty:
            f.write("No tandem duplicates detected\n")
        else:
            f.write(f"Total pairs: {len(tandem_df)}\n")
            for _, row in tandem_df.iterrows():
                f.write(f"  {row['seqid']}\t{row['gene1']} <-> {row['gene2']}\t"
                        f"distance={row['distance_kb']}kb\n")
    return path


def run_cluster_analysis(gff_file, resistance_table, pharma_table, output_dir,
                         window_kb, max_intervening, bgc_window_kb, tandem_window_kb,
                         bgc_min_core=BGC_MIN_CORE,
                         bgc_core_plus_acc=BGC_CORE_PLUS_ACCESSORY):
    ensure_dir(output_dir)
    print(f"[INFO] Parsing GFF3: {gff_file}")
    genes, gene_order = parse_gff3(gff_file)
    print(f"[INFO] GFF genes: {len(genes)}")
    tx_to_gene = parse_transcript_to_gene(gff_file)
    print(f"[INFO] Transcript to gene mappings: {len(tx_to_gene)}")
    print(f"[INFO] Sample tx to gene: {list(tx_to_gene.items())[:5]}")

    resistance_df = load_annotation_table(resistance_table) if resistance_table else pd.DataFrame()
    pharma_df = load_annotation_table(pharma_table) if pharma_table else pd.DataFrame()
    print(f"[INFO] Resistance table rows: {len(resistance_df)}")
    print(f"[INFO] Pharma table rows: {len(pharma_df)}")

    df = build_gene_table(genes, gene_order, resistance_df, pharma_df, tx_to_gene)
    df.to_csv(os.path.join(output_dir, 'gene_table.csv'), index=False)

    resistance_clusters_raw = find_clusters(df, 'is_resistance', window_kb,
                                            max_intervening, min_genes=2)
    resistance_clusters = summarize_clusters(resistance_clusters_raw, 'RC')
    resistance_clusters.to_csv(os.path.join(output_dir, 'resistance_clusters.csv'), index=False)
    print(f"[INFO] Resistance clusters: {len(resistance_clusters)}")

    bgc_clusters_raw = find_bgc_candidates(df, window_kb=bgc_window_kb,
                                           min_core=bgc_min_core,
                                           core_plus_acc=bgc_core_plus_acc)
    bgc_clusters = summarize_bgc(bgc_clusters_raw, 'BGC')
    bgc_clusters.to_csv(os.path.join(output_dir, 'biosynthesis_clusters.csv'), index=False)
    print(f"[INFO] Biosynthesis clusters (core-validated): {len(bgc_clusters)}")

    tandem_df = find_tandem_duplicates(df, 'is_resistance', window_kb=tandem_window_kb)
    tandem_df.to_csv(os.path.join(output_dir, 'tandem_duplicates.csv'), index=False)
    print(f"[INFO] Tandem duplicate pairs: {len(tandem_df)}")

    stats = {
        'GFF file': gff_file,
        'Resistance table': resistance_table or 'N/A',
        'Pharma table': pharma_table or 'N/A',
        'Cluster window (kb)': window_kb,
        'Max intervening genes': max_intervening,
        'BGC window (kb)': bgc_window_kb,
        'BGC min core genes': bgc_min_core,
        'BGC core plus accessory threshold': bgc_core_plus_acc,
        'Tandem window (kb)': tandem_window_kb,
        'Total genes': len(df),
        'Scaffolds': df['seqid'].nunique(),
        'Resistance genes': int(df['is_resistance'].sum()),
        'Biosynthesis core genes': int(df['is_biosynthesis_core'].sum()),
        'Biosynthesis accessory genes': int(df['is_biosynthesis_accessory'].sum()),
        'Resistance clusters': len(resistance_clusters),
        'Biosynthesis clusters': len(bgc_clusters),
        'Tandem duplicate pairs': len(tandem_df),
    }
    summary_path = write_summary(output_dir, stats, resistance_clusters,
                                 bgc_clusters, tandem_df, df)
    prefix = os.path.join(output_dir, 'cluster')
    plot_cluster_map(df, resistance_clusters, bgc_clusters, prefix)
    plot_cluster_sizes(resistance_clusters, bgc_clusters, prefix)

    return {
        'summary': summary_path,
        'gene_table': os.path.join(output_dir, 'gene_table.csv'),
        'resistance_clusters': os.path.join(output_dir, 'resistance_clusters.csv'),
        'biosynthesis_clusters': os.path.join(output_dir, 'biosynthesis_clusters.csv'),
        'tandem_duplicates': os.path.join(output_dir, 'tandem_duplicates.csv'),
        'n_resistance_clusters': len(resistance_clusters),
        'n_bgc_clusters': len(bgc_clusters),
        'n_tandem': len(tandem_df),
    }


def validate_gff(gff_file):
    errors = []
    warnings = []
    gene_ids = set()
    gene_re = re.compile(r'(?:^|;)ID=([^;]+)')
    parent_re = re.compile(r'(?:^|;)Parent=([^;]+)')
    with open_maybe_gzip(gff_file) as f:
        for lineno, line in enumerate(f, 1):
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                errors.append(f"Line {lineno}: expected 9 fields, got {len(fields)}")
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = fields
            try:
                int(start)
                int(end)
            except ValueError:
                errors.append(f"Line {lineno}: invalid coordinates")
                continue
            if ftype.lower() == 'gene':
                m = gene_re.search(attrs)
                if not m:
                    warnings.append(f"Line {lineno}: gene without ID")
                else:
                    gid = m.group(1)
                    if gid in gene_ids:
                        warnings.append(f"Line {lineno}: duplicate gene ID {gid}")
                    gene_ids.add(gid)
            if ftype.lower() in ('mrna', 'cds', 'exon', 'transcript'):
                if not parent_re.search(attrs):
                    warnings.append(f"Line {lineno}: {ftype} without Parent")
    return errors, warnings, len(gene_ids)


class ClusterApp:
    def __init__(self, root):
        self.root = root
        self.root.title("GFF3 Cluster Analyzer v1.3")
        self.root.geometry("980x860")
        self.root.configure(bg=BG_MAIN)
        self._configure_style()
        self._build_ui()

    def _configure_style(self):
        style = ttk.Style()
        style.theme_use('clam')
        style.configure('.', background=BG_MAIN, foreground=FG_MAIN,
                        fieldbackground=BG_INPUT, bordercolor=BORDER)
        style.configure('TFrame', background=BG_MAIN)
        style.configure('TLabel', background=BG_MAIN, foreground=FG_MAIN)
        style.configure('Title.TLabel', font=('Segoe UI', 16, 'bold'),
                        foreground=ACCENT, background=BG_MAIN)
        style.configure('Sub.TLabel', font=('Segoe UI', 10),
                        foreground=FG_MUTED, background=BG_MAIN)
        style.configure('TButton', background=ACCENT, foreground='white',
                        font=('Segoe UI', 10, 'bold'), padding=6,
                        borderwidth=0, relief='flat')
        style.map('TButton',
                  background=[('active', ACCENT_ACTIVE), ('pressed', ACCENT_DARK)],
                  foreground=[('active', 'white')])
        style.configure('Run.TButton', background=ACCENT_DARK, foreground='white',
                        font=('Segoe UI', 12, 'bold'), padding=10,
                        borderwidth=0, relief='flat')
        style.map('Run.TButton',
                  background=[('active', ACCENT_ACTIVE), ('pressed', ACCENT_DARK)])
        style.configure('TNotebook', background=BG_MAIN, bordercolor=BORDER)
        style.configure('TNotebook.Tab', background=BG_CARD, foreground=FG_MUTED,
                        padding=[16, 8], font=('Segoe UI', 10, 'bold'))
        style.map('TNotebook.Tab',
                  background=[('selected', ACCENT), ('active', BG_INPUT)],
                  foreground=[('selected', 'white'), ('active', FG_MAIN)])
        style.configure('TEntry', fieldbackground=BG_INPUT, foreground=FG_MAIN,
                        insertcolor=FG_MAIN, bordercolor=BORDER)
        style.configure('TLabelframe', background=BG_MAIN, foreground=ACCENT,
                        bordercolor=BORDER)
        style.configure('TLabelframe.Label', background=BG_MAIN, foreground=ACCENT,
                        font=('Segoe UI', 11, 'bold'))
        style.configure('TScrollbar', background=BG_CARD, troughcolor=BG_INPUT,
                        bordercolor=BORDER, arrowcolor=ACCENT)

    def _build_ui(self):
        header = ttk.Frame(self.root)
        header.pack(fill=tk.X, padx=20, pady=(15, 5))
        ttk.Label(header, text="GFF3 Cluster Analyzer",
                  style='Title.TLabel').pack(anchor='w')
        ttk.Label(header, text="Preprocessing, resistance clusters and biosynthesis clusters",
                  style='Sub.TLabel').pack(anchor='w')

        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=15, pady=10)

        self._build_main_tab()
        self._build_validation_tab()
        self._build_about_tab()

    def _row(self, parent, row, label, attr, browse='file'):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky='e', padx=8, pady=6)
        entry = ttk.Entry(parent, width=58)
        entry.grid(row=row, column=1, padx=8, pady=6, sticky='we')
        setattr(self, attr, entry)
        cmd = (lambda: self._browse_dir(entry)) if browse == 'dir' else (lambda: self._browse_file(entry))
        ttk.Button(parent, text="Browse", command=cmd).grid(row=row, column=2, padx=8, pady=6)
        return entry

    def _browse_file(self, entry):
        p = filedialog.askopenfilename()
        if p:
            entry.delete(0, tk.END)
            entry.insert(0, p)

    def _browse_dir(self, entry):
        p = filedialog.askdirectory()
        if p:
            entry.delete(0, tk.END)
            entry.insert(0, p)

    def _build_main_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Cluster Analysis")

        form = ttk.LabelFrame(tab, text="Input Files")
        form.pack(fill=tk.X, padx=15, pady=10)

        self._row(form, 0, "GFF3 file:", 'gff_entry')
        self._row(form, 1, "Resistance genes table (.xlsx):", 'res_entry')
        self._row(form, 2, "Pharma genes table (.xlsx):", 'pharma_entry')

        params = ttk.LabelFrame(tab, text="Cluster Parameters")
        params.pack(fill=tk.X, padx=15, pady=10)

        ttk.Label(params, text="Resistance window (kb):").grid(row=0, column=0, sticky='e', padx=8, pady=6)
        self.window_entry = ttk.Entry(params, width=12)
        self.window_entry.insert(0, str(CLUSTER_WINDOW_KB))
        self.window_entry.grid(row=0, column=1, sticky='w', padx=8, pady=6)

        ttk.Label(params, text="Max intervening genes:").grid(row=0, column=2, sticky='e', padx=8, pady=6)
        self.interv_entry = ttk.Entry(params, width=12)
        self.interv_entry.insert(0, str(MAX_INTERVENING_GENES))
        self.interv_entry.grid(row=0, column=3, sticky='w', padx=8, pady=6)

        ttk.Label(params, text="BGC window (kb):").grid(row=1, column=0, sticky='e', padx=8, pady=6)
        self.bgc_entry = ttk.Entry(params, width=12)
        self.bgc_entry.insert(0, str(BGC_WINDOW_KB))
        self.bgc_entry.grid(row=1, column=1, sticky='w', padx=8, pady=6)

        ttk.Label(params, text="BGC min core genes:").grid(row=1, column=2, sticky='e', padx=8, pady=6)
        self.bgc_core_entry = ttk.Entry(params, width=12)
        self.bgc_core_entry.insert(0, str(BGC_MIN_CORE))
        self.bgc_core_entry.grid(row=1, column=3, sticky='w', padx=8, pady=6)

        ttk.Label(params, text="BGC core+accessory threshold:").grid(row=2, column=0, sticky='e', padx=8, pady=6)
        self.bgc_acc_entry = ttk.Entry(params, width=12)
        self.bgc_acc_entry.insert(0, str(BGC_CORE_PLUS_ACCESSORY))
        self.bgc_acc_entry.grid(row=2, column=1, sticky='w', padx=8, pady=6)

        ttk.Label(params, text="Tandem window (kb):").grid(row=2, column=2, sticky='e', padx=8, pady=6)
        self.tandem_entry = ttk.Entry(params, width=12)
        self.tandem_entry.insert(0, str(TANDEM_WINDOW_KB))
        self.tandem_entry.grid(row=2, column=3, sticky='w', padx=8, pady=6)

        out_frame = ttk.LabelFrame(tab, text="Output")
        out_frame.pack(fill=tk.X, padx=15, pady=10)
        ttk.Label(out_frame, text="Output folder:").grid(row=0, column=0, sticky='e', padx=8, pady=6)
        self.out_entry = ttk.Entry(out_frame, width=58)
        self.out_entry.insert(0, "cluster_analysis")
        self.out_entry.grid(row=0, column=1, padx=8, pady=6, sticky='we')
        ttk.Button(out_frame, text="Browse",
                   command=lambda: self._browse_dir(self.out_entry)).grid(
            row=0, column=2, padx=8, pady=6)

        run_frame = ttk.Frame(tab)
        run_frame.pack(fill=tk.X, padx=15, pady=10)
        self.run_button = ttk.Button(run_frame, text="Run Cluster Analysis",
                                     style='Run.TButton',
                                     command=self._run_analysis)
        self.run_button.pack(pady=5)

        log_frame = ttk.LabelFrame(tab, text="Log")
        log_frame.pack(fill=tk.BOTH, expand=True, padx=15, pady=10)
        self.log = scrolledtext.ScrolledText(log_frame, height=14, bg=LOG_BG,
                                             fg=LOG_FG, insertbackground=FG_MAIN,
                                             font=('Consolas', 9), relief='flat',
                                             borderwidth=1)
        self.log.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def _build_validation_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="GFF Validation")

        form = ttk.LabelFrame(tab, text="Validate GFF3")
        form.pack(fill=tk.X, padx=15, pady=10)
        self._row(form, 0, "GFF3 file:", 'val_entry')

        self.validate_button = ttk.Button(form, text="Validate", style='Run.TButton',
                                          command=self._run_validation)
        self.validate_button.grid(row=1, column=0, columnspan=3, pady=10)

        out_frame = ttk.LabelFrame(tab, text="Validation Output")
        out_frame.pack(fill=tk.BOTH, expand=True, padx=15, pady=10)
        self.val_log = scrolledtext.ScrolledText(out_frame, height=20, bg=LOG_BG,
                                                 fg=LOG_FG, insertbackground=FG_MAIN,
                                                 font=('Consolas', 9), relief='flat',
                                                 borderwidth=1)
        self.val_log.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def _build_about_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="About")

        text = (
            "GFF3 Cluster Analyzer v1.3\n\n"
            "Steps performed:\n"
            "  1. Parses the GFF3 annotation and extracts gene coordinates.\n"
            "  2. Builds transcript-to-gene mapping from mRNA features.\n"
            "  3. Loads resistance and pharmaceutical gene tables from annoMining.\n"
            "  4. Maps annoMining transcript IDs to GFF3 gene IDs.\n"
            "  5. Detects resistance gene clusters using a sliding window.\n"
            "  6. Detects biosynthesis clusters with core vs accessory validation.\n"
            "  7. Identifies tandem duplicate pairs among resistance genes.\n"
            "  8. Produces CSV tables and publication-ready figures.\n\n"
            "Biosynthesis validation logic:\n"
            "  A cluster is reported if it contains at least N core genes (default 2),\n"
            "  OR at least 1 core gene plus enough accessory genes to reach the\n"
            "  core+accessory threshold (default 4).\n\n"
            "Core PFAMs:\n"
            "  " + ", ".join(sorted(BIOSYNTHESIS_CORE_PFAMS.keys())) + "\n\n"
            "Accessory PFAMs:\n"
            "  " + ", ".join(sorted(BIOSYNTHESIS_ACCESSORY_PFAMS.keys())) + "\n\n"
            "Default parameters:\n"
            f"  Resistance cluster window: {CLUSTER_WINDOW_KB} kb\n"
            f"  Max intervening genes: {MAX_INTERVENING_GENES}\n"
            f"  BGC window: {BGC_WINDOW_KB} kb\n"
            f"  BGC min core: {BGC_MIN_CORE}\n"
            f"  BGC core+accessory: {BGC_CORE_PLUS_ACCESSORY}\n"
            f"  Tandem window: {TANDEM_WINDOW_KB} kb\n\n"
            "Requirements:\n"
            "  Python 3.8+, pandas, numpy, matplotlib\n"
        )
        label = ttk.Label(tab, text=text, justify='left', font=('Segoe UI', 10))
        label.pack(anchor='nw', padx=20, pady=20)

    def _log(self, widget, msg):
        widget.insert(tk.END, msg + "\n")
        widget.see(tk.END)

    def _log_safe(self, widget, msg):
        self.root.after(0, lambda: self._log(widget, msg))

    def _show_info(self, title, message):
        self.root.after(0, lambda: messagebox.showinfo(title, message))

    def _show_error(self, title, message):
        self.root.after(0, lambda: messagebox.showerror(title, message))

    def _set_button_state(self, button, state):
        self.root.after(0, lambda: button.configure(state=state))

    def _run_thread(self, target, log_widget, success_msg, trigger_button=None):
        if trigger_button is not None:
            self._set_button_state(trigger_button, 'disabled')

        def wrapper():
            try:
                self._log_safe(log_widget, "Starting analysis...")
                target()
                self._log_safe(log_widget, f"Completed: {success_msg}")
                self._show_info("Success", success_msg)
            except Exception as e:
                self._log_safe(log_widget, f"Error: {e}")
                self._log_safe(log_widget, traceback.format_exc())
                self._show_error("Error", str(e))
            finally:
                if trigger_button is not None:
                    self._set_button_state(trigger_button, 'normal')

        threading.Thread(target=wrapper, daemon=True).start()

    def _run_analysis(self):
        gff = self.gff_entry.get().strip()
        res = self.res_entry.get().strip()
        pharma = self.pharma_entry.get().strip()
        outdir = self.out_entry.get().strip() or 'cluster_analysis'
        if not gff or not os.path.exists(gff):
            messagebox.showerror("Error", "Provide a valid GFF3 file")
            return
        try:
            window_kb = float(self.window_entry.get().strip() or CLUSTER_WINDOW_KB)
            max_interv = int(self.interv_entry.get().strip() or MAX_INTERVENING_GENES)
            bgc_kb = float(self.bgc_entry.get().strip() or BGC_WINDOW_KB)
            bgc_core = int(self.bgc_core_entry.get().strip() or BGC_MIN_CORE)
            bgc_acc = int(self.bgc_acc_entry.get().strip() or BGC_CORE_PLUS_ACCESSORY)
            tandem_kb = float(self.tandem_entry.get().strip() or TANDEM_WINDOW_KB)
        except ValueError:
            messagebox.showerror("Error", "Parameters must be numeric")
            return

        def task():
            result = run_cluster_analysis(gff, res, pharma, outdir,
                                          window_kb, max_interv, bgc_kb, tandem_kb,
                                          bgc_min_core=bgc_core,
                                          bgc_core_plus_acc=bgc_acc)
            self._log_safe(self.log, f"Summary: {result['summary']}")
            self._log_safe(self.log, f"Gene table: {result['gene_table']}")
            self._log_safe(self.log, f"Resistance clusters: {result['n_resistance_clusters']}")
            self._log_safe(self.log, f"Biosynthesis clusters: {result['n_bgc_clusters']}")
            self._log_safe(self.log, f"Tandem duplicate pairs: {result['n_tandem']}")

        self._run_thread(task, self.log, f"Cluster analysis done: {outdir}",
                         trigger_button=self.run_button)

    def _run_validation(self):
        gff = self.val_entry.get().strip()
        if not gff or not os.path.exists(gff):
            messagebox.showerror("Error", "Provide a valid GFF3 file")
            return

        def task():
            errors, warnings, n_genes = validate_gff(gff)
            self._log_safe(self.val_log, f"Total genes: {n_genes}")
            self._log_safe(self.val_log, f"Errors: {len(errors)}")
            self._log_safe(self.val_log, f"Warnings: {len(warnings)}")
            for e in errors[:50]:
                self._log_safe(self.val_log, f"ERROR: {e}")
            for w in warnings[:50]:
                self._log_safe(self.val_log, f"WARN: {w}")

        self._run_thread(task, self.val_log, "Validation done",
                         trigger_button=self.validate_button)


def main():
    root = tk.Tk()
    ClusterApp(root)
    root.mainloop()


if __name__ == '__main__':
    main()
