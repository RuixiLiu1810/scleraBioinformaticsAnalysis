import os
import re
import time
from pathlib import Path
from typing import Optional, Tuple, List, Dict

import numpy as np
import pandas as pd
import scanpy as sc

# =========================
# 0) 你只需要改这些配置
# =========================

# (A) 你的数据根目录：以该目录下的“每个GSE一个子文件夹”为准
# 例：E:\20260120\raw  -> raw\GSE154447\... raw\GSE162733\...
DATASET_ROOT = Path(r"E:\20260120\raw")

# (B) 输出根目录
OUT_ROOT = Path(r"E:\20260120\scanpy_out")

# (C) 若你目前只有 family.soft.gz，而没有 matrix/object：
#     设为 True 会尝试根据 GSE 号从 GEO FTP 下载 supplementary 到 DATASET_ROOT/GSE*****/
#     需要联网，且仅下载 h5ad / 10x h5 / mtx三件套
ENABLE_FTP_DOWNLOAD_IF_MISSING = False

# (D) 处理参数（通用默认）
MIN_GENES = 200
MAX_PCT_MT = 20.0
N_TOP_GENES = 3000
N_PCS = 50
LEIDEN_RES = 0.8
MARKER_TOPN = 100

RANDOM_STATE = 0


# =========================
# 1) 文件识别与加载
# =========================

def is_gse_folder(p: Path) -> bool:
    return p.is_dir() and re.match(r"^GSE\d+$", p.name, flags=re.I) is not None

def find_h5ad(gse_dir: Path) -> Optional[Path]:
    files = list(gse_dir.glob("*.h5ad"))
    return files[0] if files else None

def find_10x_h5(gse_dir: Path) -> Optional[Path]:
    # 注意：.h5 不一定都是10x，但我们先尝试 read_10x_h5
    files = [p for p in gse_dir.glob("*.h5") if p.is_file()]
    return files[0] if files else None

def has_10x_mtx_trio(gse_dir: Path) -> bool:
    has_mtx = any((gse_dir / f).exists() for f in ["matrix.mtx", "matrix.mtx.gz"])
    has_bar = any((gse_dir / f).exists() for f in ["barcodes.tsv", "barcodes.tsv.gz"])
    has_feat = any((gse_dir / f).exists() for f in ["features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"])
    return has_mtx and has_bar and has_feat

def find_10x_mtx_dir(gse_dir: Path) -> Optional[Path]:
    # root
    if has_10x_mtx_trio(gse_dir):
        return gse_dir

    # common subfolders
    for sub in ["filtered_feature_bc_matrix", "raw_feature_bc_matrix", "filtered_gene_bc_matrices", "raw_gene_bc_matrices"]:
        subdir = gse_dir / sub
        if subdir.exists() and subdir.is_dir():
            # direct
            if has_10x_mtx_trio(subdir):
                return subdir
            # one-level deeper (e.g. mm10/)
            for deeper in subdir.glob("*"):
                if deeper.is_dir() and has_10x_mtx_trio(deeper):
                    return deeper
    return None

def load_anndata_from_folder(gse: str, gse_dir: Path) -> Tuple[sc.AnnData, str]:
    """
    Load priority: h5ad > 10x h5 > 10x mtx
    """
    h5ad = find_h5ad(gse_dir)
    if h5ad:
        adata = sc.read_h5ad(h5ad)
        return adata, f"h5ad:{h5ad.name}"

    h5 = find_10x_h5(gse_dir)
    if h5:
        try:
            adata = sc.read_10x_h5(h5)
            return adata, f"10x_h5:{h5.name}"
        except Exception:
            # 如果不是10x h5，就继续尝试 mtx
            pass

    mtx_dir = find_10x_mtx_dir(gse_dir)
    if mtx_dir:
        adata = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols", cache=False)
        return adata, f"10x_mtx:{mtx_dir}"

    raise FileNotFoundError(
        f"[{gse}] No readable matrix/object found in {gse_dir}. "
        f"Need .h5ad or 10x .h5 or (matrix.mtx + barcodes.tsv + features.tsv)."
    )


# =========================
# 2) 可选：FTP 自动下载 supplementary（仅当缺文件时）
# =========================

def ftp_download_supplementary_for_gse(gse: str, target_dir: Path) -> Dict[str, str]:
    """
    Very small downloader:
    - list ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE{xxx}nnn/{gse}/suppl/
    - pick best (h5ad > 10x h5 > mtx trio)
    - download into target_dir
    Return dict with keys: status, mode, downloaded, error
    """
    from ftplib import FTP, error_perm

    GEO_FTP_HOST = "ftp.ncbi.nlm.nih.gov"
    GEO_SERIES_ROOT = "/geo/series"

    def bucket(gse_):
        num = int(re.match(r"^GSE(\d+)$", gse_, re.I).group(1))
        return f"GSE{num//1000}nnn"

    def list_dir(ftp, path):
        try:
            ftp.cwd(path)
            return ftp.nlst()
        except error_perm:
            return []

    def pick(files: List[str]) -> Tuple[List[str], str]:
        # h5ad
        for f in files:
            if re.search(r"\.h5ad(\.gz)?$|\.h5ad\.zip$", f, re.I):
                return [f], "h5ad"
        # 10x h5
        for f in files:
            if re.search(r"filtered_feature_bc_matrix\.h5(\.gz)?$", f, re.I):
                return [f], "10x_h5"
        for f in files:
            if re.search(r"\.h5(\.gz)?$", f, re.I):
                return [f], "10x_h5"
        # mtx trio
        def find_one(pat):
            for f in files:
                if re.search(pat, f, re.I):
                    return f
            return None
        mtx = find_one(r"matrix\.mtx(\.gz)?$")
        bar = find_one(r"barcodes\.tsv(\.gz)?$")
        feat = find_one(r"(features|genes)\.tsv(\.gz)?$")
        if mtx and bar and feat:
            return [mtx, bar, feat], "10x_mtx_trio"
        return [], "none"

    def retr(ftp, remote_path, local_path):
        local_path.parent.mkdir(parents=True, exist_ok=True)
        with open(local_path, "wb") as f:
            ftp.retrbinary(f"RETR {remote_path}", f.write, blocksize=1024*1024)

    info = {"status": "OK", "mode": "", "downloaded": "", "error": ""}

    try:
        ftp = FTP(GEO_FTP_HOST, timeout=60)
        ftp.login()

        supdir = f"{GEO_SERIES_ROOT}/{bucket(gse)}/{gse}/suppl"
        files = list_dir(ftp, supdir)
        if not files:
            raise FileNotFoundError(f"No supplementary files at {supdir}")

        selected, mode = pick(files)
        if not selected:
            raise FileNotFoundError("Supplementary exists but no usable matrix/object detected (h5ad/h5/mtx trio).")

        # download
        target_dir.mkdir(parents=True, exist_ok=True)
        for fn in selected:
            remote = f"{supdir}/{fn}"
            local = target_dir / fn
            if local.exists() and local.stat().st_size > 0:
                continue
            retr(ftp, remote, local)

        info["mode"] = mode
        info["downloaded"] = " | ".join(selected)

        try:
            ftp.quit()
        except Exception:
            pass

    except Exception as e:
        info["status"] = "FAIL"
        info["error"] = str(e)

    return info


# =========================
# 3) Scanpy 标准处理流程
# =========================

def ensure_unique_names(adata: sc.AnnData):
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

def add_basic_qc(adata: sc.AnnData):
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

def basic_filter(adata: sc.AnnData, min_genes=200, min_cells=3, max_pct_mt=20.0) -> sc.AnnData:
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs["n_genes_by_counts"] >= min_genes].copy()
    if "pct_counts_mt" in adata.obs.columns:
        adata = adata[adata.obs["pct_counts_mt"] <= max_pct_mt].copy()
    return adata

def run_standard_scanpy(adata: sc.AnnData):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, subset=True, flavor="seurat_v3")

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(N_PCS, max(2, adata.n_vars - 1)), svd_solver="arpack", random_state=RANDOM_STATE)

    sc.pp.neighbors(adata, n_pcs=min(N_PCS, adata.obsm["X_pca"].shape[1]))
    sc.tl.umap(adata, random_state=RANDOM_STATE)
    sc.tl.leiden(adata, resolution=LEIDEN_RES, key_added="leiden")

def rank_markers(adata: sc.AnnData, groupby="leiden", n_genes=100) -> pd.DataFrame:
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", n_genes=n_genes)
    rg = adata.uns["rank_genes_groups"]
    groups = rg["names"].dtype.names

    rows = []
    for g in groups:
        names = rg["names"][g]
        pvals = rg["pvals_adj"][g] if "pvals_adj" in rg else rg["pvals"][g]
        scores = rg["scores"][g] if "scores" in rg else [np.nan] * len(names)
        logfc = rg["logfoldchanges"][g] if "logfoldchanges" in rg else [np.nan] * len(names)

        for i in range(len(names)):
            rows.append(
                {
                    "group": g,
                    "rank": i + 1,
                    "gene": names[i],
                    "logfoldchanges": float(logfc[i]) if logfc[i] is not None else np.nan,
                    "score": float(scores[i]) if scores[i] is not None else np.nan,
                    "pvals_adj": float(pvals[i]) if pvals[i] is not None else np.nan,
                }
            )
    return pd.DataFrame(rows)


# =========================
# 4) 主流程：扫描文件夹 -> 处理所有GSE
# =========================

def main():
    # output dirs
    raw_root = DATASET_ROOT
    data_root = OUT_ROOT / "data"
    qc_root = OUT_ROOT / "qc"
    markers_root = OUT_ROOT / "results" / "markers"
    reports_root = OUT_ROOT / "results" / "reports"
    for p in [data_root, qc_root, markers_root, reports_root]:
        p.mkdir(parents=True, exist_ok=True)

    # detect dataset folders
    gse_dirs = [p for p in raw_root.iterdir() if is_gse_folder(p)]
    gse_dirs = sorted(gse_dirs, key=lambda x: x.name)

    summary = []
    print(f"Found GSE folders: {len(gse_dirs)} under {raw_root}")

    for gse_dir in gse_dirs:
        gse = gse_dir.name.upper()
        t0 = time.time()
        status = "OK"
        loader = ""
        err = ""
        downloaded = ""
        download_mode = ""

        try:
            # If missing matrix/object, optionally FTP download
            try:
                adata, loader = load_anndata_from_folder(gse, gse_dir)
            except FileNotFoundError as e:
                if ENABLE_FTP_DOWNLOAD_IF_MISSING:
                    info = ftp_download_supplementary_for_gse(gse, gse_dir)
                    if info["status"] != "OK":
                        raise FileNotFoundError(f"{e} | FTP download failed: {info['error']}")
                    downloaded = info["downloaded"]
                    download_mode = info["mode"]
                    # try load again
                    adata, loader = load_anndata_from_folder(gse, gse_dir)
                else:
                    raise

            ensure_unique_names(adata)

            # QC
            add_basic_qc(adata)
            (qc_root / gse).mkdir(parents=True, exist_ok=True)
            qc_cols = [c for c in ["total_counts", "n_genes_by_counts", "pct_counts_mt"] if c in adata.obs.columns]
            adata.obs[qc_cols].to_csv(qc_root / gse / "qc_metrics.csv")

            # Filter
            adata = basic_filter(adata, min_genes=MIN_GENES, min_cells=3, max_pct_mt=MAX_PCT_MT)

            # Pipeline
            run_standard_scanpy(adata)

            # Markers (cluster-level)
            mk = rank_markers(adata, groupby="leiden", n_genes=MARKER_TOPN)
            mk.to_csv(markers_root / f"{gse}_markers.csv", index=False, encoding="utf-8-sig")

            # Save processed
            (data_root / gse).mkdir(parents=True, exist_ok=True)
            adata.write_h5ad(data_root / gse / "processed.h5ad")

        except Exception as e:
            status = "FAIL"
            err = str(e)

        sec = round(time.time() - t0, 2)
        summary.append({
            "GSE": gse,
            "status": status,
            "seconds": sec,
            "loader": loader,
            "download_mode": download_mode,
            "downloaded": downloaded,
            "error": err
        })

        print(f"[{gse}] {status} ({sec}s) {loader}")
        if status == "FAIL":
            print(f"  -> {err}")

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(reports_root / "dataset_summary.csv", index=False, encoding="utf-8-sig")
    print(f"Done. Summary saved to: {reports_root / 'dataset_summary.csv'}")


if __name__ == "__main__":
    main()
