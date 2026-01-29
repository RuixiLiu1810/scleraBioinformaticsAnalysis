import os
import re
import shutil
import pandas as pd

# -----------------------------
# Config
# -----------------------------
CSV_PATH = r"E:\20260120\GEO_dataset_evaluation_UPDATED.csv"
SRC_ROOT_UNCLASSIFIED = r"E:\20260120"     # 未分类文件所在根目录（按 Group 名分文件夹）
DST_ROOT = r"E:\valid"                      # 分类后的目标根目录

# 你本地的 GSE 文件命名可能有多种：GSE22818_family.soft.gz / GSE22818_*.gz / GSE22818*
# 这里用正则：以 GSExxxxx 开头即可（更鲁棒）
GSE_FILE_PREFIX_RE = re.compile(r"^(GSE\d+)", re.IGNORECASE)

def safe_mkdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def find_files_for_gse(group_dir: str, gse: str) -> list[str]:
    """
    在 group_dir 内查找属于该 gse 的所有文件：
    - 文件名以 GSE#### 开头
    - 不递归子目录（如需递归可改为 os.walk）
    """
    hits = []
    if not os.path.isdir(group_dir):
        return hits

    for fn in os.listdir(group_dir):
        m = GSE_FILE_PREFIX_RE.match(fn)
        if m and m.group(1).upper() == gse.upper():
            hits.append(os.path.join(group_dir, fn))
    return hits

def move_files(files: list[str], dst_dir: str) -> None:
    safe_mkdir(dst_dir)
    for fp in files:
        try:
            shutil.move(fp, os.path.join(dst_dir, os.path.basename(fp)))
        except Exception as e:
            print(f"[WARN] Move failed: {fp} -> {dst_dir}. Reason: {e}")

def main():
    if not os.path.isfile(CSV_PATH):
        raise FileNotFoundError(f"CSV not found: {CSV_PATH}")

    df = pd.read_csv(CSV_PATH, encoding="utf-8-sig")

    required_cols = {"Group", "GSE", "Suitability_Level"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}. Got: {list(df.columns)}")

    # Normalize
    df["Group"] = df["Group"].astype(str).str.strip()
    df["GSE"] = df["GSE"].astype(str).str.strip()
    df["Suitability_Level"] = df["Suitability_Level"].astype(str).str.strip()

    # Process each group separately
    for group_name, sub in df.groupby("Group"):
        group_src_dir = os.path.join(SRC_ROOT_UNCLASSIFIED, group_name)
        group_unclassified_dir = os.path.join(SRC_ROOT_UNCLASSIFIED, group_name)  # 未分类保留在此目录

        if not os.path.isdir(group_src_dir):
            print(f"[WARN] Group source folder not found, skip: {group_src_dir}")
            continue

        # Classify per Suitability_Level
        for level, sub2 in sub.groupby("Suitability_Level"):
            # 目标：E:\valid\group名\suitability类型\
            dst_dir = os.path.join(DST_ROOT, group_name, level)
            safe_mkdir(dst_dir)

            for gse in sub2["GSE"].unique():
                if not gse or gse.lower() == "nan":
                    continue

                gse_files = find_files_for_gse(group_src_dir, gse)
                if not gse_files:
                    # 没找到就不动，留在 E:\20260120\group名（未分类目录）
                    print(f"[INFO] No files found for {gse} in {group_src_dir}; keep unclassified in {group_unclassified_dir}")
                    continue

                move_files(gse_files, dst_dir)
                print(f"[OK] {group_name} | {level} | {gse}: moved {len(gse_files)} file(s) -> {dst_dir}")

        # 未被移动的文件将自然留在 E:\20260120\group名 中，无需额外操作
        print(f"[DONE] Group finished: {group_name}. Remaining files in {group_unclassified_dir} are unclassified.")

if __name__ == "__main__":
    main()
