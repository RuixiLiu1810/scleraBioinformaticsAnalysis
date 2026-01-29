# Scripts/GEO_fetch.py
from __future__ import annotations

import os
import csv
import gzip
import io
import re
import time
import ftplib
from dataclasses import dataclass, asdict
from typing import List, Dict, Optional
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

from Bio import Entrez


# =========================
# Config
# =========================

@dataclass
class GEOConfig:
    query: str
    out_csv: str
    download_dir: str = "downloads"
    download_types: List[str] = None   # ["matrix", "soft", "raw"]
    retmax: int = 200
    email: str = ""
    api_key: Optional[str] = None
    overwrite: bool = False
    sleep: float = 0.3

    def __post_init__(self):
        if self.download_types is None:
            self.download_types = ["matrix", "soft"]


# =========================
# Metadata structure
# =========================

@dataclass
class GSEMeta:
    GSE: str
    Title: str = ""
    Summary: str = ""
    Overall_design: str = ""
    PubMed_ID: str = ""
    Organism: str = ""
    Platforms: str = ""
    N_Samples: str = ""
    Soft_URL: str = ""
    Matrix_URL: str = ""
    Suppl_URL: str = ""


# =========================
# Main Fetcher
# =========================

class GEOFetcher:
    def __init__(self, config: GEOConfig):
        self.cfg = config

        if not self.cfg.email:
            raise ValueError("NCBI Entrez requires email in config.email")

        Entrez.email = self.cfg.email
        if self.cfg.api_key:
            Entrez.api_key = self.cfg.api_key

    # ---------- Public ----------
    def run(self) -> List[GSEMeta]:
        ids = self._gds_search()
        docs = self._gds_esummary(ids)
        gses = self._extract_gse(docs)

        print(f"[INFO] Found {len(gses)} GSE accessions")

        results = []
        for i, gse in enumerate(gses, 1):
            print(f"[{i}/{len(gses)}] Processing {gse}")
            meta = self._process_gse(gse)
            results.append(meta)

        self._write_csv(results)
        return results

    # ---------- GEO search ----------
    def _gds_search(self) -> List[str]:
        handle = Entrez.esearch(
            db="gds",
            term=self.cfg.query,
            retmax=self.cfg.retmax
        )
        rec = Entrez.read(handle)
        handle.close()
        return rec.get("IdList", [])

    def _gds_esummary(self, ids: List[str]) -> List[dict]:
        if not ids:
            return []

        handle = Entrez.esummary(db="gds", id=",".join(ids))
        rec = Entrez.read(handle)
        handle.close()
        return rec

    def _extract_gse(self, docs: List[dict]) -> List[str]:
        out = []
        for doc in docs:
            for item in doc.get("Item", []):
                if item.attributes.get("Name") == "Accession":
                    if str(item).startswith("GSE"):
                        out.append(str(item))
        return sorted(set(out))

    # ---------- Per GSE ----------
    def _process_gse(self, gse: str) -> GSEMeta:
        paths = self._geo_paths(gse)

        meta = GSEMeta(
            GSE=gse,
            Soft_URL=paths["soft"],
            Matrix_URL=paths["matrix"],
            Suppl_URL=paths["suppl"],
        )

        # parse soft for metadata
        try:
            soft_bytes = self._http_get(paths["soft"])
            kv = self._parse_soft(soft_bytes)
            meta = self._normalize_meta(gse, kv, paths)
        except Exception as e:
            print(f"[WARN] metadata failed for {gse}: {e}")

        # downloads
        if "soft" in self.cfg.download_types:
            self._download(paths["soft"], f"{gse}/soft/{gse}_family.soft.gz")
        if "matrix" in self.cfg.download_types:
            self._download(paths["matrix"], f"{gse}/matrix/{gse}_series_matrix.txt.gz")
        if "raw" in self.cfg.download_types:
            self._download_raw(gse, paths["suppl"])

        time.sleep(self.cfg.sleep)
        return meta

    # ---------- Paths ----------
    def _geo_paths(self, gse: str) -> dict:
        num = gse[3:]
        prefix = num[:-3] if len(num) > 3 else "0"
        bucket = f"GSE{prefix}nnn"
        base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{gse}"
        return {
            "soft": f"{base}/soft/{gse}_family.soft.gz",
            "matrix": f"{base}/matrix/{gse}_series_matrix.txt.gz",
            "suppl": f"{base}/suppl"
        }

    # ---------- HTTP ----------
    def _http_get(self, url: str) -> bytes:
        req = Request(url, headers={"User-Agent": "GEOFetcher"})
        with urlopen(req, timeout=60) as r:
            return r.read()

    def _download(self, url: str, rel_path: str):
        path = os.path.join(self.cfg.download_dir, rel_path)
        os.makedirs(os.path.dirname(path), exist_ok=True)

        if os.path.exists(path) and not self.cfg.overwrite:
            return

        try:
            data = self._http_get(url)
            with open(path, "wb") as f:
                f.write(data)
        except Exception as e:
            print(f"[WARN] download failed {url}: {e}")

    # ---------- SOFT parsing ----------
    def _parse_soft(self, gz_bytes: bytes) -> Dict[str, List[str]]:
        out = {}
        with gzip.GzipFile(fileobj=io.BytesIO(gz_bytes)) as gf:
            for raw in gf:
                line = raw.decode("utf-8", errors="ignore").strip()
                if line.startswith("^SAMPLE"):
                    break
                if line.startswith("!Series_") and "=" in line:
                    k, v = line.split("=", 1)
                    out.setdefault(k.strip(), []).append(v.strip())
        return out

    def _normalize_meta(self, gse: str, kv: dict, paths: dict) -> GSEMeta:
        def join(k): return " ".join(kv.get(k, [])).strip()

        return GSEMeta(
            GSE=gse,
            Title=join("!Series_title"),
            Summary=join("!Series_summary"),
            Overall_design=join("!Series_overall_design"),
            PubMed_ID=join("!Series_pubmed_id"),
            Organism=join("!Series_organism_ch1"),
            Platforms=",".join(kv.get("!Series_platform_id", [])),
            N_Samples=str(len(kv.get("!Series_sample_id", []))),
            Soft_URL=paths["soft"],
            Matrix_URL=paths["matrix"],
            Suppl_URL=paths["suppl"],
        )

    # ---------- RAW tar ----------
    def _download_raw(self, gse: str, suppl_url: str):
        preferred = f"{gse}_RAW.tar"
        url = f"{suppl_url}/{preferred}"
        self._download(url, f"{gse}/suppl/{preferred}")

    # ---------- CSV ----------
    def _write_csv(self, rows: List[GSEMeta]):
        os.makedirs(os.path.dirname(self.cfg.out_csv) or ".", exist_ok=True)
        with open(self.cfg.out_csv, "w", encoding="utf-8-sig", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=asdict(rows[0]).keys())
            writer.writeheader()
            for r in rows:
                writer.writerow(asdict(r))
