"""
Generate an interactive, self-contained HTML web-summary report for a GeneExt run.

Layout inspired by the CellRanger web_summary:
  - Stat cards row  (genes extended, median / mean / max extension, peak counts)
  - Extension-length distribution chart  (histogram, interactive)
  - Peak-coverage distribution chart     (log10 overlay, genic vs non-genic)
  - Mapping-stats comparison table       (only when --estimate was used)

No extra Python dependencies beyond numpy + pandas (already required by the
pipeline).  Charts are rendered by Chart.js loaded from CDN; the rest of the
page is fully functional even without internet access.
"""

from __future__ import annotations

import base64
import datetime
import json
import math
import os

import numpy as np
import pandas as pd


# ─────────────────────────────────────────────────────────────────────────────
# Data-collection helpers
# ─────────────────────────────────────────────────────────────────────────────

def _read_bed_col(path: str, col: int) -> list[float]:
    """Return column *col* (0-indexed) of a BED file as a list of floats."""
    if not os.path.exists(path):
        return []
    try:
        df = pd.read_csv(path, sep="\t", header=None, comment="#")
        if df.shape[1] <= col:
            return []
        return pd.to_numeric(df.iloc[:, col], errors="coerce").dropna().tolist()
    except Exception:
        return []


def _count_bed_lines(path: str) -> int:
    """Return the number of data rows in a BED file (ignores comment lines)."""
    if not os.path.exists(path):
        return 0
    try:
        df = pd.read_csv(path, sep="\t", header=None, comment="#")
        return len(df)
    except Exception:
        return 0


def _read_bed_widths(path: str) -> list[int]:
    """Return peak widths (end - start, bp) from BED columns 1 and 2 (0-indexed)."""
    if not os.path.exists(path):
        return []
    try:
        df = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[1, 2])
        widths = (df.iloc[:, 1] - df.iloc[:, 0]).tolist()
        return [int(w) for w in widths if w > 0]
    except Exception:
        return []


def _histogram(values: list[float], n_bins: int = 40) -> tuple[list[str], list[int]]:
    """Return (mid-point labels, counts) for a linear-scale histogram."""
    if not values:
        return [], []
    counts, edges = np.histogram(values, bins=n_bins)
    labels = [f"{(edges[i] + edges[i + 1]) / 2:.0f}" for i in range(len(edges) - 1)]
    return labels, counts.tolist()


def _log10_histogram(
    genic: list[float],
    noov: list[float],
    n_bins: int = 50,
) -> tuple[list[str], list[int], list[int]]:
    """
    Return unified log10 histogram bins for genic and non-overlapping peaks.
    Both series share the same x-axis edges so they overlay cleanly.
    """
    all_pos = [v for v in genic + noov if v > 0]
    if not all_pos:
        return [], [], []
    log_all = np.log10(all_pos)
    _, edges = np.histogram(log_all, bins=n_bins)

    def _bin(vals: list[float]) -> list[int]:
        pos = [v for v in vals if v > 0]
        if not pos:
            return [0] * (len(edges) - 1)
        c, _ = np.histogram(np.log10(pos), bins=edges)
        return c.tolist()

    labels = [f"{(edges[i] + edges[i + 1]) / 2:.3f}" for i in range(len(edges) - 1)]
    return labels, _bin(genic), _bin(noov)


def _count_bam_reads(bam_path: str) -> int:
    """Count alignment records in a BAM file using samtools view -c."""
    import subprocess
    try:
        r = subprocess.run(
            ["samtools", "view", "-c", bam_path],
            capture_output=True, text=True, timeout=120,
        )
        if r.returncode == 0:
            return int(r.stdout.strip())
    except Exception:
        pass
    return 0


def _get_reads_info(tempdir: str) -> dict:
    """
    Return the number of reads used for peak calling and whether the BAM was subsampled.

    Priority:
      1. subsampled.bam  — used when --subsamplebam was set
      2. plus.bam + minus.bam — the strand-split BAMs actually fed to MACS2
      3. MACS2 xls "total tags" — fallback if BAM files are unavailable
    """
    subsampled_bam = os.path.join(tempdir, "subsampled.bam")
    subsampled = os.path.exists(subsampled_bam)

    # Option 1: subsampled.bam
    if subsampled:
        n = _count_bam_reads(subsampled_bam)
        if n > 0:
            return {"n_reads": n, "subsampled": True}

    # Option 2: plus.bam + minus.bam
    total = 0
    for fname in ("plus.bam", "minus.bam"):
        p = os.path.join(tempdir, fname)
        if os.path.exists(p):
            total += _count_bam_reads(p)
    if total > 0:
        return {"n_reads": total, "subsampled": subsampled}

    # Option 3: MACS2 xls fallback
    xls_total = 0
    for fname in ("plus_peaks.xls", "minus_peaks.xls"):
        path = os.path.join(tempdir, fname)
        if not os.path.exists(path):
            continue
        try:
            with open(path) as fh:
                for line in fh:
                    if line.startswith("# total tags in treatment:"):
                        xls_total += int(line.split(":")[1].strip())
                        break
        except Exception:
            pass
    return {"n_reads": xls_total, "subsampled": subsampled}


def _parse_mapping_stats(path: str) -> list[dict]:
    """
    Parse mapping_stats.txt (written by run_estimate).
    Returns a list of dicts — one per annotation (before / after).
    """
    if not os.path.exists(path):
        return []
    results: list[dict] = []
    cur: dict = {}

    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            # Header line looks like: "/path/to/file.gtf:"
            if line.endswith(":") and "/" in line:
                if cur:
                    results.append(cur)
                cur = {"label": line.rstrip(":")}
            elif line.startswith("Total reads"):
                cur["total"] = int(line.split(": ")[1].split()[0])
            elif line.startswith("Mapped reads"):
                parts = line.split("(total: ")
                cur["mapped"] = int(parts[0].split(": ")[1].strip())
                cur["mapped_pct"] = float(parts[1].split("%")[0].strip())
            elif line.startswith("Genic reads"):
                parts = line.split("(total: ")
                cur["genic"] = int(parts[0].split(": ")[1].strip())
                cur["genic_pct"] = float(parts[1].split("%")[0].strip())
            elif line.startswith("Orphan peaks"):
                parts = line.split("(total: ")
                cur["orphan"] = int(parts[0].split(": ")[1].strip())
                cur["orphan_pct"] = float(parts[1].split("%")[0].strip())
            elif line.startswith("Intergenic reads"):
                parts = line.split("(total: ")
                cur["intergenic"] = int(parts[0].split(": ")[1].strip())
                cur["intergenic_pct"] = float(parts[1].split("%")[0].strip())

    if cur:
        results.append(cur)
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Main public function
# ─────────────────────────────────────────────────────────────────────────────

def generate_html_report(
    tempdir: str,
    outputfile: str,
    genefile: str,
    infmt: str,
    coverage_percentile: float,
    count_threshold=None,
    do_estimate: bool = False,
    n_genes: int = 0,
    run_args: str = "",
) -> str:
    """
    Build the web-summary HTML and write it next to the output annotation.

    Parameters
    ----------
    tempdir            : directory containing intermediate pipeline files
    outputfile         : final annotation path (report lands at outputfile + '.Report.html')
    genefile           : original input annotation path
    infmt              : 'gtf' | 'gff' | 'bed'
    coverage_percentile: percentile used for peak filtering  (0 = disabled)
    count_threshold    : actual coverage cutoff value (may be None)
    do_estimate        : whether mapping estimation was run
    n_genes            : total gene count in the original annotation
    run_args           : command-line string for display in the report
    """
    # ── Extension lengths + table ─────────────────────────────────────────────
    ext_path = os.path.join(tempdir, "extensions.tsv")
    ext_lengths: list[float] = []
    ext_table: list[dict] = []
    if os.path.exists(ext_path):
        try:
            df_ext = pd.read_csv(ext_path, sep="\t", header=None)
            ext_lengths = (
                pd.to_numeric(df_ext.iloc[:, -1], errors="coerce").dropna().tolist()
            )
            if df_ext.shape[1] >= 3:
                for _, row in df_ext.iterrows():
                    ext_val = pd.to_numeric(row.iloc[2], errors="coerce")
                    ext_table.append({
                        "gene": str(row.iloc[0]),
                        "peak": str(row.iloc[1]),
                        "ext":  int(ext_val) if not pd.isna(ext_val) else 0,
                    })
        except Exception:
            pass

    # ── Logo (embedded as base64 for self-contained HTML) ────────────────────
    logo_b64 = ""
    logo_path = os.path.join(os.path.dirname(__file__), "..", "img", "logo.png")
    try:
        with open(logo_path, "rb") as _f:
            logo_b64 = base64.b64encode(_f.read()).decode()
    except Exception:
        pass

    n_extended  = len(ext_lengths)
    pct_extended = round(n_extended / n_genes * 100, 1) if n_genes > 0 else 0.0
    median_ext   = round(float(np.median(ext_lengths)), 1) if ext_lengths else 0.0
    mean_ext     = round(float(np.mean(ext_lengths)),   1) if ext_lengths else 0.0
    max_ext      = round(float(np.max(ext_lengths)),    1) if ext_lengths else 0.0

    ext_labels, ext_counts = _histogram(ext_lengths, n_bins=40)

    # ── Peak coverages ───────────────────────────────────────────────────────
    genic_cov = _read_bed_col(os.path.join(tempdir, "genic_peaks.bed"),   6)
    noov_cov  = _read_bed_col(os.path.join(tempdir, "allpeaks_noov.bed"), 6)

    cov_labels, cov_genic, cov_noov = _log10_histogram(genic_cov, noov_cov, n_bins=50)

    log10_thr: float | None = None
    if count_threshold is not None:
        try:
            v = float(count_threshold)
            if v > 0:
                log10_thr = round(math.log10(v), 4)
        except Exception:
            pass

    # ── Orphan peaks (prefer merged clusters) ────────────────────────────────
    orphan_bed_text = ""
    for _fname in ("orphan_merged.bed", "orphan.bed"):
        _p = os.path.join(tempdir, _fname)
        if os.path.exists(_p):
            try:
                with open(_p) as _f:
                    orphan_bed_text = _f.read()
            except Exception:
                pass
            break
    n_orphan_merged = orphan_bed_text.count("\n") if orphan_bed_text else 0

    # ── Reads (from MACS2 xls) ───────────────────────────────────────────────
    reads_info = _get_reads_info(tempdir)

    # ── Mapping stats ────────────────────────────────────────────────────────
    mapping_stats = (
        _parse_mapping_stats(os.path.join(tempdir, "mapping_stats.txt"))
        if do_estimate
        else []
    )

    # ── Assemble payload ─────────────────────────────────────────────────────
    payload = {
        "summary": {
            "n_genes":        n_genes,
            "n_extended":     n_extended,
            "pct_extended":   pct_extended,
            "median_ext":     median_ext,
            "mean_ext":       mean_ext,
            "max_ext":        max_ext,
            "n_genic_peaks":  len(genic_cov),
            "n_noov_peaks":   len(noov_cov),
            "n_orphan_peaks": n_orphan_merged,
            "cov_percentile": coverage_percentile,
            "cov_threshold":  count_threshold,
            "n_reads":        reads_info["n_reads"],
            "subsampled":     reads_info["subsampled"],
            "output_file":    os.path.basename(outputfile),
            "input_file":     os.path.basename(genefile),
            "run_date":       datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "run_args":       run_args,
        },
        "ext_hist":  {"labels": ext_labels, "counts": ext_counts},
        "cov_hist":  {
            "labels":         cov_labels,
            "counts_genic":   cov_genic,
            "counts_noov":    cov_noov,
            "log10_threshold": log10_thr,
        },
        "mapping_stats": mapping_stats,
        "ext_table": ext_table,
        "orphan_bed": orphan_bed_text,
    }

    html = _render_html(payload, logo_b64=logo_b64)
    report_path = outputfile + ".Report.html"
    with open(report_path, "w") as fh:
        fh.write(html)
    return report_path


# ─────────────────────────────────────────────────────────────────────────────
# HTML template
# ─────────────────────────────────────────────────────────────────────────────

def _render_html(data: dict, logo_b64: str = "") -> str:
    payload_json = json.dumps(data, indent=None, separators=(",", ":"))
    s = data["summary"]
    logo_html = (
        f'<img src="data:image/png;base64,{logo_b64}"'
        f' style="height:34px;width:auto;display:block">'
        if logo_b64 else
        '<div class="logo-icon">G</div>'
    )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>GeneExt — {s['output_file']}</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;
      background:#f0f2f5;color:#1a2332;min-height:100vh}}

/* ── header ─────────────────────────────────────────────────────────────── */
header{{background:linear-gradient(135deg,#0f3460 0%,#16213e 60%,#0a1628 100%);
        color:#fff;padding:16px 36px;display:flex;align-items:stretch;
        justify-content:space-between;box-shadow:0 2px 12px rgba(0,0,0,.35)}}
.logo{{display:flex;align-items:center}}
.logo img{{height:100%;max-height:64px;min-height:40px;width:auto;display:block;
           filter:brightness(0) invert(1)}}
.logo-icon{{width:50px;height:50px;background:linear-gradient(135deg,#00b4d8,#48cae4);
            border-radius:12px;display:flex;align-items:center;justify-content:center;
            font-size:26px;font-weight:900;color:#0f3460;flex-shrink:0}}
.run-meta{{text-align:right;font-size:.72rem;opacity:.65;line-height:1.8;
           display:flex;flex-direction:column;justify-content:center}}

/* ── main container ──────────────────────────────────────────────────────── */
main{{max-width:1200px;margin:28px auto;padding:0 24px}}

/* ── section title ───────────────────────────────────────────────────────── */
.section-title{{font-size:.7rem;font-weight:700;letter-spacing:1.2px;
                text-transform:uppercase;color:#5a7194;margin-bottom:12px;
                padding-bottom:6px;border-bottom:2px solid #dce3ee}}

/* ── stat cards ──────────────────────────────────────────────────────────── */
.card-group-label{{font-size:.65rem;font-weight:700;letter-spacing:1px;
                   text-transform:uppercase;color:#8695a8;margin:20px 0 8px;
                   padding-left:2px}}
.card-group-label:first-child{{margin-top:0}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));
        gap:14px;margin-bottom:6px}}
.card{{background:#fff;border-radius:12px;padding:18px 20px;
       box-shadow:0 1px 4px rgba(0,0,0,.08);border-top:3px solid transparent;
       transition:box-shadow .15s}}
.card:hover{{box-shadow:0 4px 16px rgba(0,0,0,.12)}}
.card.accent-teal  {{border-top-color:#00b4d8}}
.card.accent-green {{border-top-color:#06d6a0}}
.card.accent-orange{{border-top-color:#ff9f43}}
.card.accent-purple{{border-top-color:#a55eea}}
.card.accent-blue  {{border-top-color:#4d9ef5}}
.card.accent-red   {{border-top-color:#ee5a6f}}
.card-value{{font-size:1.9rem;font-weight:800;line-height:1;color:#1a2332}}
.card-value span{{font-size:1rem;font-weight:600;color:#8695a8;margin-left:3px}}
.card-label{{font-size:.72rem;font-weight:600;color:#8695a8;
             text-transform:uppercase;letter-spacing:.6px;margin-top:6px}}

/* ── charts ──────────────────────────────────────────────────────────────── */
.charts{{display:grid;grid-template-columns:1fr 1fr;gap:20px;margin-bottom:28px}}
@media(max-width:780px){{.charts{{grid-template-columns:1fr}}}}
.chart-card{{background:#fff;border-radius:12px;padding:20px 22px;
             box-shadow:0 1px 4px rgba(0,0,0,.08)}}
.chart-card h3{{font-size:.85rem;font-weight:700;color:#1a2332;margin-bottom:4px}}
.chart-card p{{font-size:.72rem;color:#8695a8;margin-bottom:14px}}
.chart-wrap{{position:relative;height:260px}}

/* ── mapping table ───────────────────────────────────────────────────────── */
.table-card{{background:#fff;border-radius:12px;padding:20px 22px;
             box-shadow:0 1px 4px rgba(0,0,0,.08);margin-bottom:28px}}
.table-card h3{{font-size:.85rem;font-weight:700;color:#1a2332;margin-bottom:12px}}
table{{width:100%;border-collapse:collapse;font-size:.82rem}}
thead tr{{background:#f4f7fb}}
th{{padding:9px 14px;text-align:left;font-weight:700;color:#5a7194;
    font-size:.72rem;text-transform:uppercase;letter-spacing:.5px}}
td{{padding:9px 14px;border-bottom:1px solid #eef1f7;color:#2d3e50}}
tr:last-child td{{border-bottom:none}}
.pct-bar{{display:flex;align-items:center;gap:8px}}
.bar-bg{{flex:1;height:7px;background:#eef1f7;border-radius:4px;overflow:hidden}}
.bar-fill{{height:100%;border-radius:4px;background:#4d9ef5}}
.bar-fill.genic{{background:#06d6a0}}
.bar-fill.intergenic{{background:#ee5a6f}}
.bar-fill.orphan{{background:#ff9f43}}
td.num{{text-align:right;font-variant-numeric:tabular-nums;font-weight:600}}
.diff-row{{background:#f9fbff;font-weight:700}}
.diff-row td{{color:#0f3460}}

/* ── args box ────────────────────────────────────────────────────────────── */
.args-box{{background:#fff;border-radius:12px;padding:16px 20px;
           box-shadow:0 1px 4px rgba(0,0,0,.08);margin-bottom:28px}}
.args-box pre{{font-size:.73rem;color:#5a7194;white-space:pre-wrap;
               word-break:break-all;line-height:1.55}}

/* ── footer ──────────────────────────────────────────────────────────────── */
footer{{text-align:center;font-size:.7rem;color:#aab4c0;padding:20px 0 32px}}

/* ── chart-unavailable banner ────────────────────────────────────────────── */
.chart-unavailable{{display:none;align-items:center;justify-content:center;
                    height:100%;font-size:.8rem;color:#aab4c0;flex-direction:column;gap:6px}}

/* ── clickable card ──────────────────────────────────────────────────────── */
.card.clickable{{cursor:pointer}}
.card.clickable:hover{{box-shadow:0 6px 20px rgba(0,0,0,.16);transform:translateY(-1px)}}

/* ── modal ───────────────────────────────────────────────────────────────── */
.modal-backdrop{{display:none;position:fixed;inset:0;background:rgba(10,22,40,.55);
                 z-index:1000;align-items:center;justify-content:center}}
.modal-backdrop.open{{display:flex}}
.modal{{background:#fff;border-radius:14px;width:min(860px,94vw);max-height:82vh;
        display:flex;flex-direction:column;box-shadow:0 12px 48px rgba(0,0,0,.28)}}
.modal-header{{display:flex;align-items:center;justify-content:space-between;
               padding:18px 24px;border-bottom:1px solid #dce3ee}}
.modal-header h2{{font-size:.95rem;font-weight:700;color:#1a2332}}
.modal-close{{border:none;background:none;font-size:1.3rem;color:#8695a8;
              cursor:pointer;line-height:1;padding:2px 6px;border-radius:6px}}
.modal-close:hover{{background:#f0f2f5;color:#1a2332}}
.modal-search{{padding:12px 24px;border-bottom:1px solid #eef1f7}}
.modal-search input{{width:100%;padding:8px 12px;border:1px solid #dce3ee;
                     border-radius:8px;font-size:.82rem;outline:none}}
.modal-search input:focus{{border-color:#00b4d8}}
.modal-body{{overflow-y:auto;padding:0}}
.modal-body table{{width:100%;border-collapse:collapse;font-size:.82rem}}
.modal-body thead tr{{background:#f4f7fb;position:sticky;top:0}}
.modal-body th{{padding:9px 16px;text-align:left;font-weight:700;color:#5a7194;
                font-size:.72rem;text-transform:uppercase;letter-spacing:.5px;
                cursor:pointer;user-select:none;white-space:nowrap}}
.modal-body th:hover{{color:#0f3460}}
.modal-body th .sort-icon{{margin-left:4px;opacity:.4}}
.modal-body th.sort-active .sort-icon{{opacity:1}}
.modal-body td{{padding:8px 16px;border-bottom:1px solid #eef1f7;color:#2d3e50}}
.modal-body tr:last-child td{{border-bottom:none}}
.modal-body tr:hover td{{background:#fafbff}}
.modal-footer{{padding:10px 24px;border-top:1px solid #eef1f7;font-size:.72rem;color:#8695a8}}
</style>
</head>
<body>

<header>
  <div class="logo">
    {logo_html}
  </div>
  <div class="run-meta">
    <div>Output: <strong>{s['output_file']}</strong></div>
    <div>Input: {s['input_file']}</div>
    <div>{s['run_date']}</div>
  </div>
</header>

<main>

  <!-- ── Summary cards ─────────────────────────────────────────────────── -->
  <p class="section-title">Run summary</p>
  <div id="summaryCards"></div>

  <!-- ── Charts ────────────────────────────────────────────────────────── -->
  <p class="section-title">Distributions</p>
  <div class="charts">
    <div class="chart-card">
      <h3>Gene-extension length distribution</h3>
      <p>Number of base-pairs each gene was extended at its 3&prime; end</p>
      <div class="chart-wrap">
        <canvas id="extChart"></canvas>
        <div class="chart-unavailable" id="extNA">
          <span>&#9888;</span>No extension data available
        </div>
      </div>
    </div>
    <div class="chart-card">
      <h3>Peak-coverage distribution (log&#x2081;&#x2080;)</h3>
      <p>
        <span style="color:#ff6384">&#9632;</span> Genic peaks &nbsp;
        <span style="color:#4d9ef5">&#9632;</span> Non-overlapping peaks
        <span id="thrLabel"></span>
      </p>
      <div class="chart-wrap">
        <canvas id="covChart"></canvas>
        <div class="chart-unavailable" id="covNA">
          <span>&#9888;</span>No coverage data available
        </div>
      </div>
    </div>
  </div>

  <!-- ── Mapping table (only if --estimate) ────────────────────────────── -->
  <div id="mappingSection"></div>

  <!-- ── Command args ──────────────────────────────────────────────────── -->
  <div id="argsSection"></div>

</main>

<!-- ── Extended genes modal ──────────────────────────────────────────────── -->
<div class="modal-backdrop" id="extModal" role="dialog" aria-modal="true">
  <div class="modal">
    <div class="modal-header">
      <h2>Extended genes</h2>
      <button class="modal-close" id="extModalClose" aria-label="Close">&times;</button>
    </div>
    <div class="modal-search">
      <input type="search" id="extSearch" placeholder="Filter by gene or peak ID&hellip;">
    </div>
    <div class="modal-body">
      <table id="extTable">
        <thead>
          <tr>
            <th data-col="gene">Gene ID <span class="sort-icon">&#9661;</span></th>
            <th data-col="peak">Peak ID <span class="sort-icon">&#9661;</span></th>
            <th data-col="ext" class="sort-active">Extension (bp) <span class="sort-icon">&#9660;</span></th>
          </tr>
        </thead>
        <tbody id="extTableBody"></tbody>
      </table>
    </div>
    <div class="modal-footer" id="extTableFooter"></div>
  </div>
</div>

<footer>Generated by GeneExt &bull; {s['run_date']}</footer>

<!-- Chart.js via CDN (requires internet); charts degrade gracefully offline -->
<script
  src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"
  crossorigin="anonymous"
  onerror="window.CHARTJS_FAILED=true"
></script>

<script>
const D = {payload_json};
const S = D.summary;

// ── Stat cards (grouped) ────────────────────────────────────────────────────
const readsLabel = S.n_reads
  ? (S.n_reads.toLocaleString() + (S.subsampled ? ' &#9432;' : ''))
  : '—';
const readsHint = S.subsampled
  ? '<div style="font-size:.62rem;color:#ff9f43;margin-top:4px">&#9888; BAM was subsampled</div>'
  : '';

const CARD_GROUPS = [
  {{
    title: 'Gene Extension',
    cards: [
      {{value: S.n_extended, suffix: '/' + S.n_genes, label: 'Genes extended',     cls: 'accent-teal',   id: 'cardGenesExtended', clickable: true}},
      {{value: S.pct_extended, suffix: '%',            label: 'Extension rate',     cls: 'accent-green'}},
      {{value: S.median_ext,   suffix: 'bp',           label: 'Median extension',   cls: 'accent-orange'}},
      {{value: S.mean_ext,     suffix: 'bp',           label: 'Mean extension',     cls: 'accent-blue'}},
      {{value: S.max_ext,      suffix: 'bp',           label: 'Max extension',      cls: 'accent-purple'}},
    ]
  }},
  {{
    title: 'Peaks',
    cards: [
      {{value: S.n_genic_peaks, suffix: '',  label: 'Genic peaks',              cls: 'accent-teal'}},
      {{value: S.n_noov_peaks,  suffix: '',  label: 'Non-overlapping peaks',    cls: 'accent-blue'}},
      ...(S.n_orphan_peaks ? [{{value: S.n_orphan_peaks, suffix: '', label: 'Orphan peak clusters', cls: 'accent-purple'}}] : []),
      {{value: S.cov_percentile || '—', suffix: S.cov_percentile ? 'th pct' : '', label: 'Coverage percentile', cls: 'accent-red'}},
    ]
  }},
  {{
    title: 'Reads',
    cards: [
      {{value: readsLabel, suffix: '', label: 'Reads used for peak calling', cls: 'accent-teal', extraHint: readsHint}},
    ]
  }},
];

(function buildCards() {{
  const wrap = document.getElementById('summaryCards');
  CARD_GROUPS.forEach(g => {{
    if (!g.cards || !g.cards.length) return;
    wrap.innerHTML += `<p class="card-group-label">${{g.title}}</p>`;
    let row = '<div class="cards">';
    g.cards.forEach(c => {{
      const extra = c.clickable ? ' clickable' : '';
      const idAttr = c.id ? ` id="${{c.id}}"` : '';
      const hint = c.clickable
        ? '<div style="font-size:.62rem;color:#00b4d8;margin-top:4px">Click to view &#8599;</div>'
        : (c.extraHint || '');
      row += `<div class="card ${{c.cls}}${{extra}}"${{idAttr}}>
        <div class="card-value">${{c.value}}<span>${{c.suffix}}</span></div>
        <div class="card-label">${{c.label}}</div>
        ${{hint}}
      </div>`;
    }});
    row += '</div>';
    wrap.innerHTML += row;
  }});
  const geCard = document.getElementById('cardGenesExtended');
  if (geCard) geCard.addEventListener('click', openExtModal);
}})();

// ── Extended genes modal ────────────────────────────────────────────────────
let extSortCol = 'ext', extSortAsc = false;
let extFilteredRows = D.ext_table ? [...D.ext_table] : [];

function renderExtTable() {{
  const q = (document.getElementById('extSearch').value || '').toLowerCase();
  extFilteredRows = (D.ext_table || []).filter(r =>
    r.gene.toLowerCase().includes(q) || r.peak.toLowerCase().includes(q)
  );
  extFilteredRows.sort((a, b) => {{
    const va = a[extSortCol], vb = b[extSortCol];
    if (typeof va === 'number') return extSortAsc ? va - vb : vb - va;
    return extSortAsc ? String(va).localeCompare(vb) : String(vb).localeCompare(va);
  }});
  const tbody = document.getElementById('extTableBody');
  tbody.innerHTML = extFilteredRows.map(r => `
    <tr>
      <td>${{r.gene}}</td>
      <td>${{r.peak}}</td>
      <td style="text-align:right;font-variant-numeric:tabular-nums;font-weight:600">${{r.ext.toLocaleString()}}</td>
    </tr>`).join('');
  document.getElementById('extTableFooter').textContent =
    extFilteredRows.length + ' of ' + (D.ext_table || []).length + ' genes';
  // update sort icons
  document.querySelectorAll('#extTable th').forEach(th => {{
    const icon = th.querySelector('.sort-icon');
    if (th.dataset.col === extSortCol) {{
      th.classList.add('sort-active');
      icon.innerHTML = extSortAsc ? '&#9650;' : '&#9660;';
    }} else {{
      th.classList.remove('sort-active');
      icon.innerHTML = '&#9661;';
    }}
  }});
}}

function openExtModal() {{
  document.getElementById('extModal').classList.add('open');
  document.getElementById('extSearch').value = '';
  extSortCol = 'ext'; extSortAsc = false;
  renderExtTable();
}}

function closeExtModal() {{
  document.getElementById('extModal').classList.remove('open');
}}

document.getElementById('extModalClose').addEventListener('click', closeExtModal);
document.getElementById('extModal').addEventListener('click', e => {{
  if (e.target === document.getElementById('extModal')) closeExtModal();
}});
document.addEventListener('keydown', e => {{
  if (e.key === 'Escape') closeExtModal();
}});
document.getElementById('extSearch').addEventListener('input', renderExtTable);
document.querySelectorAll('#extTable th').forEach(th => {{
  th.addEventListener('click', () => {{
    if (extSortCol === th.dataset.col) extSortAsc = !extSortAsc;
    else {{ extSortCol = th.dataset.col; extSortAsc = th.dataset.col !== 'ext'; }}
    renderExtTable();
  }});
}});

// ── Charts ─────────────────────────────────────────────────────────────────
document.addEventListener('DOMContentLoaded', function() {{
  if (window.CHARTJS_FAILED || typeof Chart === 'undefined') {{
    document.getElementById('extNA').style.display = 'flex';
    document.getElementById('covNA').style.display = 'flex';
    return;
  }}

  Chart.defaults.font.family =
    "-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif";
  Chart.defaults.font.size = 11;
  Chart.defaults.color = '#5a7194';

  const EH = D.ext_hist;
  if (EH.counts && EH.counts.length) {{
    new Chart(document.getElementById('extChart'), {{
      type: 'bar',
      data: {{
        labels: EH.labels,
        datasets: [{{
          label: 'Genes',
          data: EH.counts,
          backgroundColor: 'rgba(0,180,216,0.65)',
          borderColor:     'rgba(0,180,216,0.9)',
          borderWidth: 1,
          borderRadius: 2,
        }}]
      }},
      options: {{
        responsive: true,
        maintainAspectRatio: false,
        plugins: {{
          legend: {{display: false}},
          tooltip: {{
            callbacks: {{
              title: ctx => `~${{ctx[0].label}} bp`,
              label: ctx => ` ${{ctx.parsed.y}} gene(s)`
            }}
          }}
        }},
        scales: {{
          x: {{
            title: {{display:true, text:'Extension length (bp)', padding:{{top:6}}}},
            ticks: {{maxTicksLimit: 8, maxRotation: 0}},
            grid: {{display: false}}
          }},
          y: {{
            title: {{display:true, text:'Number of genes'}},
            beginAtZero: true,
            grid: {{color:'#eef1f7'}}
          }}
        }}
      }}
    }});
  }} else {{
    document.getElementById('extNA').style.display = 'flex';
  }}

  const CH = D.cov_hist;
  if (CH.labels && CH.labels.length) {{
    // vertical annotation line for coverage threshold
    const thrPlugin = {{
      id: 'thrLine',
      afterDraw(chart) {{
        const thr = CH.log10_threshold;
        if (thr == null) return;
        const xAxis = chart.scales.x;
        const yAxis = chart.scales.y;
        const labels = chart.data.labels;
        // find nearest label index
        let nearest = 0, minD = Infinity;
        labels.forEach((l, i) => {{
          const d = Math.abs(parseFloat(l) - thr);
          if (d < minD) {{ minD = d; nearest = i; }}
        }});
        const x = xAxis.getPixelForValue(nearest);
        const ctx2 = chart.ctx;
        ctx2.save();
        ctx2.setLineDash([5, 4]);
        ctx2.strokeStyle = '#ee5a6f';
        ctx2.lineWidth = 1.5;
        ctx2.beginPath();
        ctx2.moveTo(x, yAxis.top);
        ctx2.lineTo(x, yAxis.bottom);
        ctx2.stroke();
        ctx2.restore();
      }}
    }};

    if (CH.log10_threshold != null) {{
      document.getElementById('thrLabel').innerHTML =
        ` &nbsp;<span style="color:#ee5a6f">&#9135;</span> threshold`;
    }}

    new Chart(document.getElementById('covChart'), {{
      type: 'bar',
      plugins: [thrPlugin],
      data: {{
        labels: CH.labels,
        datasets: [
          {{
            label: 'Genic peaks',
            data: CH.counts_genic,
            backgroundColor: 'rgba(255,99,132,0.5)',
            borderColor:     'rgba(255,99,132,0.8)',
            borderWidth: 1,
            borderRadius: 2,
          }},
          {{
            label: 'Non-overlapping peaks',
            data: CH.counts_noov,
            backgroundColor: 'rgba(77,158,245,0.45)',
            borderColor:     'rgba(77,158,245,0.8)',
            borderWidth: 1,
            borderRadius: 2,
          }}
        ]
      }},
      options: {{
        responsive: true,
        maintainAspectRatio: false,
        plugins: {{
          legend: {{display: false}},
          tooltip: {{
            callbacks: {{
              title: ctx => `log10(cov) ≈ ${{ctx[0].label}}`
            }}
          }}
        }},
        scales: {{
          x: {{
            title: {{display:true, text:'log₁₀(normalized coverage)', padding:{{top:6}}}},
            ticks: {{maxTicksLimit: 8, maxRotation: 0}},
            grid: {{display: false}}
          }},
          y: {{
            title: {{display:true, text:'Number of peaks'}},
            beginAtZero: true,
            grid: {{color:'#eef1f7'}}
          }}
        }}
      }}
    }});
  }} else {{
    document.getElementById('covNA').style.display = 'flex';
  }}
}});

// ── Mapping stats table ────────────────────────────────────────────────────
(function buildMappingTable() {{
  const rows = D.mapping_stats;
  if (!rows || rows.length === 0) return;

  const cols = ['Total reads','Mapped','Genic','Orphan peaks','Intergenic'];
  const keys = ['total','mapped','genic','orphan','intergenic'];
  const pcts  = [null,'mapped_pct','genic_pct','orphan_pct','intergenic_pct'];
  const cls   = [null,'',          'genic',    'orphan',     'intergenic'];

  let html = `<div class="table-card">
    <h3>Mapping statistics</h3>
    <table>
      <thead><tr><th>Annotation</th>${{cols.map(c=>`<th>${{c}}</th>`).join('')}}</tr></thead>
      <tbody>`;

  rows.forEach(r => {{
    html += `<tr><td>${{r.label.split('/').pop()}}</td>`;
    keys.forEach((k, i) => {{
      const v   = r[k]       ?? '—';
      const pct = pcts[i] ? (r[pcts[i]] ?? null) : null;
      if (pct !== null) {{
        html += `<td><div class="pct-bar">
          <span>${{v.toLocaleString()}}</span>
          <div class="bar-bg"><div class="bar-fill ${{cls[i]}}" style="width:${{Math.min(pct,100)}}%"></div></div>
          <span style="min-width:42px;text-align:right">${{pct}}%</span>
        </div></td>`;
      }} else {{
        html += `<td class="num">${{typeof v==='number'?v.toLocaleString():v}}</td>`;
      }}
    }});
    html += `</tr>`;
  }});

  // diff row
  if (rows.length === 2 && rows[0].intergenic_pct != null && rows[1].intergenic_pct != null) {{
    const diff = Math.abs(rows[0].intergenic_pct - rows[1].intergenic_pct).toFixed(2);
    html += `<tr class="diff-row">
      <td colspan="${{cols.length + 1}}">
        &#x25B2; Intergenic read proportion reduced by ${{diff}}%
        (from ${{rows[0].intergenic_pct}}% → ${{rows[1].intergenic_pct}}%)
      </td></tr>`;
  }}

  html += `</tbody></table></div>`;
  document.getElementById('mappingSection').innerHTML = html;
}})();

// ── Command-line args ──────────────────────────────────────────────────────
(function buildArgs() {{
  if (!S.run_args) return;
  document.getElementById('argsSection').innerHTML = `
    <p class="section-title">Command</p>
    <div class="args-box"><pre>${{S.run_args}}</pre></div>`;
}})();
</script>
</body>
</html>"""
