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
    outputfile         : final annotation path (report lands at outputfile + '.web_summary.html')
    genefile           : original input annotation path
    infmt              : 'gtf' | 'gff' | 'bed'
    coverage_percentile: percentile used for peak filtering  (0 = disabled)
    count_threshold    : actual coverage cutoff value (may be None)
    do_estimate        : whether mapping estimation was run
    n_genes            : total gene count in the original annotation
    run_args           : command-line string for display in the report
    """
    # ── Extension lengths ────────────────────────────────────────────────────
    ext_path = os.path.join(tempdir, "extensions.tsv")
    ext_lengths: list[float] = []
    if os.path.exists(ext_path):
        try:
            df_ext = pd.read_csv(ext_path, sep="\t", header=None)
            ext_lengths = (
                pd.to_numeric(df_ext.iloc[:, -1], errors="coerce").dropna().tolist()
            )
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
            "cov_percentile": coverage_percentile,
            "cov_threshold":  count_threshold,
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
    }

    html = _render_html(payload)
    report_path = outputfile + ".web_summary.html"
    with open(report_path, "w") as fh:
        fh.write(html)
    return report_path


# ─────────────────────────────────────────────────────────────────────────────
# HTML template
# ─────────────────────────────────────────────────────────────────────────────

def _render_html(data: dict) -> str:
    payload_json = json.dumps(data, indent=None, separators=(",", ":"))
    s = data["summary"]
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
        color:#fff;padding:22px 36px 18px;display:flex;align-items:center;
        justify-content:space-between;box-shadow:0 2px 12px rgba(0,0,0,.35)}}
.logo{{display:flex;align-items:center;gap:14px}}
.logo-icon{{width:42px;height:42px;background:linear-gradient(135deg,#00b4d8,#48cae4);
            border-radius:10px;display:flex;align-items:center;justify-content:center;
            font-size:22px;font-weight:900;color:#0f3460;flex-shrink:0}}
.logo-text h1{{font-size:1.45rem;font-weight:700;letter-spacing:.4px}}
.logo-text p{{font-size:.75rem;opacity:.7;margin-top:2px}}
.run-meta{{text-align:right;font-size:.72rem;opacity:.65;line-height:1.7}}

/* ── main container ──────────────────────────────────────────────────────── */
main{{max-width:1200px;margin:28px auto;padding:0 24px}}

/* ── section title ───────────────────────────────────────────────────────── */
.section-title{{font-size:.7rem;font-weight:700;letter-spacing:1.2px;
                text-transform:uppercase;color:#5a7194;margin-bottom:12px;
                padding-bottom:6px;border-bottom:2px solid #dce3ee}}

/* ── stat cards ──────────────────────────────────────────────────────────── */
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));
        gap:14px;margin-bottom:28px}}
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
</style>
</head>
<body>

<header>
  <div class="logo">
    <div class="logo-icon">G</div>
    <div class="logo-text">
      <h1>GeneExt</h1>
      <p>3&prime; gene-annotation extension pipeline</p>
    </div>
  </div>
  <div class="run-meta">
    <div>Output: <strong>{s['output_file']}</strong></div>
    <div>Input:  {s['input_file']}</div>
    <div>{s['run_date']}</div>
  </div>
</header>

<main>

  <!-- ── Summary cards ─────────────────────────────────────────────────── -->
  <p class="section-title">Run summary</p>
  <div class="cards" id="cards"></div>

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

// ── Stat cards ─────────────────────────────────────────────────────────────
const CARDS = [
  {{value: S.n_extended,   suffix: '/' + S.n_genes, label: 'Genes extended',           cls: 'accent-teal'}},
  {{value: S.pct_extended, suffix: '%',              label: 'Extension rate',            cls: 'accent-green'}},
  {{value: S.median_ext,   suffix: 'bp',             label: 'Median extension',          cls: 'accent-orange'}},
  {{value: S.mean_ext,     suffix: 'bp',             label: 'Mean extension',            cls: 'accent-blue'}},
  {{value: S.max_ext,      suffix: 'bp',             label: 'Max extension',             cls: 'accent-purple'}},
  {{value: S.n_genic_peaks,suffix: '',               label: 'Genic peaks',               cls: 'accent-teal'}},
  {{value: S.n_noov_peaks, suffix: '',               label: 'Non-overlapping peaks',     cls: 'accent-blue'}},
  {{value: S.cov_percentile || '—', suffix: S.cov_percentile ? 'th pct' : '',
    label: 'Coverage percentile', cls: 'accent-red'}},
];

(function buildCards() {{
  const wrap = document.getElementById('cards');
  CARDS.forEach(c => {{
    wrap.innerHTML += `
      <div class="card ${{c.cls}}">
        <div class="card-value">${{c.value}}<span>${{c.suffix}}</span></div>
        <div class="card-label">${{c.label}}</div>
      </div>`;
  }});
}})();

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
