"""HTML report generation for raccoon commands."""
from __future__ import annotations

import base64
import csv
import os
import platform
import sys
from datetime import datetime
from typing import Iterable, Optional, Dict, Any, List

from Bio import SeqIO
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot

from .reconstruction_functions import load_tree, ensure_node_label


def _svg_data_uri(path: str) -> str:
    try:
        with open(path, "rb") as handle:
            data = handle.read()
        encoded = base64.b64encode(data).decode("ascii")
        return f"data:image/svg+xml;base64,{encoded}"
    except Exception:
        return ""


def _logo_html(data_uri: str, css_class: str, alt_text: str) -> str:
    if data_uri:
        return f'<img class="{css_class}" src="{data_uri}" alt="{alt_text}" />'
    return f'<div class="{css_class}" aria-label="{alt_text}"></div>'


def _write_html(outfile: str, title: str, summary_html: str, plots_html: List[str]) -> None:
    plots = "\n".join(plots_html)
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    raccoon_logo = _svg_data_uri(os.path.join(base_dir, "docs", "raccoon_logo.svg"))
    artic_logo = _svg_data_uri(os.path.join(base_dir, "docs", "artic-logo-small.svg"))
    raccoon_logo_html = _logo_html(raccoon_logo, "logo", "Raccoon logo")
    artic_logo_html = _logo_html(artic_logo, "logo-small", "ARTIC Network logo")
    generated_stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    html = f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <title>{title}</title>
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
    <link rel=\"stylesheet\" href=\"https://cdn.datatables.net/1.13.8/css/jquery.dataTables.min.css\">
    <style>
        body {{ font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; font-weight: 300; margin: 20px; }}
    h1 {{ margin-bottom: 0.2rem; }}
    .meta {{ color: #666; margin-bottom: 1rem; }}
    .section {{ margin: 1.5rem 0; }}
        .header {{ display: flex; align-items: center; justify-content: space-between; gap: 12px; }}
        .header-left {{ display: flex; align-items: center; gap: 12px; }}
        .logo {{ height: 48px; width: auto; }}
        .logo-small {{ height: 48px; width: auto; }}
    .card {{ border: 1px solid #ddd; border-radius: 8px; padding: 12px 16px; margin: 10px 0; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 12px; }}
        h1, h2, h3, h4 {{ font-weight: 600; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #eee; padding: 6px 8px; text-align: left; }}
    th {{ background: #e6e1ee; }}
        tbody tr.filtered {{ background: #f2d6d6 !important; }}
        details {{ border: 1px solid #eee; border-radius: 6px; padding: 8px 10px; margin: 6px 0; }}
        summary {{ cursor: pointer; font-weight: 600; }}
        .toc a {{ text-decoration: none; color: #2f5f6b; }}
        .footer {{ margin-top: 2rem; padding-top: 1rem; border-top: 1px solid #eee; color: #666; font-size: 0.9rem; display: flex; justify-content: space-between; align-items: center; }}
        .metadata-block {{ margin: 10px 0; }}
        .metadata-title {{ font-weight: 600; margin: 6px 0; }}
        .banner {{ background: #f6f3f9; border: 1px solid #e2dbea; padding: 10px 12px; border-radius: 8px; }}
        .caption {{ color: #666; font-size: 0.9rem; margin-top: 6px; }}
        .table-wrap {{ overflow-x: auto; }}
        @media print {{
            .dataTables_filter, .dataTables_length, .dataTables_info, .dataTables_paginate {{ display: none !important; }}
            details summary {{ list-style: none; }}
            details {{ border: none; padding: 0; }}
        }}
  </style>
    <script src=\"https://code.jquery.com/jquery-3.7.1.min.js\"></script>
    <script src=\"https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js\"></script>
</head>
<body>
    <div class=\"header\">
        <div class=\"header-left\">
            {raccoon_logo_html}
            <div>
                <h1>{title}</h1>
                <div class=\"meta\">Generated {generated_stamp}</div>
            </div>
        </div>
        {artic_logo_html}
    </div>
    <div class=\"section\">{summary_html}</div>
    <div class=\"section\">{plots}</div>
    <div class=\"footer\">
        <div><strong>Raccoon</strong> — Rigorous Alignment Curation: Cleanup Of Outliers and Noise</div>
        <div>Author: Áine O'Toole</div>
    </div>
        <script>
            $(document).ready(function () {{
                $('table.datatable').DataTable({{
                    paging: false,
                    info: false,
                    autoWidth: false
                }});
            }});
        </script>
</body>
</html>"""
    with open(outfile, "w") as handle:
        handle.write(html)


def _n_content(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    return seq.count("N") / len(seq)


def _apply_plot_style(fig: go.Figure) -> None:
    fig.update_layout(
        font=dict(family="Helvetica Neue", size=12, color="#111"),
        colorway=["#4BA3A8", "#7A6BB1", "#D08BA8"],
        plot_bgcolor="white",
        paper_bgcolor="white",
        showlegend=True,
        hoverlabel=dict(font=dict(color="#ffffff"), bordercolor="#ffffff"),
    )
    fig.update_xaxes(showline=True, linecolor="black", linewidth=1, gridcolor="rgba(0,0,0,0.05)")
    fig.update_yaxes(showline=True, linecolor="black", linewidth=1, gridcolor="rgba(0,0,0,0.05)")


def _plot_div(fig: go.Figure) -> str:
    return plot(fig, include_plotlyjs="cdn", output_type="div", config={"displayModeBar": False})


def _safe_mean(values: Iterable[float]) -> float:
    values = list(values)
    return sum(values) / len(values) if values else 0.0


def _safe_min(values: Iterable[int]) -> int:
    values = list(values)
    return min(values) if values else 0


def _safe_max(values: Iterable[int]) -> int:
    values = list(values)
    return max(values) if values else 0


def generate_combine_report(
    outdir: str,
    output_fasta: str,
    input_fastas: Iterable[str],
    metadata_paths: Optional[Iterable[str]] = None,
    metadata_id_field: str = "id",
    metadata_location_field: str = "location",
    metadata_date_field: str = "date",
    header_separator: str = "|",
    min_length: Optional[int] = None,
    max_n_content: Optional[float] = None,
) -> str:
    records_summary = []
    ids_by_file = {}
    seq_details_by_file: Dict[str, List[Dict[str, Any]]] = {}
    filtered_rows: List[Dict[str, Any]] = []
    for path in input_fastas:
        lengths = []
        n_contents = []
        ids = []
        seq_details = []
        for rec in SeqIO.parse(path, "fasta"):
            seq = str(rec.seq)
            seq_len = len(seq)
            n_content = _n_content(seq)
            lengths.append(seq_len)
            n_contents.append(n_content)
            ids.append(rec.id)
            seq_details.append({
                "id": rec.id,
                "length": seq_len,
                "n_content": n_content,
            })
            reasons = []
            if min_length is not None and seq_len < min_length:
                reasons.append(f"length < {min_length}")
            if max_n_content is not None and n_content > max_n_content:
                reasons.append(f"N content > {max_n_content}")
            if reasons:
                filtered_rows.append({
                    "file": os.path.basename(path),
                    "id": rec.id,
                    "length": seq_len,
                    "n_content": round(n_content, 4),
                    "reason": "; ".join(reasons),
                })
        records_summary.append({
            "file": os.path.basename(path),
            "sequences": len(lengths),
            "len_min": _safe_min(lengths),
            "len_max": _safe_max(lengths),
            "len_mean": round(_safe_mean(lengths), 2),
            "n_mean": round(_safe_mean(n_contents), 4),
        })
        ids_by_file[os.path.basename(path)] = ids
        seq_details_by_file[os.path.basename(path)] = seq_details

    metadata_summary = "No metadata used."
    metadata_tables_html = ""
    metadata_locations = pd.Series(dtype=str)
    metadata_dates = pd.Series(dtype="datetime64[ns]")
    if metadata_paths:
        try:
            frames = []
            for path in metadata_paths:
                frame = pd.read_csv(path)
                frames.append(frame)
                cols = list(frame.columns)
                header_html = "".join([f"<th>{col}</th>" for col in cols])
                rows_html = ""
                for _, row in frame.iterrows():
                    row_cells = "".join([f"<td>{row.get(col, '')}</td>" for col in cols])
                    rows_html += f"<tr>{row_cells}</tr>"
                metadata_tables_html += f"""
                    <div class=\"metadata-block\">
                        <div class=\"metadata-title\">{os.path.basename(path)} ({len(frame)} rows)</div>
                        <table>
                            <thead><tr>{header_html}</tr></thead>
                        </table>
                        <details>
                            <summary>Show rows</summary>
                            <div class=\"table-wrap\"><table class=\"datatable\">
                                <thead><tr>{header_html}</tr></thead>
                                <tbody>{rows_html}</tbody>
                            </table></div>
                        </details>
                    </div>
                """
            meta = pd.concat(frames, ignore_index=True)
            metadata_locations = meta.get(metadata_location_field, pd.Series(dtype=str)).dropna().astype(str)
            metadata_dates = pd.to_datetime(meta.get(metadata_date_field, pd.Series(dtype=str)), errors="coerce")
            date_min = metadata_dates.min()
            date_max = metadata_dates.max()
            date_range = "n/a"
            if pd.notna(date_min) and pd.notna(date_max):
                date_range = f"{date_min.date().isoformat()} → {date_max.date().isoformat()}"
            metadata_summary = (
                f"Metadata files: {', '.join([os.path.basename(p) for p in metadata_paths])}. "
                f"ID field: {metadata_id_field}. "
                f"Location field: {metadata_location_field} (unique: {metadata_locations.nunique()}). "
                f"Date field: {metadata_date_field} (range: {date_range})."
            )
        except Exception:
            metadata_summary = "Metadata provided, but summary could not be parsed."

    header_stats = {"locations": 0, "dates": ""}
    try:
        if metadata_paths and not metadata_locations.empty:
            header_stats["locations"] = int(metadata_locations.nunique())
            if not metadata_dates.empty:
                date_min = metadata_dates.min()
                date_max = metadata_dates.max()
                if pd.notna(date_min) and pd.notna(date_max):
                    header_stats["dates"] = f"{date_min.date().isoformat()} → {date_max.date().isoformat()}"
        elif output_fasta:
            locs = set()
            dates = []
            for rec in SeqIO.parse(output_fasta, "fasta"):
                parts = rec.id.split(header_separator)
                if len(parts) >= 2 and parts[1]:
                    locs.add(parts[1])
                if len(parts) >= 3:
                    try:
                        parsed = pd.to_datetime(parts[2])
                        if pd.notna(parsed):
                            dates.append(parsed)
                    except Exception:
                        pass
            header_stats["locations"] = len(locs)
            if dates:
                dmin = min(dates)
                dmax = max(dates)
                header_stats["dates"] = f"{dmin.date().isoformat()} → {dmax.date().isoformat()}"
    except Exception:
        pass

    dataset_plot_html = ""
    date_location_records = []
    if metadata_paths:
        try:
            frames = []
            for path in metadata_paths:
                frames.append(pd.read_csv(path))
            meta = pd.concat(frames, ignore_index=True)
            if metadata_date_field in meta.columns and metadata_location_field in meta.columns:
                dates = pd.to_datetime(meta[metadata_date_field], errors="coerce")
                locations = meta[metadata_location_field].astype(str)
                ids = meta.get(metadata_id_field, pd.Series(dtype=str)).astype(str)
                for date_value, location, seq_id in zip(dates, locations, ids):
                    if pd.notna(date_value) and pd.notna(location) and location:
                        date_location_records.append({
                            "date": date_value,
                            "location": location,
                            "id": seq_id,
                        })
        except Exception:
            date_location_records = []
    elif output_fasta:
        try:
            for rec in SeqIO.parse(output_fasta, "fasta"):
                parts = rec.id.split(header_separator)
                if len(parts) >= 3:
                    date_value = pd.to_datetime(parts[2], errors="coerce")
                    location = parts[1]
                    if pd.notna(date_value) and location:
                        date_location_records.append({
                            "date": date_value,
                            "location": location,
                            "id": parts[0],
                        })
        except Exception:
            date_location_records = []

    if date_location_records:
        df = pd.DataFrame(date_location_records)
        fig = go.Figure(data=[
            go.Scatter(
                x=df["date"],
                y=df["location"],
                mode="markers",
                marker=dict(size=8, opacity=0.7),
                customdata=df.get("id"),
                hovertemplate="Date: %{x}<br>Location: %{y}<br>ID: %{customdata}<extra></extra>",
            )
        ])
        fig.update_layout(
            title="Sampling dates by location",
            xaxis_title="Date",
            yaxis_title="Location",
        )
        _apply_plot_style(fig)
        dataset_plot_html = _plot_div(fig)

    filter_summary = []
    if min_length is not None:
        filter_summary.append(f"Minimum length: {min_length}")
    if max_n_content is not None:
        filter_summary.append(f"Maximum N content: {max_n_content}")
    filter_summary_text = ", ".join(filter_summary) if filter_summary else "No filters applied."

    filtered_lookup = {(row["file"], row["id"]) for row in filtered_rows}

    total_sequences = sum(row["sequences"] for row in records_summary)
    filtered_count = len(filtered_rows)

    length_plot_html = ""
    length_groups: Dict[int, List[str]] = {}
    for path in input_fastas:
        for rec in SeqIO.parse(path, "fasta"):
            length_groups.setdefault(len(rec.seq), []).append(rec.id)
    if length_groups:
        lengths_sorted = sorted(length_groups.keys())
        counts = [len(length_groups[l]) for l in lengths_sorted]
        ids_by_length = [", ".join(length_groups[l]) for l in lengths_sorted]
        fig = go.Figure(data=[
            go.Bar(
                x=lengths_sorted,
                y=counts,
                marker=dict(color="#4BA3A8"),
                customdata=ids_by_length,
                hovertemplate="Length: %{x}<br>Count: %{y}<br>IDs: %{customdata}<extra></extra>",
            )
        ])
        fig.update_layout(
            xaxis_title="Sequence length (bp)",
            yaxis_title="Count",
            showlegend=False,
        )
        if min_length is not None:
            fig.add_vline(
                x=min_length,
                line_width=2,
                line_dash="dash",
                line_color="#c77c8a",
            )
            if counts:
                fig.add_trace(go.Scatter(
                    x=[min_length],
                    y=[max(counts)],
                    mode="markers",
                    marker=dict(opacity=0),
                    hovertemplate="min length: %{x}<extra></extra>",
                    showlegend=False,
                ))
        _apply_plot_style(fig)
        length_plot_html = _plot_div(fig)

    provenance_inputs = ", ".join([os.path.basename(p) for p in input_fastas])
    provenance_metadata = ", ".join([os.path.basename(p) for p in metadata_paths]) if metadata_paths else "None"
    provenance_output = os.path.basename(output_fasta) if output_fasta else "stdout"
    generated_stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    try:
        from raccoon import __version__ as raccoon_version
    except Exception:
        raccoon_version = "unknown"

    final_rows = []
    if output_fasta:
        try:
            for rec in SeqIO.parse(output_fasta, "fasta"):
                seq = str(rec.seq)
                final_rows.append({
                    "id": rec.id,
                    "length": len(seq),
                    "n_content": round(_n_content(seq), 4),
                })
        except Exception:
            final_rows = []

    cmd_parts = ["raccoon", "combine"]
    cmd_parts.extend([os.path.basename(p) for p in input_fastas])
    if output_fasta:
        cmd_parts.extend(["--output", os.path.basename(output_fasta)])
    if metadata_paths:
        cmd_parts.append("--metadata")
        cmd_parts.extend([os.path.basename(p) for p in metadata_paths])
    if min_length is not None:
        cmd_parts.extend(["--min-length", str(min_length)])
    if max_n_content is not None:
        cmd_parts.extend(["--max-n-content", str(max_n_content)])
    cmd_line = " ".join(cmd_parts)

    summary_html = f"""
        <div class=\"grid\">
            <div class=\"card\">
                <h3>Executive summary</h3>
                <p>Total sequences: {total_sequences}</p>
                <p>Filtered: {filtered_count}</p>
                <p>Unique locations: {header_stats.get("locations", 0)}</p>
                <p>Date range: {header_stats.get("dates", "")}</p>
            </div>
            <div class=\"card\">
                <h3>Filters</h3>
                <p>Command: <code>{cmd_line}</code></p>
                <p>Filters: {filter_summary_text}</p>
            </div>
        </div>
        <div class=\"toc card\">
            <h3>Table of contents</h3>
            <ol>
            <li><a href=\"#inputs\">Inputs & sequences</a></li>
            <li><a href=\"#lengths\">Sequence length description</a></li>
            <li><a href=\"#metadata\">Metadata</a></li>
            <li><a href=\"#dataset\">Dataset description</a></li>
                <li><a href=\"#final\">Final IDs</a></li>
                <li><a href=\"#report-metadata\">Report metadata</a></li>
            </ol>
        </div>
        <div id=\"inputs\" class=\"card\">
            <h3>1. Inputs & sequences</h3>
            <div class=\"table-wrap\"><table class=\"datatable\">
                <thead><tr><th>File</th><th>Seqs</th><th>Length (min)</th><th>Length (max)</th><th>Length (mean)</th><th>N mean</th></tr></thead>
                <tbody>
    """
    for row in records_summary:
        summary_html += (
            f"<tr><td>{row['file']}</td><td>{row['sequences']}</td>"
            f"<td>{row['len_min']}</td><td>{row['len_max']}</td>"
            f"<td>{row['len_mean']}</td><td>{row['n_mean']}</td></tr>"
        )
    summary_html += """
                </tbody>
            </table></div>
    """
    summary_html += """
            <details>
                <summary>Show sequence details</summary>
    """
    for fname, seq_details in seq_details_by_file.items():
        summary_html += f"<details><summary>{fname} ({len(seq_details)} sequences)</summary>"
        summary_html += "<div class=\"table-wrap\"><table class=\"datatable\"><thead><tr><th>ID</th><th>Length</th><th>N content</th><th>Status</th></tr></thead><tbody>"
        for detail in seq_details:
            filtered = (fname, detail["id"]) in filtered_lookup
            row_class = "filtered" if filtered else ""
            status = "filtered" if filtered else "kept"
            summary_html += (
                f"<tr class=\"{row_class}\"><td>{detail['id']}</td><td>{detail['length']}</td>"
                f"<td>{round(detail['n_content'], 4)}</td><td>{status}</td></tr>"
            )
        summary_html += "</tbody></table></div></details>"
    summary_html += f"""
            </details>
        </div>
        <div id=\"lengths\" class=\"card\">
            <h3>2. Sequence length description</h3>
            {length_plot_html if length_plot_html else "<p>No length data available for plotting.</p>"}
            <p class=\"caption\">Distribution of sequence lengths across all input FASTAs.</p>
        </div>
        <div id=\"metadata\" class=\"card\">
            <h3>3. Metadata</h3>
            <p>{metadata_summary}</p>
            {metadata_tables_html if metadata_tables_html else "<p>No metadata tables available.</p>"}
        </div>
        <div id=\"dataset\" class=\"card\">
            <h3>4. Dataset description</h3>
            <p>Unique locations: {header_stats.get("locations", 0)}</p>
            <p>Date range: {header_stats.get("dates", "")}</p>
            {dataset_plot_html if dataset_plot_html else "<p>No date/location data available for plotting.</p>"}
            <p class=\"caption\">Sampling timeline by location.</p>
        </div>
        <div id=\"final\" class=\"card\">
            <h3>5. Final dataset</h3>
            <div class=\"table-wrap\"><table class=\"datatable\">
                <thead><tr><th>ID</th><th>Length</th><th>N content</th></tr></thead>
                <tbody>
    """
    for row in final_rows:
        summary_html += f"<tr><td>{row['id']}</td><td>{row['length']}</td><td>{row['n_content']}</td></tr>"
    summary_html += f"""
                </tbody>
            </table></div>
        </div>
        <div class=\"card\">
            <h3>Datafiles </h3>
            <p>Inputs: {provenance_inputs}</p>
            <p>Metadata: {provenance_metadata}</p>
            <p>Output: {provenance_output}</p>
        </div>
        <div id=\"report-metadata\" class=\"card\">
            <h3>Report metadata</h3>
            <p>Generated: {generated_stamp}</p>
            <p>Raccoon version: {raccoon_version}</p>
            <p>Python: {sys.version.split()[0]}</p>
            <p>Platform: {platform.system()} {platform.release()}</p>
        </div>
    """

    plots_html = []

    outpath = os.path.join(outdir, "combine_report.html")
    _write_html(outpath, "Raccoon combine report", summary_html, plots_html)
    return outpath


def generate_alignment_report(outdir: str, alignment_path: str, mask_file: Optional[str] = None) -> str:
    lengths = []
    n_contents = []
    completeness = []
    for rec in SeqIO.parse(alignment_path, "fasta"):
        seq = str(rec.seq)
        lengths.append(len(seq))
        n_contents.append(_n_content(seq))
        if seq:
            valid = sum(1 for c in seq.upper() if c not in ["N", "-"])
            completeness.append(valid / len(seq))

    aln_len = _safe_max(lengths)
    summary_html = f"""
    <div class=\"grid\">
      <div class=\"card\">
        <h3>Alignment summary</h3>
        <p>Sequences: {len(lengths)}</p>
        <p>Alignment length: {aln_len}</p>
        <p>Mean N content: {round(_safe_mean(n_contents), 4)}</p>
        <p>Mean completeness: {round(_safe_mean(completeness), 4)}</p>
      </div>
    </div>
    """

    plots_html = []
    if n_contents:
        fig = go.Figure(data=[go.Histogram(x=n_contents, nbinsx=20)])
        fig.update_layout(title="N content distribution", xaxis_title="N content", yaxis_title="Count")
        _apply_plot_style(fig)
        plots_html.append(_plot_div(fig))

    if mask_file and os.path.exists(mask_file):
        site_rows = []
        with open(mask_file, "r") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                notes = row.get("note", "").split(";")
                for note in notes:
                    site_rows.append({"site": int(row["Name"]), "note": note})
        if site_rows:
            df = pd.DataFrame(site_rows)
            fig = go.Figure()
            def _pretty_label(value: str) -> str:
                return value.replace("_", " ").title()

            for note, subset in df.groupby("note"):
                pretty = _pretty_label(note)
                fig.add_trace(go.Scatter(
                    x=subset["site"],
                    y=[pretty] * len(subset),
                    mode="markers",
                    name=pretty,
                    marker=dict(size=8),
                ))
            fig.update_layout(title="Flagged sites by category", xaxis_title="Site", yaxis_title="Category")
            fig.update_xaxes(range=[0, aln_len])
            _apply_plot_style(fig)
            plots_html.append(_plot_div(fig))

    outpath = os.path.join(outdir, "alignment_report.html")
    _write_html(outpath, "Raccoon alignment report", summary_html, plots_html)
    return outpath


def generate_phylo_report(outdir: str, treefile: str, flags_csv: Optional[str] = None, tree_format: str = "auto") -> str:
    my_tree = load_tree(treefile, tree_format=tree_format)
    tip_names = []
    tip_heights = []
    for node in my_tree.Objects:
        if node.branchType == 'leaf':
            label = ensure_node_label(node)
            if label:
                tip_names.append(label)
                tip_heights.append(node.height)

    summary_html = f"""
    <div class=\"grid\">
      <div class=\"card\">
        <h3>Tree summary</h3>
        <p>Tips: {len(tip_names)}</p>
        <p>Tree height: {getattr(my_tree, 'treeHeight', 'n/a')}</p>
        <p>Y span: {getattr(my_tree, 'ySpan', 'n/a')}</p>
      </div>
    </div>
    """

    plots_html = []
    if tip_heights:
        fig = go.Figure(data=[go.Histogram(x=tip_heights, nbinsx=20)])
        fig.update_layout(title="Root-to-tip distances", xaxis_title="Distance", yaxis_title="Count")
        _apply_plot_style(fig)
        plots_html.append(_plot_div(fig))

    if flags_csv and os.path.exists(flags_csv):
        flags = pd.read_csv(flags_csv)
        if not flags.empty and "mutation_type" in flags:
            counts = flags["mutation_type"].value_counts().reset_index()
            counts.columns = ["mutation_type", "count"]
            fig = go.Figure(data=[go.Bar(x=counts["mutation_type"], y=counts["count"])])
            fig.update_layout(title="Flagged mutation types", xaxis_title="Type", yaxis_title="Count")
            _apply_plot_style(fig)
            plots_html.append(_plot_div(fig))

    outpath = os.path.join(outdir, "phylo_report.html")
    _write_html(outpath, "Raccoon phylo report", summary_html, plots_html)
    return outpath
