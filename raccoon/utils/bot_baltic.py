"""Plotly-based interactive tree rendering utilities."""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import plotly.graph_objects as go
import plotly.io as pio

from .reconstruction_functions import load_tree, ensure_node_label


def _scale_values(values: List[Optional[float]], factor: float) -> List[Optional[float]]:
    return [None if v is None else v * factor for v in values]


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


def _plot_div(fig: go.Figure, div_id: Optional[str] = None) -> str:
    return pio.to_html(
        fig,
        include_plotlyjs="cdn",
        full_html=False,
        config={"displayModeBar": False},
        auto_play=False,
        div_id=div_id,
    )


def _format_traits(traits: Dict[str, Any]) -> str:
    lines = []
    for key in sorted(traits):
        value = traits.get(key)
        if value is None or value == "":
            continue
        lines.append(f"{key}: {value}")
    return "<br>".join(lines)


def _parse_label_traits(label: str) -> Dict[str, Any]:
    traits: Dict[str, Any] = {}
    if not label:
        return traits
    parts = [p.strip() for p in label.split("|")]
    if len(parts) >= 2 and parts[1]:
        traits.setdefault("location", parts[1])
    if len(parts) >= 3 and parts[2]:
        traits.setdefault("date", parts[2])
    return traits


def build_tree_plot(treefile: str, tree_format: str = "auto") -> str:
    try:
        tree = load_tree(treefile, tree_format=tree_format)
    except Exception:
        return ""

    try:
        tree.drawTree()
    except Exception:
        return ""

    line_x: List[Optional[float]] = []
    line_y: List[Optional[float]] = []

    try:
        internal_nodes = tree.getInternal()
    except Exception:
        internal_nodes = []

    for node in internal_nodes:
        if not getattr(node, "children", None):
            continue
        child_ys = []
        for child in node.children:
            if child.x is None or child.y is None or node.x is None:
                continue
            line_x.extend([node.x, child.x, None])
            line_y.extend([child.y, child.y, None])
            child_ys.append(child.y)
        if child_ys and node.x is not None:
            line_x.extend([node.x, node.x, None])
            line_y.extend([min(child_ys), max(child_ys), None])

    tip_x: List[float] = []
    tip_y: List[float] = []
    tip_hover: List[str] = []
    tip_traits: List[Dict[str, Any]] = []

    node_x: List[float] = []
    node_y: List[float] = []
    node_hover: List[str] = []

    for obj in getattr(tree, "Objects", []):
        label = ensure_node_label(obj) or getattr(obj, "name", "") or ""
        traits = dict(getattr(obj, "traits", {}) or {})
        traits.setdefault("label", label)
        traits.setdefault("branch_type", getattr(obj, "branchType", ""))
        traits.update(_parse_label_traits(label))

        if getattr(obj, "is_leaflike", lambda: False)():
            if obj.x is None or obj.y is None:
                continue
            tip_x.append(obj.x)
            tip_y.append(obj.y)
            tip_traits.append(traits)
            tip_hover.append(_format_traits(traits))
        elif getattr(obj, "is_node", lambda: False)():
            if obj.x is None or obj.y is None:
                continue
            node_x.append(obj.x)
            node_y.append(obj.y)
            node_hover.append(_format_traits(traits))

    color_keys = set()
    for traits in tip_traits:
        for key, value in traits.items():
            if value is not None and value != "":
                color_keys.add(key)
    preferred = ["location", "date", "branch_type"]
    color_keys = preferred + sorted([k for k in color_keys if k not in preferred])
    default_key = color_keys[0] if color_keys else "branch_type"

    palette = [
        "#4BA3A8",
        "#7A6BB1",
        "#D08BA8",
        "#5A9BD4",
        "#B9A44C",
        "#8BC59A",
        "#E26D5C",
        "#6C8EA5",
    ]

    color_maps: Dict[str, Dict[str, str]] = {}
    for key in color_keys:
        values = sorted({str(t.get(key, "")) for t in tip_traits if t.get(key, "") != ""})
        mapping = {}
        for idx, value in enumerate(values):
            mapping[value] = palette[idx % len(palette)]
        color_maps[key] = mapping

    def _colors_for(key: str) -> List[str]:
        mapping = color_maps.get(key, {})
        return [mapping.get(str(t.get(key, "")), "#4d4d4d") for t in tip_traits]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=line_x,
        y=line_y,
        mode="lines",
        line=dict(color="#222", width=1),
        hoverinfo="skip",
        showlegend=False,
        name="Branches",
    ))
    fig.add_trace(go.Scatter(
        x=tip_x,
        y=tip_y,
        mode="markers",
        marker=dict(size=6, color=_colors_for(default_key)),
        hovertemplate="%{text}<extra></extra>",
        text=tip_hover,
        name="Tips",
    ))
    fig.add_trace(go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers",
        marker=dict(size=4, color="#999999"),
        hovertemplate="%{text}<extra></extra>",
        text=node_hover,
        name="Nodes",
    ))

    scale_values = [round(1.0 + (i * 0.01), 2) for i in range(101)]
    y_min = min(tip_y) if tip_y else 0
    y_max = max(tip_y) if tip_y else 0
    base_height = max(500, int((y_max if tip_y else 0) * 15)) if tip_y else 500

    slider_steps = []
    for scale in scale_values:
        slider_steps.append({
            "label": f"{scale:.2f}x",
            "method": "update",
            "args": [
                {
                    "y": [
                        _scale_values(line_y, scale),
                        _scale_values(tip_y, scale),
                        _scale_values(node_y, scale),
                    ]
                },
                {
                    "yaxis": {"range": [y_min * scale, y_max * scale]},
                    "height": int(base_height * scale),
                },
            ],
        })

    color_buttons = []
    for key in color_keys:
        color_buttons.append({
            "label": key,
            "method": "restyle",
            "args": [{"marker.color": [_colors_for(key)]}, [1]],
        })

    fig.update_layout(
        xaxis_title="Branch length",
        yaxis_title="Tips",
        showlegend=False,
        height=base_height,
        yaxis=dict(range=[y_min, y_max] if tip_y else None),
        sliders=[{
            "active": 0,
            "pad": {"t": 30},
            "currentvalue": {"prefix": "Y scale: "},
            "steps": slider_steps,
        }],
        updatemenus=[{
            "buttons": color_buttons,
            "direction": "down",
            "showactive": True,
            "x": 0.02,
            "y": 1.12,
            "xanchor": "left",
            "yanchor": "top",
        }],
        margin=dict(l=60, r=20, t=40, b=40),
    )
    _apply_plot_style(fig)
    return _plot_div(fig)
