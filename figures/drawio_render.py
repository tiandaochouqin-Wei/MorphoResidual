#!/usr/bin/env python3
"""Render a (self-authored) .drawio schematic to SVG/PDF/PNG via matplotlib, matching
the drawio geometry. Used for QA + paper output of Fig1. Handles rect / rounded rect /
ellipse / text / straight-arrow cells."""
import re, sys, html
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Ellipse, FancyArrowPatch
import mrstyle as S

SRC = sys.argv[1] if len(sys.argv) > 1 else "Fig1_concept.drawio"
OUT = sys.argv[2] if len(sys.argv) > 2 else "Fig1_concept"

root = ET.parse(SRC).getroot()
model = root.find(".//mxGraphModel")
W = float(model.get("pageWidth", 1080)); H = float(model.get("pageHeight", 380))


def sty(s):
    d = {}
    for tok in (s or "").split(";"):
        if "=" in tok:
            k, v = tok.split("=", 1); d[k] = v
        elif tok:
            d[tok] = True
    return d


def txt(v):
    if not v:
        return ""
    v = v.replace("<br>", "\n").replace("<br/>", "\n")
    v = re.sub(r"</?[bi]>", "", v)
    v = html.unescape(v)
    return v.strip()


fig, ax = plt.subplots(figsize=(W / 100, H / 100))
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off"); ax.invert_yaxis()

cells = model.findall(".//mxCell")
# vertices
for c in cells:
    if c.get("vertex") != "1":
        continue
    g = c.find("mxGeometry")
    x, y = float(g.get("x", 0)), float(g.get("y", 0))
    w, h = float(g.get("width", 0)), float(g.get("height", 0))
    st = sty(c.get("style", "")); label = txt(c.get("value", ""))
    fc = "#" + st["fillColor"].lstrip("#") if st.get("fillColor") and st["fillColor"] != "none" else "none"
    ec = "#" + st["strokeColor"].lstrip("#") if st.get("strokeColor") and st["strokeColor"] != "none" else "none"
    lw = float(st.get("strokeWidth", 1)) * 0.8
    fs = float(st.get("fontSize", 10)) * 0.92
    fcol = "#" + st.get("fontColor", "333333").lstrip("#")
    bold = st.get("fontStyle") in ("1", "3")
    italic = st.get("fontStyle") in ("2", "3")
    cx, cy = x + w / 2, y + h / 2
    if "ellipse" in st:
        ax.add_patch(Ellipse((cx, cy), w, h, facecolor=fc, edgecolor=ec, lw=lw, zorder=3))
    elif "text" in st:
        pass  # text-only cell, drawn below
    elif st.get("rounded") == "1":
        pad = min(w, h) * 0.16
        ax.add_patch(FancyBboxPatch((x + pad, y + pad), w - 2 * pad, h - 2 * pad,
                     boxstyle=f"round,pad={pad},rounding_size={pad*0.9}",
                     facecolor=fc, edgecolor=ec, lw=lw, zorder=2, mutation_aspect=1))
    else:
        ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec, lw=lw, zorder=2))
    if label:
        al = st.get("align", "center")
        ha = {"left": "left", "right": "right", "center": "center"}.get(al, "center")
        tx = x + 2 if ha == "left" else (x + w - 2 if ha == "right" else cx)
        va = st.get("verticalAlign", "middle")
        ty = y + 2 if va == "top" else cy
        vae = "top" if va == "top" else "center"
        ax.text(tx, ty, label, ha=ha, va=vae, fontsize=fs, color=fcol,
                fontweight="bold" if bold else "normal",
                fontstyle="italic" if italic else "normal", zorder=6,
                linespacing=1.15)
# edges (straight arrows)
for c in cells:
    if c.get("edge") != "1":
        continue
    g = c.find("mxGeometry")
    sp = g.find("mxPoint[@as='sourcePoint']"); tp = g.find("mxPoint[@as='targetPoint']")
    if sp is None or tp is None:
        continue
    st = sty(c.get("style", ""))
    ec = "#" + st.get("strokeColor", "666666").lstrip("#")
    lw = float(st.get("strokeWidth", 1)) * 0.9
    x1, y1 = float(sp.get("x")), float(sp.get("y"))
    x2, y2 = float(tp.get("x")), float(tp.get("y"))
    if st.get("endArrow") == "none":
        ax.plot([x1, x2], [y1, y2], color=ec, lw=lw, zorder=1)
    else:
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                     mutation_scale=9, lw=lw, color=ec, shrinkA=0, shrinkB=0, zorder=4))

fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
S.save_pub(fig, OUT)
print(f"rendered {SRC} -> {OUT}.{{svg,pdf,png,tiff}}")
