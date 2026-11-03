#!/usr/bin/env python3
"""Export Figure 1, the study flow chart, to PNG and EPS.

    python3 -I tools/export_fig1.py [--svg PATH] [--outdir DIR]

Figure 1 is the one figure R does not draw. The other four come out of
07_figures.R; the flow chart is written by make_fig1.py beside this script, as
live-text SVG, so its counts stay text rather than pixels and can be read and
checked without opening an image editor.

That left it reachable only from inside the Word figures document, so a reader
who cloned this repository got four figures of five. The SVG is now kept here
beside the others, and this script renders it to the same 174 mm width and
300 dpi as the rest, plus the EPS that Springer asks for vector artwork.

Requires cairosvg, and pdftops (poppler-utils) for the EPS.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile

W_MM = 174.0          # the journal's full page width; every figure matches it
DPI = 300
# cairosvg sizes a PDF in CSS pixels at 96 dpi, not in points. Passing points
# produced a 130.5 mm EPS instead of 174 mm - short by exactly 72/96 - and with
# it type of 4.9 pt against the journal's 8 pt floor.
PX_PER_MM = 96 / 25.4


def main():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    default = os.path.join(here, "expected_2004_2025", "figures", "fig1_study_flow.svg")
    ap = argparse.ArgumentParser()
    ap.add_argument("--svg", default=default)
    ap.add_argument("--outdir")
    a = ap.parse_args()

    if not os.path.exists(a.svg):
        sys.exit("no such SVG: %s" % a.svg)
    outdir = a.outdir or os.path.dirname(a.svg)
    stem = os.path.join(outdir, os.path.splitext(os.path.basename(a.svg))[0])

    try:
        import cairosvg
    except ImportError:
        sys.exit("cairosvg is needed: pip install cairosvg")

    svg = open(a.svg, encoding="utf-8").read()
    vb = re.search(r'viewBox="\s*[\d.-]+\s+[\d.-]+\s+([\d.]+)\s+([\d.]+)', svg)
    if not vb:
        sys.exit("the SVG has no viewBox, so its aspect ratio is unknown")
    vw, vh = float(vb.group(1)), float(vb.group(2))

    w_px = round(W_MM / 25.4 * DPI)
    cairosvg.svg2png(bytestring=svg.encode("utf-8"), write_to=stem + ".png",
                     output_width=w_px, output_height=round(w_px * vh / vw),
                     background_color="white")

    # Pillow is used only to stamp the resolution, which cairosvg does not write.
    # Without it the PNG is the right number of pixels but claims no physical
    # size, and a checker reading it back cannot tell 174 mm from any other.
    try:
        from PIL import Image
        Image.open(stem + ".png").save(stem + ".png", dpi=(DPI, DPI))
    except ImportError:
        print("Pillow not installed: the PNG carries no dpi tag", file=sys.stderr)

    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
        pdf = tmp.name
    try:
        cairosvg.svg2pdf(bytestring=svg.encode("utf-8"), write_to=pdf,
                         output_width=W_MM * PX_PER_MM,
                         output_height=W_MM * PX_PER_MM * vh / vw)
        subprocess.run(["pdftops", "-eps", pdf, stem + ".eps"], check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        print("EPS not written (%s); the PNG is there" % exc, file=sys.stderr)
    finally:
        os.unlink(pdf)

    for ext in (".png", ".eps"):
        if os.path.exists(stem + ext):
            print("  %-34s %9d bytes" % (os.path.basename(stem + ext),
                                         os.path.getsize(stem + ext)))


if __name__ == "__main__":
    main()
