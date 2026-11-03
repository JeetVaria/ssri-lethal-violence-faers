#!/usr/bin/env python3
"""Draw Figure 1, the study flow chart, to SVG.

    python3 -I tools/make_fig1.py [--out PATH] [--check]

Figure 1 used to be a hand-positioned SVG: 21 rectangles and 44 text nodes at
absolute coordinates, in a 1900 x 1010 frame with every label at font-size 19.
At the journal's 174 mm full text width that works out to 4.9 pt type,
against a floor of 8 pt, and there was no way to fix it by editing numbers in
place - the labels have to grow by a third while the frame stays 174 mm, and
the old arrangement had no room.

Why the layout changed. At 8 pt the four outcome panels and the two database
boxes need 1,807 units of box width between them; the frame is 1,900, which
leaves 93 for two margins and three gaps. Flanking each spine with its outcome
panel, as the old drawing did, does not fit at any padding. So each outcome
panel keeps its own frame beside the spine it belongs to, tinted against the
spine's solid fill and joined to it by a dashed link - the same two-tone
pairing the old figure made, laid out across the frame rather than stacked,
with the boxes measured to the type rather than to fixed coordinates.

Nothing about the content changed: the same 46 labels, the same counts, the
same five-colour palette, "Post hoc" still italic.

Box widths and heights are measured from the text with real Arial metrics
(Liberation Sans is metric-compatible and is what a Linux box has), so a label
cannot silently overflow its box the way it would if the widths were guessed.
--check re-measures the result and reports the smallest type and the tightest
box in the finished drawing.
"""
import argparse
import math
import os
import sys

MM = 174.0                     # the journal's full text width
PT_PER_MM = 72.0 / 25.4
W = 1900                       # frame width in drawing units
K = MM * PT_PER_MM / W         # points per unit at final size
FLOOR_PT = 8.0                 # the journal's minimum type size
MAX_H_MM = 234.0               # the journal's maximum figure height

FS = 31                        # body type        -> 8.05 pt
FSH = 33                       # header type      -> 8.57 pt
LH = 40                        # body line height
LHH = 42                       # header line height
PAD_X = 24                     # horizontal padding inside a box
PAD_Y = 18                     # vertical padding inside a box

M = 60                         # outer margin
GUT = 90                       # gutter between the two arms
COLW = (W - 2 * M - GUT) // 2  # 845: the envelope each arm gets
CXL = M + COLW / 2             # left arm axis
CXR = W - M - COLW / 2         # right arm axis
CXM = W / 2                    # the spine above and below the split

INK = "#1A1A1A"
WHITE = "#FFFFFF"
GREEN, GREEN_L = "#4E8A3C", "#E8F3E2"
RED, RED_L = "#E8342A", "#FBE4DC"
BLUE, BLUE_L = "#2E6F95", "#8DBDD8"
BLUE_LL = "#DCEAF5"
FAMILY = "Arial, Helvetica, sans-serif"

FONTS = {
    (0, 0): "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    (1, 0): "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
    (0, 1): "/usr/share/fonts/truetype/liberation/LiberationSans-Italic.ttf",
}


def _metric(text, size, bold=False, italic=False):
    from PIL import ImageFont
    path = FONTS[(1 if bold else 0, 1 if italic else 0)]
    return ImageFont.truetype(path, 400).getlength(text) / 400.0 * size


# ---------------------------------------------------------------- text pieces

class Line:
    """One line of a box: plain text, or a run list for mixed styling."""

    def __init__(self, runs, size, bold=False):
        self.runs = [(runs, False)] if isinstance(runs, str) else runs
        self.size = size
        self.bold = bold

    @property
    def plain(self):
        return "".join(t for t, _ in self.runs)

    def width(self):
        return sum(_metric(t, self.size, self.bold, it) for t, it in self.runs)

    def svg(self, cx, y):
        if len(self.runs) == 1 and not self.runs[0][1]:
            body = esc(self.runs[0][0])
        else:
            body = "".join(
                ('<tspan font-style="italic">%s</tspan>' % esc(t)) if it else esc(t)
                for t, it in self.runs)
        # A line split into runs cannot be centred with text-anchor="middle":
        # the renderer measures the tspans separately and the space at the run
        # boundary comes out doubled. The width is already known, so a mixed
        # line is anchored at its own left edge instead, and keeps its space.
        if len(self.runs) > 1:
            pos = 'x="%g" text-anchor="start" xml:space="preserve"' % (
                cx - self.width() / 2)
        else:
            pos = 'x="%g" text-anchor="middle"' % cx
        return ('<text %s y="%g" font-family="%s" font-size="%d" '
                'font-weight="%d" fill="#000000">%s</text>'
                % (pos, y, FAMILY, self.size, 700 if self.bold else 400, body))


def H(t):
    # One weight throughout the drawing: the database and analysis headers are
    # regular, like every other label.
    return Line(t, FSH, bold=False)


def Hr(t):
    return Line(t, FSH, bold=False)


def B(t):
    return Line(t, FS)


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


# ---------------------------------------------------------------------- boxes

class Box:
    """A rounded rectangle of stacked lines, optionally with a tinted panel."""

    def __init__(self, head, body=None, stroke=INK, fill=WHITE,
                 panel_fill=None, min_w=0):
        self.head = head or []
        self.body = body or []
        self.stroke, self.fill, self.panel_fill = stroke, fill, panel_fill
        self.min_w = min_w
        self.x = self.y = 0

    def head_h(self):
        return len(self.head) * LHH + (2 * PAD_Y if self.head else 0)

    def body_h(self):
        return len(self.body) * LH + (2 * PAD_Y if self.body else 0)

    @property
    def w(self):
        widest = max([l.width() for l in self.head + self.body] or [0])
        # ceil, not round: rounding down puts the box a fraction of a unit
        # inside its own longest label, which --check then reports as overflow
        return max(self.min_w, math.ceil(widest + 2 * PAD_X))

    @property
    def h(self):
        return self.head_h() + self.body_h()

    def place(self, cx, top):
        self.x, self.y = cx - self.w / 2, top
        return self

    @property
    def cx(self):
        return self.x + self.w / 2

    @property
    def cy(self):
        return self.y + self.h / 2

    @property
    def bottom(self):
        return self.y + self.h

    def svg(self):
        out = []
        if self.panel_fill and self.body:
            # one outer outline, then the tinted panel below the header band,
            # so the pair reads as a single box rather than two stacked ones
            out.append('<rect x="%g" y="%g" rx="6" width="%g" height="%g" '
                       'fill="%s" stroke="none"/>'
                       % (self.x, self.y, self.w, self.h, self.panel_fill))
            out.append('<path d="M%g %gh%gv%gh%gZ" fill="%s"/>'
                       % (self.x, self.y + self.head_h(), self.w,
                          -(self.head_h() - 6), -self.w, self.fill))
            out.append('<rect x="%g" y="%g" rx="6" width="%g" height="%g" '
                       'fill="none" stroke="%s" stroke-width="2.5"/>'
                       % (self.x, self.y, self.w, self.h, self.stroke))
            out.append('<line x1="%g" y1="%g" x2="%g" y2="%g" stroke="%s" '
                       'stroke-width="2"/>'
                       % (self.x, self.y + self.head_h(), self.x + self.w,
                          self.y + self.head_h(), self.stroke))
        else:
            out.append('<rect x="%g" y="%g" rx="6" width="%g" height="%g" '
                       'fill="%s" stroke="%s" stroke-width="2.5"/>'
                       % (self.x, self.y, self.w, self.h, self.fill, self.stroke))
        y = self.y + PAD_Y
        for l in self.head:
            y += LHH
            out.append(l.svg(self.cx, y - LHH * 0.26))
        if self.body:
            # the header's bottom padding and the body's top padding, but only
            # when there is a header above to be padded away from
            if self.head:
                y += 2 * PAD_Y
            for l in self.body:
                y += LH
                out.append(l.svg(self.cx, y - LH * 0.27))
        return "".join(out)


def arrow(x1, y1, x2, y2):
    return ('<line x1="%g" y1="%g" x2="%g" y2="%g" stroke="#000000" '
            'stroke-width="3" marker-end="url(#ah)"/>' % (x1, y1, x2, y2))


def dashed(d, colour):
    return ('<path d="%s" stroke="%s" stroke-width="2.5" '
            'stroke-dasharray="6 5" fill="none"/>' % (d, colour))


def split(x_from, y_from, xs, y_to, mid=None):
    """One trunk down from (x_from, y_from), then a branch to each x in xs.

    `mid` is where the branch bar sits. Left to the midpoint it can land level
    with the processing sidebar, and the bar then runs behind those boxes.
    """
    mid = (y_from + y_to) / 2 if mid is None else mid
    out = ['<path d="M%g %gV%g" stroke="#000000" stroke-width="3" fill="none"/>'
           % (x_from, y_from, mid),
           '<path d="M%g %gH%g" stroke="#000000" stroke-width="3" fill="none"/>'
           % (min(xs), mid, max(xs))]
    for x in xs:
        out.append(arrow(x, mid, x, y_to))
    return "".join(out)


def join(xs, y_from, x_to, y_to):
    mid = (y_from + y_to) / 2
    out = []
    for x in xs:
        out.append('<path d="M%g %gV%g" stroke="#000000" stroke-width="3" '
                   'fill="none"/>' % (x, y_from, mid))
    out.append('<path d="M%g %gH%g" stroke="#000000" stroke-width="3" '
               'fill="none"/>' % (min(xs), mid, max(xs)))
    out.append(arrow(x_to, mid, x_to, y_to))
    return "".join(out)


# ----------------------------------------------------------------------- data

OUTCOMES = ["Completed suicide (n = %s)", "Homicide (n = %s)",
            "Composite H-S group (n = %s)", "Homicidality group (n = %s)",
            "Suicidality group (n = %s)"]

COUNTS = {
    "full":      ["77,347", "1,703", "140", "52,489", "224,487"],
    "mono":      ["16,769", "661", "64", "21,983", "76,676"],
    "ssri":      ["4,430", "334", "47", "3,669", "15,258"],
    "ssri_mono": ["1,115", "138", "20", "1,567", "5,047"],
}
N = {"raw": "23,992,742", "full": "19,687,141", "mono": "11,273,546",
     "ssri": "159,634", "ssri_mono": "65,992",
     "dedup": "4,299,374", "posthoc": "6,227"}


def outcome_lines(key):
    return [B(t % n) for t, n in zip(OUTCOMES, COUNTS[key])]


def harrow(x1, y, x2):
    return arrow(x1, y, x2, y)


def elbow(x1, y1, x2, y2):
    """Down from (x1, y1), across, then an arrow down into (x2, y2)."""
    if abs(x1 - x2) < 1:
        return arrow(x1, y1, x2, y2)
    ym = (y1 + y2) / 2.0
    return ('<path d="M%g %gV%gH%g" stroke="#000000" stroke-width="3" '
            'fill="none"/>' % (x1, y1, ym, x2)) + arrow(x2, ym, x2, y2)


def build():
    el = []

    # --- the four columns ---------------------------------------------------
    # Outcome panels flank the two spines, as in the original drawing. At 8 pt
    # that only fits because the two long right-hand headers wrap: with
    # "MONOPHARMACY DATABASE" and "SSRI MONOPHARMACY GROUP" each on one line
    # the right spine needs 509 units and the row overruns the frame, so they
    # are set over two lines and the spine comes down to 379.
    arms = []
    for keys, lhead, rhead in (
            (("full", "mono"),
             ["FULL DATABASE", "(n = %s)" % N["full"]],
             ["MONOPHARMACY", "DATABASE", "(n = %s)" % N["mono"]]),
            (("ssri", "ssri_mono"),
             ["SSRI GROUP", "(n = %s)" % N["ssri"]],
             ["SSRI MONOPHARMACY", "GROUP", "(n = %s)" % N["ssri_mono"]])):
        red = keys[0] in ("full", "mono")
        stroke, solid, tint = ((RED, RED, RED_L) if red else (BLUE, BLUE, BLUE_LL))
        arms.append({
            "outL": Box([], outcome_lines(keys[0]), stroke=stroke, fill=tint),
            "spL":  Box([H(t) for t in lhead], stroke=stroke, fill=solid),
            "spR":  Box([H(t) for t in rhead], stroke=stroke, fill=solid),
            "outR": Box([], outcome_lines(keys[1]), stroke=stroke, fill=tint),
            "dash": stroke,
        })

    # one width per column, so the two rows line up
    cols = {}
    for k in ("outL", "spL", "spR", "outR"):
        cols[k] = max(a[k].w for a in arms)
        for a in arms:
            a[k].min_w = cols[k]

    # whatever the columns leave over is split between three gaps and two
    # margins, rather than fixed and hoped for
    spare = W - sum(cols.values())
    gap = max(30, min(60, spare / 5.0))
    marg = (spare - 3 * gap) / 2.0
    if marg < 20:
        sys.exit("the four columns leave only %.0f units for margins" % marg)
    xs, run = {}, marg
    for k in ("outL", "spL", "spR", "outR"):
        xs[k] = run + cols[k] / 2.0
        run += cols[k] + gap

    # --- processing steps, three either side of DATA PROCESSING ------------
    # Six steps in one column on the right made the band six boxes tall while
    # the intake beside it was two, so the middle of the band was empty and
    # the first row sat ~250 units lower than it needed to. Split three and
    # three - parsing and deduplication left, standardization and review
    # right - the band is three boxes tall and nothing is stranded.
    steps = [("Data Parsing", False), ("Automated Deduplication", False),
             ("Cases Removed = %s" % N["dedup"], True),
             ("Standardization (MedDRA)", False), (None, False),
             ("Cases Removed = %s" % N["posthoc"], True)]
    side = []
    for label, is_cut in steps:
        line = (Line([("Post hoc ", True), ("review", False)], FS)
                if label is None else B(label))
        side.append(Box([line], stroke=GREEN if is_cut else INK,
                        fill=GREEN_L if is_cut else WHITE))
    sw = max(b.w for b in side)
    for b in side:
        b.min_w = sw
    left_col, right_col = side[:3], side[3:]

    # --- intake, on the centre line -----------------------------------------
    faers = Box([Hr("FAERS DATABASE"), Hr("(2004 \u2013 2025)")],
                stroke=GREEN, fill=GREEN_L)
    raw = Box([Hr("RAW CASE COUNT"), Hr("(n = %s)" % N["raw"])],
              stroke=GREEN, fill=GREEN_L)
    proc = Box([Hr("DATA PROCESSING")], stroke=INK, fill=WHITE)

    faers.place(CXM, M)
    raw.place(faers.x - 44 - raw.w / 2.0, M)
    el += [dashed("M%g %gH%g" % (raw.x + raw.w, raw.cy, faers.x), GREEN)]

    col_top = faers.bottom + 50
    col_h = sum(b.h for b in left_col) + 8 * (len(left_col) - 1)
    proc.place(CXM, col_top + (col_h - proc.h) / 2.0)
    el += [arrow(faers.cx, faers.bottom, faers.cx, proc.y)]

    for col, cx0 in ((left_col, marg + sw / 2.0),
                     (right_col, W - marg - sw / 2.0)):
        y = col_top
        for b in col:
            b.place(cx0, y)
            y += b.h + 8
    side_bottom = col_top + col_h

    # the rails sit midway between the step column and DATA PROCESSING, not
    # hard against the box: hugging it put the left rail on top of the arrow
    # running down into the full database
    for col, rail, edge in (
            (left_col, (left_col[0].x + left_col[0].w + proc.x) / 2.0,
             lambda b: b.x + b.w),
            (right_col, (proc.x + proc.w + right_col[0].x) / 2.0,
             lambda b: b.x)):
        el += [dashed("M%g %gH%g" % (proc.x if col is left_col
                                     else proc.x + proc.w, proc.cy, rail), BLUE),
               dashed("M%g %gV%g" % (rail, col[0].cy, col[-1].cy), BLUE)]
        el += [dashed("M%g %gH%g" % (rail, b.cy, edge(b)), BLUE) for b in col]

    # --- lay the two rows out -----------------------------------------------
    top = max(proc.bottom, side_bottom) + 66
    row_h = max(a[k].h for a in arms for k in ("outL", "spL", "spR", "outR"))
    for i, a in enumerate(arms):
        y0 = top + i * (row_h + 58)
        for k in ("outL", "spL", "spR", "outR"):
            b = a[k]
            b.place(xs[k], y0 + (row_h - b.h) / 2.0)   # centred in the band

    first, second = arms
    # into FULL DATABASE from DATA PROCESSING
    el += [elbow(proc.cx, proc.bottom, first["spL"].cx, first["spL"].y)]
    # FULL DATABASE -> MONOPHARMACY DATABASE: a restriction of it, as before
    el += [harrow(first["spL"].x + first["spL"].w, first["spL"].cy,
                  first["spR"].x)]
    # each spine down to its SSRI box
    for k in ("spL", "spR"):
        el += [arrow(first[k].cx, first[k].bottom, second[k].cx, second[k].y)]
    # the outcome panels hang off their spine on a dashed link
    for a in arms:
        el += [dashed("M%g %gH%g" % (a["outL"].x + a["outL"].w, a["spL"].cy,
                                     a["spL"].x), a["dash"]),
               dashed("M%g %gH%g" % (a["spR"].x + a["spR"].w, a["spR"].cy,
                                     a["outR"].x), a["dash"])]
    for a in arms:
        el += [a[k].svg() for k in ("outL", "spL", "spR", "outR")]

    bottom = max(a[k].bottom for a in arms for k in ("outL", "outR"))

    # --- analysis -----------------------------------------------------------
    da = Box([H("DISPROPORTIONALITY"), H("ANALYSIS")],
             [B("ROR, PRR, EBGM, & IC")], stroke=INK, fill=WHITE,
             panel_fill=WHITE)
    tto = Box([H("TTO analysis")], stroke=INK, fill=WHITE)
    cl = Box([H("Case-list analysis")], stroke=INK, fill=WHITE)
    aw = max(da.w, tto.w, cl.w)
    da.min_w = tto.min_w = cl.min_w = aw
    da.place(CXM, bottom + 74)
    tto.place(CXM, da.bottom + 10)
    cl.place(CXM, tto.bottom + 10)
    el += [join([second["spL"].cx, second["spR"].cx],
                max(second["spL"].bottom, second["spR"].bottom), CXM, da.y)]
    el += [da.svg(), tto.svg(), cl.svg()]
    for side_x, out in ((da.x, -26), (da.x + da.w, 26)):
        el += [dashed("M%g %gH%gV%gH%g"
                      % (side_x, da.bottom - 12, side_x + out, cl.cy, side_x),
                      BLUE)]

    el += [raw.svg(), faers.svg(), proc.svg()] + [b.svg() for b in side]

    height = int(cl.bottom + M)
    head = ('<svg width="%d" height="%d" viewBox="0 0 %d %d" '
            'xmlns="http://www.w3.org/2000/svg" '
            'xmlns:xlink="http://www.w3.org/1999/xlink" overflow="hidden">'
            '<defs><marker id="ah" refX="9" refY="5" markerWidth="6" '
            'markerHeight="6" markerUnits="strokeWidth" '
            'orient="auto-start-reverse" viewBox="0 0 10 10" '
            'overflow="hidden"><path d="M0 0 10 5 0 10Z"/></marker></defs>'
            '<rect x="0" y="0" width="%d" height="%d" fill="#FFFFFF"/>'
            % (W, height, W, height, W, height))
    boxes = [raw, faers, proc] + side + \
            [a[k] for a in arms for k in ("outL", "spL", "spR", "outR")] + \
            [da, tto, cl]
    return head + "".join(el) + "</svg>", height, boxes


def report(height, boxes):
    print("  frame        %d x %d units -> %.1f x %.1f mm"
          % (W, height, MM, MM * height / W))
    if MM * height / W > MAX_H_MM:
        sys.exit("  the figure is %.0f mm tall, over the %.0f mm limit"
                 % (MM * height / W, MAX_H_MM))
    print("  type         %d / %d units -> %.2f / %.2f pt   (floor %.1f)"
          % (FS, FSH, FS * K, FSH * K, FLOOR_PT))
    if min(FS, FSH) * K < FLOOR_PT:
        sys.exit("  type is under the %.1f pt floor" % FLOOR_PT)
    print("  rule         2.0-3.0 units -> %.2f-%.2f mm   (floor 0.10 mm)"
          % (2.0 * K / PT_PER_MM, 3.0 * K / PT_PER_MM))
    worst, lines = None, 0
    for b in boxes:
        for l in b.head + b.body:
            lines += 1
            slack = b.w - 2 * PAD_X - l.width()
            if worst is None or slack < worst[0]:
                worst = (slack, l.plain, b.w)
    print("  %d labels, tightest fit %.0f units spare: %r in a %g-unit box"
          % (lines, worst[0], worst[1], worst[2]))
    if worst[0] < 0:
        sys.exit("  a label overflows its box")


def main():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(
        here, "expected_2004_2025", "figures", "fig1_study_flow.svg"))
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    svg, height, boxes = build()
    report(height, boxes)
    if a.check:
        print("  check only: nothing written")
        return
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(svg)
    print("  wrote %s (%d bytes)" % (a.out, len(svg)))


if __name__ == "__main__":
    main()
