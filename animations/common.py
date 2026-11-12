"""Shared helpers for the 量子力学 manim scenes (Manim Community + LaTeX).

公式统一用 LaTeX 排版：
  MT(r"...")  纯数学公式（MathTex）
  CT(r"...")  中文 + 行内公式（Tex + ctex 模板，xelatex 编译）
依赖：TeX 发行版（TinyTeX 即可）+ ctex、standalone、preview、dvisvgm 等宏包。
"""
import numpy as np
from manim import *  # noqa: F401,F403

# ---------------------------------------------------------------- palette
C_PSI = BLUE_C      # Re ψ
C_IM = PURPLE_B     # Im ψ
C_PROB = YELLOW     # |ψ|^2
C_POT = GREY_B      # 势能
C_HL = ORANGE       # 强调
C_E = GREEN_C       # 能级

MACROS = (r"\newcommand{\ket}[1]{\left|#1\right\rangle}\newcommand{\bra}[1]{\left\langle#1\right|}"
          r"\newcommand{\ip}[2]{\left\langle#1\middle|#2\right\rangle}"
          r"\newcommand{\expval}[1]{\left\langle#1\right\rangle}\newcommand{\comm}[2]{\left[#1,#2\right]}"
          r"\newcommand{\vv}[1]{\boldsymbol{#1}}\newcommand{\dd}{\mathrm{d}}\newcommand{\ii}{\mathrm{i}}\newcommand{\ee}{\mathrm{e}}")
CTEX = TexTemplate(tex_compiler="xelatex", output_format=".xdv")
CTEX.add_to_preamble(r"\usepackage[UTF8]{ctex}" "\n" r"\usepackage{amsmath,amssymb,bm,xcolor}")
CTEX.add_to_preamble(MACROS)
CTEX.add_to_preamble(r"\setlength{\textwidth}{60cm}\setlength{\hsize}{60cm}\setlength{\linewidth}{60cm}")
MTEX = TexTemplate()
MTEX.add_to_preamble(r"\usepackage{amsmath,amssymb,bm}" "\n" + MACROS)


def MT(tex: str, size: int = 32, color=WHITE, **kw) -> MathTex:
    return MathTex(tex, font_size=size, color=color, tex_template=MTEX, **kw)


def CT(tex: str, size: int = 30, color=WHITE, **kw) -> Tex:
    """Chinese text with inline $math$ (xelatex + ctex)."""
    return Tex(tex, tex_template=CTEX, font_size=size, color=color, **kw)


def fitw(m, w):
    if m.width > w:
        m.scale_to_fit_width(w)
    return m


def title_bar(s: str, size=34, w=13.4):
    t = fitw(CT(s, size), w).to_edge(UP, buff=0.22)
    bg = BackgroundRectangle(t, fill_opacity=0.85, buff=0.08)
    return VGroup(bg, t)


def panel(*lines, size=28, buff=0.22, color=WHITE):
    mobs = [CT(l, size, color) if isinstance(l, str) else l for l in lines]
    return VGroup(*mobs).arrange(DOWN, aligned_edge=LEFT, buff=buff)


def boxed(m, color=GREY_D):
    r = SurroundingRectangle(m, color=color, buff=0.2, corner_radius=0.1)
    r.set_fill(BLACK, opacity=0.75)
    return VGroup(r, m)


def curve(ax, xs, ys, color=WHITE, sw=3.0):
    """Polyline through sampled data in Axes coordinates."""
    pts = [ax.c2p(x, y) for x, y in zip(xs, ys)]
    m = VMobject(color=color, stroke_width=sw)
    m.set_points_as_corners(pts)
    return m


def filled_curve(ax, xs, ys, color=YELLOW, opacity=0.35, sw=2.5, y0=0.0):
    pts = [ax.c2p(xs[0], y0)] + [ax.c2p(x, y) for x, y in zip(xs, ys)] + [ax.c2p(xs[-1], y0)]
    m = Polygon(*pts, color=color, stroke_width=0, fill_opacity=opacity)
    line = curve(ax, xs, ys, color, sw)
    return VGroup(m, line)
