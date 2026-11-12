"""第 5–6 章（微扰理论 / 散射）—— manim 动画源码（LaTeX 排版，自包含）。

渲染：tools/render_videos.sh ch5_6_scenes.py [Scene ...]
依赖：manim (Community)、numpy、scipy、TeX（xelatex + ctex）。
"""
import numpy as np
from manim import *  # noqa: F401,F403
from scipy.special import spherical_jn, spherical_yn

# ------------------------------------------------------------------ TeX helpers
CTEX = TexTemplate(tex_compiler="xelatex", output_format=".xdv")
CTEX.add_to_preamble(r"\usepackage[UTF8]{ctex}" "\n" r"\usepackage{amsmath,amssymb,bm,xcolor}")
CTEX.add_to_preamble(r"\setlength{\textwidth}{60cm}\setlength{\hsize}{60cm}\setlength{\linewidth}{60cm}")
MTEX = TexTemplate()
MTEX.add_to_preamble(r"\usepackage{amsmath,amssymb,bm}")

C_A = YELLOW
C_B = TEAL
C_HL = ORANGE
C_UP = RED_C
C_DN = BLUE_C


def MT(tex, size=32, color=WHITE, **kw):
    return MathTex(tex, font_size=size, color=color, tex_template=MTEX, **kw)


def CT(tex, size=28, color=WHITE, **kw):
    return Tex(tex, tex_template=CTEX, font_size=size, color=color, **kw)


def fitw(m, w):
    if m.width > w:
        m.scale_to_fit_width(w)
    return m


def title_bar(s, size=32, w=13.4):
    t = fitw(CT(s, size), w).to_edge(UP, buff=0.22)
    bg = BackgroundRectangle(t, fill_opacity=0.85, buff=0.08)
    return VGroup(bg, t)


def panel(*lines, size=26, buff=0.2, color=WHITE, w=None):
    mobs = [CT(l, size, color) if isinstance(l, str) else l for l in lines]
    g = VGroup(*mobs).arrange(DOWN, aligned_edge=LEFT, buff=buff)
    if w:
        fitw(g, w)
    return g


def boxed(m, color=GREY_D, buff=0.2):
    r = SurroundingRectangle(m, color=color, buff=buff, corner_radius=0.1)
    r.set_fill(BLACK, opacity=0.8)
    return VGroup(r, m)


def seq_cmap(v):
    """0..1 -> RGB (black-blue-orange-white)."""
    stops = np.array([[0, 0, 0, 0], [0.3, 25, 50, 150], [0.62, 230, 120, 40], [1.0, 255, 250, 210]], float)
    v = np.clip(v, 0, 1)
    out = np.zeros(v.shape + (3,))
    for c in range(3):
        out[..., c] = np.interp(v, stops[:, 0], stops[:, c + 1])
    return out


def div_cmap(v):
    """-1..1 -> RGB (blue - black - orange)."""
    v = np.clip(v, -1, 1)
    out = np.zeros(v.shape + (3,))
    pos = np.clip(v, 0, 1)[..., None] * np.array([255, 170, 60])
    neg = np.clip(-v, 0, 1)[..., None] * np.array([70, 150, 255])
    out = pos + neg
    return out


def to_image(rgb, height):
    a = np.concatenate([rgb.astype(np.uint8), 255 * np.ones(rgb.shape[:2] + (1,), np.uint8)], axis=2)
    im = ImageMobject(a)
    im.set_resampling_algorithm(RESAMPLING_ALGORITHMS["bilinear"])
    im.height = height
    return im


def level(y, x0, x1, color=WHITE, sw=4):
    return Line([x0, y, 0], [x1, y, 0], color=color, stroke_width=sw)


# ==================================================================
# 第 5 章 微扰理论
# ==================================================================
class C5_PerturbMatrix(Scene):
    """非简并微扰 = 在 H0 本征基里把矩阵对角化的级数解。"""

    def construct(self):
        title = title_bar(r"定态微扰：在 $\hat H^{(0)}$ 的本征基里，把矩阵 $\hat H^{(0)}+\lambda\hat H'$ 对角化")
        self.play(FadeIn(title))
        E0 = np.array([1.0, 2.2, 3.6])
        Hp = np.array([[0.30, 0.55, 0.25], [0.55, -0.20, 0.45], [0.25, 0.45, 0.10]])

        mat = MT(r"H=\begin{pmatrix}E_1^{(0)}+\lambda H'_{11}&\lambda H'_{12}&\lambda H'_{13}\\"
                 r"\lambda H'_{21}&E_2^{(0)}+\lambda H'_{22}&\lambda H'_{23}\\"
                 r"\lambda H'_{31}&\lambda H'_{32}&E_3^{(0)}+\lambda H'_{33}\end{pmatrix}", 26)
        mat.to_edge(LEFT, buff=0.4).shift(UP * 1.7)
        fitw(mat, 6.3)
        notes = panel(
            r"对角元 $\Rightarrow$ 一级修正：$E_n^{(1)}=H'_{nn}$",
            r"非对角元 $\Rightarrow$ 二级修正：$E_n^{(2)}=\sum_m'\dfrac{|H'_{mn}|^2}{E_n^{(0)}-E_m^{(0)}}$",
            r"判据：$\left|\dfrac{H'_{mn}}{E_n^{(0)}-E_m^{(0)}}\right|\ll1$",
            size=24, w=6.2)
        notes.next_to(mat, DOWN, buff=0.45, aligned_edge=LEFT)
        self.play(Write(mat), run_time=2)
        self.play(FadeIn(notes[0]))

        ax = Axes(x_range=[0, 1.0, 0.25], y_range=[0.4, 4.4, 1], x_length=6.0, y_length=5.4,
                  axis_config={"include_tip": False, "font_size": 20}).to_edge(RIGHT, buff=0.45).shift(DOWN * 0.4)
        xl = MT(r"\lambda", 26).next_to(ax.x_axis, DOWN, buff=0.15).align_to(ax.x_axis, RIGHT)
        yl = CT(r"能量 $E$", 22).next_to(ax.y_axis, UP, buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))
        cols = [C_DN, C_A, C_UP]

        def exact(l):
            return np.linalg.eigvalsh(np.diag(E0) + l * Hp)

        def second(n, l):
            s = sum(Hp[m, n] ** 2 / (E0[n] - E0[m]) for m in range(3) if m != n)
            return E0[n] + l * Hp[n, n] + l * l * s

        ex = VGroup(*[ax.plot(lambda l, n=n: exact(l)[n], x_range=[0, 1, 0.01], color=cols[n], stroke_width=4) for n in range(3)])
        first = VGroup(*[DashedLine(ax.c2p(0, E0[n]), ax.c2p(1, E0[n] + Hp[n, n]), color=cols[n], stroke_width=2.5, dash_length=0.08) for n in range(3)])
        sec = VGroup(*[DashedVMobject(ax.plot(lambda l, n=n: second(n, l), x_range=[0, 1, 0.01], color=WHITE, stroke_width=2), num_dashes=40) for n in range(3)])
        labs = VGroup(*[MT(f"E_{n+1}", 24, cols[n]).next_to(ax.c2p(0, E0[n]), LEFT, buff=0.1) for n in range(3)])
        self.play(FadeIn(labs), Create(first), run_time=1.5)
        leg1 = CT(r"虚线：一级 $E^{(0)}+\lambda H'_{nn}$", 20).next_to(yl, RIGHT, buff=0.5)
        self.play(FadeIn(leg1))
        self.wait(0.5)
        self.play(FadeIn(notes[1]), Create(sec), run_time=1.6)
        leg2 = CT(r"白虚线：二级；实线：精确本征值", 20).next_to(leg1, DOWN, buff=0.08, aligned_edge=LEFT)
        self.play(FadeIn(leg2))
        self.play(LaggedStart(*[Create(c) for c in ex], lag_ratio=0.2), run_time=2.0)
        lam = ValueTracker(0.0)
        dots = VGroup(*[always_redraw(lambda n=n: Dot(ax.c2p(lam.get_value(), exact(lam.get_value())[n]), color=cols[n], radius=0.07)) for n in range(3)])
        lv = always_redraw(lambda: MT(r"\lambda=%.2f" % lam.get_value(), 26, C_HL).move_to(ax.c2p(0.8, 0.62)))
        self.add(dots, lv)
        self.play(lam.animate.set_value(1.0), run_time=4, rate_func=linear)
        self.play(FadeIn(notes[2]))
        msg = boxed(CT(r"二级修正：下能级被推下、上能级被推上——\textbf{能级排斥}；$\lambda$ 越大，级数越不够用", 24, C_HL))
        fitw(msg, 12.8).to_edge(DOWN, buff=0.2)
        self.play(FadeIn(msg))
        self.wait(2.5)


class C5_AvoidedCrossing(Scene):
    """两能级耦合：交叉变为回避交叉，简并点处混合 50/50。"""

    def construct(self):
        title = title_bar(r"两能级模型：耦合 $V$ 把“交叉”变成“回避交叉”，简并处必须先在子空间里对角化")
        self.play(FadeIn(title))
        ax = Axes(x_range=[-3, 3, 1], y_range=[-3.4, 3.4, 1], x_length=6.6, y_length=5.6,
                  axis_config={"include_tip": False, "include_ticks": False}).to_edge(LEFT, buff=0.5).shift(DOWN * 0.35)
        xl = CT(r"失谐 $\varepsilon=(E_a-E_b)/2$", 22).next_to(ax.x_axis, DOWN, buff=0.15).align_to(ax, RIGHT)
        yl = CT(r"能量", 22).next_to(ax.y_axis, UP, buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))
        da = DashedLine(ax.c2p(-3, -3), ax.c2p(3, 3), color=C_A, stroke_width=2.5)
        db = DashedLine(ax.c2p(-3, 3), ax.c2p(3, -3), color=C_B, stroke_width=2.5)
        la = MT(r"|a\rangle", 26, C_A).next_to(ax.c2p(3, 3), LEFT, buff=0.1)
        lb = MT(r"|b\rangle", 26, C_B).next_to(ax.c2p(3, -3), LEFT, buff=0.1).shift(DOWN * 0.15)
        self.play(Create(da), Create(db), FadeIn(la), FadeIn(lb))
        H = MT(r"H=\begin{pmatrix}\varepsilon & V\\ V & -\varepsilon\end{pmatrix},\quad E_\pm=\pm\sqrt{\varepsilon^2+V^2}", 30)
        H.to_edge(RIGHT, buff=0.4).shift(UP * 2.0)
        fitw(H, 6.0)
        self.play(Write(H))
        V = ValueTracker(0.0)
        up = always_redraw(lambda: ax.plot(lambda e: np.sqrt(e * e + V.get_value() ** 2), x_range=[-3, 3, 0.02], color=WHITE, stroke_width=4))
        dn = always_redraw(lambda: ax.plot(lambda e: -np.sqrt(e * e + V.get_value() ** 2), x_range=[-3, 3, 0.02], color=WHITE, stroke_width=4))
        self.add(up, dn)
        vlab = always_redraw(lambda: MT(r"V=%.2f" % V.get_value(), 28, C_HL).next_to(H, DOWN, buff=0.3))
        self.add(vlab)
        self.play(V.animate.set_value(0.8), run_time=2.5)
        gap = DoubleArrow(ax.c2p(0, -0.8), ax.c2p(0, 0.8), buff=0, color=C_HL, stroke_width=3, tip_length=0.15)
        gl = MT(r"2|V|", 26, C_HL).next_to(gap, RIGHT, buff=0.1)
        self.play(GrowFromCenter(gap), FadeIn(gl))
        # composition bars for lower state
        eps = ValueTracker(-2.8)

        def lower_vec(e):
            w, v = np.linalg.eigh(np.array([[e, V.get_value()], [V.get_value(), -e]]))
            return v[:, 0]

        dot = always_redraw(lambda: Dot(ax.c2p(eps.get_value(), -np.sqrt(eps.get_value() ** 2 + V.get_value() ** 2)), color=C_HL, radius=0.09))
        base = np.array([2.6, -2.2, 0])
        frame = VGroup(Line(base + LEFT * 0.2, base + RIGHT * 3.4, color=GREY_B))
        bl = CT(r"下能级本征态的成分 $|c_a|^2,\ |c_b|^2$", 22).next_to(base + RIGHT * 1.6 + UP * 2.25, UP, buff=0.05)

        def bars():
            c = lower_vec(eps.get_value()) ** 2
            g = VGroup()
            for i, (p, col, nm) in enumerate(zip(c, [C_A, C_B], ["a", "b"])):
                r = Rectangle(width=1.0, height=max(2.0 * p, 0.001), fill_color=col, fill_opacity=0.85, stroke_width=0)
                r.move_to(base + RIGHT * (0.7 + 1.8 * i), aligned_edge=DOWN)
                g.add(r, MT(r"%.2f" % p, 22).next_to(r, UP, buff=0.05), MT(r"|%s\rangle" % nm, 24, col).next_to(base + RIGHT * (0.7 + 1.8 * i), DOWN, buff=0.12))
            return g
        bb = always_redraw(bars)
        self.play(FadeIn(frame), FadeIn(bl), FadeIn(dot), FadeIn(bb))
        self.play(eps.animate.set_value(0.0), run_time=3, rate_func=smooth)
        note = boxed(CT(r"简并点 $\varepsilon=0$：零级态是 $\tfrac{1}{\sqrt2}(|a\rangle\pm|b\rangle)$，由久期方程决定，不是 $|a\rangle$ 或 $|b\rangle$", 22, C_HL))
        fitw(note, 12.8).to_edge(DOWN, buff=0.15)
        self.play(FadeIn(note))
        self.wait(1.5)
        self.play(eps.animate.set_value(2.8), run_time=3, rate_func=smooth)
        note2 = boxed(CT(r"缓慢扫过：态从 $|b\rangle$ 连续变成 $|a\rangle$；$|\varepsilon|\gg|V|$ 时才能用非简并公式", 22, C_HL))
        fitw(note2, 12.8).to_edge(DOWN, buff=0.15)
        self.play(ReplacementTransform(note, note2))
        self.wait(2.5)


class C5_Stark(Scene):
    """氢原子 n=2 一级斯塔克效应：四重简并分裂为三条，零级态是 2s 与 2p0 的叠加。"""

    def construct(self):
        title = title_bar(r"一级斯塔克效应：$n=2$ 的四重简并在电场中分裂为三，零级态 $\tfrac{1}{\sqrt2}(\psi_{200}\mp\psi_{210})$")
        self.play(FadeIn(title))
        F = ValueTracker(0.0)
        x0, x1 = -6.4, -2.4
        y2 = 1.2
        d = 1.15
        base1 = level(-2.9, x0, x1, GREY_B)
        l1 = MT(r"n=1\ (1\text{s})", 24, GREY_B).next_to(base1, DOWN, buff=0.1)
        lup = always_redraw(lambda: level(y2 + d * F.get_value(), x0, x1, C_UP))
        lmid = always_redraw(lambda: level(y2, x0, x1, WHITE))
        ldn = always_redraw(lambda: level(y2 - d * F.get_value(), x0, x1, C_DN))
        self.play(Create(base1), FadeIn(l1), FadeIn(lup), FadeIn(lmid), FadeIn(ldn))
        tag0 = CT(r"$n=2$：$\psi_{200},\psi_{210},\psi_{211},\psi_{21,-1}$（四重简并）", 22).next_to(lmid, UP, buff=0.12).align_to(lmid, LEFT)
        fitw(tag0, 5.6)
        self.play(FadeIn(tag0))
        mat = MT(r"H'=e\mathcal{E}z:\ \begin{pmatrix}0&-3e\mathcal{E}a_0&0&0\\-3e\mathcal{E}a_0&0&0&0\\0&0&0&0\\0&0&0&0\end{pmatrix}", 26)
        mat.to_edge(RIGHT, buff=0.4).shift(UP * 1.9)
        fitw(mat, 6.4)
        rule = CT(r"只有 $\Delta l=\pm1,\ \Delta m=0$ 的矩阵元非零（宇称 + $\varphi$ 积分）", 22).next_to(mat, DOWN, buff=0.2)
        fitw(rule, 6.0).to_edge(RIGHT, buff=0.3)
        fitw(rule, 6.4)
        self.play(Write(mat), FadeIn(rule))
        self.wait(0.8)
        self.play(FadeOut(tag0), F.animate.set_value(1.0), run_time=3)
        tu = MT(r"E_2^{(0)}+3e\mathcal{E}a_0", 24, C_UP).next_to(lup, RIGHT, buff=0.15)
        tm = CT(r"$E_2^{(0)}$（$\psi_{21,\pm1}$ 仍二重简并）", 22).next_to(lmid, RIGHT, buff=0.15)
        td = MT(r"E_2^{(0)}-3e\mathcal{E}a_0", 24, C_DN).next_to(ldn, RIGHT, buff=0.15)
        fitw(tm, 3.6)
        self.play(FadeIn(tu), FadeIn(tm), FadeIn(td))
        arrs = VGroup(*[Arrow([xx, yy, 0], [xx, -2.9, 0], buff=0.02, color=c, stroke_width=3, max_tip_length_to_length_ratio=0.06)
                        for xx, yy, c in [(-5.6, y2 + d, C_UP), (-4.4, y2, WHITE), (-3.2, y2 - d, C_DN)]])
        self.play(LaggedStart(*[GrowArrow(a) for a in arrs], lag_ratio=0.3))
        sp = CT(r"一条莱曼 $\alpha$ 线 $\to$ 三条", 22, C_HL).next_to(base1, UP, buff=1.0).align_to(base1, LEFT).shift(RIGHT * 0.1)
        self.play(FadeIn(sp))
        # densities
        n = 220
        L = 14.0
        xs = np.linspace(-L, L, n)
        X, Z = np.meshgrid(xs, xs[::-1])
        R = np.sqrt(X ** 2 + Z ** 2)
        up = ((2 - R - Z) * np.exp(-R / 2)) ** 2 / 2
        dn = ((2 - R + Z) * np.exp(-R / 2)) ** 2 / 2
        m = max(up.max(), dn.max())
        im_up = to_image(seq_cmap((up / m) ** 0.45), 2.5)
        im_dn = to_image(seq_cmap((dn / m) ** 0.45), 2.5)
        im_up.move_to([2.4, -1.45, 0])
        im_dn.move_to([5.4, -1.45, 0])
        cu = MT(r"\tfrac{1}{\sqrt2}(\psi_{200}-\psi_{210})", 22, C_UP).next_to(im_up, DOWN, buff=0.1)
        cd = MT(r"\tfrac{1}{\sqrt2}(\psi_{200}+\psi_{210})", 22, C_DN).next_to(im_dn, DOWN, buff=0.1)
        Ear = Arrow(im_up.get_left() + LEFT * 0.55 + DOWN * 1.0, im_up.get_left() + LEFT * 0.55 + UP * 1.0, buff=0, color=YELLOW, stroke_width=4)
        El = MT(r"\boldsymbol{\mathcal E}", 26, YELLOW).next_to(Ear, LEFT, buff=0.08)
        self.play(FadeOut(rule), FadeIn(im_up), FadeIn(im_dn), FadeIn(cu), FadeIn(cd), GrowArrow(Ear), FadeIn(El))
        pu = CT(r"电子云偏 $+z$：能量升高", 20, C_UP).next_to(im_up, UP, buff=0.08)
        pd = CT(r"电子云偏 $-z$：能量降低", 20, C_DN).next_to(im_dn, UP, buff=0.08)
        self.play(FadeIn(pu), FadeIn(pd))
        msg = boxed(CT(r"原子像带有大小 $3ea_0$ 的\textbf{永久}电偶极矩：能移 $\propto\mathcal E$（一级）；而基态 1s 不简并，只有 $\propto\mathcal E^2$ 的二级效应", 22, C_HL))
        fitw(msg, 13.0).to_edge(DOWN, buff=0.12)
        self.play(FadeIn(msg))
        self.wait(3)


class C5_GoldenRule(Scene):
    """sin^2 跃迁概率随时间：峰高 ∝ t^2，宽 ∝ 1/t，面积 ∝ t ⇒ 跃迁速率为常数（费米黄金定则）。"""

    def construct(self):
        title = title_bar(r"常微扰下的跃迁概率：峰高 $\propto t^2$、宽度 $\propto 1/t$、面积 $\propto t$ $\Rightarrow$ 费米黄金定则")
        self.play(FadeIn(title))
        f = MT(r"W_{k\to m}(t)=\frac{|H'_{mk}|^2}{\hbar^2}\cdot\frac{\sin^2(\omega_{mk}t/2)}{(\omega_{mk}/2)^2}", 30).next_to(title, DOWN, buff=0.25)
        self.play(Write(f))
        axA = Axes(x_range=[-10, 10, 5], y_range=[0, 40, 10], x_length=5.8, y_length=3.7,
                   axis_config={"include_tip": False, "font_size": 18}).move_to([-3.4, -1.0, 0])
        axB = Axes(x_range=[-10, 10, 5], y_range=[0, 7, 1], x_length=5.8, y_length=3.7,
                   axis_config={"include_tip": False, "font_size": 18}).move_to([3.4, -1.0, 0])
        for ax in (axA, axB):
            ax.add_coordinates(font_size=18)
        la = CT(r"$\dfrac{W}{|H'|^2/\hbar^2}$ 对 $\omega_{mk}$", 22).next_to(axA, UP, buff=0.1)
        lb = CT(r"除以 $t$ 之后：$\dfrac{W/t}{|H'|^2/\hbar^2}\ \to\ 2\pi\delta(\omega_{mk})$", 22).next_to(axB, UP, buff=0.1)
        xa = MT(r"\omega_{mk}", 22).next_to(axA.x_axis, DOWN, buff=0.35)
        xb = MT(r"\omega_{mk}", 22).next_to(axB.x_axis, DOWN, buff=0.35)
        self.play(Create(axA), Create(axB), FadeIn(la), FadeIn(lb), FadeIn(xa), FadeIn(xb))
        t = ValueTracker(1.0)

        def F(w, tt):
            w = np.where(np.abs(w) < 1e-6, 1e-6, w)
            return np.sin(w * tt / 2) ** 2 / (w / 2) ** 2

        cA = always_redraw(lambda: axA.plot(lambda w: F(w, t.get_value()), x_range=[-10, 10, 0.02], color=C_A, stroke_width=3))
        cB = always_redraw(lambda: axB.plot(lambda w: F(w, t.get_value()) / t.get_value(), x_range=[-10, 10, 0.02], color=C_B, stroke_width=3))
        aB = always_redraw(lambda: axB.get_area(axB.plot(lambda w: F(w, t.get_value()) / t.get_value(), x_range=[-10, 10, 0.02]), x_range=[-10, 10], color=C_B, opacity=0.3))
        tl = always_redraw(lambda: MT(r"t=%.1f" % t.get_value(), 30, C_HL).move_to([0, 0.95, 0]))
        self.add(aB, cA, cB, tl)
        self.play(t.animate.set_value(6.0), run_time=6, rate_func=linear)
        wid = DoubleArrow(axA.c2p(-2 * PI / 6, 38), axA.c2p(2 * PI / 6, 38), buff=0, color=C_HL, stroke_width=3, tip_length=0.12)
        wl = MT(r"\Delta\omega\sim\frac{4\pi}{t}", 24, C_HL).next_to(wid, RIGHT, buff=0.1)
        ar = CT(r"阴影面积恒为 $2\pi$", 22, C_B).move_to(axB.c2p(5.5, 4.5))
        self.play(GrowFromCenter(wid), FadeIn(wl), FadeIn(ar))
        self.wait(1)
        res = boxed(VGroup(
            CT(r"末态连续分布（态密度 $\rho$）时，对主峰积分：$W\approx\dfrac{2\pi t}{\hbar}|H'_{mk}|^2\rho(m)$", 24),
            CT(r"$\Rightarrow$ 单位时间跃迁概率 $w=\dfrac{2\pi}{\hbar}|H'_{mk}|^2\rho(m)$ 与 $t$ 无关；$\ \Delta E\,\Delta t\sim\hbar$ 也藏在峰宽里", 24, C_HL),
        ).arrange(DOWN, buff=0.15))
        fitw(res, 13.0).to_edge(DOWN, buff=0.12)
        self.play(FadeIn(res))
        self.wait(3)


class C5_Helium(Scene):
    """氦原子基态：一级微扰 = 变分曲线在 Z=2 处的值；变分取 Z=27/16 更低更准。"""

    def construct(self):
        title = title_bar(r"氦原子基态：一级微扰与变分法其实在同一条曲线 $\overline H(Z)=Z^2-\tfrac{27}{8}Z$ 上")
        self.play(FadeIn(title))
        ax = Axes(x_range=[1.0, 2.4, 0.2], y_range=[-3.0, -2.0, 0.2], x_length=7.2, y_length=5.4,
                  axis_config={"include_tip": False, "font_size": 18}).to_edge(LEFT, buff=0.7).shift(DOWN * 0.4)
        ax.add_coordinates(font_size=18)
        xl = CT(r"尝试参数 $Z$（有效核电荷）", 22).next_to(ax.x_axis, DOWN, buff=0.45)
        yl = CT(r"$\overline H$ / $\dfrac{e^2}{4\pi\varepsilon_0a_0}$", 22).next_to(ax.y_axis, UP, buff=0.12).shift(RIGHT * 0.3)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))
        g = ax.plot(lambda z: z * z - 27 / 8 * z, x_range=[1.0, 2.4, 0.01], color=C_A, stroke_width=4)
        terms = panel(r"$\overline H(Z)=\underbrace{Z^2}_{\text{动能}}\underbrace{-4Z}_{\text{核吸引}}+\underbrace{\tfrac58Z}_{\text{电子互斥}}$", size=26)
        terms.to_edge(RIGHT, buff=0.4).shift(UP * 2.0)
        fitw(terms, 5.6)
        self.play(Create(g), FadeIn(terms), run_time=2)
        z = ValueTracker(2.3)
        dot = always_redraw(lambda: Dot(ax.c2p(z.get_value(), z.get_value() ** 2 - 27 / 8 * z.get_value()), color=C_HL, radius=0.08))
        self.add(dot)
        self.play(z.animate.set_value(2.0), run_time=1.2)
        p2 = Dot(ax.c2p(2.0, -2.75), color=C_UP, radius=0.09)
        t2 = CT(r"$Z=2$：$-2.75$ = 一级微扰结果", 22, C_UP).next_to(p2, RIGHT, buff=0.15).shift(UP * 0.15)
        self.play(FadeIn(p2), FadeIn(t2))
        self.play(z.animate.set_value(27 / 16), run_time=2)
        pm = Dot(ax.c2p(27 / 16, -(27 / 16) ** 2), color=C_B, radius=0.09)
        tm = CT(r"$Z=\tfrac{27}{16}$：$-2.85$（变分极小）", 22, C_B).move_to(ax.c2p(1.69, -2.955))
        ex = DashedLine(ax.c2p(1.0, -2.904), ax.c2p(2.4, -2.904), color=WHITE, stroke_width=2)
        el = CT(r"实验 $-2.904$", 20).next_to(ex, UP, buff=0.05).align_to(ex, RIGHT)
        self.play(FadeIn(pm), FadeIn(tm), Create(ex), FadeIn(el))
        bars = panel(
            r"零级（忽略互斥）：$-4$",
            r"一级微扰：$-2.75$（误差约 5\%）",
            r"变分 $Z=27/16$：$-2.85$（误差约 2\%）",
            r"实验：$-2.904$（约 $-79.0$ eV）",
            size=24, w=5.6)
        bars.next_to(terms, DOWN, buff=0.45).to_edge(RIGHT, buff=0.35)
        self.play(FadeIn(bars))
        msg = boxed(CT(r"变分上界：任意尝试函数 $\overline H\ge E_0$；一个电子“看到”的核电荷被另一个电子屏蔽成 $Z\approx1.69$", 22, C_HL))
        fitw(msg, 13.0).to_edge(DOWN, buff=0.12)
        self.play(FadeIn(msg))
        self.wait(3)


# ==================================================================
# 第 6 章 散射
# ==================================================================
class C6_ScatterWave(Scene):
    """入射平面波 + 出射球面波；散射振幅 f(θ) 决定各方向的粒子流。"""

    def construct(self):
        title = title_bar(r"散射的定态图像：$\psi\to e^{\mathrm ikz}+f(\theta)\dfrac{e^{\mathrm ikr}}{r}$（取实部，随时间演化）")
        self.play(FadeIn(title))
        n = 300
        L = 12.0
        xs = np.linspace(-L, L, n)
        Zg, Xg = np.meshgrid(xs, xs[::-1])  # z horizontal, x vertical
        R = np.sqrt(Zg ** 2 + Xg ** 2) + 1e-9
        k = 2 * np.pi / 2.4
        cos_t = Zg / R
        a = 1.0
        f = (-1.6 * (1 + 0.6 * cos_t)) * np.exp(0.5j)
        inc = np.exp(1j * k * Zg)
        sca = f * np.exp(1j * k * R) / R
        mask = R < a
        damp = np.clip((R - a) / 0.6, 0, 1)
        sca = sca * damp
        tot = (inc + sca)
        tot[mask] = 0
        sca[mask] = 0
        H = 5.9
        phase = ValueTracker(0.0)
        mode = {"which": tot, "scale": 1.6}

        def mk():
            w = mode["which"] * np.exp(-1j * phase.get_value())
            im = to_image(div_cmap(np.real(w) / mode["scale"]), H)
            im.move_to([-3.0, -0.55, 0])
            return im
        img = always_redraw(mk)
        self.add(img)
        circ = Circle(radius=a / L * H / 2, color=WHITE, stroke_width=2).move_to([-3.0, -0.55, 0])
        zl = Arrow([-6.2, -3.75, 0], [-4.6, -3.75, 0], buff=0, color=WHITE, stroke_width=3)
        zt = CT(r"入射方向 $z$", 20).next_to(zl, RIGHT, buff=0.1)
        self.add(circ, zl, zt)
        txt = panel(r"入射：平面波 $e^{\mathrm ikz}$，粒子流 $N=v$",
                    r"散射：外向球面波 $f(\theta,\varphi)\dfrac{e^{\mathrm ikr}}{r}$",
                    r"径向流 $J_r=\dfrac{v}{r^2}|f|^2$，穿过 $\mathrm dS=r^2\mathrm d\Omega$",
                    size=24, w=6.0)
        txt.to_edge(RIGHT, buff=0.35).shift(UP * 1.4)
        self.play(FadeIn(txt[0]))
        self.play(phase.animate.set_value(4 * PI), run_time=5, rate_func=linear)
        lab = CT(r"只看散射波", 24, C_HL).move_to([-3.0, 2.75, 0])
        mode["which"] = sca
        mode["scale"] = 0.7
        self.play(FadeIn(txt[1]), FadeIn(lab))
        self.play(phase.animate.set_value(8 * PI), run_time=5, rate_func=linear)
        # detector cone
        c0 = np.array([-3.0, -0.55, 0])
        th = 50 * DEGREES
        r1 = 2.7
        cone = Polygon(c0, c0 + r1 * np.array([np.cos(th - 0.12), np.sin(th - 0.12), 0]), c0 + r1 * np.array([np.cos(th + 0.12), np.sin(th + 0.12), 0]),
                       color=C_A, fill_opacity=0.35, stroke_width=2)
        arc = Arc(radius=0.9, start_angle=0, angle=th, arc_center=c0, color=C_A)
        tl = MT(r"\theta", 26, C_A).move_to(c0 + 1.15 * np.array([np.cos(th / 2), np.sin(th / 2), 0]))
        dl = MT(r"\mathrm d\Omega", 26, C_A).move_to(c0 + (r1 + 0.35) * np.array([np.cos(th), np.sin(th), 0]))
        self.play(FadeIn(cone), Create(arc), FadeIn(tl), FadeIn(dl), FadeIn(txt[2]))
        res = boxed(VGroup(
            MT(r"\mathrm dn=v|f(\theta,\varphi)|^2\,\mathrm d\Omega=q(\theta,\varphi)\,N\,\mathrm d\Omega", 30),
            MT(r"\Rightarrow\ q(\theta,\varphi)=|f(\theta,\varphi)|^2,\qquad Q=\int q\,\mathrm d\Omega", 30, C_HL)).arrange(DOWN, buff=0.2))
        fitw(res, 6.0).next_to(txt, DOWN, buff=0.4).to_edge(RIGHT, buff=0.3)
        self.play(FadeIn(res))
        self.play(phase.animate.set_value(12 * PI), run_time=5, rate_func=linear)
        self.wait(1)


class C6_PhaseShift(Scene):
    """s 波相移：势阱把径向波函数“拉进来”，外部只表现为相位移动 δ0。"""

    def construct(self):
        title = title_bar(r"分波法的核心：势场只让每个分波在远处多出一个相移 $\delta_l$（以 s 波、方势阱为例）")
        self.play(FadeIn(title))
        ax = Axes(x_range=[0, 4.0, 1], y_range=[-1.8, 1.8, 1], x_length=8.2, y_length=4.4,
                  axis_config={"include_tip": False, "include_ticks": False}).to_edge(LEFT, buff=0.4).shift(DOWN * 0.2)
        a = 1.0
        k = 1.6
        well = Rectangle(width=ax.c2p(a, 0)[0] - ax.c2p(0, 0)[0], height=4.4, fill_color=BLUE_E, fill_opacity=0.35, stroke_width=0)
        well.move_to(ax.c2p(a / 2, 0))
        wl = CT(r"势阱 $r<a$", 20, BLUE_B).next_to(well, UP, buff=0.05)
        rl = MT(r"r", 26).next_to(ax.x_axis, RIGHT, buff=0.1)
        al = MT(r"a", 24).next_to(ax.c2p(a, -1.8), DOWN, buff=0.05)
        self.play(Create(ax), FadeIn(well), FadeIn(wl), FadeIn(rl), FadeIn(al))
        k0a = ValueTracker(0.0)

        def kp():
            return np.sqrt(k ** 2 + (k0a.get_value() / a) ** 2)

        def delta():
            q = kp()
            phi = np.arctan2(k * np.sin(q * a), q * np.cos(q * a))
            phi += 2 * np.pi * np.round((q * a - phi) / (2 * np.pi))
            return phi - k * a

        def u(r):
            q = kp()
            d = delta()
            A = np.sin(k * a + d) / np.sin(q * a) if abs(np.sin(q * a)) > 1e-6 else 1.0
            return np.where(r < a, A * np.sin(q * r), np.sin(k * r + d))

        free = ax.plot(lambda r: np.sin(k * r), x_range=[0, 4, 0.01], color=GREY_B, stroke_width=2.5)
        free = DashedVMobject(free, num_dashes=60)
        fl = CT(r"无势场：$u_0=\sin kr$", 20, GREY_B).next_to(ax, UP, buff=0.3).align_to(ax, RIGHT)
        self.play(Create(free), FadeIn(fl))
        curve = always_redraw(lambda: ax.plot(lambda r: float(u(np.array(r))), x_range=[0.001, 4, 0.01], color=C_A, stroke_width=4, use_smoothing=False))
        cl = CT(r"有势场：$u_0=\sin(kr+\delta_0)$（$r>a$）", 20, C_A).next_to(fl, DOWN, buff=0.06, aligned_edge=RIGHT)
        self.play(Create(curve), FadeIn(cl))
        eqs = panel(r"$k'=\sqrt{k^2+k_0^2},\quad k_0^2=\dfrac{2m|U_0|}{\hbar^2}$",
                    r"$r=a$ 处 $u'/u$ 连续：",
                    MT(r"k\cot(ka+\delta_0)=k'\cot k'a", 28),
                    MT(r"Q_0=\frac{4\pi}{k^2}\sin^2\delta_0", 30, C_HL),
                    size=24, w=4.8)
        eqs.to_edge(RIGHT, buff=0.35).shift(UP * 1.3)
        self.play(FadeIn(eqs))
        dval = always_redraw(lambda: MT(r"k_0a=%.2f,\ \ \delta_0=%.2f\ \mathrm{rad}" % (k0a.get_value(), delta()), 28, C_HL).next_to(eqs, DOWN, buff=0.35))
        # bar for sin^2 delta
        bx = np.array([3.4, -1.75, 0])
        frame = Rectangle(width=3.0, height=0.35, color=GREY_B).move_to(bx + RIGHT * 1.5)
        bar = always_redraw(lambda: Rectangle(width=max(3.0 * np.sin(delta()) ** 2, 0.001), height=0.35,
                                              fill_color=C_HL, fill_opacity=0.9, stroke_width=0).move_to(bx, aligned_edge=LEFT))
        bl = MT(r"\sin^2\delta_0", 24).next_to(frame, UP, buff=0.08)
        self.add(dval)
        self.play(FadeIn(frame), FadeIn(bar), FadeIn(bl))
        self.play(k0a.animate.set_value(1.35), run_time=4, rate_func=smooth)
        self.wait(0.5)
        self.play(k0a.animate.set_value(2.6), run_time=4, rate_func=smooth)
        msg = boxed(CT(r"内部越“深”，波被拉得越紧，外部相位越超前；$\delta_0$ 与 $\delta_0\pm\pi$ 物理上不可区分，截面只认 $\sin^2\delta_0$", 22, C_HL))
        fitw(msg, 13.0).to_edge(DOWN, buff=0.12)
        self.play(FadeIn(msg))
        self.play(k0a.animate.set_value(0.6), run_time=4, rate_func=smooth)
        self.wait(1.5)


class C6_HardSphere(Scene):
    """硬球散射：各分波截面随 ka 变化；低能 4πa²，高能 2πa²。"""

    def construct(self):
        title = title_bar(r"硬球散射的分波截面：只有 $l\lesssim ka$ 的分波被散射；低能 $Q\to4\pi a^2$（经典值的 4 倍）")
        self.play(FadeIn(title))
        lmax = 9

        def Ql(ka):
            ls = np.arange(lmax + 1)
            j = spherical_jn(ls, ka)
            y = spherical_yn(ls, ka)
            s2 = j ** 2 / (j ** 2 + y ** 2)
            return 4 * (2 * ls + 1) * s2 / ka ** 2  # in units of pi a^2

        ax = Axes(x_range=[-0.5, lmax + 0.5, 1], y_range=[0, 4.4, 1], x_length=6.4, y_length=4.6,
                  axis_config={"include_tip": False, "font_size": 18}).to_edge(LEFT, buff=0.6).shift(DOWN * 0.3)
        ax.y_axis.add_numbers(font_size=18)
        for l in range(lmax + 1):
            ax.add(MT(str(l), 20).next_to(ax.c2p(l, 0), DOWN, buff=0.12))
        xl = CT(r"分波 $l$", 22).next_to(ax.x_axis, DOWN, buff=0.45)
        yl = CT(r"$Q_l/\pi a^2$", 22).next_to(ax.y_axis, UP, buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))
        ka = ValueTracker(0.15)

        def bars():
            q = Ql(ka.get_value())
            g = VGroup()
            for l in range(lmax + 1):
                h = min(q[l], 4.4)
                r = Rectangle(width=0.42, height=max(ax.c2p(0, h)[1] - ax.c2p(0, 0)[1], 0.001), fill_color=interpolate_color(C_A, C_UP, l / lmax), fill_opacity=0.9, stroke_width=0)
                r.move_to(ax.c2p(l, 0), aligned_edge=DOWN)
                g.add(r)
            return g
        bb = always_redraw(bars)
        kl = always_redraw(lambda: MT(r"ka=%.2f" % ka.get_value(), 30, C_HL).next_to(ax, UP, buff=0.05).align_to(ax, RIGHT))
        lline = always_redraw(lambda: DashedLine(ax.c2p(min(ka.get_value(), lmax + 0.4), 0), ax.c2p(min(ka.get_value(), lmax + 0.4), 4.3), color=WHITE, stroke_width=2))
        self.add(bb, kl, lline)
        ax2 = Axes(x_range=[0, 8, 2], y_range=[0, 4.5, 1], x_length=5.0, y_length=3.4,
                   axis_config={"include_tip": False, "font_size": 18}).to_edge(RIGHT, buff=0.5).shift(DOWN * 0.9)
        ax2.add_coordinates(font_size=18)
        x2 = MT(r"ka", 24).next_to(ax2.x_axis, DOWN, buff=0.4)
        y2 = CT(r"总截面 $Q/\pi a^2$", 22).next_to(ax2, UP, buff=0.1)
        tot = ax2.plot(lambda x: float(np.sum(Ql(x))) if x > 0.05 else 4.0, x_range=[0.05, 8, 0.02], color=C_B, stroke_width=3)
        h4 = DashedLine(ax2.c2p(0, 4), ax2.c2p(8, 4), color=GREY_B, stroke_width=1.5)
        h2 = DashedLine(ax2.c2p(0, 2), ax2.c2p(8, 2), color=GREY_B, stroke_width=1.5)
        h1 = DashedLine(ax2.c2p(0, 1), ax2.c2p(8, 1), color=C_UP, stroke_width=1.5)
        l1 = CT(r"经典 $\pi a^2$", 18, C_UP).next_to(ax2.c2p(8, 1), UP, buff=0.05).align_to(ax2.c2p(8, 1), RIGHT)
        self.play(Create(ax2), FadeIn(x2), FadeIn(y2), Create(h4), Create(h2), Create(h1), FadeIn(l1))
        dot = always_redraw(lambda: Dot(ax2.c2p(ka.get_value(), float(np.sum(Ql(ka.get_value())))), color=C_HL, radius=0.07))
        self.play(Create(tot), FadeIn(dot), run_time=1.5)
        rule = CT(r"$\tan\delta_l=\dfrac{j_l(ka)}{n_l(ka)}$，$\ Q_l=\dfrac{4\pi}{k^2}(2l+1)\sin^2\delta_l$", 22).next_to(y2, UP, buff=0.3)
        fitw(rule, 5.2)
        self.play(FadeIn(rule))
        self.wait(0.6)
        self.play(ka.animate.set_value(7.5), run_time=9, rate_func=linear)
        msg = boxed(CT(r"低能：只有 s 波，$Q=4\pi a^2$（衍射，整个球面都在散射）；高能：$Q\to2\pi a^2$ = 几何遮挡 $\pi a^2$ + 衍射“影子” $\pi a^2$", 22, C_HL))
        fitw(msg, 13.0).to_edge(DOWN, buff=0.12)
        self.play(FadeIn(msg))
        self.wait(3)


class C6_Born(Scene):
    """玻恩近似：散射振幅 = 势能对动量转移 K 的傅里叶变换。"""

    def construct(self):
        title = title_bar(r"玻恩近似：$f(\theta)=-\dfrac{m}{2\pi\hbar^2}\displaystyle\int U(\boldsymbol r)e^{-\mathrm i\boldsymbol K\cdot\boldsymbol r}\mathrm d^3r$，$\boldsymbol K=\boldsymbol k'-\boldsymbol k$")
        self.play(FadeIn(title))
        c0 = np.array([-4.6, -0.3, 0])
        R = 2.2
        circ = Circle(radius=R, color=GREY_B, stroke_width=1.5).move_to(c0)
        kv = Arrow(c0 - np.array([R, 0, 0]), c0, buff=0, color=C_A, stroke_width=5)
        kl = MT(r"\boldsymbol k", 28, C_A).next_to(kv, DOWN, buff=0.1)
        th = ValueTracker(40 * DEGREES)
        kp = always_redraw(lambda: Arrow(c0, c0 + R * np.array([np.cos(th.get_value()), np.sin(th.get_value()), 0]), buff=0, color=C_B, stroke_width=5))
        kpl = always_redraw(lambda: MT(r"\boldsymbol k'", 28, C_B).move_to(c0 + (R + 0.35) * np.array([np.cos(th.get_value()), np.sin(th.get_value()), 0])))
        K = always_redraw(lambda: Arrow(c0 + R * np.array([1, 0, 0]), c0 + R * np.array([np.cos(th.get_value()), np.sin(th.get_value()), 0]), buff=0, color=C_HL, stroke_width=5))
        kk = always_redraw(lambda: Line(c0, c0 + R * np.array([1, 0, 0]), color=C_A, stroke_width=2).set_opacity(0.5))
        Kl = always_redraw(lambda: MT(r"\boldsymbol K", 28, C_HL).next_to(K, RIGHT, buff=0.08))
        arc = always_redraw(lambda: Arc(radius=0.55, start_angle=0, angle=th.get_value(), arc_center=c0, color=WHITE))
        tl = always_redraw(lambda: MT(r"\theta", 24).move_to(c0 + 0.8 * np.array([np.cos(th.get_value() / 2), np.sin(th.get_value() / 2), 0])))
        self.play(Create(circ), GrowArrow(kv), FadeIn(kl))
        self.add(kk, kp, kpl, K, Kl, arc, tl)
        Kf = CT(r"弹性：$|\boldsymbol k'|=|\boldsymbol k|$，$\ K=2k\sin\dfrac{\theta}{2}$", 26).next_to(circ, DOWN, buff=0.35)
        self.play(FadeIn(Kf))
        self.play(th.animate.set_value(150 * DEGREES), run_time=3)
        self.play(th.animate.set_value(60 * DEGREES), run_time=2)
        # right: screened Coulomb
        ax = Axes(x_range=[0, 180, 30], y_range=[0, 4.3, 1], x_length=5.8, y_length=3.3,
                  axis_config={"include_tip": False, "font_size": 18}).to_edge(RIGHT, buff=0.5).to_edge(DOWN, buff=1.25)
        ax.x_axis.add_numbers(font_size=18)
        for v in [0, -1, -2, -3, -4]:
            ax.add(MT(r"10^{%d}" % v, 18).next_to(ax.c2p(0, v + 4), LEFT, buff=0.1))
        xl = CT(r"散射角 $\theta$ / $^\circ$", 20).next_to(ax.x_axis, DOWN, buff=0.4)
        yl = CT(r"$q(\theta)/q(0)$（对数）", 20).next_to(ax, UP, buff=0.1).align_to(ax, LEFT)
        U = panel(r"屏蔽库仑势 $U=-\dfrac{ZZ'e^2}{4\pi\varepsilon_0r}e^{-r/a}$",
                  MT(r"\tilde U(K)\propto\frac{1}{K^2+a^{-2}}\ \Rightarrow\ q\propto\frac{1}{\left(4k^2\sin^2\frac{\theta}{2}+a^{-2}\right)^2}", 28),
                  size=22, w=5.8)
        U.next_to(yl, UP, buff=0.15).align_to(ax, LEFT)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl), FadeIn(U))
        ka = ValueTracker(0.5)

        def lq(tdeg, kav):
            s = np.sin(np.radians(tdeg) / 2)
            return np.log10(1.0 / (4 * kav ** 2 * s ** 2 + 1) ** 2) + 4
        cur = always_redraw(lambda: ax.plot(lambda t: max(lq(t, ka.get_value()), 0), x_range=[0.5, 180, 0.5], color=C_A, stroke_width=3.5))
        kal = always_redraw(lambda: MT(r"ka=%.1f" % ka.get_value(), 28, C_A).move_to(ax.c2p(150, 3.7)))
        self.play(Create(cur), FadeIn(kal))
        ruth = always_redraw(lambda: DashedVMobject(ax.plot(lambda t: max(min(np.log10(1.0 / (4 * ka.get_value() ** 2 * np.sin(np.radians(t) / 2) ** 2) ** 2) + 4, 4.3), 0), x_range=[2, 180, 0.5], color=C_UP, stroke_width=2.5), num_dashes=50))
        rl = CT(r"卢瑟福 $\propto\csc^4\frac{\theta}{2}$", 20, C_UP).move_to(ax.c2p(140, 2.9))
        self.play(FadeIn(ruth), FadeIn(rl))
        self.play(ka.animate.set_value(4.0), run_time=4)
        msg = boxed(CT(r"势越“短程”，傅里叶谱越宽，散射越各向同性；$Ka\gg1$（大角、高能）时屏蔽无关紧要，回到卢瑟福公式", 22, C_HL))
        fitw(msg, 13.0).to_edge(DOWN, buff=0.12)
        self.play(FadeIn(msg))
        self.wait(2.5)
