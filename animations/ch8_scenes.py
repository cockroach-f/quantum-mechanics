"""第 8 章 全同粒子 —— manim 动画（Manim Community + LaTeX）。

渲染：manim -qm ch8_scenes.py C8_ExchangeHole  （输出改名为 docs/assets/videos/ch8_*.mp4）
"""
import numpy as np
from common import *  # noqa: F401,F403

C_D = GREY_A     # 可区分
C_S = BLUE_C     # 对称
C_A = RED_C      # 反对称


def psi(n, x):
    return np.sqrt(2.0) * np.sin(n * np.pi * x)   # a = 1


def heat_rgb(z):
    """z in [0,1] -> RGB uint8 (black -> purple -> orange -> yellow)."""
    stops = np.array([[0, 0, 0], [70, 20, 110], [200, 60, 60], [250, 170, 40], [255, 250, 200]], float)
    t = np.clip(z, 0, 1) * (len(stops) - 1)
    i = np.minimum(t.astype(int), len(stops) - 2)
    f = (t - i)[..., None]
    return (stops[i] * (1 - f) + stops[i + 1] * f).astype(np.uint8)


def density_image(kind, N=220):
    x = (np.arange(N) + 0.5) / N
    X1, X2 = np.meshgrid(x, x)           # 列 -> x1，行 -> x2
    a = psi(1, X1) * psi(2, X2)
    b = psi(2, X1) * psi(1, X2)
    if kind == "D":
        rho = a ** 2
    elif kind == "S":
        rho = 0.5 * (a + b) ** 2
    else:
        rho = 0.5 * (a - b) ** 2
    z = rho / 4.8                          # 三种情形共用同一色标（最大值约 4.74）
    img = heat_rgb(np.sqrt(z))[::-1]       # 行 0 在上方 -> 翻转使 x2 向上增大
    return ImageMobject(img)


class C8_ExchangeHole(Scene):
    """一维无限深势阱中两粒子（n=1 与 n=2）：可区分 / 对称 / 反对称 的联合概率密度。"""

    def construct(self):
        title = title_bar(r"两粒子处于 $\psi_1,\psi_2$：联合概率密度 $|\Psi(x_1,x_2)|^2$")
        self.play(FadeIn(title))
        labels = [
            ("D", r"可区分：$\psi_1(x_1)\psi_2(x_2)$", C_D),
            ("S", r"对称：$\tfrac{1}{\sqrt2}(\psi_1\psi_2+\psi_2\psi_1)$", C_S),
            ("A", r"反对称：$\tfrac{1}{\sqrt2}(\psi_1\psi_2-\psi_2\psi_1)$", C_A),
        ]
        panels = Group()
        for k, lab, col in labels:
            im = density_image(k)
            im.height = 3.0
            frame = Square(3.0, color=col, stroke_width=3).move_to(im)
            diag = DashedLine(frame.get_corner(DL), frame.get_corner(UR), color=WHITE, stroke_width=2, dash_length=0.08)
            cap = fitw(CT(lab, 24, col), 3.9).next_to(frame, UP, buff=0.15)
            ax1 = MT(r"x_1", 22).next_to(frame, DOWN, buff=0.08)
            ax2 = MT(r"x_2", 22).next_to(frame, LEFT, buff=0.08)
            panels.add(Group(im, frame, diag, cap, ax1, ax2))
        panels.arrange(RIGHT, buff=0.6).next_to(title, DOWN, buff=0.35)
        for p in panels:
            self.play(FadeIn(p[0]), Create(p[1]), Write(p[3]), FadeIn(p[4]), FadeIn(p[5]), run_time=1.1)
        self.wait(0.5)
        self.play(*[Create(p[2]) for p in panels], run_time=1.2)
        note = CT(r"虚线为对角线 $x_1=x_2$（两粒子在同一位置）", 26).next_to(panels, DOWN, buff=0.3)
        self.play(FadeIn(note))
        self.wait(1)
        hole = SurroundingRectangle(panels[2][2], color=C_A, buff=0.02)
        bunch = SurroundingRectangle(panels[1][2], color=C_S, buff=0.02)
        t1 = CT(r"反对称：对角线上处处为 0 —— \textbf{交换空穴}", 26, C_A)
        t2 = CT(r"对称：概率向对角线\textbf{聚集}", 26, C_S)
        VGroup(t2, t1).arrange(DOWN, aligned_edge=LEFT, buff=0.12).next_to(panels, DOWN, buff=0.25)
        self.play(FadeOut(note), Create(hole), FadeIn(t1))
        self.wait(1)
        self.play(Create(bunch), FadeIn(t2))
        self.wait(2)

        # ---- 第二幕：定量 —— 对角线上的密度与平均间距
        self.play(FadeOut(Group(panels, hole, bunch, t1, t2)))
        ax = Axes(x_range=[0, 1, 0.25], y_range=[0, 5, 1], x_length=6.2, y_length=3.8,
                  axis_config={"include_tip": False, "font_size": 22}).to_edge(LEFT, buff=0.6).shift(DOWN * 0.4)
        xl = CT(r"$x_1=x_2=x$（单位 $a$）", 22).next_to(ax.x_axis, DOWN, buff=0.3)
        yl = MT(r"|\Psi(x,x)|^2", 24).next_to(ax.y_axis, UP, buff=0.12)
        xs = np.linspace(0, 1, 301)
        d = (psi(1, xs) * psi(2, xs)) ** 2
        cD = curve(ax, xs, d, C_D, 3)
        cS = curve(ax, xs, 2 * d, C_S, 4)
        cA = curve(ax, xs, 0 * xs, C_A, 5)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))
        self.play(Create(cD), run_time=1.2)
        self.play(Create(cS), run_time=1.2)
        self.play(Create(cA), run_time=1.2)
        leg = VGroup(
            CT(r"可区分", 24, C_D), CT(r"对称：恰为 2 倍", 24, C_S), CT(r"反对称：恒为 0", 24, C_A)
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.1).next_to(ax, UP, buff=0.05).shift(RIGHT * 0.8)
        self.play(FadeIn(leg))
        tab = VGroup(
            MT(r"\langle (x_1-x_2)^2\rangle=\langle x^2\rangle_1+\langle x^2\rangle_2-2\langle x\rangle_1\langle x\rangle_2\mp 2|x_{12}|^2", 24),
            MT(r"x_{12}=\int_0^a x\,\psi_1\psi_2\,\mathrm dx=-\frac{16a}{9\pi^2}", 26),
            CT(r"可区分：$\sqrt{\langle (x_1-x_2)^2\rangle}=0.32\,a$", 26, C_D),
            CT(r"对称：$\sqrt{\langle (x_1-x_2)^2\rangle}=0.20\,a$（更近）", 26, C_S),
            CT(r"反对称：$\sqrt{\langle (x_1-x_2)^2\rangle}=0.41\,a$（更远）", 26, C_A),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.22)
        fitw(tab, 6.2).next_to(ax, RIGHT, buff=0.35).shift(UP * 0.2)
        for m in tab:
            self.play(Write(m), run_time=0.9)
        concl = boxed(CT(r"这不是“力”，而是波函数对称性带来的统计关联：交换力", 26, C_HL)).to_edge(DOWN, buff=0.2)
        self.play(FadeIn(concl))
        self.wait(2.5)


class C8_Filling(Scene):
    """一维势阱中逐个加入粒子：零自旋玻色子 vs 自旋 1/2 费米子。"""

    def construct(self):
        title = title_bar(r"往势阱里逐个放粒子：玻色子 vs 费米子（$E_n=n^2E_1$）")
        self.play(FadeIn(title))
        ymap = lambda n: -2.9 + 0.2 * n * n   # n^2 比例
        def ladder(cx, col):
            g = VGroup()
            for n in range(1, 6):
                ln = Line([cx - 2.0, ymap(n), 0], [cx + 2.0, ymap(n), 0], color=C_E, stroke_width=3)
                lab = MT(rf"n={n}", 20, GREY_B).next_to(ln, LEFT, buff=0.12)
                g.add(VGroup(ln, lab))
            return g
        L = ladder(-3.3, C_S)
        R = ladder(3.5, C_A)
        hL = CT(r"玻色子（自旋 0）", 28, C_S).move_to([-3.3, 2.6, 0])
        hR = CT(r"费米子（自旋 $\tfrac12$）", 28, C_A).move_to([3.5, 2.6, 0])
        self.play(FadeIn(L), FadeIn(R), FadeIn(hL), FadeIn(hR))
        cntL = MT(r"E=0", 30, C_S).move_to([-3.3, 1.95, 0])
        cntR = MT(r"E=0", 30, C_A).move_to([3.5, 1.95, 0])
        nlab = MT(r"N=0", 30).move_to([0.1, 1.95, 0])
        self.play(FadeIn(cntL), FadeIn(cntR), FadeIn(nlab))
        bos = VGroup()
        EB, EF = 0, 0
        for N in range(1, 7):
            # boson: all into n=1
            d = Dot(radius=0.11, color=C_S).move_to([-3.3 - 1.25 + 0.5 * (N - 1), ymap(1) + 0.16, 0])
            # fermion: level n = ceil(N/2), spin up then down
            n = (N + 1) // 2
            up = (N % 2 == 1)
            x = 3.5 + (-0.35 if up else 0.35)
            y0 = ymap(n)
            arr = Arrow([x, y0 - 0.02 if up else y0 + 0.42, 0], [x, y0 + 0.42 if up else y0 - 0.02, 0],
                        buff=0, color=C_A, stroke_width=5, max_tip_length_to_length_ratio=0.35)
            EB += 1
            EF += n * n
            self.play(FadeIn(d, shift=DOWN * 0.3), GrowArrow(arr),
                      Transform(nlab, MT(rf"N={N}", 30).move_to(nlab)),
                      Transform(cntL, MT(rf"E={EB}E_1", 30, C_S).move_to(cntL)),
                      Transform(cntR, MT(rf"E={EF}E_1", 30, C_A).move_to(cntR)), run_time=0.9)
            if N == 3:
                warn = CT(r"$n=1$ 已有 $\uparrow\downarrow$，第 3 个费米子只能上 $n=2$", 24, C_HL).next_to(R, DOWN, buff=0.15)
                self.play(FadeIn(warn))
                self.wait(1)
                self.play(FadeOut(warn))
            self.wait(0.3)
        s1 = CT(r"玻色子：对称波函数不限制占据，可以\textbf{全部挤在基态}", 26, C_S)
        s2 = CT(r"费米子：斯莱特行列式两行相同即为 0 —— \textbf{泡利不相容}，逐级向上填", 26, C_A)
        VGroup(s1, s2).arrange(DOWN, buff=0.12).to_edge(DOWN, buff=0.15)
        bg = BackgroundRectangle(VGroup(s1, s2), fill_opacity=0.85, buff=0.08)
        self.play(FadeIn(bg), FadeIn(s1), FadeIn(s2))
        self.wait(2.5)


def ket(a, b, size=34, color=WHITE):
    ar = {1: r"\uparrow", -1: r"\downarrow"}
    return MT(rf"\ket{{{ar[a]}{ar[b]}}}", size, color)


class C8_SingletTriplet(Scene):
    """两电子自旋：四个乘积态 -> 交换算符 -> 三重态 + 单态。"""

    def construct(self):
        title = title_bar(r"两个电子的自旋：$2\times2=4$ 个态如何重新组合")
        self.play(FadeIn(title))
        prods = VGroup(ket(1, 1), ket(1, -1), ket(-1, 1), ket(-1, -1)).arrange(DOWN, buff=0.45).shift(LEFT * 4.6 + DOWN * 0.3)
        hdr = CT(r"乘积态 $\chi_{\alpha}(s_{1z})\chi_{\beta}(s_{2z})$", 24).next_to(prods, UP, buff=0.35)
        self.play(FadeIn(hdr), LaggedStart(*[FadeIn(p) for p in prods], lag_ratio=0.2))
        note = CT(r"约定 $\ket{\uparrow\downarrow}=\chi_+(s_{1z})\chi_-(s_{2z})$：左边是电子 1", 22, GREY_B).to_edge(DOWN, buff=0.3)
        self.play(FadeIn(note))
        self.wait(0.8)
        # 交换算符作用
        P = MT(r"\hat P_{12}", 34, C_HL).move_to(prods.get_right() + RIGHT * 1.15 + UP * 1.0 + UP * prods.height / 2)
        self.play(Write(P))
        res = VGroup(ket(1, 1), ket(-1, 1), ket(1, -1), ket(-1, -1))
        for r, p in zip(res, prods):
            r.next_to(p, RIGHT, buff=2.0)
        arrows = VGroup(*[Arrow(p.get_right(), r.get_left(), buff=0.15, stroke_width=3, color=GREY_B) for p, r in zip(prods, res)])
        self.play(LaggedStart(*[GrowArrow(a) for a in arrows], lag_ratio=0.15), LaggedStart(*[FadeIn(r) for r in res], lag_ratio=0.15))
        c1 = CT(r"不变", 22, C_S).next_to(res[0], RIGHT)
        c4 = CT(r"不变", 22, C_S).next_to(res[3], RIGHT)
        c23 = CT(r"互换：都不是 $\hat P_{12}$ 的本征态", 22, C_A).next_to(VGroup(res[1], res[2]), RIGHT)
        self.play(FadeIn(c1), FadeIn(c4), FadeIn(c23))
        self.wait(1.5)
        self.play(FadeOut(VGroup(P, res, arrows, c1, c4, c23, note)))
        # 组合
        T = VGroup(
            MT(r"\ket{1,1}=\ket{\uparrow\uparrow}", 32, C_S),
            MT(r"\ket{1,0}=\tfrac{1}{\sqrt2}\big(\ket{\uparrow\downarrow}+\ket{\downarrow\uparrow}\big)", 32, C_S),
            MT(r"\ket{1,-1}=\ket{\downarrow\downarrow}", 32, C_S),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
        Sg = MT(r"\ket{0,0}=\tfrac{1}{\sqrt2}\big(\ket{\uparrow\downarrow}-\ket{\downarrow\uparrow}\big)", 32, C_A)
        right = VGroup(T, Sg).arrange(DOWN, aligned_edge=LEFT, buff=0.6).shift(RIGHT * 1.3 + DOWN * 0.1)
        self.play(TransformFromCopy(prods[0], T[0]))
        self.play(TransformFromCopy(VGroup(prods[1], prods[2]), T[1]))
        self.play(TransformFromCopy(prods[3], T[2]))
        self.play(TransformFromCopy(VGroup(prods[1], prods[2]), Sg))
        bT = Brace(T, RIGHT, color=C_S)
        bS = Brace(Sg, RIGHT, color=C_A)
        lT = VGroup(CT(r"三重态 $S=1$", 24, C_S), MT(r"\hat S^2=2\hbar^2,\ \hat P_{12}=+1", 24, C_S)).arrange(DOWN, buff=0.08).next_to(bT, RIGHT, buff=0.1)
        lS = VGroup(CT(r"单态 $S=0$", 24, C_A), MT(r"\hat S^2=0,\ \hat P_{12}=-1", 24, C_A)).arrange(DOWN, buff=0.08).next_to(bS, RIGHT, buff=0.1)
        grp = VGroup(right, bT, bS, lT, lS)
        if grp.get_right()[0] > 6.9:
            grp.shift(LEFT * (grp.get_right()[0] - 6.9))
        self.play(GrowFromCenter(bT), FadeIn(lT))
        self.play(GrowFromCenter(bS), FadeIn(lS))
        self.wait(1.5)
        key = boxed(VGroup(
            MT(r"\hat{\boldsymbol S}_1\cdot\hat{\boldsymbol S}_2=\tfrac12\big[\hat S^2-\tfrac32\hbar^2\big]=\begin{cases}+\tfrac14\hbar^2,&S=1\\-\tfrac34\hbar^2,&S=0\end{cases}", 28),
            CT(r"总波函数反对称：$\Phi=\psi_{\mathrm S}(\boldsymbol r_1,\boldsymbol r_2)\chi_{\mathrm A}$ 或 $\psi_{\mathrm A}(\boldsymbol r_1,\boldsymbol r_2)\chi_{\mathrm S}$", 26, C_HL),
        ).arrange(DOWN, buff=0.18)).to_edge(DOWN, buff=0.12)
        self.play(FadeOut(hdr), grp.animate.shift(UP * 0.55), prods.animate.shift(UP * 0.55))
        self.play(FadeIn(key))
        self.wait(3)


class C8_HeliumLevels(Scene):
    """氦原子：零级能级 -> 库仑能 K -> 交换能 ±J，仲氦（单态）与正氦（三重态）。"""

    def construct(self):
        title = title_bar(r"氦原子能级：库仑能 $K$ 整体抬高，交换能 $J$ 把单态与三重态劈开")
        self.play(FadeIn(title))
        xs = [-4.6, -0.9, 2.9]
        heads = [CT(r"无相互作用", 24), CT(r"仲氦（单态 $\chi_{\mathrm A}$）", 24, C_S), CT(r"正氦（三重态 $\chi_{\mathrm S}$）", 24, C_A)]
        for h, x in zip(heads, xs):
            h.move_to([x, 2.55, 0])
        self.play(*[FadeIn(h) for h in heads])
        def lvl(x, y, col, lab, side=UP):
            ln = Line([x - 1.25, y, 0], [x + 1.25, y, 0], color=col, stroke_width=4)
            t = MT(lab, 24, col).next_to(ln, side, buff=0.08)
            return VGroup(ln, t)
        g0 = lvl(xs[0], -1.75, WHITE, r"2\epsilon_1")
        e0 = lvl(xs[0], 0.8, WHITE, r"\epsilon_1+\epsilon_2")
        self.play(Create(g0), Create(e0))
        self.wait(0.5)
        g1 = lvl(xs[1], -0.95, C_S, r"2\epsilon_1+K")
        dg = DashedLine(g0[0].get_right(), g1[0].get_left(), color=GREY_B)
        self.play(Create(dg), Create(g1))
        n1 = CT(r"基态：两电子同处 1s，空间对称 $\Rightarrow$ 只能是单态", 22, GREY_A).next_to(g1, DOWN, buff=0.15).shift(RIGHT * 1.3)
        self.play(FadeIn(n1))
        mid = 1.45
        Kl = DashedLine([xs[0] + 1.25, 0.8, 0], [xs[1] - 1.25, mid, 0], color=GREY_B)
        p1 = lvl(xs[1], mid + 0.55, C_S, r"\epsilon_1+\epsilon_2+K+J")
        o1 = lvl(xs[2], mid - 0.55, C_A, r"\epsilon_1+\epsilon_2+K-J", DOWN)
        kref = DashedLine([xs[1] - 1.25, mid, 0], [xs[2] + 1.25, mid, 0], color=GREY_D, dash_length=0.06)
        kt = MT(r"\epsilon_1+\epsilon_2+K", 20, GREY_B).move_to([(xs[1] + xs[2]) / 2, mid + 0.22, 0])
        self.play(Create(Kl), Create(kref), FadeIn(kt))
        self.play(Create(p1), Create(o1))
        bJ1 = DoubleArrow([xs[2] + 1.45, mid, 0], [xs[2] + 1.45, mid - 0.55, 0], buff=0, color=C_HL, tip_length=0.12)
        bJ2 = DoubleArrow([xs[1] - 1.45, mid, 0], [xs[1] - 1.45, mid + 0.55, 0], buff=0, color=C_HL, tip_length=0.12)
        self.play(GrowFromCenter(bJ1), GrowFromCenter(bJ2), FadeIn(MT("J", 22, C_HL).next_to(bJ2, LEFT, buff=0.05)))
        deg = VGroup(CT(r"简并度 1", 20, C_S).next_to(p1, DOWN, buff=0.08),
                     CT(r"简并度 3", 20, C_A).next_to(o1, UP, buff=0.08).shift(LEFT * 0.0))
        self.play(FadeIn(deg))
        self.wait(1)
        box = boxed(VGroup(
            MT(r"K=\iint\frac{e^2|\psi_n(\boldsymbol r_1)|^2|\psi_m(\boldsymbol r_2)|^2}{4\pi\varepsilon_0 r_{12}}\,\mathrm dV_1\mathrm dV_2,\quad J=\iint\frac{e^2\psi_n^*(\boldsymbol r_1)\psi_m^*(\boldsymbol r_2)\psi_m(\boldsymbol r_1)\psi_n(\boldsymbol r_2)}{4\pi\varepsilon_0 r_{12}}\,\mathrm dV_1\mathrm dV_2", 22),
            CT(r"$J>0$：空间反对称（三重态）带交换空穴，电子相距更远，排斥能更小 $\Rightarrow$ 正氦更低", 22, C_HL),
        ).arrange(DOWN, buff=0.14))
        fitw(box, 13.2).to_edge(DOWN, buff=0.12)
        self.play(FadeOut(n1), FadeIn(box))
        self.wait(3)
