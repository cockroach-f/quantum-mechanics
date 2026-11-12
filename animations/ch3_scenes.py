"""第 3 章 量子力学中的力学量 —— manim 动画源码（LaTeX 排版）."""
from common import *  # noqa: F401,F403

S3, S6, S30 = np.sqrt(3), np.sqrt(6), np.sqrt(30)
RNL = {
    (1, 0): lambda r: 2 * np.exp(-r),
    (2, 0): lambda r: (1 / np.sqrt(2)) * (1 - r / 2) * np.exp(-r / 2),
    (2, 1): lambda r: r * np.exp(-r / 2) / (2 * S6),
    (3, 0): lambda r: 2 / (3 * S3) * (1 - 2 * r / 3 + 2 * r ** 2 / 27) * np.exp(-r / 3),
    (3, 1): lambda r: 8 / (27 * S6) * r * (1 - r / 6) * np.exp(-r / 3),
    (3, 2): lambda r: 4 / (81 * S30) * r ** 2 * np.exp(-r / 3),
}
Y2 = {  # |Y_lm|^2 as function of cos(theta)
    (0, 0): lambda c: np.full_like(c, 1 / (4 * np.pi)),
    (1, 0): lambda c: 3 / (4 * np.pi) * c ** 2,
    (1, 1): lambda c: 3 / (8 * np.pi) * (1 - c ** 2),
    (2, 0): lambda c: 5 / (16 * np.pi) * (3 * c ** 2 - 1) ** 2,
    (2, 1): lambda c: 15 / (8 * np.pi) * (1 - c ** 2) * c ** 2,
    (2, 2): lambda c: 15 / (32 * np.pi) * (1 - c ** 2) ** 2,
}


# --------------------------------------------------------------------------
# 3.3 氢原子径向概率密度
# --------------------------------------------------------------------------
class C3_HydrogenRadial(Scene):
    def construct(self):
        title = title_bar(r"氢原子径向概率密度 $W_{nl}(r)=r^2R_{nl}^2(r)$：$n_r=n-l-1$ 个节点")
        self.play(FadeIn(title))
        ax = Axes(x_range=[0, 30, 5], y_range=[0, 0.6, 0.1], x_length=8.4, y_length=5.0,
                  axis_config={"include_tip": False, "font_size": 22},
                  x_axis_config={"numbers_to_include": [5, 10, 15, 20, 25, 30]},
                  y_axis_config={"numbers_to_include": [0.1, 0.2, 0.3, 0.4, 0.5]}).shift(LEFT * 2.2 + DOWN * 0.6)
        xl = MT(r"r/a_0", 26).next_to(ax.x_axis, DOWN, buff=0.4)
        self.play(Create(ax), FadeIn(xl))
        r = np.linspace(0, 30, 600)
        cols = {(1, 0): YELLOW, (2, 0): BLUE_C, (2, 1): GREEN_C, (3, 0): RED_C, (3, 1): PURPLE_B, (3, 2): ORANGE}
        info = {
            (1, 0): r"1s：最概然半径 $r=a_0$（玻尔半径）",
            (2, 0): r"2s：1 个节点",
            (2, 1): r"2p：无节点，峰在 $4a_0$",
            (3, 0): r"3s：2 个节点",
            (3, 1): r"3p：1 个节点",
            (3, 2): r"3d：无节点，峰在 $9a_0$",
        }
        rows = VGroup()
        prev = None
        for key in [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]:
            W = r ** 2 * RNL[key](r) ** 2
            c = curve(ax, r, W, cols[key], 3.5)
            t = CT(info[key], 24, cols[key])
            rows.add(t)
            rows.arrange(DOWN, aligned_edge=LEFT, buff=0.25).to_corner(UR, buff=0.4).shift(DOWN * 0.8)
            anims = [Create(c), FadeIn(t)]
            if prev is not None:
                anims.append(prev.animate.set_stroke(opacity=0.35))
            self.play(*anims, run_time=1.4)
            prev = c
            self.wait(0.4)
        k = boxed(panel(r"$l=n-1$ 时峰在 $r=n^2a_0$——恰是玻尔轨道半径",
                        r"$n$ 越大，电子离核越远；$l$ 越大，节点越少", size=22))
        fitw(k, 5.4)
        k.next_to(rows, DOWN, buff=0.4).to_edge(RIGHT, buff=0.25)
        self.play(FadeIn(k))
        self.wait(2.5)


def density_image(n, l, m, R=None, N=220):
    R = R or {1: 5, 2: 14, 3: 26}[n]
    g = np.linspace(-R, R, N)
    X, Z = np.meshgrid(g, -g)
    r = np.sqrt(X ** 2 + Z ** 2) + 1e-9
    rho = RNL[(n, l)](r) ** 2 * Y2[(l, abs(m))](Z / r)
    v = (rho / rho.max()) ** 0.45
    stops = np.array([[0, 0, 0], [60, 20, 110], [200, 60, 120], [255, 160, 40], [255, 250, 200]], float)
    idx = v * (len(stops) - 1)
    i0 = np.clip(np.floor(idx).astype(int), 0, len(stops) - 2)
    f = (idx - i0)[..., None]
    rgb = stops[i0] * (1 - f) + stops[i0 + 1] * f
    return ImageMobject(rgb.astype(np.uint8))


# --------------------------------------------------------------------------
# 3.3 氢原子电子云：|ψ_nlm|^2 在 xz 平面的截面
# --------------------------------------------------------------------------
class C3_HydrogenOrbitals(Scene):
    def construct(self):
        title = title_bar(r"氢原子电子云 $\lvert\psi_{nlm}\rvert^2=R_{nl}^2\lvert Y_{lm}\rvert^2$（$xz$ 截面，$z$ 轴竖直，绕 $z$ 轴旋转对称）")
        self.play(FadeIn(title))
        states = [(1, 0, 0), (2, 0, 0), (2, 1, 0), (2, 1, 1), (3, 2, 0), (3, 2, 1), (3, 2, 2), (3, 1, 0)]
        names = ["1s", "2s", "2p", "2p", "3d", "3d", "3d", "3p"]
        cells = Group()
        for (n, l, m), nm in zip(states, names):
            im = density_image(n, l, m)
            im.height = 2.35
            lab = MT(rf"\psi_{{{n}{l}{'' if m == 0 else ('\\pm' + str(m))}{'0' if m == 0 else ''}}}", 26).next_to(im, DOWN, buff=0.08)
            lab2 = CT(rf"{nm}，$m={'0' if m == 0 else r'\pm' + str(m)}$", 20, GREY_B).next_to(lab, DOWN, buff=0.04)
            cells.add(Group(im, lab, lab2))
        cells.arrange_in_grid(rows=2, cols=4, buff=(0.35, 0.25)).next_to(title, DOWN, buff=0.25)
        for i in range(0, 8, 2):
            self.play(FadeIn(cells[i]), FadeIn(cells[i + 1]), run_time=1.0)
            self.wait(0.6)
        note = CT(r"$m$ 不同，电子云形状不同但能量相同（$E_n$ 只与 $n$ 有关）：能级 $E_n$ 的简并度为 $n^2$", 24, C_HL).to_edge(DOWN, buff=0.15)
        fitw(note, 13)
        self.play(FadeIn(note))
        self.wait(3)


# --------------------------------------------------------------------------
# 3.2 角动量的矢量模型（l = 2）
# --------------------------------------------------------------------------
class C3_AngularMomentum(Scene):
    def construct(self):
        title = title_bar(r"角动量矢量模型（$l=2$）：$\lvert\vv L\rvert=\sqrt{l(l+1)}\,\hbar$，$L_z=m\hbar$ 只能取 $2l+1$ 个值")
        self.play(FadeIn(title))
        u = 1.05
        O = LEFT * 2.6 + DOWN * 0.6
        Lmag = np.sqrt(6)
        zax = Arrow(O + DOWN * 3.0, O + UP * 3.1, buff=0, stroke_width=3, color=GREY_B)
        zl = MT(r"z", 28, GREY_B).next_to(zax, UP, buff=0.05)
        circ = Circle(radius=Lmag * u, color=GREY_D, stroke_width=1.5).move_to(O)
        self.play(GrowArrow(zax), FadeIn(zl), Create(circ))
        vecs = VGroup(); labs = VGroup(); cones = VGroup()
        for m in range(-2, 3):
            x = np.sqrt(6 - m * m)
            tip = O + np.array([x * u, m * u, 0])
            v = Arrow(O, tip, buff=0, color=C_HL, stroke_width=4, max_tip_length_to_length_ratio=0.12)
            dl = DashedLine(O + UP * m * u + LEFT * x * u, O + UP * m * u + RIGHT * x * u, color=BLUE_C, stroke_width=1.5)
            ell = Ellipse(width=2 * x * u, height=0.35 * x * u / 2.4, color=BLUE_C, stroke_width=2).move_to(O + UP * m * u)
            lab = MT(rf"m={m:+d}" if m else "m=0", 24, BLUE_C).next_to(O + UP * m * u + RIGHT * x * u, RIGHT, buff=0.12)
            vecs.add(v); labs.add(lab); cones.add(VGroup(dl, ell))
        self.play(LaggedStart(*[GrowArrow(v) for v in vecs], lag_ratio=0.3), LaggedStart(*[FadeIn(l) for l in labs], lag_ratio=0.3), run_time=2.5)
        txt = boxed(panel(r"$\hat L^2Y_{lm}=l(l+1)\hbar^2Y_{lm}$，$\hat L_zY_{lm}=m\hbar Y_{lm}$",
                          r"$l=2$：$\lvert\vv L\rvert=\sqrt6\,\hbar\approx2.45\hbar$，而 $L_z^{\max}=2\hbar$",
                          r"$\Rightarrow$ $\vv L$ 永远不能完全指向 $z$ 轴", size=24))
        txt.to_corner(UR, buff=0.3).shift(DOWN * 0.9)
        self.play(FadeIn(txt))
        self.play(FadeIn(cones))
        # 进动：m=1 的矢量绕 z 轴转
        m = 1; x = np.sqrt(5)
        ang = ValueTracker(0)
        pv = always_redraw(lambda: Arrow(O, O + np.array([x * u * np.cos(ang.get_value()), m * u + 0.35 * x * u / 4.8 * np.sin(ang.get_value()), 0]),
                                         buff=0, color=YELLOW, stroke_width=5, max_tip_length_to_length_ratio=0.12))
        self.play(vecs.animate.set_opacity(0.25), FadeIn(pv))
        t2 = CT(r"$L_x,L_y$ 与 $L_z$ 不对易：$L_z$ 确定时 $L_x,L_y$ 不确定，$\langle L_x\rangle=\langle L_y\rangle=0$", 24, YELLOW)
        t2b = CT(r"可以想象 $\vv L$ 在圆锥面上“进动”，但它只是帮助记忆的图像", 24, YELLOW)
        g2 = VGroup(t2, t2b).arrange(DOWN, aligned_edge=LEFT, buff=0.2).next_to(txt, DOWN, buff=0.4)
        fitw(g2, 6.6)
        self.play(FadeIn(g2))
        self.play(ang.animate.set_value(2 * TAU), run_time=6, rate_func=linear)
        self.wait(1.5)


# --------------------------------------------------------------------------
# 3.7 不确定关系：高斯波包在坐标空间与动量空间的宽度此消彼长
# --------------------------------------------------------------------------
class C3_Uncertainty(Scene):
    def construct(self):
        title = title_bar(r"不确定关系：坐标分布越窄，动量分布越宽，$\Delta x\,\Delta p\geqslant\hbar/2$（高斯态取等号）")
        self.play(FadeIn(title))
        axL = Axes(x_range=[-4, 4, 1], y_range=[0, 1.6, 0.5], x_length=5.6, y_length=3.6,
                   axis_config={"include_tip": False, "include_ticks": False}).shift(LEFT * 3.4 + DOWN * 0.1)
        axR = Axes(x_range=[-4, 4, 1], y_range=[0, 1.6, 0.5], x_length=5.6, y_length=3.6,
                   axis_config={"include_tip": False, "include_ticks": False}).shift(RIGHT * 3.4 + DOWN * 0.1)
        lL = CT(r"坐标空间 $\lvert\psi(x)\rvert^2$", 26, C_PROB).next_to(axL, UP, buff=0.15)
        lR = CT(r"动量空间 $\lvert c(p)\rvert^2$", 26, GREEN_C).next_to(axR, UP, buff=0.15)
        xL = MT("x", 26).next_to(axL.x_axis, RIGHT, buff=0.05)
        xR = MT("p", 26).next_to(axR.x_axis, RIGHT, buff=0.05)
        self.play(Create(axL), Create(axR), FadeIn(lL), FadeIn(lR), FadeIn(xL), FadeIn(xR))
        sig = ValueTracker(1.0)
        g = np.linspace(-4, 4, 400)

        def gauss(s):
            return np.exp(-g ** 2 / (2 * s ** 2)) / np.sqrt(2 * np.pi * s ** 2)
        L = always_redraw(lambda: filled_curve(axL, g, gauss(sig.get_value()), C_PROB, 0.35, 3))
        R = always_redraw(lambda: filled_curve(axR, g, gauss(0.5 / sig.get_value()), GREEN_C, 0.35, 3))
        self.play(FadeIn(L), FadeIn(R))
        dxn = DecimalNumber(1.0, num_decimal_places=2, font_size=30)
        dpn = DecimalNumber(0.5, num_decimal_places=2, font_size=30)
        prn = DecimalNumber(0.5, num_decimal_places=2, font_size=30, color=C_HL)
        row = VGroup(MT(r"\Delta x=", 30), dxn, MT(r",\quad \Delta p=", 30), dpn, MT(r"\hbar,\quad \Delta x\,\Delta p=", 30), prn,
                     MT(r"\hbar", 30)).arrange(RIGHT, buff=0.08).shift(DOWN * 2.75)
        dxn.add_updater(lambda m: m.set_value(sig.get_value()))
        dpn.add_updater(lambda m: m.set_value(0.5 / sig.get_value()))
        prn.add_updater(lambda m: m.set_value(0.5))
        self.add(row)
        self.play(sig.animate.set_value(0.35), run_time=3)
        self.wait(0.5)
        self.play(sig.animate.set_value(1.8), run_time=4)
        self.wait(0.5)
        self.play(sig.animate.set_value(0.9), run_time=2)
        k = boxed(panel(r"$\psi(x)$ 与 $c(p)$ 互为傅里叶变换：窄脉冲 $\leftrightarrow$ 宽频谱（信号与系统的时宽–带宽积）",
                        r"根源：$[\hat x,\hat p]=\ii\hbar\neq0$ $\Rightarrow$ $\Delta A\,\Delta B\geqslant\frac12\lvert\langle[\hat A,\hat B]\rangle\rvert$", size=22))
        fitw(k, 12.8)
        k.to_edge(DOWN, buff=0.2)
        self.play(row.animate.next_to(title, DOWN, buff=0.2), FadeIn(k))
        self.wait(3)
