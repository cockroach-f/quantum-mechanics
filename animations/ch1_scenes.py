"""第 1 章 绪论 —— manim 动画源码（LaTeX 排版）."""
from common import *  # noqa: F401,F403

HC_K = 14388.0  # hc/k_B, 单位 μm·K


def planck(lam, T):
    return lam ** -5 / np.expm1(HC_K / (lam * T))


# --------------------------------------------------------------------------
# 1.2 黑体辐射：普朗克公式 vs 瑞利–金斯 / 维恩，维恩位移
# --------------------------------------------------------------------------
class C1_Blackbody(Scene):
    def construct(self):
        title = title_bar(r"黑体辐射：普朗克公式同时拟合长波与短波，峰值波长 $\lambda_m\propto 1/T$")
        self.play(FadeIn(title))
        ax = Axes(x_range=[0, 3, 0.5], y_range=[0, 1.2, 0.2], x_length=8.2, y_length=5.0,
                  axis_config={"include_tip": False, "font_size": 22},
                  x_axis_config={"numbers_to_include": [0.5, 1, 1.5, 2, 2.5, 3]}).shift(LEFT * 2.2 + DOWN * 0.55)
        xl = CT(r"波长 $\lambda/\mu\mathrm{m}$", 24).next_to(ax.x_axis, DOWN, buff=0.35)
        yl = CT(r"能量密度 $\rho_\lambda$（相对值）", 24).rotate(PI / 2).next_to(ax.y_axis, LEFT, buff=0.2)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))
        norm = planck(HC_K / (4.965 * 6000), 6000)
        lam = np.linspace(0.08, 3, 400)
        T = ValueTracker(6000)

        def pl():
            t = T.get_value()
            y = np.clip(planck(lam, t) / norm, 0, 1.25)
            return curve(ax, lam, y, YELLOW, 4)

        def rj():
            t = T.get_value()
            l2 = lam[lam > 0.35]
            y = np.clip(t / HC_K * l2 ** -4 / norm, 0, 1.2)
            return curve(ax, l2[y < 1.2], y[y < 1.2], RED_C, 2.5)

        def wien():
            t = T.get_value()
            return curve(ax, lam, np.clip(lam ** -5 * np.exp(-HC_K / (lam * t)) / norm, 0, 1.2), BLUE_C, 2.5)

        P = always_redraw(pl)
        legend = VGroup(
            CT(r"普朗克公式（与实验一致）", 24, YELLOW),
            CT(r"瑞利–金斯公式：短波发散（紫外灾难）", 24, RED_C),
            CT(r"维恩公式：长波偏低", 24, BLUE_C),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.2).to_corner(UR, buff=0.4).shift(DOWN * 0.8)
        self.play(Create(P), FadeIn(legend[0]), run_time=1.5)
        R = always_redraw(rj)
        self.play(Create(R), FadeIn(legend[1]))
        W = always_redraw(wien)
        self.play(Create(W), FadeIn(legend[2]))
        self.wait(1.0)
        f1 = boxed(VGroup(
            MT(r"\rho_\nu\,\dd\nu=\frac{8\pi h\nu^3}{c^3}\,\frac{1}{\ee^{h\nu/k_BT}-1}\,\dd\nu", 30),
            CT(r"关键假设：能量只能是 $h\nu$ 的整数倍", 22, C_HL)).arrange(DOWN, buff=0.18))
        f1.next_to(legend, DOWN, buff=0.45)
        self.play(FadeIn(f1))
        self.wait(1.2)
        self.play(FadeOut(R), FadeOut(W), FadeOut(legend[1:]))
        peak = always_redraw(lambda: Dot(ax.c2p(2898 / T.get_value(), planck(2898 / T.get_value(), T.get_value()) / norm), color=C_HL))
        vline = always_redraw(lambda: DashedLine(ax.c2p(2898 / T.get_value(), 0), peak.get_center(), color=C_HL))
        tnum = DecimalNumber(6000, num_decimal_places=0, font_size=26, color=C_HL)
        lnum = DecimalNumber(0.48, num_decimal_places=2, font_size=26, color=C_HL)
        lab = VGroup(MT(r"T=", 26, C_HL), tnum, MT(r"\,\mathrm{K},\quad \lambda_m=", 26, C_HL), lnum, MT(r"\,\mu\mathrm{m}", 26, C_HL)).arrange(RIGHT, buff=0.08)
        lab.next_to(f1, DOWN, buff=0.35)
        tnum.add_updater(lambda m: m.set_value(T.get_value()))
        lnum.add_updater(lambda m: m.set_value(2898 / T.get_value()))
        self.play(FadeIn(peak), Create(vline), FadeIn(lab))
        self.play(T.animate.set_value(3500), run_time=4, rate_func=smooth)
        self.play(T.animate.set_value(5000), run_time=2.5)
        wd = CT(r"维恩位移定律：$\lambda_m T=b\approx2.90\times10^{-3}\,\mathrm{m\cdot K}$", 26, WHITE)
        fitw(wd, 5.0)
        wd.next_to(lab, DOWN, buff=0.3).to_edge(RIGHT, buff=0.3)
        self.play(Write(wd))
        self.wait(2)


# --------------------------------------------------------------------------
# 1.2 光电效应：频率决定能否发射，强度只决定电子数目
# --------------------------------------------------------------------------
class C1_Photoelectric(Scene):
    def construct(self):
        title = title_bar(r"光电效应：一个光子把全部能量 $h\nu$ 交给一个电子，$E_k=h\nu-W_0$")
        self.play(FadeIn(title))
        metal = Rectangle(width=0.5, height=3.6, color=GREY_B, fill_opacity=0.6).move_to(LEFT * 4.6 + DOWN * 1.2)
        mlab = CT(r"金属（逸出功 $W_0$）", 22, GREY_B).next_to(metal, DOWN, buff=0.15)
        self.add(metal, mlab)
        ax = Axes(x_range=[0, 10, 2], y_range=[-1, 5, 1], x_length=5.6, y_length=4.0,
                  axis_config={"include_tip": False}).shift(RIGHT * 3.3 + DOWN * 0.6)
        xl = MT(r"\nu", 28).next_to(ax.x_axis, RIGHT, buff=0.1)
        yl = MT(r"E_k", 28).next_to(ax.y_axis, UP, buff=0.1)
        nu0 = 4.0
        line = ax.plot(lambda v: max(0.0, 0.8 * (v - nu0)), x_range=[0, 10], color=YELLOW)
        ext = DashedLine(ax.c2p(nu0, 0), ax.c2p(0, -0.8 * nu0 / 1.0 * 0.25), color=GREY_B)
        n0 = MT(r"\nu_0=W_0/h", 24, C_HL).next_to(ax.c2p(nu0, 0), DOWN, buff=0.15)
        slope = CT(r"斜率 $=h$（与金属无关）", 22, YELLOW).next_to(ax.c2p(8.5, 3.6), UP, buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))
        self.play(Create(line), FadeIn(n0), FadeIn(slope), Create(ext))

        def shoot(nu, n_photons, color, label):
            txt = CT(label, 24, color).move_to(LEFT * 2.6 + UP * 2.2)
            self.play(FadeIn(txt), run_time=0.4)
            dot = Dot(ax.c2p(nu, max(0, 0.8 * (nu - nu0))), color=color, radius=0.1)
            self.play(FadeIn(dot), run_time=0.3)
            anims = []
            ek = max(0.0, 0.8 * (nu - nu0))
            for k in range(n_photons):
                y = -2.4 + 2.6 * (k + 0.5) / n_photons
                ph = FunctionGraph(lambda x: 0.12 * np.sin(14 * x), x_range=[0, 1.2], color=color).move_to(LEFT * 2.4 + UP * (y - 0.6))
                anims.append(Succession(FadeIn(ph, run_time=0.1),
                                        ph.animate(run_time=0.8).move_to(metal.get_right() + RIGHT * 0.65 + UP * (y - 0.6)),
                                        FadeOut(ph, run_time=0.1)))
            self.play(LaggedStart(*anims, lag_ratio=0.15))
            if ek > 0:
                es = VGroup(*[Dot(metal.get_right() + UP * (-2.4 + 2.6 * (k + 0.5) / n_photons - 0.6), radius=0.07, color=BLUE_B)
                              for k in range(n_photons)])
                self.play(FadeIn(es, run_time=0.2))
                self.play(es.animate.shift(RIGHT * (0.6 + 0.9 * ek) + UP * 0.0), run_time=1.0)
                self.play(FadeOut(es), run_time=0.3)
            else:
                no = CT(r"无电子逸出", 26, RED).next_to(metal, RIGHT, buff=0.4)
                self.play(FadeIn(no))
                self.play(FadeOut(no))
            self.play(FadeOut(txt), run_time=0.3)
            return dot

        d1 = shoot(2.5, 3, RED_C, r"红光 $\nu<\nu_0$：弱光")
        d2 = shoot(2.5, 9, RED_C, r"红光 $\nu<\nu_0$：强光也不行")
        d3 = shoot(6.5, 3, BLUE_C, r"紫光 $\nu>\nu_0$：弱光也立即发射")
        d4 = shoot(8.5, 3, PURPLE_B, r"频率更高：电子动能更大")
        concl = boxed(panel(r"$\bullet$ 有无光电子只看频率，不看光强",
                            r"$\bullet$ 光强 $\uparrow$ $\Rightarrow$ 光子数 $\uparrow$ $\Rightarrow$ 光电流 $\uparrow$",
                            r"$\bullet$ 无需积累时间：一个光子、一次吸收", size=24))
        fitw(concl, 6.0)
        concl.to_corner(UL, buff=0.3).shift(DOWN * 0.75)
        self.play(FadeOut(VGroup(d1, d2)), FadeIn(concl))
        self.wait(2.5)


# --------------------------------------------------------------------------
# 1.4 电子双缝：单个电子落点随机，大量电子累积出干涉条纹
# --------------------------------------------------------------------------
class C1_DoubleSlit(Scene):
    def construct(self):
        title = title_bar(r"电子双缝实验：每个电子随机落点，累积起来却是干涉条纹——概率波")
        self.play(FadeIn(title))
        rng = np.random.default_rng(7)
        # 屏上强度分布 I(y) ∝ cos^2(π d y/λL) sinc^2(π a y/λL)
        def I(y):
            return np.cos(3.2 * y) ** 2 * np.sinc(0.55 * y) ** 2
        ys = np.linspace(-3, 3, 2000)
        cdf = np.cumsum(I(ys)); cdf /= cdf[-1]
        # 光源、双缝、屏
        src = Dot(LEFT * 6.4 + DOWN * 0.4, color=BLUE_B, radius=0.12)
        slab = CT(r"电子枪", 22).next_to(src, DOWN)
        wall = VGroup(Line(LEFT * 4.7 + UP * 2.6, LEFT * 4.7 + UP * 0.05), Line(LEFT * 4.7 + DOWN * 0.25, LEFT * 4.7 + DOWN * 0.55),
                      Line(LEFT * 4.7 + DOWN * 0.85, LEFT * 4.7 + DOWN * 3.4)).set_stroke(GREY_B, 5)
        screen = Rectangle(width=2.2, height=6.0, color=GREY_C, fill_color=BLACK, fill_opacity=1).move_to(LEFT * 1.9 + DOWN * 0.4)
        self.play(FadeIn(src), FadeIn(slab), Create(wall), Create(screen))
        counter = Integer(0, font_size=34, color=C_HL)
        clab = CT(r"电子数 $N=$", 28).next_to(screen, RIGHT, buff=2.2).shift(UP * 2.2)
        counter.next_to(clab, RIGHT)
        self.add(clab, counter)

        def sample(n):
            u = rng.random(n)
            yy = np.interp(u, cdf, ys)
            xx = rng.uniform(-0.95, 0.95, n)
            return yy, xx
        dots_all = VGroup()
        # 单个电子飞行
        for _ in range(4):
            y, x = sample(1)
            e = Dot(src.get_center(), radius=0.06, color=BLUE_B)
            tgt = screen.get_center() + np.array([x[0], y[0] * 0.95, 0])
            self.play(MoveAlongPath(e, VMobject().set_points_smoothly([src.get_center(), LEFT * 4.7 + DOWN * 0.4, tgt])), run_time=0.7)
            d = Dot(tgt, radius=0.025, color=WHITE)
            self.remove(e); self.add(d); dots_all.add(d)
            counter.set_value(len(dots_all))
        note = CT(r"单个电子：落点完全随机，像粒子", 26, BLUE_B).next_to(clab, DOWN, buff=0.5, aligned_edge=LEFT)
        self.play(FadeIn(note))
        for n_target, rt in [(60, 1.5), (400, 2.0), (2500, 2.5)]:
            n = n_target - len(dots_all)
            y, x = sample(n)
            g = VGroup(*[Dot(screen.get_center() + np.array([x[i], y[i] * 0.95, 0]), radius=0.022, color=WHITE) for i in range(n)])
            self.play(LaggedStart(*[FadeIn(dd) for dd in g], lag_ratio=1.0 / max(n, 1) * 2), counter.animate.set_value(n_target), run_time=rt)
            dots_all.add(*g)
        note2 = CT(r"大量电子：亮暗相间的干涉条纹，像波", 26, C_PROB).next_to(note, DOWN, buff=0.3, aligned_edge=LEFT)
        self.play(FadeIn(note2))
        # 强度曲线（转 90°）
        prof = VMobject(color=C_PROB, stroke_width=3)
        pts = [screen.get_right() + RIGHT * (0.15 + 1.6 * I(y)) + UP * (y * 0.95) for y in np.linspace(-3, 3, 300)]
        prof.set_points_as_corners(pts)
        self.play(Create(prof), run_time=1.5)
        key = boxed(panel(r"$\lvert\psi(\vv r)\rvert^2\,\dd\tau$：在 $\dd\tau$ 内找到粒子的概率",
                          r"$\psi=\psi_1+\psi_2\ \Rightarrow\ \lvert\psi\rvert^2=\lvert\psi_1\rvert^2+\lvert\psi_2\rvert^2+2\Re(\psi_1^*\psi_2)$",
                          r"干涉项来自概率幅相加，而不是概率相加", size=22))
        fitw(key, 5.6)
        key.to_corner(DR, buff=0.3)
        self.play(FadeIn(key))
        self.wait(2.5)
