"""第 2 章 波函数和薛定谔方程 —— manim 动画源码（LaTeX 排版）."""
from math import factorial

from numpy.polynomial.hermite import hermval

from common import *  # noqa: F401,F403


def osc_psi(n, xi):
    c = np.zeros(n + 1); c[n] = 1
    return hermval(xi, c) * np.exp(-xi ** 2 / 2) / np.sqrt(2.0 ** n * factorial(n) * np.sqrt(np.pi))


# --------------------------------------------------------------------------
# 2.1/2.3 自由粒子波包：群速度运动 + 扩散
# --------------------------------------------------------------------------
class C2_WavePacket(Scene):
    def construct(self):
        title = title_bar(r"自由粒子高斯波包：$\lvert\Psi\rvert^2$ 以群速度 $v_g=\hbar k_0/m$ 前进，同时不断变宽")
        self.play(FadeIn(title))
        N, L = 2048, 120.0
        x = np.linspace(-L / 2, L / 2, N, endpoint=False)
        k = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
        s0, k0, x0 = 1.5, 3.0, -18.0
        psi0 = (2 * np.pi * s0 ** 2) ** -0.25 * np.exp(-(x - x0) ** 2 / (4 * s0 ** 2) + 1j * k0 * x)
        phi0 = np.fft.fft(psi0)

        def psi_t(t):
            return np.fft.ifft(phi0 * np.exp(-0.5j * k ** 2 * t))

        ax = Axes(x_range=[-25, 25, 5], y_range=[-0.6, 0.6, 0.2], x_length=12, y_length=4.6,
                  axis_config={"include_tip": False, "include_ticks": False}).shift(DOWN * 0.7)
        xl = MT("x", 28).next_to(ax.x_axis, RIGHT, buff=0.1)
        self.play(Create(ax), FadeIn(xl))
        m = (x > -25) & (x < 25)
        xs = x[m]
        t = ValueTracker(0.0)

        def draw():
            p = psi_t(t.get_value())[m]
            re = curve(ax, xs, 0.9 * p.real, C_PSI, 2)
            pr = filled_curve(ax, xs, 2.0 * np.abs(p) ** 2, C_PROB, 0.35, 3)
            return VGroup(re, pr)
        g = always_redraw(draw)
        leg = VGroup(CT(r"$\Re\Psi$（相速度传播的载波）", 24, C_PSI), CT(r"$\lvert\Psi\rvert^2$（概率密度，包络）", 24, C_PROB)
                     ).arrange(RIGHT, buff=0.8).next_to(title, DOWN, buff=0.25)
        self.play(FadeIn(g), FadeIn(leg))
        tn = DecimalNumber(0, num_decimal_places=1, font_size=26)
        sn = DecimalNumber(s0, num_decimal_places=2, font_size=26)
        info = VGroup(MT(r"t=", 26), tn, MT(r",\quad \sigma(t)=\sigma_0\sqrt{1+(\hbar t/2m\sigma_0^2)^2}=", 26), sn).arrange(RIGHT, buff=0.08).to_corner(DL, buff=0.35)
        tn.add_updater(lambda m: m.set_value(t.get_value()))
        sn.add_updater(lambda m: m.set_value(s0 * np.sqrt(1 + (t.get_value() / (2 * s0 ** 2)) ** 2)))
        self.add(info)
        self.play(t.animate.set_value(12.0), run_time=9, rate_func=linear)
        note = boxed(panel(r"$\bullet$ 波包 = 许多平面波 $\ee^{\ii(kx-\omega t)}$ 的叠加（傅里叶积分）",
                           r"$\bullet$ 色散关系 $\omega=\hbar k^2/2m$ 非线性 $\Rightarrow$ 各分量相速度不同 $\Rightarrow$ 波包扩散",
                           r"$\bullet$ 所以“电子就是波包”的说法不成立：波包会散开，电子不会", size=22))
        note.to_corner(UR, buff=0.3).shift(DOWN * 1.1)
        fitw(note, 7.4)
        self.play(FadeIn(note))
        self.wait(3)


# --------------------------------------------------------------------------
# 2.6 无限深方势阱：本征态 + 叠加态的时间演化
# --------------------------------------------------------------------------
class C2_InfiniteWell(Scene):
    def construct(self):
        title = title_bar(r"一维无限深方势阱：定态是驻波，$E_n\propto n^2$；两个定态叠加后 $\lvert\Psi\rvert^2$ 来回振荡")
        self.play(FadeIn(title))
        ax = Axes(x_range=[0, 1, 0.5], y_range=[0, 18, 4], x_length=4.5, y_length=5.8,
                  axis_config={"include_tip": False, "include_ticks": False}).shift(LEFT * 4.0 + DOWN * 0.55)
        walls = VGroup(Line(ax.c2p(0, 0), ax.c2p(0, 18)), Line(ax.c2p(1, 0), ax.c2p(1, 18))).set_stroke(C_POT, 5)
        self.play(Create(walls), Create(ax.x_axis))
        xs = np.linspace(0, 1, 200)
        levels = VGroup()
        for n in range(1, 5):
            E = n ** 2  # 以 E_1 为单位
            lv = DashedLine(ax.c2p(0, E), ax.c2p(1, E), color=C_E, stroke_width=1.5)
            wf = curve(ax, xs, E + 1.3 * np.sin(n * np.pi * xs), C_PSI, 3)
            lab = MT(rf"n={n},\ E_{n}={n * n}E_1", 22, C_E).next_to(ax.c2p(1, E), RIGHT, buff=0.15)
            levels.add(VGroup(lv, wf, lab))
        self.play(LaggedStart(*[FadeIn(l) for l in levels], lag_ratio=0.5), run_time=3)
        f = boxed(panel(MT(r"\psi_n(x)=\sqrt{\tfrac{2}{a}}\sin\frac{n\pi x}{a},\quad E_n=\frac{n^2\pi^2\hbar^2}{2ma^2}", 28),
                        r"$\psi_n$ 有 $n-1$ 个节点；基态能量 $E_1>0$（零点能）", size=24))
        f.to_corner(UR, buff=0.3).shift(DOWN * 0.9)
        self.play(FadeIn(f))
        self.wait(1.5)
        # 叠加态
        ax2 = Axes(x_range=[0, 1, 0.5], y_range=[0, 3.2, 1], x_length=5.6, y_length=3.0,
                   axis_config={"include_tip": False, "include_ticks": False}).shift(RIGHT * 2.8 + DOWN * 1.7)
        w2 = VGroup(Line(ax2.c2p(0, 0), ax2.c2p(0, 3.2)), Line(ax2.c2p(1, 0), ax2.c2p(1, 3.2))).set_stroke(C_POT, 4)
        sup = CT(r"$\Psi=\frac{1}{\sqrt2}\left[\psi_1\ee^{-\ii E_1t/\hbar}+\psi_2\ee^{-\ii E_2t/\hbar}\right]$", 26).next_to(ax2, UP, buff=0.25)
        self.play(Create(ax2.x_axis), Create(w2), Write(sup))
        t = ValueTracker(0.0)
        s1, s2 = np.sqrt(2) * np.sin(np.pi * xs), np.sqrt(2) * np.sin(2 * np.pi * xs)

        def dens():
            ph = t.get_value()  # 相位差 (E2-E1)t/ħ
            rho = 0.5 * (s1 ** 2 + s2 ** 2 + 2 * s1 * s2 * np.cos(ph))
            xm = 0.5 - 16 / (9 * np.pi ** 2) * np.cos(ph)
            d = Dot(ax2.c2p(xm, 0), color=C_HL, radius=0.09)
            return VGroup(filled_curve(ax2, xs, rho, C_PROB, 0.4, 3), d)
        D = always_redraw(dens)
        lab = CT(r"橙点：$\langle x\rangle=\frac a2-\frac{16a}{9\pi^2}\cos\frac{(E_2-E_1)t}{\hbar}$", 22, C_HL).next_to(ax2, DOWN, buff=0.2)
        self.play(FadeIn(D), FadeIn(lab))
        self.play(t.animate.set_value(4 * PI), run_time=8, rate_func=linear)
        nt = CT(r"振荡角频率 $\omega_{21}=(E_2-E_1)/\hbar$——正是玻尔频率条件", 26, C_HL)
        fitw(nt, 6.6).move_to(f)
        self.play(FadeOut(f), FadeIn(nt))
        self.play(t.animate.set_value(6 * PI), run_time=4, rate_func=linear)
        self.wait(1.5)


# --------------------------------------------------------------------------
# 2.7 线性谐振子：等间距能级、波函数、经典对应
# --------------------------------------------------------------------------
class C2_Oscillator(Scene):
    def construct(self):
        title = title_bar(r"线性谐振子：能级等间距 $E_n=(n+\tfrac12)\hbar\omega$，波函数可渗入经典禁区")
        self.play(FadeIn(title))
        ax = Axes(x_range=[-5, 5, 1], y_range=[0, 6.2, 1], x_length=7.2, y_length=5.8,
                  axis_config={"include_tip": False, "include_ticks": False}).shift(LEFT * 2.6 + DOWN * 0.55)
        pot = ax.plot(lambda x: 0.5 * x * x, x_range=[-3.5, 3.5], color=C_POT, stroke_width=3)
        xl = MT(r"\xi=\alpha x", 26).next_to(ax.x_axis, RIGHT, buff=0.1)
        self.play(Create(ax.x_axis), Create(pot), FadeIn(xl))
        xi = np.linspace(-5, 5, 400)
        states = VGroup()
        for n in range(5):
            E = n + 0.5
            a = np.sqrt(2 * E)
            lv = Line(ax.c2p(-a, E), ax.c2p(a, E), color=C_E, stroke_width=1.5)
            wf = curve(ax, xi, E + 0.75 * osc_psi(n, xi), C_PSI, 3)
            lab = MT(rf"n={n}", 22, C_E).next_to(ax.c2p(4.6, E), RIGHT, buff=0.05)
            states.add(VGroup(lv, wf, lab))
        self.play(LaggedStart(*[FadeIn(s) for s in states], lag_ratio=0.5), run_time=3.5)
        br = BraceBetweenPoints(ax.c2p(-4.0, 1.5), ax.c2p(-4.0, 0.5), direction=LEFT)
        bl = MT(r"\hbar\omega", 26, C_HL).next_to(br, LEFT, buff=0.1)
        self.play(FadeIn(br), FadeIn(bl))
        f = boxed(panel(MT(r"\psi_n(\xi)=N_n H_n(\xi)\,\ee^{-\xi^2/2},\ \ \alpha=\sqrt{m\omega/\hbar}", 26),
                        r"$\bullet$ 零点能 $E_0=\tfrac12\hbar\omega$",
                        r"$\bullet$ 宇称 $(-1)^n$，$n$ 个节点",
                        r"$\bullet$ 经典转折点 $\xi=\pm\sqrt{2n+1}$ 之外 $\psi\neq0$", size=22))
        f.to_corner(UR, buff=0.3).shift(DOWN * 0.9)
        self.play(FadeIn(f))
        self.wait(2)
        # n = 10 与经典概率对比
        self.play(FadeOut(states), FadeOut(br), FadeOut(bl))
        n = 10
        E = n + 0.5
        A = np.sqrt(2 * E)
        lv = Line(ax.c2p(-A, E * 0.5), ax.c2p(A, E * 0.5), color=C_E, stroke_width=1)
        rho = osc_psi(n, xi) ** 2
        q = filled_curve(ax, xi, 1.0 + 9.0 * rho, C_PROB, 0.35, 2.5, y0=1.0)
        xc = np.linspace(-A + 0.02, A - 0.02, 300)
        cl = curve(ax, xc, 1.0 + 9.0 / (np.pi * np.sqrt(A ** 2 - xc ** 2)), RED_C, 3)
        t1 = CT(r"$n=10$：量子 $\lvert\psi_{10}\rvert^2$（黄）与经典 $w(\xi)\propto 1/\sqrt{\xi_0^2-\xi^2}$（红）", 24).next_to(f, DOWN, buff=0.35)
        fitw(t1, 6.0).align_to(f, RIGHT)
        self.play(FadeIn(q), FadeIn(t1))
        self.play(Create(cl))
        t2 = CT(r"$n$ 越大，量子分布的平均越接近经典——对应原理", 24, C_HL).next_to(t1, DOWN, buff=0.3)
        fitw(t2, 6.0).align_to(f, RIGHT)
        self.play(FadeIn(t2))
        self.wait(2.5)


# --------------------------------------------------------------------------
# 2.8 势垒穿透：波包撞上方势垒（分步傅里叶法数值求解）
# --------------------------------------------------------------------------
class C2_Tunneling(Scene):
    def construct(self):
        title = title_bar(r"隧道效应：$E<U_0$ 的波包撞上方势垒，一部分反射、一部分穿透")
        self.play(FadeIn(title))
        N, L = 4096, 400.0
        x = np.linspace(-L / 2, L / 2, N, endpoint=False)
        dx = L / N
        k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
        U0, a = 1.0, 2.0
        V = np.where((x > 0) & (x < a), U0, 0.0)
        k0, s0, x0 = 1.2, 6.0, -40.0
        psi = (2 * np.pi * s0 ** 2) ** -0.25 * np.exp(-(x - x0) ** 2 / (4 * s0 ** 2) + 1j * k0 * x)
        dt, steps, every = 0.05, 1500, 10
        half = np.exp(-0.5j * V * dt)
        kin = np.exp(-0.5j * k ** 2 * dt)
        frames = [psi.copy()]
        for i in range(steps):
            psi = half * np.fft.ifft(kin * np.fft.fft(half * psi))
            if (i + 1) % every == 0:
                frames.append(psi.copy())
        Tprob = float(np.sum(np.abs(frames[-1][x > a]) ** 2) * dx)
        m = (x > -70) & (x < 70)
        xs = x[m]
        ax = Axes(x_range=[-70, 70, 10], y_range=[-0.25, 0.45, 0.1], x_length=12.6, y_length=4.6,
                  axis_config={"include_tip": False, "include_ticks": False}).shift(DOWN * 0.9)
        bar = Rectangle(width=ax.x_length * a / 140, height=ax.y_length * 0.4 / 0.7,
                        color=C_POT, fill_opacity=0.25)
        bar.move_to(ax.c2p(a / 2, 0.2))
        blab = CT(r"势垒 $U_0$", 22, C_POT).next_to(bar, UP, buff=0.05)
        Eline = DashedLine(ax.c2p(-70, 0.72 * 0.4), ax.c2p(70, 0.72 * 0.4), color=C_E, stroke_width=1.5)
        Elab = CT(r"粒子平均能量 $E\approx0.72U_0$", 20, C_E).next_to(Eline.get_start(), UP, buff=0.05).shift(RIGHT * 1.6)
        self.play(Create(ax.x_axis), FadeIn(bar), FadeIn(blab), Create(Eline), FadeIn(Elab))
        idx = ValueTracker(0)

        def draw():
            p = frames[int(round(idx.get_value()))][m]
            return VGroup(curve(ax, xs, 0.9 * p.real, C_PSI, 1.6), filled_curve(ax, xs, 1.6 * np.abs(p) ** 2, C_PROB, 0.4, 2.6))
        g = always_redraw(draw)
        self.play(FadeIn(g))
        self.play(idx.animate.set_value(len(frames) - 1), run_time=12, rate_func=linear)
        res = boxed(panel(rf"数值结果：透射概率 $\approx{Tprob:.2f}$，反射概率 $\approx{1 - Tprob:.2f}$",
                          MT(r"D\approx D_0\,\ee^{-\frac{2}{\hbar}\sqrt{2m(U_0-E)}\,a}", 28),
                          r"$D$ 对势垒宽度 $a$ 和 $\sqrt{m(U_0-E)}$ 指数敏感", size=22))
        res.to_corner(UR, buff=0.25).shift(DOWN * 0.9)
        self.play(FadeIn(res))
        self.wait(3)
