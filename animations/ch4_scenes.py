"""第 4 章 态和力学量的表象 —— manim 动画源码（LaTeX 排版）."""
from common import *  # noqa: F401,F403


# --------------------------------------------------------------------------
# 4.1/4.4 表象 = 基矢；幺正变换 = 换一组正交归一基，态矢量本身不变
# --------------------------------------------------------------------------
class C4_BasisChange(Scene):
    def construct(self):
        title = title_bar(r"表象就是坐标系：同一个态矢量 $\ket{\psi}$，换一组基，分量变了，长度不变")
        self.play(FadeIn(title))
        O = LEFT * 3.0 + DOWN * 0.7
        u = 2.4
        psi = np.array([0.8, 0.6, 0])
        grid = NumberPlane(x_range=[-1.4, 1.4, 0.5], y_range=[-1.2, 1.2, 0.5], x_length=2.8 * u, y_length=2.4 * u,
                           background_line_style={"stroke_color": BLUE_E, "stroke_width": 1, "stroke_opacity": 0.4},
                           axis_config={"stroke_opacity": 0}).move_to(O)
        self.add(grid)
        v = Arrow(O, O + psi * u, buff=0, color=C_PROB, stroke_width=6, max_tip_length_to_length_ratio=0.1)
        vl = MT(r"\ket{\psi}", 34, C_PROB).next_to(v.get_end(), UR, buff=0.05)
        self.play(GrowArrow(v), FadeIn(vl))
        th = ValueTracker(0.0)

        def basis():
            t = th.get_value()
            e1 = np.array([np.cos(t), np.sin(t), 0]); e2 = np.array([-np.sin(t), np.cos(t), 0])
            a1 = Arrow(O, O + e1 * u, buff=0, color=BLUE_C, stroke_width=4, max_tip_length_to_length_ratio=0.1)
            a2 = Arrow(O, O + e2 * u, buff=0, color=GREEN_C, stroke_width=4, max_tip_length_to_length_ratio=0.1)
            c1, c2 = psi @ e1, psi @ e2
            p1 = DashedLine(O + psi * u, O + c1 * e1 * u, color=BLUE_C, stroke_width=2)
            p2 = DashedLine(O + psi * u, O + c2 * e2 * u, color=GREEN_C, stroke_width=2)
            d1 = Dot(O + c1 * e1 * u, color=BLUE_C); d2 = Dot(O + c2 * e2 * u, color=GREEN_C)
            return VGroup(a1, a2, p1, p2, d1, d2)
        B = always_redraw(basis)
        lb = always_redraw(lambda: VGroup(
            MT(r"\ket{1}", 28, BLUE_C).next_to(O + u * np.array([np.cos(th.get_value()), np.sin(th.get_value()), 0]), DR, buff=0.05),
            MT(r"\ket{2}", 28, GREEN_C).next_to(O + u * np.array([-np.sin(th.get_value()), np.cos(th.get_value()), 0]), UL, buff=0.05)))
        self.play(FadeIn(B), FadeIn(lb))
        c1n = DecimalNumber(0.8, num_decimal_places=2, font_size=30, color=BLUE_C)
        c2n = DecimalNumber(0.6, num_decimal_places=2, font_size=30, color=GREEN_C)
        nn = DecimalNumber(1.0, num_decimal_places=2, font_size=30, color=C_HL)
        col = VGroup(MT(r"c_1=\ip{1}{\psi}=", 30, BLUE_C), c1n).arrange(RIGHT, buff=0.1)
        col2 = VGroup(MT(r"c_2=\ip{2}{\psi}=", 30, GREEN_C), c2n).arrange(RIGHT, buff=0.1)
        col3 = VGroup(MT(r"\lvert c_1\rvert^2+\lvert c_2\rvert^2=", 30, C_HL), nn).arrange(RIGHT, buff=0.1)
        panel_ = VGroup(col, col2, col3).arrange(DOWN, aligned_edge=LEFT, buff=0.3).to_corner(UR, buff=0.5).shift(DOWN * 0.9 + LEFT * 0.3)
        c1n.add_updater(lambda m: m.set_value(psi @ np.array([np.cos(th.get_value()), np.sin(th.get_value()), 0])))
        c2n.add_updater(lambda m: m.set_value(psi @ np.array([-np.sin(th.get_value()), np.cos(th.get_value()), 0])))
        self.play(FadeIn(panel_))
        self.play(th.animate.set_value(0.65), run_time=3)
        self.play(th.animate.set_value(-0.5), run_time=3)
        self.play(th.animate.set_value(np.arctan2(0.6, 0.8)), run_time=2.5)
        tx = CT(r"基矢转到与 $\ket\psi$ 同向：$\ket\psi$ 成为该表象的“本征态”，$c=(1,0)^{\mathrm T}$", 24, C_HL)
        k = boxed(panel(r"$\bullet$ 新旧分量：$c'=S^\dagger c$，$S_{mn}=\ip{m}{n'}$ 幺正：$S^\dagger S=I$",
                        r"$\bullet$ 算符：$F'=S^\dagger FS$，本征值、迹、$\langle F\rangle$ 都不变",
                        r"$\bullet$ 坐标表象、动量表象、能量表象 = 选不同的基", size=24))
        g = VGroup(tx, k).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
        fitw(g, 6.3)
        g.next_to(panel_, DOWN, buff=0.45).to_edge(RIGHT, buff=0.3)
        self.play(FadeIn(g))
        self.wait(3)


# --------------------------------------------------------------------------
# 4.6 占有数表象：升降算符在等间距能级间“爬梯子”
# --------------------------------------------------------------------------
class C4_Ladder(Scene):
    def construct(self):
        title = title_bar(r"占有数表象：$\hat a^\dagger\ket n=\sqrt{n+1}\ket{n+1}$，$\hat a\ket n=\sqrt n\ket{n-1}$，$\hat N=\hat a^\dagger\hat a$")
        self.play(FadeIn(title))
        X0 = LEFT * 4.2
        ys = [-2.9 + 1.25 * n for n in range(5)]
        lv = VGroup(); lb = VGroup()
        for n, y in enumerate(ys):
            l = Line(X0 + LEFT * 1.4 + UP * y, X0 + RIGHT * 1.4 + UP * y, color=C_E, stroke_width=3)
            t = MT(rf"\ket{{{n}}}", 28).next_to(l, LEFT, buff=0.15)
            e = MT(rf"E_{n}=({n}+\tfrac12)\hbar\omega" if n else r"E_0=\tfrac12\hbar\omega", 24, C_E).next_to(l, RIGHT, buff=0.15)
            lv.add(l); lb.add(VGroup(t, e))
        self.play(LaggedStart(*[Create(l) for l in lv], lag_ratio=0.2), LaggedStart(*[FadeIn(b) for b in lb], lag_ratio=0.2))
        ball = Dot(X0 + UP * ys[0], color=C_HL, radius=0.15)
        self.play(FadeIn(ball))
        for n in range(3):
            ar = CurvedArrow(X0 + RIGHT * 0.5 + UP * ys[n], X0 + RIGHT * 0.5 + UP * ys[n + 1], angle=-PI / 2, color=YELLOW)
            tl = MT(rf"\hat a^\dagger:\ \sqrt{{{n + 1}}}", 24, YELLOW).next_to(ar, RIGHT, buff=0.05).shift(LEFT * 0.2)
            self.play(Create(ar), FadeIn(tl), ball.animate.move_to(X0 + UP * ys[n + 1]), run_time=0.9)
            self.play(FadeOut(ar), FadeOut(tl), run_time=0.3)
        for n in range(3, 0, -1):
            ar = CurvedArrow(X0 + LEFT * 0.5 + UP * ys[n], X0 + LEFT * 0.5 + UP * ys[n - 1], angle=-PI / 2, color=BLUE_C)
            tl = MT(rf"\hat a:\ \sqrt{{{n}}}", 24, BLUE_C).next_to(ar, LEFT, buff=0.05).shift(RIGHT * 0.2)
            self.play(Create(ar), FadeIn(tl), ball.animate.move_to(X0 + UP * ys[n - 1]), run_time=0.9)
            self.play(FadeOut(ar), FadeOut(tl), run_time=0.3)
        z = CT(r"$\hat a\ket0=0$：梯子到底了（零点能不可去除）", 24, RED_C).next_to(lv[0], DOWN, buff=0.2).shift(RIGHT * 0.9)
        self.play(FadeIn(z))
        M = MT(r"a^\dagger=\begin{pmatrix}0&0&0&0&\cdots\\1&0&0&0&\cdots\\0&\sqrt2&0&0&\cdots\\0&0&\sqrt3&0&\cdots\\\vdots&&&\ddots&\end{pmatrix}", 30)
        Mh = MT(r"H=\hbar\omega\left(\hat N+\tfrac12\right)=\hbar\omega\,\mathrm{diag}\left(\tfrac12,\tfrac32,\tfrac52,\dots\right)", 28)
        Mx = MT(r"\hat x=\sqrt{\frac{\hbar}{2m\omega}}\left(\hat a+\hat a^\dagger\right),\quad \hat p=\ii\sqrt{\frac{m\omega\hbar}{2}}\left(\hat a^\dagger-\hat a\right)", 28)
        Mc = MT(r"[\hat a,\hat a^\dagger]=1", 32, C_HL)
        g = VGroup(M, Mh, Mx, Mc).arrange(DOWN, buff=0.35)
        fitw(g, 6.4)
        g.to_corner(UR, buff=0.4).shift(DOWN * 0.9)
        hd = CT(r"能量表象中的矩阵：只有次对角线非零", 24, GREY_B).next_to(g, UP, buff=0.15)
        self.play(FadeIn(hd), Write(M))
        self.play(FadeIn(Mh))
        self.play(FadeIn(Mx), FadeIn(Mc))
        self.wait(3)


# --------------------------------------------------------------------------
# 4.1 坐标表象 vs 动量表象：高斯波包是一对傅里叶变换，宽度此消彼长
# --------------------------------------------------------------------------
class C4_FourierPair(Scene):
    def construct(self):
        title = title_bar(r"同一个态，两种表象：$\psi(x)$ 与 $c(p)$ 互为傅里叶变换（取 $\hbar=1$）")
        self.play(FadeIn(title))
        kw = dict(x_length=5.6, y_length=3.1, tips=False,
                  axis_config={"stroke_width": 2, "include_ticks": False})
        axl = Axes(x_range=[-6, 6, 1], y_range=[-1.0, 1.6, 0.5], **kw).move_to(LEFT * 3.35 + UP * 0.25)
        axr = Axes(x_range=[-6, 6, 1], y_range=[-1.0, 1.6, 0.5], **kw).move_to(RIGHT * 3.35 + UP * 0.25)
        hl = CT(r"坐标表象 $\ip{x}{\psi}=\psi(x)$", 26, C_PSI).next_to(axl, UP, buff=0.12)
        hr = CT(r"动量表象 $\ip{p}{\psi}=c(p)$", 26, C_HL).next_to(axr, UP, buff=0.12)
        xl = MT("x", 26).next_to(axl.x_axis, RIGHT, buff=0.08)
        pl = MT("p", 26).next_to(axr.x_axis, RIGHT, buff=0.08)
        self.play(Create(axl), Create(axr), FadeIn(hl), FadeIn(hr), FadeIn(xl), FadeIn(pl))
        sig = ValueTracker(1.0)   # Δx
        p0 = ValueTracker(0.0)
        xs = np.linspace(-6, 6, 400)

        def left():
            s, k = sig.get_value(), p0.get_value()
            env = (2 * np.pi * s ** 2) ** -0.25 * np.exp(-xs ** 2 / (4 * s ** 2))
            re = env * np.cos(k * xs)
            rho = env ** 2
            return VGroup(filled_curve(axl, xs, rho, C_PROB, 0.3),
                          curve(axl, xs, re, C_PSI, 2.5))

        def right():
            s, k = sig.get_value(), p0.get_value()
            sp = 1 / (2 * s)
            rho = (2 * np.pi * sp ** 2) ** -0.5 * np.exp(-(xs - k) ** 2 / (2 * sp ** 2))
            return filled_curve(axr, xs, rho, C_HL, 0.35)
        L = always_redraw(left)
        R = always_redraw(right)
        leg = VGroup(CT(r"蓝：$\mathrm{Re}\,\psi(x)$　黄：$\lvert\psi(x)\rvert^2$", 22),
                     CT(r"橙：$\lvert c(p)\rvert^2$", 22)).arrange(RIGHT, buff=2.6).next_to(VGroup(axl, axr), DOWN, buff=0.05)
        self.play(Create(L), Create(R), FadeIn(leg))
        dx = DecimalNumber(1.0, num_decimal_places=2, font_size=30, color=C_PROB)
        dp = DecimalNumber(0.5, num_decimal_places=2, font_size=30, color=C_HL)
        pr = DecimalNumber(0.5, num_decimal_places=2, font_size=30, color=C_E)
        dx.add_updater(lambda m: m.set_value(sig.get_value()))
        dp.add_updater(lambda m: m.set_value(0.5 / sig.get_value()))
        pr.add_updater(lambda m: m.set_value(0.5))
        row = VGroup(VGroup(MT(r"\Delta x=", 30, C_PROB), dx).arrange(RIGHT, buff=0.1),
                     VGroup(MT(r"\Delta p=", 30, C_HL), dp).arrange(RIGHT, buff=0.1),
                     VGroup(MT(r"\Delta x\,\Delta p=", 30, C_E), pr, MT(r"=\tfrac{\hbar}{2}", 30, C_E)).arrange(RIGHT, buff=0.1)).arrange(RIGHT, buff=0.7)
        f1 = MT(r"c(p)=\frac{1}{\sqrt{2\pi\hbar}}\int_{-\infty}^{\infty}\psi(x)\,\ee^{-\ii px/\hbar}\,\dd x,\qquad"
                r"\psi(x)=\frac{1}{\sqrt{2\pi\hbar}}\int_{-\infty}^{\infty}c(p)\,\ee^{\ii px/\hbar}\,\dd p", 28)
        fitw(f1, 12.6)
        g = VGroup(row, f1).arrange(DOWN, buff=0.22).to_edge(DOWN, buff=0.2)
        self.play(FadeIn(g))
        self.wait(1)
        self.play(sig.animate.set_value(0.45), run_time=3)
        self.play(sig.animate.set_value(1.8), run_time=3.5)
        self.play(sig.animate.set_value(1.0), run_time=2)
        note = CT(r"乘相因子 $\ee^{\ii p_0x/\hbar}$：$\lvert\psi\rvert^2$ 不变，$c(p)$ 整体平移到 $p_0$", 26, C_HL)
        note.next_to(title, DOWN, buff=0.12)
        self.play(FadeIn(note))
        self.play(p0.animate.set_value(2.5), run_time=3)
        self.wait(1)
        self.play(p0.animate.set_value(-2.0), run_time=3)
        self.wait(1)
        end = CT(r"$\int\lvert\psi\rvert^2\dd x=\int\lvert c\rvert^2\dd p=1$（帕塞瓦尔定理）：换表象不改变态，只改变“看它的角度”", 26, C_E)
        fitw(end, 12.8)
        end.move_to(note)
        self.play(FadeOut(note), FadeIn(end))
        self.wait(3)


# --------------------------------------------------------------------------
# 4.2/4.3 算符的矩阵：F_mn = <m|F|n>；第 n 列 = F|n> 的分量；厄米 ⇔ 关于对角线共轭对称
# --------------------------------------------------------------------------
class C4_OperatorMatrix(Scene):
    def construct(self):
        title = title_bar(r"算符 = 矩阵：$F_{mn}=\langle m|\hat F|n\rangle$　以谐振子能量表象中的 $\hat x$ 为例")
        self.play(FadeIn(title))
        N = 5
        ent = [["0"] * N for _ in range(N)]
        for n in range(N - 1):
            v = "1" if n == 0 else rf"\sqrt{{{n + 1}}}"
            ent[n][n + 1] = v
            ent[n + 1][n] = v
        M = Matrix(ent, element_to_mobject=lambda s: MT(s, 32), h_buff=0.88, v_buff=0.72,
                   bracket_h_buff=0.2)
        pre = MT(r"x=\sqrt{\frac{\hbar}{2m\omega}}", 28)
        M.move_to(LEFT * 2.25 + DOWN * 0.45)
        rows = VGroup(*[MT(rf"\bra{{{m}}}", 26, GREY_B).next_to(M.get_rows()[m], LEFT, buff=0.35).set_x(M.get_left()[0] - 0.4) for m in range(N)])
        rows.shift(LEFT * 0.0)
        cols = VGroup(*[MT(rf"\ket{{{n}}}", 26, GREY_B).next_to(M.get_columns()[n], UP, buff=0.35) for n in range(N)])
        cols.set_y(M.get_top()[1] + 0.3)
        pre.next_to(rows, LEFT, buff=0.25)
        self.play(FadeIn(pre), Create(M.get_brackets()), FadeIn(cols), FadeIn(rows))
        self.play(LaggedStart(*[FadeIn(e) for e in M.get_entries()], lag_ratio=0.03), run_time=1.5)
        RX = 3.65
        t1 = CT(r"由 $\hat x\propto\hat a+\hat a^\dagger$：", 26).move_to(RIGHT * RX + UP * 2.2)
        t1b = MT(r"\langle m|\hat x|n\rangle=\sqrt{\tfrac{\hbar}{2m\omega}}\big(\sqrt n\,\delta_{m,n-1}+\sqrt{n+1}\,\delta_{m,n+1}\big)", 26)
        fitw(t1b, 5.7).next_to(t1, DOWN, buff=0.15)
        self.play(FadeIn(t1), FadeIn(t1b))
        # 第 n=1 列
        colbox = SurroundingRectangle(M.get_columns()[1], color=YELLOW, buff=0.12)
        t2 = CT(r"第 $n$ 列 = $\hat x\ket n$ 在各基矢上的分量", 26, YELLOW)
        t2b = MT(r"\hat x\ket1=\sqrt{\tfrac{\hbar}{2m\omega}}\left(\ket0+\sqrt2\,\ket2\right)", 28, YELLOW)
        fitw(VGroup(t2, t2b).arrange(DOWN, buff=0.15), 5.7).move_to(RIGHT * RX + UP * 0.6)
        self.play(Create(colbox), FadeIn(t2), FadeIn(t2b))
        self.wait(1.5)
        # 厄米：镜像对称
        diag = DashedLine(M.get_entries()[0].get_center() + UL * 0.3, M.get_entries()[N * N - 1].get_center() + DR * 0.3, color=C_HL)
        a = M.get_entries()[0 * N + 1]; b = M.get_entries()[1 * N + 0]
        c = M.get_entries()[2 * N + 3]; d = M.get_entries()[3 * N + 2]
        pairs = VGroup(*[SurroundingRectangle(e, color=C_HL, buff=0.08) for e in (a, b, c, d)])
        t3 = CT(r"厄米：$F_{mn}=F_{nm}^{*}$，关于对角线共轭对称", 26, C_HL)
        t3b = CT(r"$\Rightarrow$ 本征值为实数，可用幺正变换对角化", 24, C_HL)
        fitw(VGroup(t3, t3b).arrange(DOWN, buff=0.12), 5.7).move_to(RIGHT * RX + DOWN * 0.9)
        self.play(FadeOut(colbox), Create(diag), Create(pairs), FadeIn(t3), FadeIn(t3b))
        self.play(Swap(a, b), Swap(c, d), run_time=1.5)
        self.wait(1)
        # 换成 H：自身表象中对角
        hent = [["0"] * N for _ in range(N)]
        for n in range(N):
            hent[n][n] = str(2 * n + 1)
        H = Matrix(hent, element_to_mobject=lambda s: MT(s, 32, C_E), h_buff=0.88, v_buff=0.72, bracket_h_buff=0.2)
        H.move_to(M)
        preH = MT(r"H=\frac{\hbar\omega}{2}", 30, C_E).move_to(pre, aligned_edge=RIGHT)
        t4 = CT(r"算符在自身表象中是对角矩阵，对角元 = 本征值", 26, C_E)
        t4b = MT(r"\langle m|\hat H|n\rangle=\left(n+\tfrac12\right)\hbar\omega\,\delta_{mn}", 28, C_E)
        fitw(VGroup(t4, t4b).arrange(DOWN, buff=0.12), 5.7).move_to(RIGHT * RX + DOWN * 2.4)
        self.play(FadeOut(pairs), FadeOut(diag), ReplacementTransform(M, H), ReplacementTransform(pre, preH))
        self.play(FadeIn(t4), FadeIn(t4b))
        self.wait(1.5)
        fin = CT(r"期望值 $\langle F\rangle=\Psi^\dagger F\Psi$，本征方程 $F\Psi=\lambda\Psi$，薛定谔方程 $\ii\hbar\dot\Psi=H\Psi$", 26, WHITE)
        fitw(fin, 12.6)
        bx = boxed(fin).to_edge(DOWN, buff=0.15)
        self.play(FadeIn(bx))
        self.wait(3)
