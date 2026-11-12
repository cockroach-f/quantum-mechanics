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
