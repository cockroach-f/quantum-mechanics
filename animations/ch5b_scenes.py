"""第 5 章 光的发射和吸收 —— 激光：粒子数反转与受激发射的雪崩（manim，LaTeX 排版）."""
from common import *  # noqa: F401,F403

C_EXC = ORANGE
C_GND = GREY_B


def bar(height, color, w=0.5):
    return Rectangle(width=max(height, 0.02), height=w, fill_color=color, fill_opacity=0.85, stroke_width=0)


class C5_Laser(Scene):
    def construct(self):
        title = title_bar(r"激光：先把上能级的原子数抽得比下能级多（粒子数反转），受激发射才能一变二、二变四")
        self.play(FadeIn(title))
        # ---------------- 第一幕：热平衡
        h1 = CT(r"热平衡：$\dfrac{N_2}{N_1}=\ee^{-h\nu/k_BT}$", 28).to_corner(UL, buff=0.35).shift(DOWN * 0.8)
        lv1, lv2 = Line(LEFT * 1.8, RIGHT * 0.2).shift(LEFT * 3.6 + DOWN * 1.8), Line(LEFT * 1.8, RIGHT * 0.2).shift(LEFT * 3.6 + UP * 0.6)
        lb = VGroup(MT("E_1", 26).next_to(lv1, LEFT), MT("E_2", 26).next_to(lv2, LEFT))
        N1 = ValueTracker(4.0); N2 = ValueTracker(0.15)
        b1 = always_redraw(lambda: bar(N1.get_value(), C_GND).next_to(lv1, RIGHT, buff=0.3).align_to(lv1, DOWN).shift(UP * 0.0))
        b2 = always_redraw(lambda: bar(N2.get_value(), C_EXC).next_to(lv2, RIGHT, buff=0.3).align_to(lv2, DOWN))
        n1l = always_redraw(lambda: MT(r"N_1", 24, C_GND).next_to(b1, RIGHT, buff=0.1))
        n2l = always_redraw(lambda: MT(r"N_2", 24, C_EXC).next_to(b2, RIGHT, buff=0.1))
        self.play(FadeIn(h1), Create(lv1), Create(lv2), FadeIn(lb), FadeIn(b1), FadeIn(b2), FadeIn(n1l), FadeIn(n2l))
        r1 = boxed(panel(r"吸收速率 $\propto B_{12}\rho N_1$，受激发射 $\propto B_{21}\rho N_2$，$B_{12}=B_{21}$",
                         r"$N_1>N_2$：吸收总是赢，光穿过介质被削弱",
                         r"可见光、$300\,\mathrm K$：$h\nu/k_BT\approx 80$，$N_2/N_1\sim10^{-35}$", size=22))
        fitw(r1, 6.6).to_corner(UR, buff=0.3).shift(DOWN * 0.9)
        self.play(FadeIn(r1))
        self.wait(2.5)
        # ---------------- 第二幕：三能级抽运
        lv3 = Line(LEFT * 1.8, RIGHT * 0.2).shift(LEFT * 3.6 + UP * 2.0)
        lb3 = MT("E_3", 26).next_to(lv3, LEFT)
        h2 = boxed(panel(r"三能级（红宝石，$694.3\,\mathrm{nm}$）：", r"$E_1\to E_3$ 抽运；$E_3\to E_2$ 快速无辐射弛豫；", r"$E_2$ 是寿命约 $3\,\mathrm{ms}$ 的亚稳态，原子在这里堆积", size=24))
        fitw(h2, 6.6).move_to(r1)
        self.play(FadeOut(r1), FadeOut(h1), Create(lv3), FadeIn(lb3), FadeIn(h2), lv2.animate.shift(UP * 0.4), lb[1].animate.shift(UP * 0.4))
        pump = Arrow(lv1.get_center() + RIGHT * 0.3 + UP * 0.05, lv3.get_center() + RIGHT * 0.3, buff=0, color=GREEN_C, stroke_width=5)
        relax = CurvedArrow(lv3.get_center() + LEFT * 0.5, lv2.get_center() + LEFT * 0.5, color=GREY_B, angle=0.6)
        pl = CT("抽运光", 20, GREEN_C).next_to(pump, LEFT, buff=0.05)
        rl = CT("快", 20, GREY_B).next_to(relax, LEFT, buff=0.05)
        self.play(GrowArrow(pump), FadeIn(pl))
        self.play(Create(relax), FadeIn(rl))
        self.play(N1.animate.set_value(1.6), N2.animate.set_value(2.5), run_time=3)
        inv = CT(r"$N_2>N_1$：粒子数反转", 28, C_EXC).next_to(h2, DOWN, buff=0.4)
        self.play(FadeIn(inv))
        self.wait(1.5)
        self.play(*[FadeOut(m) for m in [lv1, lv2, lv3, lb, lb3, b1, b2, n1l, n2l, pump, relax, pl, rl, h2, inv]])
        # ---------------- 第三幕：谐振腔中的雪崩
        mL = Rectangle(width=0.15, height=2.6, fill_color=GREY_A, fill_opacity=1, stroke_width=0).move_to(LEFT * 6.0 + DOWN * 0.6)
        mR = Rectangle(width=0.15, height=2.6, fill_color=GREY_A, fill_opacity=0.5, stroke_width=0).move_to(RIGHT * 6.0 + DOWN * 0.6)
        ml = VGroup(CT("全反射镜", 20).next_to(mL, DOWN), CT("部分透射镜", 20).next_to(mR, DOWN))
        xs = np.linspace(-4.8, 4.8, 9)
        atoms = VGroup(*[Dot([x, -0.6 + (0.5 if i % 2 else -0.5), 0], radius=0.13, color=C_EXC) for i, x in enumerate(xs)])
        self.play(FadeIn(mL), FadeIn(mR), FadeIn(ml), FadeIn(atoms))
        cnt = ValueTracker(1)
        cl = always_redraw(lambda: CT(rf"腔内同模光子数 $n=2^{{{int(round(np.log2(max(cnt.get_value(), 1))))}}}$", 26, YELLOW).to_corner(UL, buff=0.35).shift(DOWN * 0.8))
        self.add(cl)
        note = CT(r"受激发射的光子与入射光子同频率、同相位、同方向、同偏振", 22).to_edge(UP, buff=0.95).shift(RIGHT * 2.2)
        fitw(note, 7.5)
        self.play(FadeIn(note))

        def packet(x, n):
            amp = min(0.12 + 0.035 * np.log2(max(n, 1)), 0.6)
            g = FunctionGraph(lambda t: amp * np.sin(10 * t), x_range=[-0.7, 0.7], color=YELLOW, stroke_width=min(2 + 0.4 * np.log2(max(n, 1)), 7))
            return g.move_to([x, -0.6, 0])
        n = 1
        direction = 1
        for sweep in range(2):
            order = range(len(xs)) if direction > 0 else range(len(xs) - 1, -1, -1)
            x_prev = -5.6 if direction > 0 else 5.6
            P = packet(x_prev, n)
            self.add(P)
            for i in order:
                x = xs[i]
                self.play(P.animate.move_to([x, -0.6, 0]), run_time=0.35, rate_func=linear)
                if atoms[i].get_color() != ManimColor(C_GND):
                    n *= 2
                    newP = packet(x, n)
                    self.play(atoms[i].animate.set_color(C_GND), Transform(P, newP), cnt.animate.set_value(n), run_time=0.3)
            xe = 5.6 if direction > 0 else -5.6
            self.play(P.animate.move_to([xe, -0.6, 0]), run_time=0.35, rate_func=linear)
            if direction > 0:
                out = packet(6.6, max(n // 8, 1)).set_opacity(0.8)
                self.play(FadeIn(out), out.animate.shift(RIGHT * 0.4), run_time=0.4)
                self.remove(out)
            self.remove(P)
            # 抽运把原子重新激发
            self.play(*[a.animate.set_color(C_EXC) for a in atoms], run_time=0.5)
            direction *= -1
        concl = boxed(panel(r"$\bullet$ 增益：每经过一个激发态原子，光子数按比例增加（$N_2>N_1$ 时净增益）",
                            r"$\bullet$ 谐振腔：只让沿轴、满足驻波条件的模式来回放大",
                            r"$\bullet$ 自发发射 / 受激发射 $=\ee^{h\nu/k_BT}-1$：没有反转和腔，可见光几乎全是自发发射", size=22))
        fitw(concl, 11.5).to_edge(DOWN, buff=0.25)
        self.play(FadeIn(concl))
        self.wait(3)
