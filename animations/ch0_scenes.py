"""预备 原子 —— manim 动画源码：汤姆孙模型与卢瑟福模型的 alpha 粒子散射对照."""
from common import *  # noqa: F401,F403


def trajectory(b, k, R, x0=-6.0, x1=6.0, n=4000, model="point"):
    """二维库仑斥力下的轨迹（RK4）。k=Z1Z2e^2/(4πε0)/(m v0^2) 的等效长度（=a/2）。
    model="point": 点电荷；model="sphere": 半径 R 的均匀带电球（外部同点电荷，内部线性力）。"""
    def acc(p):
        r = np.hypot(*p)
        if model == "sphere" and r < R:
            return k * p / R ** 3
        return k * p / max(r, 1e-6) ** 3
    p = np.array([x0, b]); v = np.array([1.0, 0.0])
    dt = (x1 - x0) * 1.6 / n
    pts = [p.copy()]
    for _ in range(n):
        k1v = acc(p); k1p = v
        k2v = acc(p + 0.5 * dt * k1p); k2p = v + 0.5 * dt * k1v
        k3v = acc(p + 0.5 * dt * k2p); k3p = v + 0.5 * dt * k2v
        k4v = acc(p + dt * k3p); k4p = v + dt * k3v
        p = p + dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
        v = v + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
        pts.append(p.copy())
        if abs(p[0]) > 6.2 or abs(p[1]) > 3.6:
            break
    return np.array(pts), np.degrees(np.arctan2(v[1], v[0]))


class C0_Rutherford(Scene):
    def construct(self):
        title = title_bar(r"同样的总电荷，两种分布：汤姆孙“布丁”只能让 $\alpha$ 粒子轻轻偏一下，点状原子核能把它弹回来")
        self.play(FadeIn(title))
        S = 0.55  # 画面缩放：模型坐标 -> 屏幕
        R = 2.2   # 原子半径（模型单位）
        k = 0.06  # = a/2，a 为对头碰最近距离
        left_c, right_c = LEFT * 3.5 + DOWN * 0.5, RIGHT * 3.5 + DOWN * 0.5
        tl = CT("汤姆孙模型：正电荷均匀铺满原子", 26).move_to(left_c + UP * 2.55)
        tr = CT("卢瑟福模型：正电荷集中在核里", 26).move_to(right_c + UP * 2.55)
        ball = Circle(radius=R * S, color=RED_B, fill_opacity=0.25, stroke_width=2).move_to(left_c)
        atom = DashedVMobject(Circle(radius=R * S, color=GREY_B, stroke_width=2), num_dashes=40).move_to(right_c)
        nuc = Dot(right_c, radius=0.07, color=RED)
        es_l = VGroup(*[Dot(left_c + S * np.array([x, y, 0]), radius=0.04, color=BLUE_B)
                        for x, y in [(-1.2, 0.6), (0.3, 1.3), (1.1, -0.4), (-0.4, -1.1), (0.6, 0.3), (-1.4, -0.5)]])
        es_r = es_l.copy().shift(right_c - left_c)
        self.play(FadeIn(tl), FadeIn(tr), FadeIn(ball), FadeIn(es_l), Create(atom), FadeIn(nuc), FadeIn(es_r))
        bs = [-1.6, -0.9, -0.45, -0.18, -0.06, 0.02, 0.18, 0.45, 0.9, 1.6]
        paths_l, paths_r, angs = VGroup(), VGroup(), []
        for b in bs:
            for model, c, grp in [("sphere", left_c, paths_l), ("point", right_c, paths_r)]:
                pts, ang = trajectory(b, k, R, x0=-5.2, x1=5.2, model=model)
                pts = pts[(np.abs(pts[:, 0]) < 5.2) & (np.abs(pts[:, 1]) < 3.4)]
                m = VMobject(stroke_width=2.2, color=YELLOW if model == "point" else ORANGE)
                m.set_points_smoothly([c + S * np.array([x, y, 0]) for x, y in pts[::25]])
                grp.add(m)
                if model == "point":
                    angs.append((b, ang))
        self.play(LaggedStart(*[Create(p) for p in paths_l], lag_ratio=0.15), run_time=4)
        nl = CT(r"最大偏转角 $\sim\dfrac{a}{R}$：金箔、$5\,\mathrm{MeV}$ 时约 $5\times10^{-4}\,\mathrm{rad}$" "\n" r"（图中把偏转放大了几百倍）", 22, ORANGE)
        fitw(nl, 6.2).move_to(left_c + DOWN * 2.45)
        self.play(FadeIn(nl))
        self.play(LaggedStart(*[Create(p) for p in paths_r], lag_ratio=0.15), run_time=5)
        nr = fitw(CT(r"$b$ 越小偏得越狠：$b=\dfrac a2\cot\dfrac\theta2$", 24, YELLOW), 6.2).move_to(right_c + DOWN * 2.45)
        self.play(FadeIn(nr))
        self.wait(1.5)
        # 第二幕：瞄准距离与散射角
        self.play(*[FadeOut(m) for m in [paths_l, paths_r, nl, nr, ball, es_l, tl, es_r, atom, tr]],
                  nuc.animate.move_to(LEFT * 2.2 + DOWN * 0.3))
        c0 = nuc.get_center()
        bt = ValueTracker(1.2)
        S2 = 0.9

        def traj():
            b = bt.get_value()
            pts, ang = trajectory(b, 0.12, R, x0=-5.0, x1=5.0, model="point")
            pts = pts[(np.abs(pts[:, 0]) < 5.0) & (np.abs(pts[:, 1]) < 3.2)]
            m = VMobject(stroke_width=3, color=YELLOW)
            m.set_points_smoothly([c0 + S2 * np.array([x, y, 0]) for x, y in pts[::20]])
            bl = DoubleArrow(c0 + S2 * np.array([-4.2, 0, 0]), c0 + S2 * np.array([-4.2, b, 0]), buff=0,
                             stroke_width=2, tip_length=0.12, color=BLUE_B)
            axis = DashedLine(c0 + S2 * np.array([-5, 0, 0]), c0 + S2 * np.array([2.6, 0, 0]), color=GREY_B, stroke_width=1.5)
            th = 2 * np.degrees(np.arctan(0.12 / b))
            lab = MT(rf"b={b / 0.24:.2f}\,a,\ \ \theta={th:.0f}^\circ", 30, YELLOW).to_corner(UR, buff=0.4).shift(DOWN * 0.8)
            return VGroup(axis, m, bl, lab)
        T = always_redraw(traj)
        self.play(FadeIn(T))
        form = boxed(panel(MT(r"b=\frac a2\cot\frac\theta2,\qquad a=\frac{Z_1Z_2e^2}{4\pi\varepsilon_0E}", 30),
                           r"金核、$E=7.68\,\mathrm{MeV}$：$a\approx30\,\mathrm{fm}$", r"$b=10\,\mathrm{fm}\Rightarrow\theta=112^\circ$；$b=1000\,\mathrm{fm}\Rightarrow\theta=1.7^\circ$", size=24))
        form.to_corner(DR, buff=0.3)
        self.play(FadeIn(form))
        self.play(bt.animate.set_value(0.08), run_time=7, rate_func=smooth)
        self.wait(0.5)
        self.play(bt.animate.set_value(0.5), run_time=3)
        end = CT(r"大角散射只来自 $b\lesssim10^{-14}\,\mathrm m$ 的极少数粒子 $\Rightarrow$ 核比原子小 $10^4$ 倍", 26, C_HL)
        fitw(end, 6.6).to_corner(DL, buff=0.3)
        self.play(FadeIn(end))
        self.wait(2)
