"""第 7 章 自旋 —— manim 动画源码（LaTeX 排版）。

渲染（720p30）：
  cd animations && manim -qm ch7_scenes.py C7_SternGerlach
输出文件复制为 docs/assets/videos/ch7_*.mp4，封面为 docs/assets/img/ch7_poster_*.jpg。
"""
from common import *  # noqa: F401,F403

C_UP = RED_C
C_DN = BLUE_C
C_B = TEAL_C


# --------------------------------------------------------------------------
# 7.1 施特恩–格拉赫实验：经典预期（连续带） vs 实验结果（两条线）
# --------------------------------------------------------------------------
class C7_SternGerlach(Scene):
    def construct(self):
        title = title_bar(r"施特恩–格拉赫实验：$s$ 态氢原子束在非均匀磁场中只分成两束")
        self.play(FadeIn(title))
        y0 = -0.35
        xo, x1, x2, xs = -6.0, -1.6, 1.4, 5.4
        L = x2 - x1
        oven = VGroup(Rectangle(width=0.8, height=0.9, color=GREY_B, fill_opacity=0.3).move_to([xo, y0, 0]),
                      MT(r"K", 26).move_to([xo, y0, 0]))
        slit = VGroup(Line([-3.7, y0 + 0.12, 0], [-3.7, y0 + 0.8, 0], stroke_width=6, color=GREY_B),
                      Line([-3.7, y0 - 0.12, 0], [-3.7, y0 - 0.8, 0], stroke_width=6, color=GREY_B))
        north = Polygon([x1, y0 + 2.2, 0], [x2, y0 + 2.2, 0], [x2, y0 + 0.75, 0], [(x1 + x2) / 2, y0 + 0.35, 0], [x1, y0 + 0.75, 0],
                        color=RED_D, fill_opacity=0.55, stroke_width=2)
        south = Polygon([x1, y0 - 2.2, 0], [x2, y0 - 2.2, 0], [x2, y0 - 0.55, 0], [(x1 + x2) / 2 + 0.35, y0 - 0.55, 0],
                        [(x1 + x2) / 2, y0 - 0.9, 0], [(x1 + x2) / 2 - 0.35, y0 - 0.55, 0], [x1, y0 - 0.55, 0],
                        color=BLUE_D, fill_opacity=0.55, stroke_width=2)
        nl = MT(r"\mathrm N", 30).move_to(north.get_center() + UP * 0.4)
        sl = MT(r"\mathrm S", 30).move_to(south.get_center() + DOWN * 0.4)
        screen = Line([xs, y0 - 2.4, 0], [xs, y0 + 2.4, 0], color=GREY_A, stroke_width=5)
        scl = CT(r"屏 P", 24, GREY_A).next_to(screen, DOWN, buff=0.1)
        zax = VGroup(Arrow([6.3, y0 - 0.6, 0], [6.3, y0 + 0.9, 0], buff=0, stroke_width=3, color=GREY_B),
                     MT(r"z", 26, GREY_B).move_to([6.6, y0 + 0.8, 0]))
        self.play(FadeIn(oven), FadeIn(slit), FadeIn(north), FadeIn(south), FadeIn(nl), FadeIn(sl),
                  Create(screen), FadeIn(scl), FadeIn(zax))
        form = boxed(panel(r"$U=-\vv\mu\cdot\vv B=-\mu B_z\cos\theta$",
                           r"$F_z=-\dfrac{\partial U}{\partial z}=\mu\,\dfrac{\partial B_z}{\partial z}\cos\theta$", size=24))
        form.to_corner(UL, buff=0.25).shift(DOWN * 0.75)
        self.play(FadeIn(form))

        def path(d):
            a = d / (L * L / 2 + L * (xs - x2))
            pts = [[xo + 0.4, y0, 0], [x1, y0, 0]]
            for x in np.linspace(x1, x2, 16)[1:]:
                pts.append([x, y0 + a * (x - x1) ** 2 / 2, 0])
            pts.append([xs, y0 + a * L * L / 2 + a * L * (xs - x2), 0])
            m = VMobject(); m.set_points_as_corners([np.array(p) for p in pts])
            return m

        D = 1.7
        # ---- classical
        cl = CT(r"经典预期：$\vv\mu$ 取向任意，$\cos\theta\in[-1,1]$ 连续 $\Rightarrow$ 屏上一条连续的带", 26, YELLOW)
        cl.to_edge(DOWN, buff=0.3)
        self.play(FadeIn(cl))
        rng = np.random.default_rng(3)
        cs = rng.uniform(-1, 1, 36)
        dots = [Dot(radius=0.05, color=YELLOW).move_to([xo + 0.4, y0, 0]) for _ in cs]
        self.play(LaggedStart(*[MoveAlongPath(d, path(D * c), rate_func=linear) for d, c in zip(dots, cs)],
                              lag_ratio=0.06, run_time=5))
        band = Rectangle(width=0.14, height=2 * D, color=YELLOW, fill_opacity=0.6, stroke_width=0).move_to([xs, y0, 0])
        self.play(FadeIn(band), *[FadeOut(d) for d in dots])
        self.wait(0.6)
        # ---- quantum
        qm = CT(r"实验结果：只有两条分立的线，对应 $\cos\theta=\pm1$——磁矩是\textbf{空间量子化}的", 26, C_HL)
        qm.to_edge(DOWN, buff=0.3)
        self.play(FadeOut(cl), FadeIn(qm), band.animate.set_opacity(0.15))
        cs2 = np.array([1, -1] * 18)
        dots2 = [Dot(radius=0.05, color=C_UP if c > 0 else C_DN).move_to([xo + 0.4, y0, 0]) for c in cs2]
        self.play(LaggedStart(*[MoveAlongPath(d, path(D * c), rate_func=linear) for d, c in zip(dots2, cs2)],
                              lag_ratio=0.06, run_time=5))
        spots = VGroup(Ellipse(width=0.22, height=0.32, color=C_UP, fill_opacity=0.9).move_to([xs, y0 + D, 0]),
                       Ellipse(width=0.22, height=0.32, color=C_DN, fill_opacity=0.9).move_to([xs, y0 - D, 0]))
        self.play(FadeIn(spots), *[FadeOut(d) for d in dots2], FadeOut(band))
        lab = VGroup(MT(r"\mu_z=-\mu_{\mathrm B}\ (S_z=+\tfrac{\hbar}{2})", 22, C_UP).next_to(spots[0], LEFT, buff=0.15),
                     MT(r"\mu_z=+\mu_{\mathrm B}\ (S_z=-\tfrac{\hbar}{2})", 22, C_DN).next_to(spots[1], LEFT, buff=0.15))
        lab[0].shift(UP * 0.35); lab[1].shift(DOWN * 0.35)
        self.play(FadeIn(lab))
        concl = boxed(panel(r"$s$ 态：$l=0$，没有轨道磁矩 $\Rightarrow$ 磁矩来自电子本身",
                            r"$\mu_{sz}=\pm\mu_{\mathrm B}$，$\mu_{\mathrm B}=\dfrac{e\hbar}{2m_{\mathrm e}}$", size=24, color=C_HL))
        concl.to_corner(UL, buff=0.25).shift(DOWN * 0.75)
        self.play(FadeOut(form), FadeIn(concl))
        self.wait(3)


# --------------------------------------------------------------------------
# 7.2 布洛赫球与自旋在磁场中的进动
# --------------------------------------------------------------------------
class C7_SpinPrecession(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=68 * DEGREES, theta=35 * DEGREES, zoom=0.95)
        R = 2.0
        title = title_bar(r"布洛赫球：自旋态 $\chi=\big(\cos\tfrac\theta2,\ \ee^{\ii\varphi}\sin\tfrac\theta2\big)^{\mathrm T}$ 对应球面上一点")
        self.add_fixed_in_frame_mobjects(title)
        self.play(FadeIn(title))
        sph = Sphere(radius=R, resolution=(18, 36), fill_opacity=0.08, stroke_width=0.4,
                     checkerboard_colors=[BLUE_E, BLUE_E], stroke_color=GREY_D)
        axes = VGroup(*[Line3D(-R * 1.25 * v, R * 1.25 * v, thickness=0.01, color=GREY_B) for v in (RIGHT, UP, OUT)])
        eq = Circle(radius=R, color=GREY_C, stroke_width=1.5)
        self.play(FadeIn(sph), FadeIn(axes), Create(eq))
        labs = [(MT(r"\ket{\uparrow}", 28, C_UP), OUT * (R + 0.4)), (MT(r"\ket{\downarrow}", 28, C_DN), IN * (R + 0.4)),
                (MT(r"\ket{+x}", 26), RIGHT * (R + 0.55)), (MT(r"\ket{+y}", 26), UP * (R + 0.55)),
                (MT(r"\ket{-x}", 22, GREY_B), LEFT * (R + 0.45))]
        for m, p in labs:
            m.move_to(p)
            self.add_fixed_orientation_mobjects(m)
        th0, ph0 = 55 * DEGREES, 0.0
        tt = ValueTracker(th0); pp = ValueTracker(ph0)

        def tip():
            t, p = tt.get_value(), pp.get_value()
            return R * np.array([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)])
        vec = always_redraw(lambda: Arrow3D(ORIGIN, tip(), color=YELLOW, thickness=0.025, resolution=8))
        self.play(FadeIn(vec))
        info = boxed(panel(r"$\langle\hat{\vv S}\rangle=\tfrac\hbar2(\sin\theta\cos\varphi,\ \sin\theta\sin\varphi,\ \cos\theta)$",
                           r"测 $S_z$：$P(\uparrow)=\cos^2\tfrac\theta2$，$P(\downarrow)=\sin^2\tfrac\theta2$", size=22))
        info.to_corner(UL, buff=0.2).shift(DOWN * 0.7)
        self.add_fixed_in_frame_mobjects(info)
        self.play(FadeIn(info))
        # 扫过几个特殊态
        self.play(tt.animate.set_value(PI / 2), run_time=1.5)
        self.play(pp.animate.set_value(PI / 2), run_time=1.5)
        self.play(tt.animate.set_value(th0), pp.animate.set_value(0), run_time=1.5)
        # 磁场沿 z
        Bz = Arrow([3.4, -0.8, 0], [3.4, 1.6, 0], buff=0, color=C_B, stroke_width=6)
        Bl = CT(r"$\vv B\parallel z$", 28, C_B).next_to(Bz, RIGHT, buff=0.15)
        self.add_fixed_in_frame_mobjects(Bz, Bl)
        self.play(FadeIn(Bz), FadeIn(Bl))
        ham = boxed(panel(r"$\hat H=\dfrac{e}{m_{\mathrm e}}B\hat S_z=\dfrac{\hbar\omega}{2}\sigma_z$，$\omega=\dfrac{eB}{m_{\mathrm e}}$",
                          r"$\chi(t)=\big(\cos\tfrac\theta2\,\ee^{-\ii\omega t/2},\ \sin\tfrac\theta2\,\ee^{\ii\omega t/2}\big)^{\mathrm T}$",
                          r"$\Rightarrow$ $\varphi(t)=\varphi_0+\omega t$，$\theta$ 不变", size=22, color=C_B))
        ham.to_corner(DR, buff=0.2)
        self.add_fixed_in_frame_mobjects(ham)
        self.play(FadeIn(ham))
        trace = always_redraw(lambda: ParametricFunction(
            lambda s: R * np.array([np.sin(tt.get_value()) * np.cos(s), np.sin(tt.get_value()) * np.sin(s), np.cos(tt.get_value())]),
            t_range=[0, max(pp.get_value(), 1e-3)], color=YELLOW, stroke_width=3))
        self.add(trace)
        self.begin_ambient_camera_rotation(rate=0.08)
        self.play(pp.animate.set_value(2 * TAU), run_time=8, rate_func=linear)
        self.stop_ambient_camera_rotation()
        note = CT(r"$\langle\hat{\vv S}\rangle$ 绕 $\vv B$ 进动（经典：$\dd\vv S/\dd t=\vv\mu\times\vv B$）；$P(\uparrow)$ 始终为 $\cos^2\tfrac\theta2$", 24, YELLOW)
        fitw(note, 7.4)
        note.to_corner(DL, buff=0.3)
        self.add_fixed_in_frame_mobjects(note)
        self.play(FadeIn(note))
        self.wait(2)
        # 磁场沿 x：自旋翻转
        self.remove(trace)
        self.play(FadeOut(Bz), FadeOut(Bl), FadeOut(ham), FadeOut(note))
        self.play(tt.animate.set_value(0.001), pp.animate.set_value(PI / 2), run_time=1.2)
        Bx = Arrow3D(ORIGIN, RIGHT * (R + 1.0), color=C_B, thickness=0.03).shift(IN * 0.0 + DOWN * 2.6)
        self.play(FadeIn(Bx))
        ham2 = boxed(panel(r"$\vv B$ 沿 $x$，初态 $\ket{\uparrow}$：",
                           r"$P_\downarrow(t)=\sin^2\dfrac{\omega t}{2}$（自旋翻转）", size=22, color=C_B))
        ham2.to_corner(DR, buff=0.2)
        self.add_fixed_in_frame_mobjects(ham2)
        self.play(FadeIn(ham2))
        # 绕 x 轴转：在 yz 平面内的大圆，theta 从 0 到 pi（phi=±pi/2）
        ang = ValueTracker(0)

        def tip2():
            a = ang.get_value()
            return R * np.array([0.0, -np.sin(a), np.cos(a)])
        vec2 = always_redraw(lambda: Arrow3D(ORIGIN, tip2(), color=YELLOW, thickness=0.025, resolution=8))
        tr2 = always_redraw(lambda: ParametricFunction(lambda s: R * np.array([0.0, -np.sin(s), np.cos(s)]),
                                                       t_range=[0, max(ang.get_value(), 1e-3)], color=YELLOW, stroke_width=3))
        self.remove(vec); self.add(vec2, tr2)
        self.play(ang.animate.set_value(TAU), run_time=6, rate_func=linear)
        self.wait(2)


# --------------------------------------------------------------------------
# 7.2 连续的施特恩–格拉赫装置：测量会改变态
# --------------------------------------------------------------------------
class C7_SequentialSG(Scene):
    def sg(self, ax, pos, col):
        r = RoundedRectangle(width=1.25, height=0.9, corner_radius=0.12, color=col, fill_opacity=0.2).move_to(pos)
        t = CT(rf"SG-${ax}$", 24, col).move_to(pos)
        return VGroup(r, t)

    def beam(self, a, b, w, col=YELLOW):
        return Arrow(a, b, buff=0.05, stroke_width=max(1.5, 7 * w), color=col, max_tip_length_to_length_ratio=0.08)

    def construct(self):
        title = title_bar(r"连续的施特恩–格拉赫实验：测量 $S_x$ 会“抹掉”$S_z$ 的信息")
        self.play(FadeIn(title))
        rows = [2.0, -0.05, -2.1]
        N = 1000
        cZ, cX = C_UP, GREEN_C
        # ---------- row 1
        y = rows[0]
        src = CT(rf"$N={N}$", 22).move_to([-6.2, y, 0])
        a1 = self.sg("z", [-4.3, y, 0], cZ)
        b0 = self.beam([-5.6, y, 0], a1.get_left(), 1)
        up = self.beam(a1.get_right() + UP * 0.2, [-2.1, y + 0.45, 0], 0.5, C_UP)
        dn = self.beam(a1.get_right() + DOWN * 0.2, [-2.6, y - 0.45, 0], 0.5, C_DN)
        blk = Line([-2.5, y - 0.7, 0], [-2.5, y - 0.25, 0], stroke_width=8, color=GREY_B)
        a2 = self.sg("z", [-1.0, y + 0.45, 0], cZ)
        up2 = self.beam(a2.get_right() + UP * 0.1, [1.2, y + 0.75, 0], 0.5, C_UP)
        n1 = MT(r"\uparrow:500", 22, C_UP).next_to(up2, RIGHT, buff=0.1)
        n2 = MT(r"\downarrow:0", 22, C_DN).move_to([1.6, y + 0.15, 0])
        r1 = VGroup(CT(r"$\Rightarrow$ 重复测量结果不变", 24, GREY_A),
                    CT(r"$\chi_{+}$ 是 $\hat S_z$ 的本征态", 24, GREY_A)).arrange(DOWN, aligned_edge=LEFT, buff=0.12).move_to([4.7, y + 0.45, 0])
        self.play(FadeIn(src), GrowArrow(b0), FadeIn(a1))
        self.play(GrowArrow(up), GrowArrow(dn), FadeIn(blk))
        self.play(FadeIn(a2), GrowArrow(up2), FadeIn(n1), FadeIn(n2))
        self.play(FadeIn(r1))
        # ---------- row 2
        y = rows[1]
        src = CT(rf"$N={N}$", 22).move_to([-6.2, y, 0])
        a1 = self.sg("z", [-4.3, y, 0], cZ)
        b0 = self.beam([-5.6, y, 0], a1.get_left(), 1)
        up = self.beam(a1.get_right() + UP * 0.2, [-2.1, y + 0.45, 0], 0.5, C_UP)
        dn = self.beam(a1.get_right() + DOWN * 0.2, [-2.6, y - 0.45, 0], 0.5, C_DN)
        blk = Line([-2.5, y - 0.7, 0], [-2.5, y - 0.25, 0], stroke_width=8, color=GREY_B)
        a2 = self.sg("x", [-1.0, y + 0.45, 0], cX)
        px = self.beam(a2.get_right() + UP * 0.1, [1.2, y + 0.75, 0], 0.25, cX)
        mx = self.beam(a2.get_right() + DOWN * 0.1, [1.2, y + 0.15, 0], 0.25, GREEN_E)
        n1 = MT(r"+x:250", 22, cX).next_to(px, RIGHT, buff=0.1)
        n2 = MT(r"-x:250", 22, GREEN_E).next_to(mx, RIGHT, buff=0.1)
        r2 = MT(r"|\langle\pm x|\!\uparrow\rangle|^2=\tfrac12", 26, GREY_A).move_to([4.6, y + 0.45, 0])
        self.play(FadeIn(src), GrowArrow(b0), FadeIn(a1), GrowArrow(up), GrowArrow(dn), FadeIn(blk))
        self.play(FadeIn(a2), GrowArrow(px), GrowArrow(mx), FadeIn(n1), FadeIn(n2))
        self.play(FadeIn(r2))
        # ---------- row 3
        y = rows[2]
        src = CT(rf"$N={N}$", 22).move_to([-6.2, y, 0])
        a1 = self.sg("z", [-4.3, y, 0], cZ)
        b0 = self.beam([-5.6, y, 0], a1.get_left(), 1)
        up = self.beam(a1.get_right() + UP * 0.2, [-2.1, y + 0.45, 0], 0.5, C_UP)
        dn = self.beam(a1.get_right() + DOWN * 0.2, [-2.6, y - 0.45, 0], 0.5, C_DN)
        blk = Line([-2.5, y - 0.7, 0], [-2.5, y - 0.25, 0], stroke_width=8, color=GREY_B)
        a2 = self.sg("x", [-1.0, y + 0.45, 0], cX)
        px = self.beam(a2.get_right() + UP * 0.1, [1.0, y + 0.75, 0], 0.25, cX)
        mx = self.beam(a2.get_right() + DOWN * 0.1, [0.4, y + 0.15, 0], 0.25, GREEN_E)
        blk2 = Line([0.5, y - 0.05, 0], [0.5, y + 0.35, 0], stroke_width=8, color=GREY_B)
        a3 = self.sg("z", [1.8, y + 0.75, 0], cZ)
        u3 = self.beam(a3.get_right() + UP * 0.12, [3.7, y + 1.05, 0], 0.125, C_UP)
        d3 = self.beam(a3.get_right() + DOWN * 0.12, [3.7, y + 0.45, 0], 0.125, C_DN)
        n1 = MT(r"\uparrow:125", 22, C_UP).next_to(u3, RIGHT, buff=0.1)
        n2 = MT(r"\downarrow:125", 22, C_DN).next_to(d3, RIGHT, buff=0.1)
        self.play(FadeIn(src), GrowArrow(b0), FadeIn(a1), GrowArrow(up), GrowArrow(dn), FadeIn(blk))
        self.play(FadeIn(a2), GrowArrow(px), GrowArrow(mx), FadeIn(blk2))
        self.play(FadeIn(a3), GrowArrow(u3), GrowArrow(d3))
        self.play(FadeIn(n1), FadeIn(n2), Flash(n2, color=C_DN))
        expl = boxed(panel(r"已被挡掉的 $\downarrow$ 又出现了！",
                           r"$\chi_{+x}=\tfrac1{\sqrt2}(\chi_{\uparrow}+\chi_{\downarrow})$：测完 $S_x$，态坍缩到 $\chi_{+x}$",
                           r"$[\hat S_x,\hat S_z]=-\ii\hbar\hat S_y\neq0$，二者不能同时确定", size=22, color=C_HL))
        expl.to_corner(DR, buff=0.2)
        self.play(FadeIn(expl))
        self.wait(3)


# --------------------------------------------------------------------------
# 7.3 塞曼效应：强磁场中 2p、1s 能级分裂，谱线一分为三
# --------------------------------------------------------------------------
class C7_Zeeman(Scene):
    def construct(self):
        title = title_bar(r"（简单）塞曼效应：$E_{nlm}=E_{nl}+\mu_{\mathrm B}B\,(m+2m_s)$，谱线一分为三")
        self.play(FadeIn(title))
        ax = Axes(x_range=[0, 1.05, 0.5], y_range=[0, 7.0, 1], x_length=5.2, y_length=6.0,
                  axis_config={"include_tip": True, "include_ticks": False, "color": GREY_B}).to_edge(LEFT, buff=0.7).shift(DOWN * 0.35)
        xl = MT(r"B", 28).next_to(ax.x_axis, RIGHT, buff=0.1)
        yl = MT(r"E", 28).next_to(ax.y_axis, UP, buff=0.05)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))
        k = 0.42
        E2p, E1s = 5.3, 1.2
        b = ValueTracker(0.0001)
        lv2p = [2, 1, 0, -1, -2]
        cols = {2: RED_B, 1: ORANGE, 0: YELLOW, -1: GREEN_B, -2: BLUE_B}
        lines = VGroup()
        for q in lv2p:
            lines.add(always_redraw(lambda q=q: Line(ax.c2p(0, E2p), ax.c2p(b.get_value(), E2p + k * q * b.get_value()),
                                                     color=cols[q], stroke_width=3)))
        for q in (1, -1):
            lines.add(always_redraw(lambda q=q: Line(ax.c2p(0, E1s), ax.c2p(b.get_value(), E1s + k * q * b.get_value()),
                                                     color=C_UP if q > 0 else C_DN, stroke_width=3)))
        l2p = MT(r"2\mathrm p", 28).next_to(ax.c2p(0, E2p), LEFT, buff=0.15)
        l1s = MT(r"1\mathrm s", 28).next_to(ax.c2p(0, E1s), LEFT, buff=0.15)
        self.add(lines)
        self.play(FadeIn(l2p), FadeIn(l1s))
        # 右侧光谱
        sp = NumberLine(x_range=[-1.6, 1.6, 0.8], length=5.0, include_ticks=False, color=GREY_B).move_to([3.7, -2.5, 0])
        spl = MT(r"\omega", 26).next_to(sp, RIGHT, buff=0.1)
        w0l = MT(r"\omega_0", 22).next_to(sp.n2p(0), DOWN, buff=0.15)
        hspec = CT(r"光谱", 24).next_to(sp, UP, buff=1.7)
        sc = 1.2

        def specline(d, col):
            return always_redraw(lambda: Line(sp.n2p(d * sc * b.get_value()), sp.n2p(d * sc * b.get_value()) + UP * 1.4,
                                              color=col, stroke_width=5))
        spec = VGroup(specline(1, RED_C), specline(0, GREEN_C), specline(-1, BLUE_C))
        self.play(Create(sp), FadeIn(spl), FadeIn(w0l), FadeIn(hspec), FadeIn(spec))
        info = boxed(panel(r"$\hat H'=\dfrac{eB}{2m_{\mathrm e}}(\hat L_z+2\hat S_z)$",
                           r"$m_s=\pm\tfrac12$：2p 分成 5 条（$m+2m_s=0$ 两重）", size=22))
        info.to_corner(UR, buff=0.25).shift(DOWN * 0.7)
        self.play(FadeIn(info))
        self.play(b.animate.set_value(1.0), run_time=4)
        # 标出 m+2m_s
        qlabs = VGroup(*[MT(rf"{q:+d}" if q else "0", 20, cols[q]).next_to(ax.c2p(1.0, E2p + k * q), RIGHT, buff=0.08) for q in lv2p])
        qlabs.add(MT(r"+1", 20, C_UP).next_to(ax.c2p(1.0, E1s + k), RIGHT, buff=0.08),
                  MT(r"-1", 20, C_DN).next_to(ax.c2p(1.0, E1s - k), RIGHT, buff=0.08))
        qh = MT(r"m+2m_s", 20, GREY_A).next_to(ax.c2p(1.0, E2p + k * 2), UP, buff=0.35)
        self.play(FadeIn(qlabs), FadeIn(qh))
        # 跃迁：Δm_s=0, Δm=0,±1
        trans = [(2, 1, RED_C, 0.30), (1, 1, GREEN_C, 0.42), (0, 1, BLUE_C, 0.54),
                 (0, -1, RED_C, 0.70), (-1, -1, GREEN_C, 0.82), (-2, -1, BLUE_C, 0.94)]
        arrs = VGroup(*[Arrow(ax.c2p(x, E2p + k * q * x), ax.c2p(x, E1s + k * s * x), buff=0, color=c, stroke_width=3,
                              max_tip_length_to_length_ratio=0.05) for q, s, c, x in trans])
        sel = CT(r"选择定则 $\Delta m=0,\pm1$，$\Delta m_s=0$ $\Rightarrow$ 6 个跃迁只有 3 种频率", 22, C_HL)
        fitw(sel, 6.2)
        sel.next_to(info, DOWN, buff=0.3).to_edge(RIGHT, buff=0.3)
        self.play(LaggedStart(*[GrowArrow(a) for a in arrs], lag_ratio=0.25), FadeIn(sel), run_time=3)
        res = MT(r"\omega=\omega_0,\quad \omega_0\pm\frac{eB}{2m_{\mathrm e}}", 30, C_HL).next_to(sel, DOWN, buff=0.35)
        self.play(Write(res))
        self.wait(3)


# --------------------------------------------------------------------------
# 7.4 角动量耦合：矢量模型 + 态数计数
# --------------------------------------------------------------------------
def proj(v):
    """简单的斜投影：3D (x, y, z) -> 屏幕，z 向上，y 朝里。"""
    return np.array([v[0] + 0.35 * v[1], v[2] + 0.18 * v[1], 0.0])


class C7_AngularCoupling(Scene):
    def construct(self):
        title = title_bar(r"两个角动量的耦合：$\hat{\vv J}=\hat{\vv L}+\hat{\vv S}$，$j=l+\tfrac12,\ l-\tfrac12$")
        self.play(FadeIn(title))
        # ---------- 矢量模型（l=1, s=1/2, j=3/2, m=1/2）
        u = 1.45
        O = np.array([-3.6, -2.4, 0])
        l, s, j, mj = 1, 0.5, 1.5, 1.5
        Lm, Sm, Jm = np.sqrt(l * (l + 1)), np.sqrt(s * (s + 1)), np.sqrt(j * (j + 1))
        cLJ = (Jm ** 2 + Lm ** 2 - Sm ** 2) / (2 * Jm * Lm)
        TJ = np.arccos(mj / Jm)
        zax = Arrow(O, O + UP * 4.3, buff=0, color=GREY_B, stroke_width=3)
        zl = MT("z", 26, GREY_B).next_to(zax, UP, buff=0.05)
        self.play(GrowArrow(zax), FadeIn(zl))
        t = ValueTracker(0)

        def vecs():
            a = t.get_value()
            n = np.array([np.sin(TJ) * np.cos(a), np.sin(TJ) * np.sin(a), np.cos(TJ)])
            e1 = np.cross(n, [0, 0, 1.0]); e1 /= np.linalg.norm(e1)
            e2 = np.cross(n, e1)
            ph = 3.2 * a
            Lv = Lm * (cLJ * n + np.sqrt(1 - cLJ ** 2) * (np.cos(ph) * e1 + np.sin(ph) * e2))
            return Jm * n, Lv
        def mk():
            Jv, Lv = vecs()
            J = Arrow(O, O + u * proj(Jv), buff=0, color=YELLOW, stroke_width=5, max_tip_length_to_length_ratio=0.1)
            Lr = Arrow(O, O + u * proj(Lv), buff=0, color=C_HL, stroke_width=4, max_tip_length_to_length_ratio=0.12)
            Sr = Arrow(O + u * proj(Lv), O + u * proj(Jv), buff=0, color=TEAL_C, stroke_width=4, max_tip_length_to_length_ratio=0.25)
            return VGroup(J, Lr, Sr)
        grp = always_redraw(mk)
        # cone ellipse for J around z
        h = Jm * np.cos(TJ); rr = Jm * np.sin(TJ)
        cone = VMobject(color=GREY_C, stroke_width=1.5)
        cone.set_points_smoothly([O + u * proj([rr * np.cos(a), rr * np.sin(a), h]) for a in np.linspace(0, TAU, 60)])
        self.play(FadeIn(grp), Create(cone))
        leg = VGroup(MT(r"\vv J", 28, YELLOW), MT(r"\vv L", 28, C_HL), MT(r"\vv S", 28, TEAL_C)).arrange(RIGHT, buff=0.5)
        leg.next_to(O, DOWN, buff=0.25).shift(RIGHT * 0.6)
        self.play(FadeIn(leg))
        txt = boxed(panel(r"$\hat{\vv L}\cdot\hat{\vv S}=\tfrac12(\hat J^2-\hat L^2-\hat S^2)$ 守恒：",
                          r"$\vv L$、$\vv S$ 绕 $\vv J$ 进动，$\vv J$ 绕 $z$ 进动",
                          r"$J^2,J_z,L^2,S^2$ 守恒；$L_z$、$S_z$ 不守恒",
                          r"$\Rightarrow$ $l,j,m_j$ 是好量子数，$m_l,m_s$ 不是", size=22))
        txt.to_corner(UR, buff=0.25).shift(DOWN * 0.75)
        self.play(FadeIn(txt))
        self.play(t.animate.set_value(TAU), run_time=9, rate_func=linear)
        self.play(FadeOut(grp), FadeOut(cone), FadeOut(leg), FadeOut(zax), FadeOut(zl), FadeOut(txt))
        # ---------- 计数：j1=2, j2=1
        h2 = CT(r"为什么 $j=j_1+j_2,\dots,|j_1-j_2|$？数一数态（$j_1=2,\ j_2=1$）", 26).next_to(title, DOWN, buff=0.25)
        self.play(FadeIn(h2))
        j1, j2 = 2, 1
        g = 0.75
        G0 = np.array([-3.2, -0.7, 0])
        dots = {}
        dg = VGroup()
        for m1 in range(-j1, j1 + 1):
            for m2 in range(-j2, j2 + 1):
                d = Dot(G0 + np.array([m1 * g, m2 * g, 0]), radius=0.11, color=GREY_B)
                dots[(m1, m2)] = d; dg.add(d)
        axm1 = MT(r"m_1", 24).next_to(G0 + RIGHT * (j1 + 0.6) * g, RIGHT, buff=0.05)
        axm2 = MT(r"m_2", 24).next_to(G0 + UP * (j2 + 0.6) * g, UP, buff=0.05)
        m1t = VGroup(*[MT(str(m), 20, GREY_A).next_to(G0 + np.array([m * g, -(j2 + 0.6) * g, 0]), DOWN, buff=0.02) for m in range(-j1, j1 + 1)])
        m2t = VGroup(*[MT(str(m), 20, GREY_A).next_to(G0 + np.array([-(j1 + 0.6) * g, m * g, 0]), LEFT, buff=0.02) for m in range(-j2, j2 + 1)])
        self.play(LaggedStart(*[FadeIn(d) for d in dg], lag_ratio=0.04), FadeIn(axm1), FadeIn(axm2), FadeIn(m1t), FadeIn(m2t))
        diag = VGroup()
        for m in range(-3, 4):
            pts = [G0 + np.array([m1 * g, (m - m1) * g, 0]) for m1 in range(-j1, j1 + 1) if abs(m - m1) <= j2]
            ln = DashedLine(pts[0] + (pts[0] - pts[-1]) * 0.0 + np.array([0.25, -0.25, 0]) * (1 if len(pts) > 0 else 0),
                            pts[-1] + np.array([-0.25, 0.25, 0]), color=GREY_D, stroke_width=1.5) if len(pts) > 1 else VGroup()
            diag.add(ln)
        cnt = VGroup(CT(r"$m=m_1+m_2$：$3,\ 2,\ 1,\ 0,\ -1,\ -2,\ -3$", 24),
                     CT(r"态数：$\ 1,\ 2,\ 3,\ 3,\ 3,\ 2,\ 1$", 24)).arrange(DOWN, aligned_edge=LEFT, buff=0.15).move_to([3.4, 1.0, 0])
        self.play(FadeIn(diag), FadeIn(cnt))
        cols = [YELLOW, C_HL, TEAL_C]
        jtxt = VGroup()
        taken = set()
        for idx, jj in enumerate([3, 2, 1]):
            pick = []
            for m in range(-jj, jj + 1):
                cands = [k for k in sorted(dots, key=lambda k: -k[0]) if k[0] + k[1] == m and k not in taken]
                pick.append(cands[0]); taken.add(cands[0])
            self.play(*[dots[k].animate.set_color(cols[idx]).scale(1.25) for k in pick], run_time=1.0)
            lab = CT(rf"$j={jj}$：$2j+1={2 * jj + 1}$ 个态", 24, cols[idx])
            jtxt.add(lab)
            jtxt.arrange(DOWN, aligned_edge=LEFT, buff=0.18).move_to([3.4, -0.9, 0])
            self.play(FadeIn(lab))
        tot = boxed(panel(r"$7+5+3=15=(2j_1+1)(2j_2+1)$",
                          r"$\sum_{j=j_{\min}}^{j_1+j_2}(2j+1)=(2j_1+1)(2j_2+1)\Rightarrow j_{\min}=|j_1-j_2|$", size=22, color=C_HL))
        tot.to_edge(DOWN, buff=0.25).shift(RIGHT * 2.2)
        self.play(FadeIn(tot))
        self.wait(3)


# --------------------------------------------------------------------------
# 7.5 精细结构：自旋–轨道耦合把 nl 能级分成 j=l±1/2 两条
# --------------------------------------------------------------------------
class C7_FineStructure(Scene):
    def construct(self):
        title = title_bar(r"精细结构：$\xi(r)\,\hat{\vv L}\cdot\hat{\vv S}$ 把 $nl$ 能级（$l>0$）分成 $j=l\pm\tfrac12$ 两条")
        self.play(FadeIn(title))
        x0, x1, x2, x3 = -6.7, -4.9, -4.1, -2.3
        E0 = 0.3
        base = Line([x0, E0, 0], [x1, E0, 0], color=WHITE, stroke_width=4)
        bl = CT(r"$n\mathrm p\ (l=1)$：6 个态", 24).next_to(base, UP, buff=0.12)
        self.play(Create(base), FadeIn(bl))
        a = 0.9
        up = Line([x2, E0 + a, 0], [x3, E0 + a, 0], color=C_HL, stroke_width=4)
        dn = Line([x2, E0 - 2 * a, 0], [x3, E0 - 2 * a, 0], color=TEAL_C, stroke_width=4)
        c1 = DashedLine([x1, E0, 0], [x2, E0 + a, 0], color=GREY_C)
        c2 = DashedLine([x1, E0, 0], [x2, E0 - 2 * a, 0], color=GREY_C)
        ul = CT(r"$j=\tfrac32$（4 个态）：$+\tfrac12\hbar^2\zeta$", 24, C_HL).next_to(up, RIGHT, buff=0.15)
        dl = CT(r"$j=\tfrac12$（2 个态）：$-\hbar^2\zeta$", 24, TEAL_C).next_to(dn, RIGHT, buff=0.15)
        cen = DashedLine([x2, E0, 0], [x3, E0, 0], color=GREY_B, stroke_width=1.5)
        self.play(Create(c1), Create(c2), Create(up), Create(dn), FadeIn(ul), FadeIn(dl))
        self.play(Create(cen))
        info = boxed(panel(r"$E^{(1)}_{nlj}=\dfrac{\hbar^2}{2}\big[j(j+1)-l(l+1)-\tfrac34\big]\,\zeta_{nl}$",
                           r"$\zeta_{nl}=\int_0^\infty R_{nl}^2\,\xi(r)\,r^2\dd r>0$",
                           r"$j=l+\tfrac12$：$+\tfrac{l}{2}\hbar^2\zeta$；$j=l-\tfrac12$：$-\tfrac{l+1}{2}\hbar^2\zeta$",
                           r"重心不变：$(2l+2)\cdot l-2l\cdot(l+1)=0$", size=22))
        info.to_corner(UR, buff=0.25).shift(DOWN * 0.75)
        self.play(FadeIn(info))
        self.wait(1)
        # 钠 D 线
        nah = CT(r"钠原子 D 双线：$3^2\mathrm P_{3/2},\,3^2\mathrm P_{1/2}\to3^2\mathrm S_{1/2}$", 24, YELLOW)
        nah.next_to(info, DOWN, buff=0.35)
        fitw(nah, 6.6)
        s_lev = Line([x2, -2.9, 0], [x3, -2.9, 0], color=WHITE, stroke_width=4)
        sl = MT(r"3^2\mathrm S_{1/2}", 24).next_to(s_lev, RIGHT, buff=0.1)
        ar1 = Arrow([x2 + 0.6, E0 + a, 0], [x2 + 0.7, -2.9, 0], buff=0, color=YELLOW, stroke_width=3)
        ar2 = Arrow([x2 + 1.3, E0 - 2 * a, 0], [x2 + 1.6, -2.9, 0], buff=0, color=GOLD, stroke_width=3)
        w1 = MT(r"589.0\,\mathrm{nm}", 20, YELLOW).next_to(ar1, LEFT, buff=0.05)
        w2 = MT(r"589.6\,\mathrm{nm}", 20, GOLD).next_to(ar2, RIGHT, buff=0.05)
        self.play(FadeIn(nah), Create(s_lev), FadeIn(sl))
        self.play(GrowArrow(ar1), GrowArrow(ar2), FadeIn(w1), FadeIn(w2))
        est = boxed(panel(r"$\Delta E=\dfrac{hc\,\Delta\lambda}{\lambda^2}\approx2.1\times10^{-3}\ \mathrm{eV}$",
                          r"类氢：$\Delta E\sim mc^2\alpha^4Z^4/n^3$（$\alpha\approx\tfrac1{137}$）", size=22, color=YELLOW))
        est.next_to(nah, DOWN, buff=0.3)
        self.play(FadeIn(est))
        self.wait(3)
