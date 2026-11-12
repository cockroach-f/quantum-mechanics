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


# --------------------------------------------------------------------------
# 3.2 / 3.3 球谐函数：|Y_lm|^2 的角分布（绕 z 轴旋转对称）与实球谐（化学中的轨道形状）
# --------------------------------------------------------------------------
def _sph_surface(rfun, scale, res=(28, 56), color=BLUE_C, sign=None):
    """rfun(theta, phi) >= 0 为极坐标半径；sign(x,y,z) 给出每个面片的符号用于着色。"""
    def f(u, v):
        r = rfun(u, v) * scale
        return np.array([r * np.sin(u) * np.cos(v), r * np.sin(u) * np.sin(v), r * np.cos(u)])
    s = Surface(f, u_range=[1e-3, PI - 1e-3], v_range=[0, TAU], resolution=res,
                fill_opacity=0.85, stroke_width=0.3, stroke_color=GREY_D, checkerboard_colors=False)
    s.set_fill(color)
    if sign is not None:
        for face in s:
            c = face.get_center()
            face.set_fill(BLUE_C if sign(*c) >= 0 else ORANGE, opacity=0.85)
    return s


class C3_SphericalHarmonics(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=68 * DEGREES, theta=-50 * DEGREES)
        title = title_bar(r"球谐函数：角向概率 $W_{lm}(\theta)=\lvert Y_{lm}\rvert^2$ 与 $\varphi$ 无关，图形绕 $z$ 轴旋转对称")
        self.add_fixed_in_frame_mobjects(title)
        axes = ThreeDAxes(x_range=[-2.6, 2.6], y_range=[-2.6, 2.6], z_range=[-2.6, 2.6],
                          x_length=5.2, y_length=5.2, z_length=5.2, axis_config={"stroke_width": 1.5, "include_tip": False})
        axes.set_color(GREY_B)
        zl = MT("z", 28, GREY_B).rotate(PI / 2, RIGHT).move_to([0, 0, 2.85])
        self.play(FadeIn(title), Create(axes), FadeIn(zl))
        Y2fun = {k: (lambda f: (lambda u, v: f(np.cos(u))))(f) for k, f in Y2.items()}
        seq = [((0, 0), r"Y_{00}=\frac{1}{\sqrt{4\pi}}", r"s 态：球面，各方向等概率"),
               ((1, 0), r"\lvert Y_{10}\rvert^2=\frac{3}{4\pi}\cos^2\theta", r"沿 $z$ 轴最大，$xy$ 平面为零"),
               ((1, 1), r"\lvert Y_{1,\pm1}\rvert^2=\frac{3}{8\pi}\sin^2\theta", r"“甜甜圈”：$z$ 轴上为零"),
               ((2, 0), r"\lvert Y_{20}\rvert^2\propto(3\cos^2\theta-1)^2", r"两个主瓣 + 赤道环，节锥 $\cos^2\theta=\tfrac13$"),
               ((2, 1), r"\lvert Y_{2,\pm1}\rvert^2\propto\sin^2\theta\cos^2\theta", r"峰在 $\theta=45^\circ,135^\circ$"),
               ((2, 2), r"\lvert Y_{2,\pm2}\rvert^2\propto\sin^4\theta", r"更扁的赤道环")]
        cur = None; lab = None
        self.begin_ambient_camera_rotation(rate=0.18)
        for (l, m), f, d in seq:
            mx = max(Y2[(l, m)](np.linspace(-1, 1, 801)))
            surf = _sph_surface(Y2fun[(l, m)], 2.3 / mx, color=[BLUE_C, GREEN_C, PURPLE_B][l])
            mstr = r"\pm" + str(m) if m else "0"
            lzs = (r"\pm\hbar" if m == 1 else rf"\pm{m}\hbar") if m else "0"
            qn = CT(rf"$l={l}$，$m={mstr}$：$L^2={l*(l+1)}\hbar^2$，$L_z={lzs}$" if l else r"$l=0$，$m=0$：$L^2=0$，$L_z=0$", 24, YELLOW)
            new_lab = VGroup(MT(f, 32), CT(d, 24, GREY_A), qn).arrange(DOWN, aligned_edge=LEFT, buff=0.18).to_corner(DL, buff=0.35)
            self.add_fixed_in_frame_mobjects(new_lab)
            if cur is None:
                self.play(Create(surf), FadeIn(new_lab), run_time=1.6)
            else:
                self.remove(lab)
                self.play(FadeOut(cur), FadeIn(surf), FadeIn(new_lab), run_time=1.3)
            cur, lab = surf, new_lab
            self.wait(1.6)
        # 实球谐：±m 叠加成实函数，颜色表示正负号
        note = CT(r"把 $Y_{l,m}$ 与 $Y_{l,-m}$ 叠加成实函数，就是化学中的 $p_x$、$d_{xy}$ 轨道（蓝 $+$，橙 $-$）", 24, C_HL).to_edge(DOWN, buff=0.2)
        fitw(note, 13)
        self.remove(lab)
        self.add_fixed_in_frame_mobjects(note)
        real = [(r"p_x\propto\frac{Y_{1,-1}-Y_{11}}{\sqrt2}\propto\sin\theta\cos\varphi",
                 lambda u, v: np.abs(np.sin(u) * np.cos(v)), lambda x, y, z: x),
                (r"d_{z^2}\propto Y_{20}\propto 3\cos^2\theta-1",
                 lambda u, v: np.abs(3 * np.cos(u) ** 2 - 1) / 2, lambda x, y, z: 2 * z * z - x * x - y * y),
                (r"d_{xy}\propto\frac{Y_{2,2}-Y_{2,-2}}{\sqrt2\,\ii}\propto\sin^2\theta\sin2\varphi",
                 lambda u, v: np.abs(np.sin(u) ** 2 * np.sin(2 * v)), lambda x, y, z: x * y)]
        for i, (f, rf, sg) in enumerate(real):
            surf = _sph_surface(rf, 2.3, sign=sg)
            lab2 = MT(f, 30).to_corner(DL, buff=0.35).shift(UP * 0.6)
            self.add_fixed_in_frame_mobjects(lab2)
            anims = [FadeOut(cur), FadeIn(surf), FadeIn(lab2)]
            if i == 0:
                anims.append(FadeIn(note))
            else:
                self.remove(lab)
            self.play(*anims, run_time=1.3)
            cur, lab = surf, lab2
            self.wait(2.0)
        self.stop_ambient_camera_rotation()
        self.wait(0.5)


# --------------------------------------------------------------------------
# 3.4 / 3.5 本征函数展开与测量：投影 → 概率 → 坍缩 → 统计
# --------------------------------------------------------------------------
class C3_Measurement(Scene):
    def construct(self):
        title = title_bar(r"按本征函数展开：$\psi=\sum_n c_n\phi_n$，测得 $E_n$ 的概率为 $\lvert c_n\rvert^2$，测量后态坍缩为 $\phi_n$")
        self.play(FadeIn(title))
        xs = np.linspace(0, 1, 300)
        phi = lambda n: np.sqrt(2) * np.sin(n * np.pi * xs)
        # ---- Part 1: ψ = √30 x(1-x) 的正弦级数部分和
        ax = Axes(x_range=[0, 1, 0.5], y_range=[-0.2, 2.0, 0.5], x_length=6.2, y_length=3.6,
                  axis_config={"include_tip": False, "font_size": 22}).shift(LEFT * 3.1 + DOWN * 0.5)
        lab = CT(r"例：势阱 $0<x<a$ 中 $\psi=\sqrt{30/a^5}\,x(a-x)$", 24).next_to(ax, UP, buff=0.15)
        psi0 = np.sqrt(30) * xs * (1 - xs)
        self.play(Create(ax), FadeIn(lab))
        target = curve(ax, xs, psi0, C_PROB, 5)
        self.play(Create(target))
        cn = lambda n: 0 if n % 2 == 0 else 8 * np.sqrt(15) / (n * np.pi) ** 3
        partial = np.zeros_like(xs)
        info = VGroup()
        prev = None
        for n in [1, 3, 5]:
            partial = partial + cn(n) * phi(n)
            c = curve(ax, xs, partial, C_PSI, 3)
            t = CT(rf"$n\le{n}$：$c_{n}={cn(n):.4f}$，$\lvert c_{n}\rvert^2={cn(n)**2:.4f}$", 24, C_PSI)
            info.add(t)
            info.arrange(DOWN, aligned_edge=LEFT, buff=0.22).to_edge(RIGHT, buff=0.5).shift(UP * 0.8)
            self.play(Create(c) if prev is None else Transform(prev, c), FadeIn(t), run_time=1.2)
            if prev is None:
                prev = c
            self.wait(0.5)
        k = boxed(panel(r"$c_n=\ip{\phi_n}{\psi}=\int\phi_n^*\psi\,\dd x$（投影）",
                        r"偶数 $n$ 的 $c_n=0$：宇称不同，重叠为零",
                        r"$\lvert c_1\rvert^2\approx0.9986$：抛物线几乎就是 $\phi_1$", size=24))
        k.next_to(info, DOWN, buff=0.4)
        fitw(k, 6.2)
        self.play(FadeIn(k))
        self.wait(2)
        self.play(*[FadeOut(m) for m in [ax, lab, target, prev, info, k]])
        # ---- Part 2: 三能级叠加态，测量与坍缩
        cs = [np.sqrt(0.5), 0.5, 0.5]
        ax2 = Axes(x_range=[0, 1, 0.5], y_range=[-2.2, 2.6, 1], x_length=5.6, y_length=3.8,
                   axis_config={"include_tip": False, "include_ticks": False}).shift(LEFT * 3.4 + DOWN * 0.4)
        st = MT(r"\psi=\tfrac{1}{\sqrt2}\phi_1+\tfrac12\phi_2+\tfrac12\phi_3", 30).next_to(ax2, UP, buff=0.15)
        wav = lambda cc: sum(c * phi(i + 1) for i, c in enumerate(cc))
        wc = curve(ax2, xs, wav(cs), C_PSI, 4)
        # bars
        base = RIGHT * 1.6 + DOWN * 2.5
        bw, H = 0.9, 3.1
        bars = VGroup(); blabs = VGroup()
        for i, c in enumerate(cs):
            b = Rectangle(width=bw, height=max(H * c * c, 1e-3), fill_color=[YELLOW, GREEN_C, RED_C][i], fill_opacity=0.8, stroke_width=0)
            b.move_to(base + RIGHT * 1.4 * i, aligned_edge=DOWN)
            bars.add(b)
            blabs.add(MT(rf"E_{i+1}", 28).next_to(base + RIGHT * 1.4 * i, DOWN, buff=0.15))
        ptxt = VGroup(*[MT(f"{c*c:.2f}", 26).next_to(bars[i], UP, buff=0.08) for i, c in enumerate(cs)])
        hd = CT(r"测量能量的概率 $\lvert c_n\rvert^2$", 26).move_to(base + RIGHT * 1.4 + UP * (H + 0.75))
        self.play(Create(ax2), FadeIn(st), Create(wc))
        self.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in bars], lag_ratio=0.3), FadeIn(blabs), FadeIn(ptxt), FadeIn(hd))
        self.wait(1)
        flash = CT(r"测量！结果：$E_2$", 32, C_HL).next_to(hd, UP, buff=0.2)
        self.play(FadeIn(flash, scale=1.3), Indicate(bars[1], color=WHITE))
        new = [0, 1, 0]
        wc2 = curve(ax2, xs, wav(new), GREEN_C, 4)
        st2 = CT(r"测量后 $\psi\to\phi_2$（坍缩）", 28, GREEN_C).move_to(st)
        newbars = VGroup(*[Rectangle(width=bw, height=max(H * c, 1e-3), fill_color=bars[i].get_fill_color(), fill_opacity=0.8, stroke_width=0)
                           .move_to(base + RIGHT * 1.4 * i, aligned_edge=DOWN) for i, c in enumerate(new)])
        newp = VGroup(*[MT(f"{c:.2f}", 26).next_to(newbars[i], UP, buff=0.08) for i, c in enumerate(new)])
        self.play(Transform(wc, wc2), Transform(st, st2), Transform(bars, newbars), Transform(ptxt, newp), run_time=1.6)
        again = CT(r"立即再测一次：必然仍得 $E_2$", 26, GREEN_C).next_to(ax2, DOWN, buff=0.25)
        self.play(FadeIn(again))
        self.wait(1.5)
        # ---- Part 3: 大量全同体系的统计
        self.play(*[FadeOut(m) for m in [wc, st, again, flash, ptxt]])
        st3 = CT(r"对 $N$ 个处于同一 $\psi$ 的体系各测一次：", 26).next_to(ax2, UP, buff=0.15).shift(UP * 1.0)
        self.play(FadeOut(ax2), FadeIn(st3))
        rng = np.random.default_rng(7)
        outcomes = rng.choice(3, size=600, p=[c * c for c in cs])
        Nt = ValueTracker(1)

        def mkbars():
            n = int(Nt.get_value())
            cnt = np.bincount(outcomes[:n], minlength=3) / n
            g = VGroup()
            for i in range(3):
                g.add(Rectangle(width=bw, height=max(H * cnt[i], 1e-3), fill_color=[YELLOW, GREEN_C, RED_C][i],
                                fill_opacity=0.8, stroke_width=0).move_to(base + RIGHT * 1.4 * i, aligned_edge=DOWN))
            return g
        hist = always_redraw(mkbars)
        ideal = VGroup(*[DashedLine(base + RIGHT * 1.4 * i + LEFT * 0.55 + UP * H * c * c, base + RIGHT * 1.4 * i + RIGHT * 0.55 + UP * H * c * c,
                                    color=WHITE, stroke_width=3) for i, c in enumerate(cs)])
        nlab = VGroup(MT(r"N=", 30), Integer(1, font_size=30)).arrange(RIGHT, buff=0.1).to_edge(LEFT, buff=1.0).shift(UP * 0.6)
        nlab[1].add_updater(lambda m: m.set_value(int(Nt.get_value())))
        self.remove(bars)
        self.add(hist, nlab)
        self.play(FadeIn(ideal))
        self.play(Nt.animate.set_value(600), run_time=5, rate_func=rate_functions.ease_in_quad)
        fin = boxed(panel(r"频率 $\to\lvert c_n\rvert^2$（虚线）",
                          r"$\langle E\rangle=\sum_n E_n\lvert c_n\rvert^2=\int\psi^*\hat H\psi\,\dd x$",
                          r"单次结果随机；统计规律确定", size=26))
        fin.next_to(nlab, DOWN, buff=0.4).to_edge(LEFT, buff=0.5)
        self.play(FadeIn(fin))
        self.wait(2.5)


# --------------------------------------------------------------------------
# 3.6 对易关系：绕 x、y 轴的转动不可交换  ->  [L_x, L_y] = iħ L_z
# --------------------------------------------------------------------------
class C3_Commutator(ThreeDScene):
    def make_book(self, center):
        b = Prism(dimensions=[1.6, 1.1, 0.35]).move_to(center)
        cols = [BLUE_C, BLUE_C, RED_C, RED_C, YELLOW, YELLOW]
        for face, c in zip(b, cols):
            face.set_fill(c, opacity=0.95).set_stroke(WHITE, 1)
        mark = Text("QM", font_size=28, color=BLACK).move_to(center + OUT * 0.18)
        return VGroup(b, mark)

    def construct(self):
        self.set_camera_orientation(phi=70 * DEGREES, theta=-90 * DEGREES)
        title = title_bar(r"对易子度量“先后次序”的影响：先绕 $x$ 再绕 $y$ 转 $\neq$ 先绕 $y$ 再绕 $x$ 转")
        self.add_fixed_in_frame_mobjects(title)
        cL, cR = np.array([-3.0, 0, 0.9]), np.array([3.0, 0, 0.9])
        A = self.make_book(cL); B = self.make_book(cR)
        axL = ThreeDAxes(x_range=[-1.3, 1.3], y_range=[-1.3, 1.3], z_range=[-1.3, 1.3], x_length=2.6, y_length=2.6, z_length=2.6,
                         axis_config={"include_tip": True, "stroke_width": 1.5}).move_to(cL).set_color(GREY_B)
        axR = axL.copy().move_to(cR)
        alabs = VGroup()
        for c in (cL, cR):
            for t, d in (("x", RIGHT), ("y", UP), ("z", OUT)):
                alabs.add(MT(t, 26, GREY_B).rotate(PI / 2, RIGHT).move_to(c + d * 1.55))
        tl = CT(r"先 $R_x(90^\circ)$，再 $R_y(90^\circ)$", 26, YELLOW).to_corner(UL, buff=0.4).shift(DOWN * 0.8)
        tr = CT(r"先 $R_y(90^\circ)$，再 $R_x(90^\circ)$", 26, YELLOW).to_corner(UR, buff=0.4).shift(DOWN * 0.8)
        self.add_fixed_in_frame_mobjects(tl, tr)
        self.play(FadeIn(title), Create(axL), Create(axR), FadeIn(alabs), FadeIn(A), FadeIn(B), FadeIn(tl), FadeIn(tr))
        self.wait(0.5)
        self.play(Rotate(A, PI / 2, axis=RIGHT, about_point=cL), Rotate(B, PI / 2, axis=UP, about_point=cR), run_time=1.6)
        self.wait(0.3)
        self.play(Rotate(A, PI / 2, axis=UP, about_point=cL), Rotate(B, PI / 2, axis=RIGHT, about_point=cR), run_time=1.6)
        res = CT(r"结果不同：转动不可交换", 30, C_HL).to_edge(DOWN, buff=0.6)
        self.add_fixed_in_frame_mobjects(res)
        self.play(FadeIn(res))
        self.wait(1)
        self.move_camera(phi=60 * DEGREES, theta=-80 * DEGREES, run_time=2)
        self.wait(1)
        k = boxed(panel(r"无穷小转动 $R_\alpha(\varepsilon)\approx1-\dfrac{\ii\varepsilon}{\hbar}\hat L_\alpha$（$\hat{\vv L}$ 是转动的生成元）",
                        r"$R_x(\varepsilon)R_y(\varepsilon)-R_y(\varepsilon)R_x(\varepsilon)=-\dfrac{\varepsilon^2}{\hbar^2}[\hat L_x,\hat L_y]\approx R_z(\varepsilon^2)-1$",
                        r"$\Rightarrow\ [\hat L_x,\hat L_y]=\ii\hbar\hat L_z$：两个分量不能同时有确定值", size=26))
        fitw(k, 12.6)
        k.to_edge(DOWN, buff=0.2)
        self.add_fixed_in_frame_mobjects(k)
        self.play(FadeOut(res), FadeIn(k))
        self.wait(3)
