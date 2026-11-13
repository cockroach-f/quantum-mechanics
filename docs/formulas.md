# 公式速查

按章整理全书最常用的公式。记号约定见[绪论](index.md)。每个公式的推导和适用条件见对应章节。第 6 章整章、以及各章标了「延伸」的公式，课件不要求。

## 常量与常用组合

| 量 | 数值 |
|---|---|
| $hc$ | $1240\ \mathrm{eV\cdot nm}$ |
| $\hbar c$ | $197.3\ \mathrm{eV\cdot nm}$ |
| $m_ec^2$ | $0.511\ \mathrm{MeV}$ |
| $e^2/4\pi\varepsilon_0$ | $1.44\ \mathrm{eV\cdot nm}$ |
| $\alpha=e^2/4\pi\varepsilon_0\hbar c$ | $1/137.04$ |
| $a_0=4\pi\varepsilon_0\hbar^2/m_ee^2$ | $0.0529\ \mathrm{nm}$ |
| 里德伯能量 $e^2/8\pi\varepsilon_0a_0$ | $13.6\ \mathrm{eV}$ |
| $\mu_B=e\hbar/2m_e$ | $9.27\times10^{-24}\ \mathrm{J/T}=5.79\times10^{-5}\ \mathrm{eV/T}$ |
| 电子康普顿波长 $h/m_ec$ | $2.43\times10^{-3}\ \mathrm{nm}$ |

## 第 1 章 绪论

| 内容 | 公式 |
|---|---|
| 普朗克公式 | $\rho_\nu\,\dd\nu=\dfrac{8\pi h\nu^3}{c^3}\dfrac{\dd\nu}{\ee^{h\nu/k_BT}-1}$ |
| 维恩位移定律 | $\lambda_mT=2.90\times10^{-3}\ \mathrm{m\cdot K}$ |
| 光电效应 | $\tfrac12m_ev_m^2=h\nu-W_0$ |
| 康普顿散射 | $\Delta\lambda=\dfrac{h}{m_0c}(1-\cos\theta)$ |
| 玻尔频率条件 | $h\nu=E_n-E_m$ |
| 玻尔–索末菲量子化 | $\oint p\,\dd q=nh$ |
| 德布罗意关系 | $E=\hbar\omega$，$\vv p=\hbar\vv k$，$\lambda=h/p$ |
| 电子德布罗意波长 | $\lambda\approx1.226/\sqrt{V}\ \mathrm{nm}$（$V$ 以伏特计） |

## 第 2 章 波函数和薛定谔方程

| 内容 | 公式 |
|---|---|
| 归一化 | $\int\lvert\Psi\rvert^2\dd\tau=1$ |
| 动量波函数 | $c(p)=\dfrac{1}{\sqrt{2\pi\hbar}}\int\psi(x)\ee^{-\ii px/\hbar}\dd x$，$\psi(x)=\dfrac{1}{\sqrt{2\pi\hbar}}\int c(p)\ee^{\ii px/\hbar}\dd p$ |
| 薛定谔方程 | $\ii\hbar\dfrac{\partial\Psi}{\partial t}=\left[-\dfrac{\hbar^2}{2m}\nabla^2+U\right]\Psi$ |
| 概率流密度 | $\vv J=\dfrac{\ii\hbar}{2m}(\Psi\nabla\Psi^*-\Psi^*\nabla\Psi)$，$\dfrac{\partial w}{\partial t}+\nabla\cdot\vv J=0$ |
| 定态 | $\Psi=\psi\,\ee^{-\ii Et/\hbar}$，$\hat H\psi=E\psi$ |
| 一般解 | $\Psi(x,t)=\sum_nc_n\psi_n(x)\ee^{-\ii E_nt/\hbar}$，$c_n=\int\psi_n^*\Psi(x,0)\dd x$ |
| 无限深方势阱 $0<x<a$ | $E_n=\dfrac{n^2\pi^2\hbar^2}{2ma^2}$，$\psi_n=\sqrt{\dfrac2a}\sin\dfrac{n\pi x}{a}$，$n=1,2,\dots$ |
| 无限深方势阱 $\lvert x\rvert<a$ | $E_n=\dfrac{n^2\pi^2\hbar^2}{8ma^2}$，$\psi_n=\dfrac{1}{\sqrt a}\sin\dfrac{n\pi(x+a)}{2a}$ |
| 线性谐振子 | $E_n=\left(n+\tfrac12\right)\hbar\omega$，$\psi_n=N_n\ee^{-\alpha^2x^2/2}H_n(\alpha x)$，$\alpha=\sqrt{m\omega/\hbar}$，$N_n=\sqrt{\dfrac{\alpha}{\sqrt\pi2^nn!}}$ |
| 厄米多项式 | $H_0=1$，$H_1=2\xi$，$H_2=4\xi^2-2$，$H_3=8\xi^3-12\xi$；$H_{n+1}=2\xi H_n-2nH_{n-1}$ |
| 谐振子递推 | $x\psi_n=\dfrac1\alpha\left[\sqrt{\tfrac n2}\psi_{n-1}+\sqrt{\tfrac{n+1}2}\psi_{n+1}\right]$ |
| 方势垒透射（$E<U_0$） | $D=\dfrac{4k_1^2k_3^2}{(k_1^2+k_3^2)^2\sinh^2k_3a+4k_1^2k_3^2}\approx D_0\ee^{-\frac2\hbar\sqrt{2m(U_0-E)}a}$ |
| 任意势垒 | $D\approx D_0\exp\left[-\dfrac2\hbar\displaystyle\int_a^b\sqrt{2m(U-E)}\,\dd x\right]$ |
| $\delta$ 势衔接条件 | $\psi'(a^+)-\psi'(a^-)=\dfrac{2mA}{\hbar^2}\psi(a)$（$U=A\delta(x-a)$） |
| $\delta$ 势垒透射 | $D=\left(1+\dfrac{mA^2}{2\hbar^2E}\right)^{-1}$ |
| $\delta$ 势阱束缚态 | $E=-\dfrac{m\gamma^2}{2\hbar^2}$，$\psi=\sqrt\kappa\,\ee^{-\kappa\lvert x\rvert}$，$\kappa=m\gamma/\hbar^2$ |

## 第 3 章 量子力学中的力学量

| 内容 | 公式 |
|---|---|
| 期望值 | $\langle F\rangle=\int\psi^*\hat F\psi\,\dd\tau=\sum_n\lambda_n\lvert c_n\rvert^2$ |
| 厄米算符 | $\int\psi^*\hat F\phi\,\dd\tau=\int(\hat F\psi)^*\phi\,\dd\tau$ |
| 动量本征函数 | $\psi_{\vv p}=(2\pi\hbar)^{-3/2}\ee^{\ii\vv p\cdot\vv r/\hbar}$，$\ip{\psi_{\vv p'}}{\psi_{\vv p}}=\delta(\vv p-\vv p')$；箱归一化 $p_x=2\pi\hbar n_x/L$ |
| 角动量算符 | $\hat L_z=-\ii\hbar\dfrac{\partial}{\partial\varphi}$，$\hat L^2=-\hbar^2\left[\dfrac{1}{\sin\theta}\dfrac{\partial}{\partial\theta}\left(\sin\theta\dfrac{\partial}{\partial\theta}\right)+\dfrac{1}{\sin^2\theta}\dfrac{\partial^2}{\partial\varphi^2}\right]$ |
| 球谐函数 | $\hat L^2Y_{lm}=l(l+1)\hbar^2Y_{lm}$，$\hat L_zY_{lm}=m\hbar Y_{lm}$，$\lvert m\rvert\leqslant l$ |
| 低阶球谐函数 | $Y_{00}=\dfrac1{\sqrt{4\pi}}$，$Y_{10}=\sqrt{\dfrac3{4\pi}}\cos\theta$，$Y_{1,\pm1}=\mp\sqrt{\dfrac{3}{8\pi}}\sin\theta\,\ee^{\pm\ii\varphi}$ |
| 径向方程 | $-\dfrac{\hbar^2}{2m}u''+\left[U+\dfrac{l(l+1)\hbar^2}{2mr^2}\right]u=Eu$，$u=rR$，$u(0)=0$ |
| 氢原子能级 | $E_n=-\dfrac{m}{2\hbar^2}\left(\dfrac{Ze^2}{4\pi\varepsilon_0}\right)^2\dfrac1{n^2}$，$n=n_r+l+1$，简并度 $n^2$ |
| 氢原子基态 | $\psi_{100}=\dfrac{1}{\sqrt{\pi a_0^3}}\ee^{-r/a_0}$，$r_{\text{最概然}}=a_0$，$\langle r\rangle=\tfrac32a_0$ |
| 低阶径向函数 | $R_{10}=2a_0^{-3/2}\ee^{-r/a_0}$，$R_{20}=\dfrac{1}{\sqrt2a_0^{3/2}}\left(1-\dfrac{r}{2a_0}\right)\ee^{-r/2a_0}$，$R_{21}=\dfrac{1}{2\sqrt6a_0^{3/2}}\dfrac r{a_0}\ee^{-r/2a_0}$ |
| 光谱与轨道磁矩 | $\nu=R_\infty c\left(\dfrac1{n'^2}-\dfrac1{n^2}\right)$，$R_\infty=-E_1/hc$；$M_z=-m\mu_{\mathrm B}$，$M_z/L_z=-\dfrac{e}{2m_e}$ |
| 径向期望值 | $\langle r\rangle_{nl}=\dfrac{a_0}{2}\left[3n^2-l(l+1)\right]$，$\left\langle\dfrac1r\right\rangle_{nl}=\dfrac{1}{n^2a_0}$ |
| 正交归一、完全 | $\ip{\phi_m}{\phi_n}=\delta_{mn}$，$\sum_n\phi_n(x)\phi_n^*(x')=\delta(x-x')$ |
| 测量与谱分解 | $P(\lambda_n)=\lvert\ip{\phi_n}{\psi}\rvert^2$，$\hat F=\sum_n\lambda_n\ket{\phi_n}\bra{\phi_n}$ |
| 宇称 | $\hat P\psi(\vv r)=\psi(-\vv r)$，本征值 $\pm1$；$Y_{lm}$ 的宇称为 $(-1)^l$ |
| 基本对易关系 | $[\hat x_\alpha,\hat p_\beta]=\ii\hbar\delta_{\alpha\beta}$，$\hat{\vv L}\times\hat{\vv L}=\ii\hbar\hat{\vv L}$，$[\hat L^2,\hat L_\alpha]=0$ |
| 对易子恒等式 | $[\hat A,\hat B\hat C]=[\hat A,\hat B]\hat C+\hat B[\hat A,\hat C]$；$[\hat x,f(\hat p)]=\ii\hbar f'(\hat p)$，$[\hat p,U(x)]=-\ii\hbar U'(x)$ |
| 不确定关系 | $\Delta F\Delta G\geqslant\tfrac12\lvert\langle[\hat F,\hat G]\rangle\rvert$，$\Delta x\Delta p_x\geqslant\hbar/2$ |
| 期望值演化 | $\dfrac{\dd\langle F\rangle}{\dd t}=\left\langle\dfrac{\partial\hat F}{\partial t}\right\rangle+\dfrac{1}{\ii\hbar}\langle[\hat F,\hat H]\rangle$ |
| 埃伦费斯特定理 | $\dfrac{\dd\langle x\rangle}{\dd t}=\dfrac{\langle p\rangle}{m}$，$\dfrac{\dd\langle p\rangle}{\dd t}=-\left\langle\dfrac{\partial U}{\partial x}\right\rangle$ |

## 第 4 章 态和力学量的表象

| 内容 | 公式 |
|---|---|
| 态在 $Q$ 表象 | $\ket\psi\to(a_1,a_2,\dots)^{\mathrm T}$，$a_n=\ip{n}{\psi}$ |
| 算符矩阵元 | $F_{mn}=\mel{m}{\hat F}{n}=\int u_m^*\hat Fu_n\,\dd\tau$，厄米：$F=F^\dagger$ |
| 期望值 | $\langle F\rangle=a^\dagger Fa=\sum_{mn}a_m^*F_{mn}a_n$ |
| 本征值问题 | $\det(F-\lambda I)=0$（久期方程），$(F-\lambda I)a=0$ |
| 薛定谔方程 | $\ii\hbar\dfrac{\dd a}{\dd t}=Ha$ |
| 幺正变换 | $S_{n\alpha}=\ip{n}{\alpha}$，$S^\dagger S=SS^\dagger=I$；$a'=S^\dagger a$，$F'=S^\dagger FS$ |
| 不变量 | 本征值、迹 $\Tr F$、行列式、$\langle F\rangle$、$\ip\phi\psi$ |
| 狄拉克符号 | $\sum_n\ket n\bra n=\hat I$，$\int\ket x\bra x\dd x=\hat I$，$\ip{x}{\psi}=\psi(x)$ |
| 坐标与动量表象 | $\ip{x}{p}=\dfrac{1}{\sqrt{2\pi\hbar}}\ee^{\ii px/\hbar}$，$\mel{x}{\hat p}{\psi}=-\ii\hbar\dfrac{\partial\psi}{\partial x}$，$\mel{p}{\hat x}{\psi}=\ii\hbar\dfrac{\partial c}{\partial p}$ |
| 升降算符 | $\hat a=\sqrt{\dfrac{m\omega}{2\hbar}}\left(\hat x+\dfrac{\ii\hat p}{m\omega}\right)$，$\hat a^\dagger=\sqrt{\dfrac{m\omega}{2\hbar}}\left(\hat x-\dfrac{\ii\hat p}{m\omega}\right)$，$[\hat a,\hat a^\dagger]=1$ |
| 占有数表象 | $\hat H=\hbar\omega\left(\hat N+\tfrac12\right)$，$\hat N=\hat a^\dagger\hat a$，$\hat N\ket n=n\ket n$ |
| 升降作用 | $\hat a\ket n=\sqrt n\ket{n-1}$，$\hat a^\dagger\ket n=\sqrt{n+1}\ket{n+1}$，$\ket n=\dfrac{(\hat a^\dagger)^n}{\sqrt{n!}}\ket0$ |
| $x,p$ 用升降算符 | $\hat x=\sqrt{\dfrac{\hbar}{2m\omega}}(\hat a+\hat a^\dagger)$，$\hat p=\ii\sqrt{\dfrac{m\omega\hbar}{2}}(\hat a^\dagger-\hat a)$ |
| 谐振子 $\ket n$ 中 | $\langle x^2\rangle=\dfrac{\hbar}{m\omega}\left(n+\tfrac12\right)$，$\langle p^2\rangle=m\omega\hbar\left(n+\tfrac12\right)$，$\Delta x\Delta p=\left(n+\tfrac12\right)\hbar$ |

## 第 5 章 微扰理论

| 内容 | 公式 |
|---|---|
| 非简并一级修正 | $E_n^{(1)}=H'_{nn}$，$\psi_n^{(1)}=\sum_{m\neq n}\dfrac{H'_{mn}}{E_n^{(0)}-E_m^{(0)}}\psi_m^{(0)}$ |
| 非简并二级能量 | $E_n^{(2)}=\sum_{m\neq n}\dfrac{\lvert H'_{mn}\rvert^2}{E_n^{(0)}-E_m^{(0)}}$（基态二级修正恒为负） |
| 适用条件 | $\left\lvert\dfrac{H'_{mn}}{E_n^{(0)}-E_m^{(0)}}\right\rvert\ll1$ |
| 简并微扰 | $\det\left[H'_{\mu\nu}-E^{(1)}\delta_{\mu\nu}\right]=0$（$\mu,\nu$ 遍历简并子空间） |
| 氢原子一级斯塔克效应 | $n=2$：$E^{(1)}=\pm3e\mathcal Ea_0,\ 0,\ 0$ |
| 变分法 | $E_0\leqslant\bar H(\lambda)=\dfrac{\int\psi^*\hat H\psi}{\int\psi^*\psi}$，$\dfrac{\partial\bar H}{\partial\lambda}=0$ |
| 氦原子基态（变分） | $Z_{\mathrm{eff}}=Z-\tfrac{5}{16}=\tfrac{27}{16}$，$E_0\approx-2\left(\tfrac{27}{16}\right)^2\times13.6\ \mathrm{eV}\approx-77.5\ \mathrm{eV}$ |
| 一级跃迁振幅 | $a_m^{(1)}(t)=\dfrac{1}{\ii\hbar}\displaystyle\int_0^tH'_{mk}(t')\ee^{\ii\omega_{mk}t'}\dd t'$，$\omega_{mk}=(E_m-E_k)/\hbar$ |
| 常微扰跃迁概率 | $W_{k\to m}=\dfrac{4\lvert H'_{mk}\rvert^2\sin^2(\omega_{mk}t/2)}{\hbar^2\omega_{mk}^2}$ |
| 费米黄金定则 | $w_{k\to m}=\dfrac{2\pi}{\hbar}\lvert H'_{mk}\rvert^2\rho(E_m)$ |
| 周期微扰 | $\hat H'=\hat F\ee^{-\ii\omega t}+\hat F^\dagger\ee^{\ii\omega t}$：$w=\dfrac{2\pi}{\hbar}\lvert F_{mk}\rvert^2\delta(E_m-E_k\mp\hbar\omega)$ |
| 能量–时间不确定关系 | $\Delta E\,\Delta t\sim\hbar$ |
| 爱因斯坦系数 | $\dfrac{A_{mk}}{B_{mk}}=\dfrac{\hbar\omega_{mk}^3}{\pi^2c^3}$，$B_{km}=B_{mk}$ |
| 电偶极跃迁 | $B_{mk}=\dfrac{\pi e^2}{3\varepsilon_0\hbar^2}\lvert\vv r_{mk}\rvert^2$，$A_{mk}=\dfrac{e^2\omega_{mk}^3}{3\pi\varepsilon_0\hbar c^3}\lvert\vv r_{mk}\rvert^2$ |
| 选择定则 | $\Delta l=\pm1$，$\Delta m=0,\pm1$（电偶极、单电子） |

## 第 6 章 散射

| 内容 | 公式 |
|---|---|
| 微分截面 | $\dd n=Nq(\theta,\varphi)\dd\Omega$，$\sigma_t=\int q\,\dd\Omega$ |
| 渐近波函数 | $\psi\xrightarrow{r\to\infty}A\left[\ee^{\ii kz}+f(\theta,\varphi)\dfrac{\ee^{\ii kr}}{r}\right]$，$q=\lvert f\rvert^2$ |
| 平面波的分波展开 | $\ee^{\ii kr\cos\theta}=\sum_l(2l+1)\ii^lj_l(kr)P_l(\cos\theta)$ |
| 分波法 | $f(\theta)=\dfrac1k\sum_l(2l+1)\ee^{\ii\delta_l}\sin\delta_l\,P_l(\cos\theta)$ |
| 总截面 | $\sigma_t=\dfrac{4\pi}{k^2}\sum_l(2l+1)\sin^2\delta_l$；光学定理 $\sigma_t=\dfrac{4\pi}{k}\Im f(0)$ |
| 需要的分波数 | $l\lesssim ka$（$a$ 为力程） |
| 低能硬球 | $\delta_0=-ka$，$\sigma_t\to4\pi a^2$ |
| 玻恩近似 | $f(\theta)=-\dfrac{2m}{\hbar^2q}\displaystyle\int_0^\infty rU(r)\sin(qr)\,\dd r$，$q=2k\sin\dfrac\theta2$ |
| 玻恩近似（一般） | $f=-\dfrac{m}{2\pi\hbar^2}\displaystyle\int U(\vv r)\ee^{-\ii\vv q\cdot\vv r}\dd\tau$，$\vv q=\vv k-\vv k_0$ |
| 卢瑟福公式 | $q(\theta)=\left(\dfrac{Z_1Z_2e^2}{4\pi\varepsilon_0}\right)^2\dfrac{1}{16E^2\sin^4(\theta/2)}$ |
| 质心系 ↔ 实验室系 | $\tan\theta_L=\dfrac{\sin\theta}{\cos\theta+m_1/m_2}$；$q_L(\theta_L)\sin\theta_L\dd\theta_L=q(\theta)\sin\theta\dd\theta$ |

## 第 7 章 自旋

| 内容 | 公式 |
|---|---|
| 自旋算符 | $\hat{\vv S}=\dfrac\hbar2\hat{\vv\sigma}$，$\hat{\vv S}\times\hat{\vv S}=\ii\hbar\hat{\vv S}$，$S^2=\tfrac34\hbar^2$，$S_z=\pm\tfrac\hbar2$ |
| 泡利矩阵 | $\sigma_x=\begin{pmatrix}0&1\\1&0\end{pmatrix}$，$\sigma_y=\begin{pmatrix}0&-\ii\\\ii&0\end{pmatrix}$，$\sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}$ |
| 泡利矩阵性质 | $\sigma_\alpha^2=1$，$\sigma_x\sigma_y=\ii\sigma_z$（轮换），$\sigma_\alpha\sigma_\beta+\sigma_\beta\sigma_\alpha=2\delta_{\alpha\beta}$ |
| 自旋波函数 | $\chi_{1/2}=\begin{pmatrix}1\\0\end{pmatrix}$，$\chi_{-1/2}=\begin{pmatrix}0\\1\end{pmatrix}$；$\Psi=\begin{pmatrix}\psi_1(\vv r,t)\\\psi_2(\vv r,t)\end{pmatrix}$ |
| 任意方向 $\vv n$ 的自旋 | $\sigma_n=\vv\sigma\cdot\vv n$ 的本征值 $\pm1$；$\chi_{\vv n,+}=\begin{pmatrix}\cos\frac\theta2\\\ee^{\ii\varphi}\sin\frac\theta2\end{pmatrix}$，$\chi_{\vv n,-}=\begin{pmatrix}\sin\frac\theta2\\-\ee^{\ii\varphi}\cos\frac\theta2\end{pmatrix}$ |
| 布洛赫球 | 两态的跃迁概率 $\lvert\chi_{\vv n}^\dagger\chi_{\vv n'}\rvert^2=\cos^2\dfrac\Theta2$（$\Theta$ 为 $\vv n,\vv n'$ 夹角）；$\langle\hat{\vv S}\rangle=\dfrac\hbar2\vv n$ |
| 自旋磁矩 | $\vv\mu_s=-\dfrac{e}{m_e}\vv S$（$g_s=2$），$\vv\mu_l=-\dfrac{e}{2m_e}\vv L$（$g_l=1$）；统一写成 $\vv\mu=-g\dfrac{\mu_B}{\hbar}\vv J$ |
| 磁场中的进动 | $\hat H=-\vv\mu_s\cdot\vv B=\mu_B\vv\sigma\cdot\vv B$；拉莫尔频率 $\omega=eB/m_e$ |
| 简单塞曼效应（强场） | $E=E_{nl}+\mu_BB(m_l+2m_s)$；谱线分裂为三条，$\Delta\omega=\dfrac{eB}{2m_e}$（$\Delta m=0,\pm1$） |
| 弱场（反常塞曼效应） | $\Delta E=g_j\mu_BBm_j$，$g_j=1+\dfrac{j(j+1)-l(l+1)+s(s+1)}{2j(j+1)}$ |
| 角动量耦合 | $\vv J=\vv J_1+\vv J_2$：$j=j_1+j_2,\dots,\lvert j_1-j_2\rvert$，$\sum_j(2j+1)=(2j_1+1)(2j_2+1)$ |
| 克莱布什–戈登（$j_2=\tfrac12$） | $\psi_{l,l\pm\frac12,m}=\pm\sqrt{\dfrac{l\pm m+\frac12}{2l+1}}\,Y_{l,m-\frac12}\chi_++\sqrt{\dfrac{l\mp m+\frac12}{2l+1}}\,Y_{l,m+\frac12}\chi_-$ |
| 自旋–轨道 | $j=l\pm\tfrac12$；$\hat{\vv L}\cdot\hat{\vv S}=\tfrac12\left(\hat J^2-\hat L^2-\hat S^2\right)$ |
| 精细结构（自旋–轨道） | $\hat H'=\xi(r)\hat{\vv L}\cdot\hat{\vv S}$，$\xi(r)=\dfrac{1}{2m_e^2c^2}\dfrac1r\dfrac{\dd U}{\dd r}$ |
| 一级修正 | $E^{(1)}_{nlj}=\dfrac{\hbar^2}2\left[j(j+1)-l(l+1)-\tfrac34\right]\zeta_{nl}$，$\zeta_{nl}=\displaystyle\int_0^\infty R_{nl}^2\xi(r)r^2\dd r$ |
| 碱金属双线 | $j=l+\tfrac12$ 上移 $\tfrac l2\hbar^2\zeta$，$j=l-\tfrac12$ 下移 $\tfrac{l+1}2\hbar^2\zeta$，差 $\Delta E=\left(l+\tfrac12\right)\hbar^2\zeta$（钠 D 线 $589.0/589.6\ \mathrm{nm}$，$\Delta E\approx2.1\times10^{-3}\ \mathrm{eV}$） |

## 第 8 章 全同粒子

| 内容 | 公式 |
|---|---|
| 交换算符 | $\hat P_{12}\Psi(q_1,q_2)=\Psi(q_2,q_1)$，本征值 $\pm1$，$[\hat P_{12},\hat H]=0$ |
| 玻色子 / 费米子 | 整数自旋：$\Psi_S$ 对称；半整数自旋：$\Psi_A$ 反对称 |
| 两粒子（无相互作用） | $\Psi_{S,A}=\dfrac{1}{\sqrt2}\left[\phi_i(q_1)\phi_j(q_2)\pm\phi_j(q_1)\phi_i(q_2)\right]$ |
| 斯莱特行列式 | $\Psi_A=\dfrac{1}{\sqrt{N!}}\det\left[\phi_i(q_j)\right]$；两个单粒子态相同 ⇒ $\Psi_A=0$（泡利原理） |
| 两电子自旋态 | 三重态 $S=1$：$\chi_{11}=\uparrow\uparrow$，$\chi_{10}=\tfrac{1}{\sqrt2}(\uparrow\downarrow+\downarrow\uparrow)$，$\chi_{1,-1}=\downarrow\downarrow$（对称）；单态 $S=0$：$\chi_{00}=\tfrac1{\sqrt2}(\uparrow\downarrow-\downarrow\uparrow)$（反对称） |
| 总自旋 | $\hat S^2=\hbar^2\left[\tfrac32+\tfrac12\hat{\vv\sigma}_1\cdot\hat{\vv\sigma}_2\right]$，$\langle\vv\sigma_1\cdot\vv\sigma_2\rangle=1$（三重态），$-3$（单态） |
| 氦原子（微扰） | $E=E_1^{(0)}+E_2^{(0)}+K\pm J$；正氦（三重态）取 $-J$，仲氦（单态）取 $+J$ |
| 库仑积分、交换积分 | $K=\iint\dfrac{e^2}{4\pi\varepsilon_0r_{12}}\lvert\psi_1(1)\rvert^2\lvert\psi_2(2)\rvert^2$，$J=\iint\dfrac{e^2}{4\pi\varepsilon_0r_{12}}\psi_1^*(1)\psi_2^*(2)\psi_1(2)\psi_2(1)$ |
| 氢分子（海特勒–伦敦） | $E_\pm=2E_{\mathrm H}+\dfrac{e^2}{4\pi\varepsilon_0R_{AB}}+\dfrac{K\pm J}{1\pm\Delta^2}$：$E_+$ 为自旋单态（空间对称），$E_-$ 为三重态（空间反对称）；$\Delta$ 为重叠积分，$K,J$ 含电子–核吸引项，见 (8.5.8–9)。成键距离上 $J<0$，单态形成稳定化学键 |
