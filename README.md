# 量子力学从零到精通

《量子力学》讲义。章节结构与周世勋《量子力学教程》（第 3 版）一致，讲第 1–8 章。

**在线阅读：** <https://cockroach-f.github.io/quantum-mechanics/>

姊妹站点：[复变函数从零到精通](https://cockroach-f.github.io/complex-functions/) · [电磁场与电磁波从零到精通](https://cockroach-f.github.io/em-fields-waves/) · [ISAR 目标识别入门与文献带读](https://cockroach-f.github.io/isar-tutorial/)

## 内容结构

| 部分 | 章节 |
|---|---|
| 绪论与知识地图 | `docs/index.md` |
| 第一部分 量子力学的诞生与基本框架 | 第 1 章 绪论 · 第 2 章 波函数和薛定谔方程 |
| 第二部分 力学量与表象 | 第 3 章 量子力学中的力学量 · 第 4 章 态和力学量的表象 |
| 第三部分 近似方法与散射 | 第 5 章 微扰理论 · 第 6 章 散射 |
| 第四部分 自旋与全同粒子 | 第 7 章 自旋 · 第 8 章 全同粒子 |
| 全课程总结 · 公式速查 | `docs/summary.md` · `docs/formulas.md` |

每章包括：编号的定义与定理框、例题、易错点、知识点归纳总结（思维导图、公式表、方法），以及 8 道带折叠答案的练习题。文中嵌入 manim 动画，源码在 `animations/`，用 `tools/render_videos.sh` 渲染。

## 本地预览

```bash
pip install -r requirements.txt
mkdocs serve
```
