// 统一的量子力学记号宏：狄拉克符号 \ket{\psi}、\bra{\phi}、\ip{\phi}{\psi}、\mel{m}{\hat A}{n}、\expval{\hat A}、\comm{A}{B} ...
const QMMACROS = {
  "\\ket": "\\left|#1\\right\\rangle",
  "\\bra": "\\left\\langle#1\\right|",
  "\\ip": "\\left\\langle#1\\middle|#2\\right\\rangle",
  "\\mel": "\\left\\langle#1\\middle|#2\\middle|#3\\right\\rangle",
  "\\expval": "\\left\\langle#1\\right\\rangle",
  "\\comm": "\\left[#1,#2\\right]",
  "\\op": "\\hat{#1}",
  "\\vv": "\\boldsymbol{#1}",
  "\\dd": "\\mathrm{d}", "\\ii": "\\mathrm{i}", "\\ee": "\\mathrm{e}",
  "\\Tr": "\\operatorname{Tr}",
  "\\Re": "\\operatorname{Re}", "\\Im": "\\operatorname{Im}"
};
document$.subscribe(() => {
  renderMathInElement(document.body, {
    delimiters: [
      { left: "$$",  right: "$$",  display: true },
      { left: "$",   right: "$",   display: false },
      { left: "\\(", right: "\\)", display: false },
      { left: "\\[", right: "\\]", display: true }
    ],
    macros: QMMACROS,
    throwOnError: false
  });
});
