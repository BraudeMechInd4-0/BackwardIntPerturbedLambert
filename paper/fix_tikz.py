#!/usr/bin/env python3
"""
paper/fix_tikz.py
Transform benchmark_results tikz files to publication-ready format for JGCD.

Read-only access to benchmark_results/.
Writes EXCLUSIVELY to paper/newImages/.
"""

from pathlib import Path
import re

# ─── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR     = Path(__file__).parent           # paper/
ROOT           = SCRIPT_DIR.parent               # AMDR root
SOURCE_FIGURES = ROOT / "benchmark_results" / "paper_figures"
SOURCE_CONV    = ROOT / "benchmark_results" / "convergence_analysis"
OUTPUT_DIR     = SCRIPT_DIR / "newImages"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ─── Safety guard ──────────────────────────────────────────────────────────────
def safe_write(path: Path, content: str) -> None:
    resolved = str(path.resolve())
    guard    = str(OUTPUT_DIR.resolve())
    assert resolved.startswith(guard), \
        f"SAFETY: refusing to write outside paper/newImages/: {path}"
    path.write_text(content, encoding="utf-8")
    print(f"  [OK] {path.name}")

# ─── Color mappings ────────────────────────────────────────────────────────────

# For convergence trajectory plots: solver → colorlet expression (with !40 opacity)
CONV_COLORLETS = {
    "T":     "cyan!85!black!40",
    "TR":    "blue!85!black!40",
    "TRB":   "blue!80!cyan!40",
    "TB":    "violet!70!blue!40",
    "MPS":   "orange!85!black!40",
    "MPSR":  "red!85!black!40",
    "MPSRB": "red!80!black!40",
    "MPSB":  "orange!70!yellow!40",
    "B":     "green!65!black!40",
    "BR":    "olive!80!black!40",
    "BRB":   "green!70!yellow!40",
}

# Darker color for the median line — matches AST-2026 hand-tuned values
CONV_MEDIAN_COLORS = {
    "T":     "cyan!85!black!80",
    "TR":    "blue!85!black!80",
    "TRB":   "blue!80!cyan",
    "TB":    "violet!70!blue",
    "MPS":   "orange!85!black",
    "MPSR":  "red!85!black!80",
    "MPSRB": "red!80!black",
    "MPSB":  "orange!70!yellow",
    "B":     "green!65!black!80",
    "BR":    "olive!80!black!80",
    "BRB":   "green!70!yellow",
}

# For bar charts: replaces the 7-color definecolor block with 8 family colorlets
BAR_COLORLETS = "\n".join([
    r"\colorlet{mycolor1}{cyan}%",
    r"\colorlet{mycolor2}{orange!85!black}%",
    r"\colorlet{mycolor3}{green!60!black}%",
    r"\colorlet{mycolor4}{violet!70!blue}%",
    r"\colorlet{mycolor5}{blue!80!cyan}%",
    r"\colorlet{mycolor6}{orange!70!yellow}%",
    r"\colorlet{mycolor7}{red!80!black}%",
    r"\colorlet{mycolor8}{green!70!yellow}%",
])

# Solver order in single-bar charts (x = 1..8)
SOLVERS_8 = ["T", "MPS", "B", "TB", "TRB", "MPSB", "MPSRB", "BRB"]

# Fill colors per solver for per-solver addplot reconstruction
SOLVER_FILLS = {
    "T":     "cyan",
    "MPS":   "orange!85!black",
    "B":     "green!60!black",
    "TB":    "violet!70!blue",
    "TRB":   "blue!80!cyan",
    "MPSB":  "orange!70!yellow",
    "MPSRB": "red!80!black",
    "BRB":   "green!70!yellow",
}

# ─── Shared transforms ─────────────────────────────────────────────────────────

def comment_title(content: str) -> str:
    """Comment out title style and title lines (handles any indentation)."""
    content = re.sub(r'(\n\s*)(title style=\{[^}]+\},?)', r'\1%\2', content)
    content = re.sub(r'(\n\s*)(title=\{[^}]*\},?)',        r'\1%\2', content)
    return content


def replace_definecolor7_block(content: str) -> str:
    """Replace the 7-line \\definecolor block with 8 family colorlet lines."""
    pattern = re.compile(
        r'(?:\\definecolor\{mycolor\d+\}\{rgb\}\{[^}]+\}%\n){2,}',
        re.MULTILINE,
    )
    replacement = BAR_COLORLETS + "\n"
    replaced, n = pattern.subn(lambda _: replacement, content, count=1)
    if n == 0:
        print("    [WARN] definecolor block not found – skipping color replacement")
    return replaced


def fix_brb_fill(content: str) -> str:
    """
    Grouped bar charts (legend-based): change BRB fill from mycolor1 → mycolor8.
    Finds the section after the LAST \\addlegendentry{MPSRB} and replaces there.
    """
    if r"\addlegendentry{BRB}" not in content:
        return content
    marker = r"\addlegendentry{MPSRB}"
    idx = content.rfind(marker)
    if idx == -1:
        return content
    idx += len(marker)
    pre, post = content[:idx], content[idx:]
    post = re.sub(r"fill=mycolor1\b", "fill=mycolor8", post)
    return pre + post


def fix_brb_distribution_color(content: str) -> str:
    """
    Distribution plots (no legend): change color=mycolor1 → color=mycolor8
    only in \\addplot blocks whose data x-values are near 8 (BRB position).
    """
    def fixer(match):
        full = match.group(0)
        opts = match.group(1)
        data = match.group(2)
        if "mycolor1" not in opts:
            return full
        # Extract first column of each data row
        x_vals = re.findall(r"^\s*([\d.]+)\s", data, re.MULTILINE)
        if x_vals and any(7.0 <= float(x) <= 9.0 for x in x_vals):
            return full.replace("color=mycolor1", "color=mycolor8", 1)
        return full

    pattern = re.compile(
        r"\\addplot \[([^\]]*)\]\n\s+table\[row sep=crcr\]\{%\n(.*?)\};",
        re.DOTALL,
    )
    return pattern.sub(fixer, content)


# ─── Batch 1: Convergence trajectory plots ─────────────────────────────────────

def process_conv(src: Path, dest: Path, solver: str) -> None:
    c = src.read_text(encoding="utf-8")

    # 1. Replace \\definecolor{mycolor1} with \\colorlet + add mycolor2 for text
    new_defs = (
        f"\\colorlet{{mycolor1}}{{{CONV_COLORLETS[solver]}}}%\n"
        f"\\definecolor{{mycolor2}}{{rgb}}{{0.12941,0.12941,0.12941}}%"
    )
    c = re.sub(
        r"\\definecolor\{mycolor1\}\{rgb\}\{[^}]+\}%",
        lambda _: new_defs,
        c,
    )

    # 2. Fix ylabel to LaTeX math notation
    ylabel_replacement = r"ylabel={Position error $\|\delta\mathbf{r}_f\|$ (km)}"
    c = re.sub(
        r"ylabel=\{[^}]*\}",
        lambda _: ylabel_replacement,
        c,
        count=1,
    )

    # 3. Comment title lines
    c = comment_title(c)

    # 4. Replace axis-text color
    c = c.replace("white!15!black", "mycolor2")

    # 5. Add thin-line global style so background traces don't overpower the median
    c = c.replace(
        'legend style={legend cell align=left, align=left, draw=mycolor2}',
        'legend style={legend cell align=left, align=left, draw=mycolor2},\nevery axis plot/.style={line width=0.2pt}',
    )

    # 6. Make median line darker so it stands out from the faded background traces
    median_color = CONV_MEDIAN_COLORS[solver]
    c = re.sub(
        r'\\addplot \[color=mycolor1, line width=2\.5pt\]',
        lambda _: f'\\addplot [color={median_color}, line width=2.5pt]',
        c,
        count=1,
    )

    safe_write(dest, c)


# ─── Batch 2a: Grouped bar charts (separate \\addplot per solver, with legend) ─

def process_grouped_bar(src: Path, dest: Path) -> None:
    c = src.read_text(encoding="utf-8")
    c = replace_definecolor7_block(c)
    c = fix_brb_fill(c)
    c = comment_title(c)
    safe_write(dest, c)


# ─── Batch 2b: Distribution box plots (separate \\addplot per solver, no legend) ─

def process_distribution(src: Path, dest: Path) -> None:
    c = src.read_text(encoding="utf-8")
    c = replace_definecolor7_block(c)
    c = fix_brb_distribution_color(c)
    c = comment_title(c)
    safe_write(dest, c)


# ─── Batch 2c: Single-addplot bar charts (mean_iterations and variants) ────────

def _parse_table(raw: str):
    """Parse 'X\\tY\\\\' rows → list of float lists."""
    rows = []
    for line in raw.strip().splitlines():
        vals = line.strip().rstrip("\\").split()
        if vals:
            try:
                rows.append([float(v) for v in vals])
            except ValueError:
                pass
    return rows


def process_mean_iterations(src: Path, dest: Path) -> None:
    c = src.read_text(encoding="utf-8")
    c = c.replace('bar shift auto,\n', '')

    # Extract main ybar data block
    m_main = re.search(
        r"\\addplot\[ybar[^\]]*\] table\[row sep=crcr\] \{%\n(.*?)\};",
        c, re.DOTALL,
    )
    if not m_main:
        print(f"    [WARN] No ybar block found in {src.name} – copying as-is")
        safe_write(dest, comment_title(c))
        return

    main_rows = _parse_table(m_main.group(1))
    data = {int(r[0]): r[1] for r in main_rows}

    # Extract error bar data block
    m_err = re.search(
        r"y error plus index=2, y error minus index=3\]\{%\n(.*?)\};",
        c, re.DOTALL,
    )
    errs = {}
    if m_err:
        for r in _parse_table(m_err.group(1)):
            if len(r) >= 4:
                errs[int(r[0])] = (r[2], r[3])

    # Build per-solver addplots
    per_solver = []
    err_data_rows = []
    for i, solver in enumerate(SOLVERS_8):
        x = i + 1
        y = data.get(x, 0.0)
        ep, em = errs.get(x, (0.0, 0.0))
        per_solver.append(
            f"\\addplot[ybar, bar width=0.8, fill={SOLVER_FILLS[solver]}, draw=black] "
            f"coordinates {{({x}, {y:.10g})}};\n"
        )
        err_data_rows.append(f"{x}\t{y:.10g}\t{ep:.10g}\t{em:.10g}\\\\")

    err_block = (
        "\\addplot [color=black, line width=1.5pt, only marks, mark=*, "
        "mark options={solid, black}, forget plot]\n"
        " plot [error bars/.cd, y dir=both, y explicit, "
        "error bar style={line width=1.5pt}, "
        "error mark options={line width=1.5pt, mark size=6.0pt, rotate=90}]\n"
        " table[row sep=crcr, y error plus index=2, y error minus index=3]{%\n"
        + "\n".join(err_data_rows) + "\n};\n"
    )

    # Find optional red reference line (norm variants have it)
    m_red = re.search(
        r"(\\addplot \[color=red, dashed.*?\];)", c[m_main.start():], re.DOTALL
    )
    red_str = m_red.group(1) + "\n" if m_red else ""

    # Find the baseline forget-plot line
    m_base = re.search(
        r"(\\addplot\[forget plot, color=white!15!black\].*?\};)",
        c[m_main.start():], re.DOTALL,
    )
    base_str = m_base.group(1) + "\n" if m_base else ""

    # Reassemble: pre-axis-block + axis header + per-solver plots + ... + \end{axis}...
    pre = c[: m_main.start()]
    m_end = re.search(r"\\end\{axis\}", c[m_main.start():])
    post = c[m_main.start() + m_end.start():]

    result = pre + "".join(per_solver) + base_str + err_block + red_str + post
    result = comment_title(result)
    safe_write(dest, result)


# ─── Batch 2d: Speedup comparison (5 bars, per-bar colors) ────────────────────

SPEEDUP_COMPARISONS = [
    (r"$\text{T}\rightarrow\text{TB}$",       "violet!70!blue",   1),
    (r"$\text{T}\rightarrow\text{TRB}$",      "blue!80!cyan",     2),
    (r"$\text{MPS}\rightarrow\text{MPSB}$",   "orange!70!yellow", 3),
    (r"$\text{MPS}\rightarrow\text{MPSRB}$",  "red!80!black",     4),
    (r"$\text{B}\rightarrow\text{BRB}$",      "green!70!yellow",  5),
]


def process_speedup(src: Path, dest: Path) -> None:
    c = src.read_text(encoding="utf-8")

    # Extract main ybar data (5 bars — MATLAB now generates all 5 after re-run)
    m_main = re.search(
        r"\\addplot\[ybar[^\]]*\] table\[row sep=crcr\] \{%\n(.*?)\};",
        c, re.DOTALL,
    )
    if not m_main:
        print(f"    [WARN] No ybar block in {src.name} – copying as-is")
        safe_write(dest, comment_title(c))
        return

    main_rows = _parse_table(m_main.group(1))
    data = {int(r[0]): r[1] for r in main_rows}

    # Extract error bars
    m_err = re.search(
        r"y error plus index=2, y error minus index=3\]\{%\n(.*?)\};",
        c, re.DOTALL,
    )
    errs = {}
    if m_err:
        for r in _parse_table(m_err.group(1)):
            if len(r) >= 4:
                errs[int(r[0])] = (r[2], r[3])

    # Update axis x-range and tick labels for 5 bars
    c2 = re.sub(r"xmax=\d+\.?\d*", "xmax=6.2", c)
    c2 = re.sub(r"xtick=\{[\d,]+\}", "xtick={1,2,3,4,5}", c2)
    # Build xticklabels string and replace using lambda (avoids backslash escape issues)
    new_xticklabels = "xticklabels={" + ",".join(f"{{{label}}}" for label, _, _ in SPEEDUP_COMPARISONS) + "}"
    c2 = re.sub(r"xticklabels=\{[^\n]*\}", lambda _: new_xticklabels, c2)
    # Extend baseline/reference lines' upper x-bound to 6.2.
    # Matches lines like "5.2\t0\\" — only exact-zero or exact-one y-values appear in
    # forget-plot/reference-line tables; bar data has fractional y-values.
    c2 = re.sub(r"\b\d+\.?\d*\t0\\\\", lambda _: "6.2\t0\\\\", c2)
    c2 = re.sub(r"\b\d+\.?\d*\t1\\\\", lambda _: "6.2\t1\\\\", c2)
    c2 = c2.replace('bar shift auto,\n', '')

    # Build per-comparison addplots
    per_comp = []
    err_data_rows = []
    for _, color, x in SPEEDUP_COMPARISONS:
        y = data.get(x, 0.0)
        ep, em = errs.get(x, (0.0, 0.0))
        per_comp.append(
            f"\\addplot[ybar, bar width=0.8, fill={color}, draw=black] "
            f"coordinates {{({x}, {y:.10g})}};\n"
        )
        err_data_rows.append(f"{x}\t{y:.10g}\t{ep:.10g}\t{em:.10g}\\\\")

    err_block = (
        "\\addplot [color=black, line width=1.5pt, only marks, mark=*, "
        "mark options={solid, black}, forget plot]\n"
        " plot [error bars/.cd, y dir=both, y explicit, "
        "error bar style={line width=1.5pt}, "
        "error mark options={line width=1.5pt, mark size=6.0pt, rotate=90}]\n"
        " table[row sep=crcr, y error plus index=2, y error minus index=3]{%\n"
        + "\n".join(err_data_rows) + "\n};\n"
    )

    m_main2 = re.search(
        r"\\addplot\[ybar[^\]]*\] table\[row sep=crcr\] \{%\n(.*?)\};",
        c2, re.DOTALL,
    )

    m_base = re.search(
        r"(\\addplot\[forget plot, color=white!15!black\].*?\};)",
        c2[m_main2.start():], re.DOTALL,
    )
    base_str = m_base.group(1) + "\n" if m_base else ""

    m_red = re.search(
        r"(\\addplot \[color=red, dashed.*?\];)", c2[m_main2.start():], re.DOTALL
    )
    red_str = m_red.group(1) + "\n" if m_red else ""

    pre = c2[: m_main2.start()]
    m_end = re.search(r"\\end\{axis\}", c2[m_main2.start():])
    post = c2[m_main2.start() + m_end.start():]

    result = pre + "".join(per_comp) + base_str + err_block + red_str + post
    result = comment_title(result)
    safe_write(dest, result)


# ─── Batch 3: Switch sensitivity (heatmaps – comment titles only) ──────────────

def process_switch(src: Path, dest: Path) -> None:
    c = src.read_text(encoding="utf-8")
    c = comment_title(c)
    safe_write(dest, c)


# ─── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    errors = 0

    # ── Batch 1: Convergence trajectory plots ──
    print("\n=== Batch 1: Convergence trajectory plots ===")
    for solver in CONV_COLORLETS:
        src = SOURCE_CONV / f"conv_{solver}.tikz"
        if src.exists():
            try:
                process_conv(src, OUTPUT_DIR / f"conv_{solver}.tikz", solver)
            except Exception as e:
                print(f"  [ERR] conv_{solver}.tikz: {e}"); errors += 1
        else:
            print(f"  [SKIP] conv_{solver}.tikz not found in source")

    # ── Batch 2a: Grouped bar charts with legend ──
    print("\n=== Batch 2a: Grouped bar charts ===")
    grouped_bar_bases = [
        "convergence_vs_min_altitude",
        "success_vs_revolutions",
        "convergence_by_regime",
    ]
    suffixes_6 = ["", "_hard", "_nonhard", "_norm", "_norm_hard", "_norm_nonhard"]
    for base in grouped_bar_bases:
        for sfx in suffixes_6:
            src = SOURCE_FIGURES / f"{base}{sfx}.tikz"
            if src.exists():
                try:
                    process_grouped_bar(src, OUTPUT_DIR / src.name)
                except Exception as e:
                    print(f"  [ERR] {src.name}: {e}"); errors += 1
            else:
                print(f"  [SKIP] {src.name}")

    # ── Batch 2b: Distribution plots (no legend) ──
    print("\n=== Batch 2b: Distribution plots ===")
    dist_bases = ["runtime_distribution", "iteration_distribution"]
    for base in dist_bases:
        for sfx in suffixes_6:
            src = SOURCE_FIGURES / f"{base}{sfx}.tikz"
            if src.exists():
                try:
                    process_distribution(src, OUTPUT_DIR / src.name)
                except Exception as e:
                    print(f"  [ERR] {src.name}: {e}"); errors += 1
            else:
                print(f"  [SKIP] {src.name}")

    # ── Batch 2c: mean_iterations (single addplot → per-solver) ──
    print("\n=== Batch 2c: Mean iterations charts ===")
    for sfx in suffixes_6:
        src = SOURCE_FIGURES / f"mean_iterations{sfx}.tikz"
        if src.exists():
            try:
                process_mean_iterations(src, OUTPUT_DIR / src.name)
            except Exception as e:
                print(f"  [ERR] {src.name}: {e}"); errors += 1
        else:
            print(f"  [SKIP] mean_iterations{sfx}.tikz")

    # ── Batch 2d: Speedup comparison ──
    print("\n=== Batch 2d: Speedup comparison ===")
    for sfx in ["", "_hard", "_nonhard"]:
        src = SOURCE_FIGURES / f"speedup_comparison{sfx}.tikz"
        if src.exists():
            try:
                process_speedup(src, OUTPUT_DIR / src.name)
            except Exception as e:
                print(f"  [ERR] {src.name}: {e}"); errors += 1
        else:
            print(f"  [SKIP] speedup_comparison{sfx}.tikz")

    # ── Batch 3: Switch sensitivity ──
    print("\n=== Batch 3: Switch sensitivity ===")
    for kind in ["conv", "runtime"]:
        for sfx in ["", "_hard", "_nonhard"]:
            src = SOURCE_FIGURES / f"switch_sensitivity_{kind}{sfx}.tikz"
            if src.exists():
                try:
                    process_switch(src, OUTPUT_DIR / src.name)
                except Exception as e:
                    print(f"  [ERR] {src.name}: {e}"); errors += 1
            else:
                print(f"  [SKIP] switch_sensitivity_{kind}{sfx}.tikz")

    # ── Summary ──
    written = list(OUTPUT_DIR.glob("*.tikz"))
    print(f"\n{'='*50}")
    print(f"Done. {len(written)} files written to {OUTPUT_DIR}")
    if errors:
        print(f"WARNING: {errors} error(s) occurred – review output above.")
    else:
        print("No errors.")


if __name__ == "__main__":
    main()
