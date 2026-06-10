from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


OUT_DIR = Path(__file__).resolve().parent
SVG_OUT = OUT_DIR / "Figure_8_Adnectin_PCSK9_evidence_chain.svg"
PNG_OUT = OUT_DIR / "Figure_8_Adnectin_PCSK9_evidence_chain.png"

rcParams["font.sans-serif"] = [
    "Microsoft YaHei",
    "SimHei",
    "Arial Unicode MS",
    "DejaVu Sans",
]
rcParams["axes.unicode_minus"] = False


INK = "#182230"
MUTED = "#475467"
BLUE = "#2563eb"
BLUE_2 = "#dbeafe"
TEAL = "#0f766e"
TEAL_2 = "#ccfbf1"
AMBER = "#b45309"
AMBER_2 = "#fef3c7"
ROSE = "#be123c"
ROSE_2 = "#ffe4e6"
GREEN = "#15803d"
GREEN_2 = "#dcfce7"
LINE = "#d0d7e2"


def draw_panel(ax, x, y, w, h, title, tag, color):
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=1.3,
        edgecolor=LINE,
        facecolor="white",
    )
    ax.add_patch(box)
    ax.add_patch(Rectangle((x, y + h - 0.055), w, 0.055, facecolor=color, edgecolor="none"))
    ax.text(x + 0.018, y + h - 0.035, tag, color="white", fontsize=13, fontweight="bold", va="center")
    ax.text(x + 0.07, y + h - 0.035, title, color="white", fontsize=12.5, fontweight="bold", va="center")


def small_box(ax, x, y, w, h, text, fc, ec=None, fs=9.5, weight="normal", color=INK, ha="center"):
    if ec is None:
        ec = fc
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.008,rounding_size=0.012",
        linewidth=1,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(box)
    ax.text(
        x + w / 2 if ha == "center" else x + 0.012,
        y + h / 2,
        text,
        ha=ha,
        va="center",
        fontsize=fs,
        color=color,
        fontweight=weight,
        linespacing=1.15,
    )


def arrow(ax, x1, y1, x2, y2, color="#64748b", rad=0):
    patch = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle="-|>",
        mutation_scale=14,
        linewidth=1.4,
        color=color,
        connectionstyle=f"arc3,rad={rad}",
    )
    ax.add_patch(patch)


def main():
    fig = plt.figure(figsize=(15.8, 8.8), dpi=220)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("#f8fafc")

    ax.text(
        0.04,
        0.965,
        "Figure 8 | PCSK9 antibody reference landscape and Adnectin-FBBP scaffold-centric evidence chain",
        fontsize=17,
        fontweight="bold",
        color=INK,
        va="top",
    )
    ax.text(
        0.04,
        0.928,
        "成熟 PCSK9 抗体药物生态作为参照；FBBP-DB 展示小型 Adnectin/FN3 binder 的模块化序列、loop/interface、柔性与开发性证据。",
        fontsize=10.5,
        color=MUTED,
        va="top",
    )

    draw_panel(ax, 0.035, 0.58, 0.30, 0.30, "PCSK9 therapeutic antibody panel", "A", BLUE)
    draw_panel(ax, 0.365, 0.58, 0.27, 0.30, "FBBP query and target split", "B", TEAL)
    draw_panel(ax, 0.665, 0.58, 0.30, 0.30, "Structure and active-loop interface", "C", AMBER)
    draw_panel(ax, 0.035, 0.10, 0.45, 0.36, "Sequence, format and developability", "D", GREEN)
    draw_panel(ax, 0.515, 0.10, 0.45, 0.36, "Added value of FBBP", "E", ROSE)

    # Panel A
    ax.text(0.055, 0.815, "Thera-SAbDab/SAbDab PCSK9 records", fontsize=10.5, fontweight="bold", color=INK)
    small_box(ax, 0.055, 0.765, 0.075, 0.045, "10\nrecords", BLUE_2, "#93c5fd", fs=10, weight="bold", color=BLUE)
    small_box(ax, 0.142, 0.765, 0.075, 0.045, "6\napproved", GREEN_2, "#86efac", fs=10, weight="bold", color=GREEN)
    small_box(ax, 0.229, 0.765, 0.075, 0.045, "4\nclinical/stop", AMBER_2, "#fcd34d", fs=9.2, weight="bold", color=AMBER)

    for i, name in enumerate(["alirocumab", "evolocumab", "ebronucimab", "ongericimab", "recaticimab", "tafolecimab"]):
        x = 0.055 + (i % 2) * 0.135
        y = 0.715 - (i // 2) * 0.045
        small_box(ax, x, y, 0.12, 0.032, name, "#f1f5f9", "#d8e0ea", fs=8.2, color=INK)
    ax.text(
        0.055,
        0.585,
        "Structure field: bococizumab 100% match;\nralpancizumab 95-98% match; others none in download table.",
        fontsize=8.8,
        color=MUTED,
        va="bottom",
    )

    # Panel B
    small_box(ax, 0.39, 0.785, 0.095, 0.05, "1,996\nFBBP records", TEAL_2, "#5eead4", fs=10, weight="bold", color=TEAL)
    small_box(ax, 0.505, 0.785, 0.095, 0.05, "8\nPCSK9 text hits", "#ecfeff", "#67e8f9", fs=10, weight="bold", color="#0e7490")
    arrow(ax, 0.485, 0.81, 0.505, 0.81, TEAL)
    small_box(ax, 0.39, 0.705, 0.095, 0.055, "4 direct\nPCSK9 binders", GREEN_2, "#86efac", fs=9.5, weight="bold", color=GREEN)
    small_box(ax, 0.505, 0.705, 0.095, 0.055, "4 ALB/HSA\ncontext", AMBER_2, "#fcd34d", fs=9.5, weight="bold", color=AMBER)
    arrow(ax, 0.552, 0.785, 0.552, 0.76, TEAL)
    arrow(ax, 0.552, 0.705, 0.487, 0.705, TEAL, rad=0.2)
    ax.text(0.385, 0.635, "Curated target fields keep\ntarget binding separate from\nhalf-life/fusion context.", fontsize=9.2, color=MUTED, va="top")

    # Panel C
    ax.text(0.69, 0.815, "1459D05-PCSK9 complex (PDB 4OV6)", fontsize=10.2, fontweight="bold", color=INK)
    small_box(ax, 0.695, 0.748, 0.085, 0.055, "FN3\nAdnectin", TEAL_2, "#5eead4", fs=8.9, weight="bold", color=TEAL)
    small_box(ax, 0.845, 0.742, 0.09, 0.065, "PCSK9", BLUE_2, "#93c5fd", fs=10.5, weight="bold", color=BLUE)
    arrow(ax, 0.78, 0.775, 0.845, 0.775, AMBER)
    labels = [("Loop 1", 68.1, 1, "#fde68a"), ("Loop 2", 0, 0, "#e5e7eb"), ("Loop 3", 524.8, 8, "#fb923c")]
    y0 = 0.700
    for i, (lab, val, residues, color) in enumerate(labels):
        y = y0 - i * 0.055
        ax.text(0.69, y + 0.014, lab, fontsize=9.2, color=INK, va="center")
        ax.add_patch(Rectangle((0.755, y), 0.16, 0.026, facecolor="#f1f5f9", edgecolor="#d8e0ea", linewidth=0.8))
        width = 0 if val == 0 else 0.16 * (val / 524.8)
        ax.add_patch(Rectangle((0.755, y), width, 0.026, facecolor=color, edgecolor="none"))
        ax.text(0.925, y + 0.013, f"{residues} res; {val:.1f} A2", fontsize=8.0, color=MUTED, va="center", ha="left")
    ax.text(0.69, 0.595, "Loop 3 dominates; Loop 1 is auxiliary; Loop 2 is non-contacting under 4.5 A.", fontsize=8.3, color=MUTED, va="bottom")

    # Panel D
    small_box(ax, 0.06, 0.352, 0.12, 0.05, "Whole mAb\n~150 kDa", BLUE_2, "#93c5fd", fs=9.5, weight="bold", color=BLUE)
    small_box(ax, 0.205, 0.352, 0.12, 0.05, "single FN3\n103 aa", TEAL_2, "#5eead4", fs=9.5, weight="bold", color=TEAL)
    small_box(ax, 0.35, 0.352, 0.10, 0.05, "tandem\n219 aa", GREEN_2, "#86efac", fs=9.5, weight="bold", color=GREEN)
    ax.text(0.06, 0.312, "1459D05 vs BMS-962476", fontsize=10.2, fontweight="bold", color=INK)
    small_box(ax, 0.06, 0.265, 0.12, 0.04, "90.3% identity", "#f1f5f9", "#d8e0ea", fs=9.2, color=INK)
    small_box(ax, 0.195, 0.265, 0.13, 0.04, "7/10 substitutions\nin active loops", AMBER_2, "#fcd34d", fs=8.5, color=AMBER, weight="bold")
    small_box(ax, 0.34, 0.265, 0.11, 0.04, "4 in Loop 3", ROSE_2, "#fecdd3", fs=9.2, color=ROSE, weight="bold")
    ax.text(0.06, 0.213, "FBBP evidence-card fields", fontsize=10.2, fontweight="bold", color=INK)
    fields = [
        "KD 1.58 / 1.3 / 0.515 nM",
        "loop flexibility profiles",
        "predicted soluble expression",
        "low oral rating, high risk flag",
    ]
    for i, field in enumerate(fields):
        small_box(ax, 0.06 + (i % 2) * 0.205, 0.170 - (i // 2) * 0.047, 0.185, 0.034, field, "#f8fafc", "#d8e0ea", fs=8.3, color=MUTED)

    # Panel E
    items = [
        ("Antibody panel", "mature clinical landscape;\nwhole-mAb reference formats", BLUE_2, BLUE),
        ("FBBP query", "target/context spli<local_path_removed>", TEAL_2, TEAL),
        ("FBBP structure", "compact scaffold with\nresidue-level loop interface", AMBER_2, AMBER),
        ("FBBP resource", "sequence, loop, flexibility,\ndevelopability and provenance", ROSE_2, ROSE),
    ]
    for i, (head, body, fc, color) in enumerate(items):
        y = 0.352 - i * 0.068
        small_box(ax, 0.545, y, 0.13, 0.046, head, fc, color, fs=8.8, weight="bold", color=color)
        ax.text(0.69, y + 0.023, body, fontsize=8.5, color=MUTED, va="center")
    ax.text(
        0.545,
        0.126,
        "Interpretation boundary: FBBP is presented as a compact, modular, loop-interpretable\nresource layer, not as a categorical replacement for approved PCSK9 antibodies.",
        fontsize=8.6,
        color=INK,
        va="bottom",
    )

    arrow(ax, 0.335, 0.73, 0.365, 0.73, "#94a3b8")
    arrow(ax, 0.635, 0.73, 0.665, 0.73, "#94a3b8")
    arrow(ax, 0.50, 0.58, 0.50, 0.46, "#94a3b8")

    ax.text(0.04, 0.055, "Abbreviations: mAb, monoclonal antibody; FN3, fibronectin type III; HSA, human serum albumin. A2 denotes square angstrom in this schematic label.", fontsize=8.5, color="#64748b")
    ax.text(0.04, 0.035, "Source layers: Thera-SAbDab/SAbDab PCSK9 therapeutic antibody fields; FBBP-DB PCSK9-Adnectin records; PDB 4OV6 loop-interface analysis.", fontsize=8.5, color="#64748b")

    fig.savefig(SVG_OUT, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(PNG_OUT, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(SVG_OUT)
    print(PNG_OUT)


if __name__ == "__main__":
    main()
