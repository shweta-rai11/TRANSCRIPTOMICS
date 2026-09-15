# Figure: female vs male RA log2FC for the 53 diagnosis-by-sex interaction genes (FDR<0.05, qnorm + batch_full model),
# coloured by per-sex pattern; PTPN22 and the 7 cell-composition-robust genes are labelled.
# Source: 11_SEX_INTERACTION_TESTING/TABLE/DEG_interaction_significant.csv
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = pd.read_csv("11_SEX_INTERACTION_TESTING/TABLE/DEG_interaction_significant.csv")
order  = ["MALE-restricted", "MAGNITUDE difference", "OPPOSITE direction", "FEMALE-restricted", "neither sex significant"]
colors = dict(zip(order, ["#2a78d6", "#eb6834", "#1baf7a", "#e87ba4", "#eda100"]))
markers= dict(zip(order, ["o", "s", "^", "D", "v"]))
robust = d.loc[d.interaction_FDR_celladj < 0.05, "gene"].tolist()

fig, ax = plt.subplots(figsize=(7.2, 6.2), dpi=200)
fig.patch.set_facecolor("#fcfcfb"); ax.set_facecolor("#fcfcfb")
lim = 1.55
ax.plot([-lim, lim], [-lim, lim], color="#c3c2b7", lw=1, ls="--", zorder=1)
ax.axhline(0, color="#c3c2b7", lw=0.8, zorder=1); ax.axvline(0, color="#c3c2b7", lw=0.8, zorder=1)
for p in order:
    s = d[d.pattern == p]
    n = len(s)
    ax.scatter(s.female_logFC, s.male_logFC, s=55, marker=markers[p], c=colors[p],
               edgecolors="#fcfcfb", linewidths=1.2, label=f"{p} (n = {n})", zorder=3)
off = {"MAP4K5": (-75, 16), "RIF1": (-10, 26), "BCLAF1": (-58, -12), "ZNF800": (-60, -12),
       "UBXN4": (8, 4), "RAB2A": (8, -10), "KDM5D": (8, -2)}
for _, r in d[d.gene.isin(robust)].iterrows():
    ax.annotate(r.gene, (r.female_logFC, r.male_logFC), xytext=off[r.gene], textcoords="offset points",
                fontsize=8.5, color="#1a1a19", zorder=4,
                arrowprops=dict(arrowstyle="-", color="#6b6b68", lw=0.6, shrinkA=0, shrinkB=3))
px, py = d.loc[d.gene=="PTPN22","female_logFC"].item(), d.loc[d.gene=="PTPN22","male_logFC"].item()
ax.scatter([px], [py], s=130, facecolors="none", edgecolors="#1a1a19", linewidths=1.2, zorder=5)
ax.annotate("PTPN22\nfemale +0.17 (FDR 0.046)\nmale +0.89 (FDR 1.2 × 10⁻⁶)\ninteraction +0.58 (FDR 0.038)",
            (px, py), xytext=(0.55, 1.05), textcoords="data", fontsize=8, color="#1a1a19", va="bottom",
            fontweight="normal", arrowprops=dict(arrowstyle="-", color="#1a1a19", lw=0.8, shrinkB=6), zorder=6)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
ax.set_xlabel("Female RA vs control, log₂ fold change", fontsize=10)
ax.set_ylabel("Male RA vs control, log₂ fold change", fontsize=10)
ax.set_title("Diagnosis × sex interaction genes (53 at FDR < 0.05)\nlabelled: PTPN22 and the 7 genes robust to cell-composition adjustment",
             fontsize=10.5, loc="left")
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
for sp in ["left", "bottom"]: ax.spines[sp].set_color("#c3c2b7")
ax.tick_params(colors="#4a4a48", labelsize=9)
ax.legend(frameon=False, fontsize=8, loc="lower right", title="Per-sex pattern", title_fontsize=8.5)
ax.text(-lim + 0.05, -lim + 0.05, "dashed line: equal effect in both sexes", fontsize=7.5, color="#6b6b68")
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(f"11_SEX_INTERACTION_TESTING/FIGURE/fig_interaction_female_vs_male_logFC.{ext}", facecolor=fig.get_facecolor())
print("saved")
