"""Draw the stage-2e figures from maps.json (needs only numpy + matplotlib, no torch).

Figures (PNG, one per file, in --outdir):
  importance.png          permutation importance per condition, L1/L2/L3
  pd_<level>.png          partial dependence P(level) vs each condition
  slice_<x>_x_<y>.png     2D slices (other conditions at the 2b reference), one panel per level
  accuracy.png            COMPASS vs Sobol accuracy curves per level (all and hard)
"""
import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LOG_AXES = {"S", "mu", "food", "km"}
LABELS = {"S": "division size S", "sigma": "parasite advantage σ", "mu": "mutation rate μ", "K": "gene types K",
          "d": "degradation d", "food": "food", "km": "membrane synthesis km", "dm": "membrane turnover dm",
          "z1": "distractor z1", "z2": "distractor z2"}
COLORS = {"L1": "#4C78A8", "L2": "#F58518", "L3": "#54A24B"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default="maps.json")
    ap.add_argument("--outdir", default="figures")
    a = ap.parse_args()
    r = json.load(open(a.json))
    os.makedirs(a.outdir, exist_ok=True)
    names, levels = r["names"], ["L1", "L2", "L3"]
    plt.rcParams.update({"font.size": 10, "figure.dpi": 110})

    # importance
    fig, ax = plt.subplots(figsize=(8, 3.6))
    w = 0.27
    xs = np.arange(len(names))
    for i, L in enumerate(levels):
        v = [r["levels"][L]["importance"][n] for n in names]
        ax.bar(xs + (i - 1) * w, v, w, label=L, color=COLORS[L])
    ax.set_xticks(xs, names)
    ax.set_ylabel("accuracy drop when shuffled")
    ax.set_title("Which conditions matter for each life level (permutation importance)")
    ax.axhline(0, color="k", lw=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "importance.png"))
    plt.close(fig)

    # partial dependence, one figure per level
    for L in levels:
        fig, axs = plt.subplots(2, 5, figsize=(13, 5), sharey=True)
        for ax, n in zip(axs.ravel(), names):
            pd_ = r["levels"][L]["partial_dependence"][n]
            ax.plot(pd_["value"], pd_["p"], color=COLORS[L])
            if n in LOG_AXES:
                ax.set_xscale("log")
            ax.set_xlabel(LABELS[n])
            ax.set_ylim(0, 1)
            ax.grid(alpha=0.3)
        axs[0, 0].set_ylabel(f"P({L})")
        axs[1, 0].set_ylabel(f"P({L})")
        fig.suptitle(f"{L}: partial dependence (averaged over random other conditions)")
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, f"pd_{L}.png"))
        plt.close(fig)

    # 2D slices
    for key in r["levels"]["L1"]["slices"]:
        fig, axs = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
        for ax, L in zip(axs, levels):
            sl = r["levels"][L]["slices"][key]
            X, Y, P = np.array(sl["x_values"]), np.array(sl["y_values"]), np.array(sl["p"])
            im = ax.pcolormesh(X, Y, P, vmin=0, vmax=1, cmap="viridis", shading="auto")
            ax.contour(X, Y, P, levels=[0.5], colors="w", linewidths=1.2)
            if sl["x"] in LOG_AXES:
                ax.set_xscale("log")
            if sl["y"] in LOG_AXES:
                ax.set_yscale("log")
            ax.set_xlabel(LABELS[sl["x"]])
            ax.set_ylabel(LABELS[sl["y"]])
            ax.set_title(L)
        fig.colorbar(im, ax=axs, label="P(level)", shrink=0.85)
        fig.suptitle(f"{sl['x']} × {sl['y']} (others at the 2b reference; white line = 50%)")
        fig.savefig(os.path.join(a.outdir, f"slice_{key}.png"), bbox_inches="tight")
        plt.close(fig)

    # accuracy curves
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
    for L in levels:
        for meth, ls in (("engine4_soft", "-"), ("sobol_1seed", "--")):
            c = np.array(r["levels"][L][meth]["curve"])
            axs[0].plot(c[:, 0], c[:, 1], ls, color=COLORS[L], label=f"{L} {'COMPASS' if ls == '-' else 'Sobol'}")
            axs[1].plot(c[:, 0], c[:, 2], ls, color=COLORS[L])
    for ax, t in zip(axs, ("all test conditions", "hard conditions")):
        ax.set_xscale("log")
        ax.set_xlabel("simulations")
        ax.set_ylabel("accuracy")
        ax.set_title(t)
        ax.grid(alpha=0.3)
    axs[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "accuracy.png"))
    plt.close(fig)
    print("figures:", sorted(os.listdir(a.outdir)))


if __name__ == "__main__":
    main()
