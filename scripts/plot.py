import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from scipy.interpolate import PPoly, splrep

matplotlib.use("pgf")
matplotlib.rcParams.update({
    "lines.linewidth": 1,
    "pgf.texsystem": "pdflatex",
    "font.family": "serif",
    "text.usetex": True,
    "pgf.rcfonts": False,
})


def main():
    def fmt(x):
        return f"${x:.2f}$"

    data = np.loadtxt("slice.csv", delimiter=",")
    input_shape = [
        np.unique(data[:, i]).shape[0] for i in range(data.shape[1]-1)]
    print(input_shape)
    data = data[np.lexsort(data[:, :-1:][:, ::-1].T)]

    time = data[:, 1].reshape(input_shape)
    taper_angle = data[:, 0].reshape(input_shape)
    filling_ratio = data[:, 2].reshape(input_shape)

    fig, ax = plt.subplots(1, 1, squeeze=True)
    ax.imshow(filling_ratio, aspect="auto", interpolation="bicubic", extent=(
        np.min(time), np.max(time), np.min(taper_angle), np.max(taper_angle)),
        origin="lower", cmap="viridis", vmin=0, vmax=1.6)

    # plt.colorbar()
    ax.set_title("{\\large Filling Ratio}\n{\\footnotesize AR=$8$, s=$0.008$}")
    ax.set_xlabel("Normalized time")
    ax.set_ylabel("Taper angle [°]")

    t_angle, t_fr = data[data[:, 1] == 1.0][:, [0, 2]].T
    sorted_index = np.argsort(t_angle)
    t_angle = t_angle[sorted_index]
    t_fr = t_fr[sorted_index]

    tck = splrep(t_angle, t_fr-1.0, s=0)
    ppoly = PPoly.from_spline(tck)
    root = ppoly.roots(extrapolate=False)
    if root:
        root = root[0]
        # The arrow for the first search direction
        # ax.arrow(x=1.0, y=0.0, dx=0, dy=root, width=0.01*mx,
        #          head_length=0.03*my, color="orange", length_includes_head=True,
        #          zorder=1000)
        # horizontal line for showing the maximum angle
        ax.axhspan(root, root+2, linestyle="dashed",
                   color="lightgrey", alpha=0.5)
        # ax.axhline(root + 2, linestyle="dashed", color="lightgrey")
        num_below = np.count_nonzero(t_angle < root + 2)

        # The contour for the isoline following
        clines = ax.contour(time[:num_below, :],
                            taper_angle[:num_below, :],
                            filling_ratio[:num_below, :],
                            levels=[1.0],
                            colors="orange",
                            )
        plt.setp(clines.collections, linewidth=3)
        plt.setp(clines.collections, zorder=10000)
        t_min = np.asarray(clines.collections[0].get_paths()[
                           0].vertices)[:, 0].min()
        # The final target mark
        ax.scatter(t_min, root+2, marker="*", s=100,
                   color="orange", zorder=1000)

    n_levels = 5
    levels = list(np.linspace(1/n_levels, 1, n_levels))
    levels = sorted(levels)

    CS = ax.contour(time, taper_angle, filling_ratio,
                    levels=levels,
                    colors="white",
                    )
    clabels = ax.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=10)
    # [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=1)) for txt in clabels]
    plt.tight_layout()
    plt.savefig("slice.pdf")
    # plt.show()


if __name__ == "__main__":
    main()
