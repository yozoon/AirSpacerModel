import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

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
    input_shape = tuple(
        np.unique(data[:, i]).shape[0] for i in range(data.shape[1]-1))
    print(input_shape)
    data = data[np.lexsort(data[:, :-1:][:, ::-1].T)]

    time = data[:, 1].reshape(input_shape)
    taper_angle = data[:, 0].reshape(input_shape)
    filling_ratio = data[:, 2].reshape(input_shape)
    n_levels = 5
    levels = list(np.linspace(1/n_levels, 1, n_levels))
    levels = sorted(levels)

    fig, ax = plt.subplots(1, 1, squeeze=True)
    ax.imshow(filling_ratio, aspect="auto", interpolation="bicubic", extent=(
        np.min(time), np.max(time), np.min(taper_angle), np.max(taper_angle)),
        origin="lower", cmap="viridis", vmin=0, vmax=1.6)
    CS = ax.contour(time, taper_angle, filling_ratio,
                    levels=levels,
                    # colors=["white"]*(n_levels-1) + ["#E7298A"],
                    colors="white",
                    )

    # plt.colorbar()
    ax.set_title("Filling Ratio")
    ax.set_xlabel("normalized time")
    ax.set_ylabel("taper angle [°]")
    mx, my = np.max(time), np.max(taper_angle)

    # for i in range(2):
    # ax.clabel(CS, CS.levels, inline=False, fmt=fmt, fontsize=10)
    clabels = ax.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=10)
    # [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=1)) for txt in clabels]
    plt.tight_layout()
    # ax.scatter(0.8, 5.5, marker="X", s=200, color="magenta")
    # ax.arrow(x=1.0, y=0.0, dx=0, dy=5.5, width=0.01*mx,
    #          head_length=0.03*my, color="purple", length_includes_head=True,
    #          clip_on=True, clip_box=ax.get_tightbbox(),)
    # ax.arrow(x=1.0, y=5.5, dx=-0.2, dy=0, width=0.01*my,
    #          head_length=0.03*mx, color="purple", length_includes_head=True,
    #          clip_on=True, clip_box=ax.get_tightbbox())
    plt.savefig("filling_ratio_slice.pdf")


if __name__ == "__main__":
    main()
