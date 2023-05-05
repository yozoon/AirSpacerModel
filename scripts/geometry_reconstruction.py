""" Example demonstrating the extraction and reconstruction of geometry """
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import viennals2d as vls
import viennaps2d as vps
from trench_features import TrenchFeatureExtractor, TrenchFeatures


def print_surface(levelset: vls.lsDomain, filename: str) -> None:
    mesh = vls.lsMesh()
    vls.lsToSurfaceMesh(levelset, mesh).apply()
    vls.lsVTKWriter(mesh, filename).apply()


def arggroup(
    data,
    tol: float = 1e-4,
    predicate: Callable[[np.ndarray], bool] | None = None,
):
    arg = np.arange(data.shape[0])
    sections = np.where(np.abs(data[1:] - data[:-1]) >= tol)

    if sections:
        groups = np.split(arg, sections[0] + 1)
    else:
        groups = [arg]

    if predicate:
        return [g for g in groups if predicate(g)]

    return groups


def features_to_stamp(
    grid: vls.hrleGrid,
    origin: list[float],
    trench_depth: float,
    trench_top_width: float,
    features: TrenchFeatures,
    min_group_size: int = 3,
) -> vls.lsDomain:
    x = origin[0] + trench_top_width * np.hstack([
        features.right_sidewall_above,
        features.right_sidewall_below,
        features.left_sidewall_below,
        features.left_sidewall_above,
    ])

    y = origin[1] + trench_depth * np.hstack([
        features.vertical_pos_above[::-1],
        features.vertical_pos_below[::-1],
    ])
    y = np.hstack([y, y[::-1]])

    # Take difference between right and left sidewall distances (from top to bottom)
    diffs = x[:len(x) // 2] - x[len(x) // 2:][::-1]

    # Replace NaNs with -1
    diffs[np.isnan(diffs)] = -1

    pinchoff = (diffs <= 0) * 1

    index_groups = arggroup(
        pinchoff,
        predicate=lambda grp: (not pinchoff[grp[0]] > 0 and len(grp) >= min_group_size),
    )

    stamp = vls.lsDomain(grid)

    for i, idx in enumerate(index_groups):
        idx = np.hstack([idx, (len(x) - idx - 1)[::-1]])
        profile: list[list[float]] = [[]]
        if i == 0:
            profile = np.vstack([
                np.hstack([x[idx[0]], x[idx], x[idx[-1]]]),
                np.hstack([y[idx[0]] + 1, y[idx], y[idx[-1]] + 1]),
            ]).T.tolist()
        else:
            profile = np.vstack([x[idx].tolist(), y[idx].tolist()]).T.tolist()

        dom = vps.psDomain()
        vps.psExtrudeProfile(psDomain=dom,
                             gridDelta=grid.getGridDelta(),
                             xExtent=2 * trench_top_width,
                             yExtent=1.0,
                             profile=profile,
                             extrusionLength=1,
                             fullExtent=True).apply()
        dom.printSurface(f"hull_mesh_{i}.vtp")
        hull = dom.getLevelSets()[0]
        vls.lsBooleanOperation(stamp, hull, vls.lsBooleanOperationEnum.UNION).apply()
    return stamp


def main() -> None:
    grid_delta = 0.2
    sticking_probability = 0.5
    trench_top_width = 4
    process_duration = 5
    trench_depth = 10
    taper_angle = 0.0

    geometry = vps.psDomain()
    vps.psMakeTrench(psDomain=geometry,
                     gridDelta=grid_delta,
                     xExtent=2 * trench_top_width,
                     yExtent=1.0,
                     baseHeight=-trench_depth,
                     trenchWidth=trench_top_width,
                     trenchHeight=trench_depth,
                     taperingAngle=taper_angle).apply()

    model = vps.SimpleDeposition(stickingProbability=sticking_probability,
                                 sourceExponent=1.0)

    mesh = vls.lsMesh()
    vls.lsToDiskMesh(geometry.getLevelSets()[0], mesh).apply()
    nodes = np.asarray(mesh.getNodes())

    process = vps.psProcess()
    process.setDomain(geometry)
    process.setProcessModel(model.getProcessModel())
    process.setNumberOfRaysPerPoint(1000)
    process.setProcessDuration(process_duration / sticking_probability)

    geometry.printSurface("initial.vtp")

    process.apply()

    geometry.printSurface("final.vtp")

    origin = [0., 0, 0.]

    extractor = TrenchFeatureExtractor(
        grid_delta=grid_delta,
        trench_top_width=trench_top_width,
        trench_depth=trench_depth,
        origin=origin,
        n_samples_below=256,
        n_samples_above=64,
    )

    features = extractor.extract(geometry)

    if features:
        grid = geometry.getLevelSets()[-1].getGrid()
        reconstructed_geometry = vls.lsDomain(grid)
        stamp = features_to_stamp(
            grid=grid,
            origin=origin,
            trench_depth=trench_depth,
            trench_top_width=trench_top_width,
            features=features,
        )

        print_surface(stamp, "stamp.vtp")

        plane_origin = origin
        plane_origin[1] += process_duration
        normal = [0, 1, 0]
        vls.lsMakeGeometry(reconstructed_geometry, vls.lsPlane(plane_origin,
                                                               normal)).apply()

        vls.lsBooleanOperation(reconstructed_geometry, stamp,
                               vls.lsBooleanOperationEnum.RELATIVE_COMPLEMENT).apply()

        print_surface(reconstructed_geometry, "reconstructed.vtp")

        # x = origin[0] + trench_top_width * np.hstack([
        #     features.right_sidewall_above,
        #     features.right_sidewall_below,
        #     features.left_sidewall_below,
        #     features.left_sidewall_above,
        # ])

        # idx = np.arange(len(x))
        # nan_mask = np.isnan(x)
        # idx = idx[~nan_mask]

        # y = origin[1] - trench_depth + trench_depth * np.hstack([
        #     features.vertical_pos_above[::-1],
        #     features.vertical_pos_below[::-1],
        # ])
        # y = np.hstack([y, y[::-1]])

        # idx = np.hstack([idx, idx[0]])
        # x = x[idx]
        # y = y[idx]

        # plt.scatter(nodes[:, 0], nodes[:, 1], marker=".", color="k")

        # mesh = vls.lsMesh()
        # vls.lsToDiskMesh(geometry.getLevelSets()[0], mesh).apply()
        # nodes = np.asarray(mesh.getNodes())

        # plt.scatter(nodes[:, 0], nodes[:, 1], marker=".", color="k")

        # # plt.plot(x, y)

        # plt.axis("equal")
        # plt.tight_layout()
        # plt.grid()
        # plt.show()


if __name__ == "__main__":
    main()
