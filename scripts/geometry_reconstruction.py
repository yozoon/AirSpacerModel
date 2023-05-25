""" Example demonstrating the extraction and reconstruction of geometry """
import logging
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import viennals2d as vls
import viennaps2d as vps
from trench_features import TrenchFeatureExtractor, TrenchFeatures


def print_surface(levelset: vls.lsDomain, filename: str) -> None:
    """ Save a surface mesh of the provided levelset """
    mesh = vls.lsMesh()
    vls.lsToSurfaceMesh(levelset, mesh).apply()
    vls.lsVTKWriter(mesh, filename).apply()


def extrude_profile(
    grid: vls.hrleGrid,
    profile: list[tuple[float, float]] | list[list[float]] | np.ndarray,
) -> vls.lsDomain:
    """ Create a levelset from the provided profile """
    substrate = vls.lsDomain(grid)
    mesh = vls.lsMesh()
    # Insert all points of the profile into the mesh and connect them by lines
    for i, pt in enumerate(profile):
        mesh.insertNextNode([*pt, 0.])
        mesh.insertNextLine([i, (i + 1) % len(profile)])

    # Now create the levelset from the mesh
    vls.lsFromSurfaceMesh(substrate, mesh).apply()

    return substrate


def arggroup(
    data,
    tol: float = 1e-4,
    predicate: Callable[[np.ndarray], bool] | None = None,
):
    """ Split the provided data into sections and return indices of data in these sections

    Splits the data at points where the absolute change between into sections and return
    indices of data in these sections consecutive values is larger than the provided
    tolerance. If predicate is provided, the list of sections is additionally filtered
    using the predicate."""
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
    """ Converts a trench feature object into a stamp (levelset) that can be subtracted
    from another levelset to reproduce the geometry described by the feature object."""
    x = origin[0] + trench_top_width * features.get_data()

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

        hull = extrude_profile(grid, profile)
        vls.lsBooleanOperation(stamp, hull, vls.lsBooleanOperationEnum.UNION).apply()
    return stamp


def main() -> None:
    grid_delta = 0.2
    sticking_probability = 0.5
    trench_top_width = 4
    process_duration = 7
    trench_depth = 8
    taper_angle = 0.0
    x_extent = 2 * trench_top_width
    y_extent = 1.0

    geometry = vps.psDomain()
    vps.psMakeTrench(psDomain=geometry,
                     gridDelta=grid_delta,
                     xExtent=x_extent,
                     yExtent=y_extent,
                     baseHeight=-trench_depth,
                     trenchWidth=trench_top_width,
                     trenchHeight=trench_depth,
                     taperingAngle=taper_angle).apply()

    model = vps.SimpleDeposition(stickingProbability=sticking_probability,
                                 sourceExponent=1.0)

    process = vps.psProcess()
    process.setDomain(geometry)
    process.setProcessModel(model)
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

    features = extractor.extract(levelset=geometry.getLevelSets()[-1])

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
        normal = [0.0, 1.0, 0.0]
        vls.lsMakeGeometry(
            reconstructed_geometry,
            vls.lsPlane(plane_origin, normal),
        ).apply()

        vls.lsBooleanOperation(
            reconstructed_geometry,
            stamp,
            vls.lsBooleanOperationEnum.RELATIVE_COMPLEMENT,
        ).apply()

        print_surface(reconstructed_geometry, "reconstructed.vtp")


if __name__ == "__main__":
    main()
