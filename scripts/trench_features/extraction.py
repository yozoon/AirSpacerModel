from typing import Optional

import numpy as np
import viennals2d as vls
from scipy.spatial import KDTree

from .features import TrenchFeatures


class TrenchFeatureExtractor:

    def __init__(
        self,
        grid_delta: float,
        trench_depth: float,
        trench_top_width: float,
        origin: tuple[float, float, float] | list[float] = (0.0, 0.0, 0.0),
        n_samples_below: int = 256,
        n_samples_above: int = 128,
    ):
        self.grid_delta = grid_delta
        self.trench_depth = trench_depth
        self.trench_top_width = trench_top_width
        self.origin = origin
        self.n_samples_below = n_samples_below
        self.n_samples_above = n_samples_above

        self.sample_locations_below = np.linspace(self.origin[1] - self.trench_depth,
                                                  self.origin[1], self.n_samples_below)
        self.sample_locations_above = np.linspace(
            self.origin[1], self.origin[1] + 2 * self.trench_top_width,
            self.n_samples_above)

    def _extract_widths(
        self,
        tree: KDTree,
        sample_locations: np.ndarray,
        nodes: np.ndarray,
    ) -> np.ndarray:
        query_results = tree.query_ball_point(
            x=sample_locations.reshape(-1, 1),
            r=4 * self.grid_delta,
            p=1,  # Manhattan Distance (faster than euclidean, but equivalent in 1D)
        )

        features = np.full_like(sample_locations, fill_value=np.nan)

        for i, (loc, neighbor_indices) in enumerate(zip(sample_locations, query_results)):
            neighbors = nodes[neighbor_indices]

            lower_neighbors = neighbors[neighbors[:, 1] < loc][::-1]
            upper_neighbors = neighbors[neighbors[:, 1] >= loc]

            # print(f"======== {loc:.3f} ========")
            # print(f"n_neigbors={len(neighbor_indices)}")
            # print(f"lower={loc - lower_neighbors[:, 1]}")
            # print(f"upper={upper_neighbors[:, 1] - loc}")

            upper_distance = 0
            upper_width = 0
            if len(upper_neighbors) > 0:
                upper_distance = upper_neighbors[0, 1] - loc
                upper_width = upper_neighbors[0, 0] - self.origin[0]

            lower_distance = 0
            lower_width = 0
            if len(lower_neighbors) > 0:
                lower_distance = loc - lower_neighbors[0, 1]
                lower_width = lower_neighbors[0, 0] - self.origin[0]

            if len(upper_neighbors) > 0 and upper_distance < self.grid_delta / 4:
                # If the vertical position of the upper point coincides with the sample
                # position up to a certain epsilon, use its width.
                features[i] = upper_width / self.trench_top_width
            elif len(lower_neighbors) > 0 and lower_distance < self.grid_delta / 4:
                # If the vertical position of the lower point coincides with the sample
                # position up to a certain epsilon, use its width.
                features[i] = lower_width / self.trench_top_width
            elif len(upper_neighbors) > 0 and len(lower_neighbors) > 0:
                # Otherwise linearly interpolate between the two widths based on the
                # offset to to sample position.
                total_distance = lower_distance + upper_distance
                features[i] = (
                    lower_distance * upper_width + upper_distance * lower_width) / \
                    total_distance / self.trench_top_width
        return features

    def _extract_features(
        self,
        nodes: np.ndarray,
        sample_locations: tuple[np.ndarray, ...],
    ) -> tuple[np.ndarray, ...]:
        # Instantiate the KDTree for neighbors lookup
        tree = KDTree(np.c_[nodes[:, 1].ravel()])

        # Extract the widths for each of the provided sample_location arrays
        return tuple(
            self._extract_widths(tree=tree, nodes=nodes, sample_locations=sl)
            for sl in sample_locations)

    def extract(
        self,
        levelset: vls.lsDomain,
        normals_threshold: float = 0.03,
    ) -> Optional[TrenchFeatures]:
        mesh = vls.lsMesh()
        vls.lsToDiskMesh(levelset, mesh).apply()
        nodes = np.asarray(mesh.getNodes())
        normals = np.asarray(mesh.getCellData().getVectorData(
            vls.lsCalculateNormalVectors.normalVectorsLabel))

        # Extract the nodes on the right sidewall and sort them according to
        # their vertical position
        right_mask = normals[:, 0] < -normals_threshold
        right_nodes = nodes[right_mask]
        sorted(right_nodes, key=lambda x: x[1])

        # Extract the nodes on the left sidewall and sort them according to
        # their vertical position
        left_mask = normals[:, 0] > normals_threshold
        left_nodes = nodes[left_mask]
        sorted(left_nodes, key=lambda x: x[1])

        # print(sample_locations_above)
        # print(sample_locations_below)

        right_sidewall_below, right_sidewall_above = self._extract_features(
            right_nodes,
            sample_locations=(self.sample_locations_below[::-1],
                              self.sample_locations_above[::-1]),
        )
        left_sidewall_below, left_sidewall_above = self._extract_features(
            left_nodes,
            sample_locations=(self.sample_locations_below, self.sample_locations_above),
        )

        return TrenchFeatures(
            vertical_pos_below=(self.sample_locations_below - self.origin[1]) /
            self.trench_depth,
            vertical_pos_above=(self.sample_locations_above - self.origin[1]) /
            self.trench_depth,
            right_sidewall_below=right_sidewall_below,
            right_sidewall_above=right_sidewall_above,
            left_sidewall_below=left_sidewall_below,
            left_sidewall_above=left_sidewall_above,
        )
