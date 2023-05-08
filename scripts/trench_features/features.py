from dataclasses import dataclass

import numpy as np


@dataclass(kw_only=True, slots=True)
class TrenchFeatures:
    # Bottom to top
    vertical_pos_below: np.ndarray
    vertical_pos_above: np.ndarray

    # The sidewall features are stored in clockwise order, starting in the top right
    # Top to bottom
    right_sidewall_above: np.ndarray
    right_sidewall_below: np.ndarray

    # Bottom to top
    left_sidewall_below: np.ndarray
    left_sidewall_above: np.ndarray

    def get_metadata(self) -> dict[str, np.ndarray]:
        return {
            "vertical_pos_below": self.vertical_pos_below,
            "vertical_pos_above": self.vertical_pos_above,
        }

    def get_data(self) -> np.ndarray:
        return np.hstack([
            self.right_sidewall_above,
            self.right_sidewall_below,
            self.left_sidewall_below,
            self.left_sidewall_above,
        ])
