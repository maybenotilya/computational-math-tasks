import numpy as np
from svd.AbstractSVD import AbstractSVD
from typing import Tuple


class SimpleSVD(AbstractSVD):
    def __call__(self, matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        ...
