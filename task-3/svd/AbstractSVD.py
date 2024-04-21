from abc import ABC, abstractmethod
import numpy as np
from typing import Tuple

# A = U @ S @ V
class AbstractSVD(ABC):
    @abstractmethod
    def __call__(self, matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        ...
