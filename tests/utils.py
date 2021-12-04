from pathlib import Path
import numpy as np


def np_load(file: Path) -> np.ndarray:
    return np.load(file)
