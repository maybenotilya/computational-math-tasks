import numpy as np
from pathlib import Path
from PIL import Image
from typing import Tuple, Dict
from svd.AbstractSVD import AbstractSVD


class Compressor:
    @staticmethod
    def _load_image_channels(filepath: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        R, G, B = Image.open(filepath).convert('RGB').split()
        return np.array(R), np.array(G), np.array(B)

    @staticmethod
    def _save_image_channels(filepath: Path, channels: Tuple[np.ndarray, np.ndarray, np.ndarray]):
        image = Image.fromarray(np.dstack(channels).astype(np.uint8), mode='RGB')
        image.save(filepath)

    @staticmethod
    def _load_matrices(filepath: Path) -> np.ndarray:
        return np.load(filepath)

    @staticmethod
    def _save_matrices(filepath: Path, **matrices) -> None:
        np.savez(filepath, **matrices)

    @staticmethod
    def _compress_channels(
            channels: Tuple[np.ndarray, np.ndarray, np.ndarray],
            svd: AbstractSVD,
            s: int
    ) -> Dict[str, np.ndarray]:
        R, G, B = channels
        matrices = {}
        matrices['UR'], matrices['SR'], matrices['VR'] = svd(R)
        matrices['UG'], matrices['SG'], matrices['VG'] = svd(G)
        matrices['UB'], matrices['SB'], matrices['VB'] = svd(B)
        for chan in ('R', 'G', 'B'):
            matrices[f'U{chan}'] = matrices[f'U{chan}'][:, :s]
            matrices[f'S{chan}'] = matrices[f'S{chan}'][:s]
            matrices[f'V{chan}'] = matrices[f'V{chan}'][:s, :]
        return matrices

    @staticmethod
    def compress(
            input_path: Path,
            output_path: Path,
            svd: AbstractSVD,
            N: int
    ) -> None:
        R, G, B = Compressor._load_image_channels(input_path)
        matrices = Compressor._compress_channels((R, G, B), svd, 3)
        Compressor._save_matrices(output_path, **matrices)

    @staticmethod
    def decompress(
            input_path: Path,
            output_path: Path,
    ) -> None:
        matrices = Compressor._load_matrices(input_path)
        R = matrices['UR'] @ np.diag(matrices['SR']) @ matrices['VR']
        G = matrices['UG'] @ np.diag(matrices['SG']) @ matrices['VG']
        B = matrices['UB'] @ np.diag(matrices['SB']) @ matrices['VB']
        image = Image.fromarray(np.dstack((R, G, B)).astype(np.uint8), mode='RGB')
        image.save(output_path)
