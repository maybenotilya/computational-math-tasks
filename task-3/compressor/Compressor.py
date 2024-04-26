import numpy as np
import os
from pathlib import Path
from PIL import Image
from typing import Tuple, Dict
from svd.AbstractSVD import AbstractSVD
from constants import INT_SIZE, FLOAT_SIZE


class Compressor:
    @staticmethod
    def _cals_s(original_size: int, N: float, shape: Tuple[int, ...]):
        n, m = shape
        size = lambda k: 3 * INT_SIZE + 3 * FLOAT_SIZE * (n * k + k + k * m)
        l, r = 0, min(n, m)
        while r - l > 1:
            mid = (l + r) // 2
            if size(mid) <= int(original_size / N):
                l = mid
            else:
                r = mid
        return l

    @staticmethod
    def _load_image_channels(filepath: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        R, G, B = Image.open(filepath).convert('RGB').split()
        return np.array(R), np.array(G), np.array(B)

    @staticmethod
    def _save_image_channels(filepath: Path, channels: Tuple[np.ndarray, np.ndarray, np.ndarray]):
        image = Image.fromarray(np.dstack(channels).astype(np.uint8), mode='RGB')
        image.save(filepath)

    @staticmethod
    def _load_matrices(filepath: Path) -> Dict[str, np.ndarray]:
        with open(filepath, 'rb') as f:
            data = bytearray(f.read())

        n = int(np.frombuffer(data, offset=0, dtype=np.uint32, count=1)[0])
        k = int(np.frombuffer(data, offset=4, dtype=np.uint32, count=1)[0])
        m = int(np.frombuffer(data, offset=8, dtype=np.uint32, count=1)[0])
        offset = 12
        matrices = dict()
        for chan in ('R', 'G', 'B'):
            matrices[f'U{chan}'] = np.frombuffer(data, offset=offset, dtype=np.float64, count=n * k).reshape((n, k))
            offset += n * k * FLOAT_SIZE

            matrices[f'S{chan}'] = np.frombuffer(data, offset=offset, dtype=np.float64, count=k)
            offset += k * FLOAT_SIZE

            matrices[f'V{chan}'] = np.frombuffer(data, offset=offset, dtype=np.float64, count=k * m).reshape((k, m))
            offset += k * m * FLOAT_SIZE
        return matrices

    @staticmethod
    def _save_matrices(filepath: Path, **matrices: np.ndarray) -> None:
        data = bytearray()
        data.extend(np.uint32(matrices['UR'].shape[0]).tobytes())
        data.extend(np.uint32(matrices['SR'].shape[0]).tobytes())
        data.extend(np.uint32(matrices['VR'].shape[1]).tobytes())

        for chan in ('R', 'G', 'B'):
            for matrix in ('U', 'S', 'V'):
                data.extend(matrices[f'{matrix}{chan}'].flatten().tobytes())

        with open(filepath, 'wb') as f:
            f.write(data)

    @staticmethod
    def _compress_channels(
            channels: Tuple[np.ndarray, np.ndarray, np.ndarray],
            svd: AbstractSVD,
            s: int
    ) -> Dict[str, np.ndarray]:
        R, G, B = channels
        matrices = dict()
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
            N: float,
    ) -> None:
        R, G, B = Compressor._load_image_channels(input_path)
        s = Compressor._cals_s(os.path.getsize(input_path), N, R.shape)
        print(f'Using matrices of rank of {s}')
        matrices = Compressor._compress_channels((R, G, B), svd, s)
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
