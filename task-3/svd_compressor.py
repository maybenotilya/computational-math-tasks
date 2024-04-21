import numpy as np

from argparser import parse_args
from svd.NumpySVD import NumpySVD
from compressor.Compressor import Compressor
from pathlib import Path

from PIL import Image

if __name__ == '__main__':
    svd = NumpySVD()
    compressor = Compressor()
    compressor.compress(Path('basov.jpeg'), Path('temp'), svd, 10)
    compressor.decompress(Path('temp.npz'), Path('babasov.png'))
