from argparser import parse_args
from svd.NumpySVD import NumpySVD
from compressor.Compressor import Compressor
from pathlib import Path


if __name__ == '__main__':
    svd = NumpySVD()
    compressor = Compressor()
    compressor.compress(Path('basov.jpeg'), Path('temp.penis'), svd, 1)
    compressor.decompress(Path('temp.penis'), Path('babasov.png'))
