import argparse
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(
        prog='SVD compressor',
        description='Utility to compress BMP images using singular value decomposition'
    )

    parser.add_argument(
        '--mode',
        required=True,
        type=str,
        choices=['compress', 'decompress'],
        help='SVD compressor processing mode'
    )
    parser.add_argument(
        '--method',
        required=False,
        type=str,
        choices=['numpy', 'simple', 'advanced'],
        help='Method of compression. Only used in compress mode'
    )
    parser.add_argument(
        '--compression',
        required=False,
        type=float,
        help='Compression ratio. Only used in compress mode'
    )
    parser.add_argument(
        '--in_file',
        required=True,
        type=Path,
        help='Input file path'
    )
    parser.add_argument(
        '--out_file',
        required=True,
        type=Path,
        help='Output file path'
    )

    return parser.parse_args()
