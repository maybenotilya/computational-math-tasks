import argparse
from pathlib import Path
from PIL import Image


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--in_file',
        required=True,
        type=Path,
        help='Input file path'
    )
    parser.add_argument(
        '--out_file',
        required=False,
        type=Path,
        help='Output file path'
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if args.out_file is None:
        args.out_file = args.in_file.with_suffix('.BMP')
    image = Image.open(args.in_file).convert('RGB')
    image.save(args.out_file)
