from argparser import parse_args
from svd.NumpySVD import NumpySVD
from compressor.Compressor import Compressor

if __name__ == '__main__':
    args = parse_args()
    compressor = Compressor()
    if args.mode == "decompress":
        compressor.decompress(args.in_file, args.out_file)
    elif args.mode == "compress":
        if args.method == "numpy":
            svd = NumpySVD()
        else:
            print("Wrong method is given")
            exit(1)
        compressor.compress(args.in_file, args.out_file, svd, args.compression)
    else:
        print("Wrong mode is given.")
        exit(1)
