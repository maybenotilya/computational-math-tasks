from argparser import parse_args
from svd.NumpySVD import NumpySVD
from svd.SimpleSVD import SimpleSVD
from svd.AdvancedSVD import AdvancedSVD
from compressor.Compressor import Compressor

if __name__ == '__main__':
    args = parse_args()
    compressor = Compressor()
    if args.mode == "decompress":
        compressor.decompress(args.in_file, args.out_file)
    elif args.mode == "compress":
        match args.method:
            case "numpy":
                svd = NumpySVD()
            case "simple":
                svd = SimpleSVD()
            case "advanced":
                svd = AdvancedSVD()
            case _:
                print("Wrong SVD method is given.")
                exit(1)

        compressor.compress(args.in_file, args.out_file, svd, args.compression)
    else:
        print("Wrong mode is given.")
        exit(1)
