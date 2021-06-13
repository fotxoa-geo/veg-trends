import geoprocess
import argparse


def run_preprocess(base_directory: str, instrument: str):
    gp = geoprocess.reflectance(base_directory=base_directory, instrument=instrument)
    gp.layer_stack(modis_product='MOD09A1')
    gp.layer_stack(modis_product='MOD11A2')


def main():
    parser = argparse.ArgumentParser(description='Run vegetation trend workflow')
    parser.add_argument('-bd', type=str, help='specify base directory')
    parser.add_argument('-ins', type=str, help="Instrument to Convolve Wavelengths", default='MODIS')
    parser.add_argument("-mode", type=str, help="set the run mode.", choices=['preprocess', 'trends', 'all'],
                        default="trends")
    args = parser.parse_args()

    if args.mode in ['preprocess']:
        run_preprocess(base_directory=args.bd, instrument=args.ins)


if __name__ == '__main__':
    main()
