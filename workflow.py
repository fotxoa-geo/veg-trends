import geoprocess
import argparse
from osgeo import _gdalconst


def run_preprocess(base_directory: str, instrument: str):
    gp = geoprocess.modis_process(base_directory=base_directory, instrument=instrument)

    # Process reflectance and temperature data
    # gp.layer_stack(modis_product='MOD09A1', search_kw='sur_refl_b', output_kw='refl', output_type=_gdalconst.GDT_Int16,
    #                no_data=-9999)
    # gp.layer_stack(modis_product='MOD11A2', search_kw='LST_', output_kw='LST', output_type=_gdalconst.GDT_UInt16,
    #                no_data=0)
    # gp.layer_stack(modis_product='MOD09A1', search_kw='sur_refl_qc_500m', output_kw='qc',
    #                output_type=_gdalconst.GDT_UInt32, no_data=9999)
    # gp.layer_stack(modis_product='MOD11A2', search_kw='QC_', output_kw='qc', output_type=_gdalconst.GDT_Byte,
    #                no_data=255)
    gp.temporal_stack(search_kw='LST_', modis_product='MOD11A2', data_type='8-bit')


def main():
    parser = argparse.ArgumentParser(description='Run vegetation trend workflow')
    parser.add_argument('-bd', type=str, help='Specify base directory')
    parser.add_argument('-ins', type=str, help="Instrument to Convolve Wavelengths", default='MODIS')
    parser.add_argument("-mode", type=str, help="set the run mode.", choices=['preprocess', 'trends', 'all'],
                        default="trends")
    args = parser.parse_args()

    if args.mode in ['preprocess']:
        run_preprocess(base_directory=args.bd, instrument=args.ins)


if __name__ == '__main__':
    main()
