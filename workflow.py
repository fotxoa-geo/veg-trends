import geoprocess


def run_clip():
    gp = geoprocess.reflectance(base_directory=r'D:\veg-trends\\', instrument='MODIS')
    print(gp.layer_name(input_hdf=r'D:\veg-trends\data\MOD09A1\MOD09A1.A2000049.h08v05.006.2015136143537.hdf'))
    print(gp.layer_stack(modis_product='MOD09A1'))

run_clip()