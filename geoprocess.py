import os
from glob import glob
import datetime
from collections import OrderedDict
from osgeo import gdal
from tqdm import tqdm


class reflectance:
    def __init__(self, base_directory: str, instrument: str):
        """ This class will process HDF files to
        :param base_directory: Data to base folder
        :param instrument: instrment being used in analysis (default = MODIS)

        These are files typically used by this class object
         shapefile: A .shp file of the AOI  boundary (CS has to be GCS = WGS84). Data will be clipped to here
         """

        self.base_directory = base_directory
        self.instrument = instrument
        self.shapefile = os.path.join(base_directory, "data", "gis", "Jornada_GCS.shp")

    def layer_name(self, input_hdf: str):
        ds = gdal.Open(input_hdf)
        if ds is None:
            print("Problem opening file %s!" % input_hdf)
        else:
            print("File %s opened fine" % input_hdf)
        
        subdatasets = ds.GetSubDatasets()
        for fname, name in subdatasets:
            print(name)
            print("\t", fname)
    
    def layer_stack(self, modis_product: str):
        hdfs = glob(os.path.join(self.base_directory, "data", modis_product, "*.hdf"))

        dates = []

        for hdf in hdfs:
            file_date = os.path.basename(hdf).split(".")[1]
            mod_grid = os.path.basename(hdf).split(".")[2]
            date = file_date[1:8]
            dates.append(date)

        # Remove duplicates of dates: Keeps the order
        dates = list(OrderedDict.fromkeys(dates))

        for date in tqdm(dates,ncols=80, total=len(dates), desc='Layer Stacking', colour='white'):
            # Get list of HDFS that satisfy the date
            hdfs_date = glob(os.path.join(self.base_directory, "data", modis_product, '*A' + date + '*.hdf'))

            
