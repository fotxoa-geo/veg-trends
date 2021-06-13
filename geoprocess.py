import os
from glob import glob
import datetime
from collections import OrderedDict

import numpy as np
import pandas as pd
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
        self.output_directory = os.path.join(base_directory, "output")

    def layer_name(self, input_hdf: str):
        ds = gdal.Open(input_hdf)
        subdatasets = ds.GetSubDatasets()
        file_templates = []
        band_names = []
        for fname, name in subdatasets:
            file_templates.append(fname)
            band_names.append(name)

        return file_templates, band_names
    
    def layer_stack(self, modis_product: str):
        hdfs = glob(os.path.join(self.base_directory, "data", modis_product, "*.hdf"))

        dates = []

        for hdf in hdfs:
            file_date = os.path.basename(hdf).split(".")[1]
            mod_grid = os.path.basename(hdf).split(".")[2]
            date = file_date[1:8]
            dates.append(date)

        # Get band info
        num_bands = self.layer_name(input_hdf=hdfs[0])

        # Remove duplicates of dates: Keeps the order
        dates = list(OrderedDict.fromkeys(dates))

        for date in tqdm(dates, ncols=80, desc='Layer Stacking', colour='white'):
            # Get list of HDFS that satisfy the date
            hdfs_date = glob(os.path.join(self.base_directory, "data", modis_product, '*A' + date + '*.hdf'))

            # Create a dataframe to store selected bands with same date
            df = pd.DataFrame(columns=num_bands[1])

            for hdf in hdfs_date:
                file_date = os.path.basename(hdf).split(".")[1]
                mod_grid = os.path.basename(hdf).split(".")[2]
                date = file_date[1:8]
                hdf_subbset = self.layer_name(input_hdf=hdf)[0]
                df.loc[len(df)] = hdf_subbset

            # Save each band as GCS vrt
            for col in df.columns:
                gdal.BuildVRT(os.path.join(self.output_directory, "scratch", col.split(" ")[1] + ".vrt"),
                              list(df[col].values))

                # Convert the VRT into a GCS WGS 1984 VRT
                gdal.Warp(os.path.join(self.output_directory, "scratch", "gcs_" + col.split(" ")[1] + ".vrt"),
                          os.path.join(self.output_directory, "scratch", col.split(" ")[1] + ".vrt"),
                          dstSRS='EPSG:4326')

            # Layer stack the extracted bands
            gcs_vrts = glob(os.path.join(self.output_directory, "scratch", 'gcs_*.vrt'))
            gdal.BuildVRT(os.path.join(self.output_directory, "scratch", 'stacked.vrt'), gcs_vrts, separate=True)

            # Save the stack to the final destination
            outras = os.path.join(self.base_directory, "data", "jornada", modis_product, 'gcs_sur_refl_' + date + '.tif')
            outtif = gdal.Warp(outras, os.path.join(self.output_directory, "scratch", 'stacked.vrt'),
                               cutlineDSName=self.shapefile,
                               cropToCutline=True, dstNodata=0)
            outtif = None









            
