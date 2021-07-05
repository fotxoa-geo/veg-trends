import os
from glob import glob
import datetime
from collections import OrderedDict
import shutil
import numpy as np
import pandas as pd
from osgeo import gdal, osr
from tqdm import tqdm
from spectral.io import envi



class modis_process:
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
    
    def layer_stack(self, modis_product: str, search_kw: str, output_kw: str, output_type, no_data: int):
        # Check to see if scratch file exists;
        scratch_directory = os.path.join(self.base_directory, "output", "scratch")

        # Clear data from previous runs
        if os.path.isdir(scratch_directory):
            shutil.rmtree(scratch_directory, ignore_errors=True)
            os.mkdir(scratch_directory)
        else:
            os.mkdir(scratch_directory)

        hdfs = glob(os.path.join(self.base_directory, "data", modis_product, "*.hdf"))

        dates = []

        for hdf in hdfs:
            file_date = os.path.basename(hdf).split(".")[1]
            mod_grid = os.path.basename(hdf).split(".")[2]
            date = file_date[1:8]
            dates.append(date)

        # Get band info
        num_bands = self.layer_name(input_hdf=hdfs[0])
        num_bands = [x for x in num_bands[1] if search_kw in x]

        # Remove duplicates of dates: Keeps the order
        dates = list(OrderedDict.fromkeys(dates))

        for date in tqdm(dates, ncols=80, desc=os.path.join('Layer Stacking ', modis_product), colour='white'):
            # Get list of HDFS that satisfy the date
            hdfs_date = glob(os.path.join(self.base_directory, "data", modis_product, '*A' + date + '*.hdf'))

            # Create a dataframe to store selected bands with same date
            df = pd.DataFrame(columns=num_bands)

            for hdf in hdfs_date:
                file_date = os.path.basename(hdf).split(".")[1]
                mod_grid = os.path.basename(hdf).split(".")[2]
                date = file_date[1:8]
                hdf_subbset = self.layer_name(input_hdf=hdf)[0]
                matching = [x for x in hdf_subbset if search_kw in x]
                df.loc[len(df)] = matching

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
            outras = os.path.join(self.base_directory, "data", "jornada", modis_product, output_kw + '_' + date + '.tif')
            outtif = gdal.Warp(outras, os.path.join(self.output_directory, "scratch", 'stacked.vrt'), dstNodata=no_data,
                               cutlineDSName=self.shapefile, cropToCutline=True, outputType=output_type)
            outtif = None

    def temporal_stack(self, search_kw: str, modis_product: str, data_type: str):
        tifs = glob(os.path.join(self.base_directory, 'data', 'jornada', modis_product, search_kw + '*.tif'))

        # Bits lengths; usually taken from the MODIS product page
        type_dict = {'8-bit': '{0:08b}', '16-bit': '{0:016b}', '32-bit': '{0:032b}'}

        # Read first tiff for col/row infomration
        init_tif = gdal.Open(tifs[0], gdal.GA_ReadOnly)
        init_tif_array = init_tif.ReadAsArray().transpose((1, 2, 0))

        # Create the temporary stack - day and nighttime temperature
        temp_stack_day = np.zeros((init_tif_array.shape[0], init_tif_array.shape[1], len(tifs)))
        temp_stack_night = np.zeros((init_tif_array.shape[0], init_tif_array.shape[1], len(tifs)))

        # Function to convert quality bits to strings
        vfunc = np.vectorize(lambda i: type_dict[data_type].format(i))

        # Data to write tif files
        originX, pixelWidth, b, originY, d, pixelHeight = init_tif.GetGeoTransform()
        driver = gdal.GetDriverByName('GTiff')

        for _tif, tif in tqdm(enumerate(tifs), ncols=80, total=len(tifs), desc='Temporal Stack: ' + modis_product):
            date = os.path.basename(tif).split(".")[0].split("_")[1]

            # Open LST - 1st band is day, 2nd band is nighttime
            ds_lst = gdal.Open(tif, gdal.GA_ReadOnly)
            lst_array = ds_lst.ReadAsArray().transpose((1, 2, 0))

            # Open qc array
            qc_path = os.path.join(self.base_directory, 'data', 'jornada', modis_product, 'qc_' + date + '.tif')
            ds_qc = gdal.Open(qc_path, gdal.GA_ReadOnly)
            qc_array = ds_qc.ReadAsArray().transpose((1, 2, 0))
            qc_bits = vfunc(qc_array)

            # Loop through pixels and perform quality check
            for _row, row in enumerate(lst_array):
                for _col, col in enumerate(row):
                    for _band, band in enumerate(col):

                        lst_val = lst_array[_row, _col, _band]
                        qc_bit = qc_bits[_row, _col, _band]

                        if 7500 <= lst_val <= 65535:
                            if qc_bit[:2] == '00':
                                if _band == 0:
                                    temp_stack_day[_row, _col, _tif] = (lst_val * 0.02) - 273.15
                                else:
                                    temp_stack_night[_row, _col, _tif] = (lst_val * 0.02) - 273.15
                            else:
                                if _band == 0:
                                    temp_stack_day[_row, _col, _tif] = 0
                                else:
                                    temp_stack_night[_row, _col, _tif] = 0

                        else:
                            if _band == 0:
                                temp_stack_day[_row, _col, _tif] = 0
                            else:
                                temp_stack_night[_row, _col, _tif] = 0

        # Save rasters
        day_raster = (os.path.join(self.output_directory, 'jornada_day_temp' + ".tif"))
        outRaster = driver.Create(day_raster, init_tif_array.shape[1], init_tif_array.shape[0], len(tifs),
                                  gdal.GDT_UInt16)
        outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))

        # Loop over all bands.
        for b in range(len(tifs)):
            outband = outRaster.GetRasterBand(b + 1)
            outband.WriteArray(temp_stack_day[:, :, b])
            outband.SetNoDataValue(0)

        # setteing srs from input tif file.
        prj = init_tif.GetProjection()
        outRasterSRS = osr.SpatialReference(wkt=prj)
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()

        night_raster = (os.path.join(self.output_directory, 'jornada_night_temp' + ".tif"))
        outRaster = driver.Create(night_raster, init_tif_array.shape[1], init_tif_array.shape[0], len(tifs),
                                  gdal.GDT_UInt16)
        outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))

        # Loop over all bands.
        for b in range(len(tifs)):
            outband = outRaster.GetRasterBand(b + 1)
            outband.WriteArray(temp_stack_night[:, :, b])
            outband.SetNoDataValue(0)

        # setteing srs from input tif file.
        prj = init_tif.GetProjection()
        outRasterSRS = osr.SpatialReference(wkt=prj)
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()


        # # Metadata for envi files
        # meta = {
        #     'lines': lst_array.shape[0],
        #     'samples': lst_array.shape[1],
        #     'bands': len(tifs),
        #     'interleave': 'bil',
        #     'data type': 4,
        #     'file_type': 'ENVI Standard',
        #     'byte order': 0,
        #     'header offset': 0,
        # }
        #
        # # Save daytime temp to .bil format
        # outDataset_lst = envi.create_image(os.path.join(self.output_directory, 'jornada_day_temp' + ".hdr"),
        #                                    meta, ext='', force=True)
        # mm = outDataset_lst.open_memmap(interleave='bip', writable=True)
        # mm[...] = temp_stack_day
        # del mm
        #
        # # Save nighttime temp to .bil format
        # outDataset_lst = envi.create_image(os.path.join(self.output_directory, 'jornada_night_temp' + ".hdr"),
        #                                    meta, ext='', force=True)
        # mm = outDataset_lst.open_memmap(interleave='bip', writable=True)
        # mm[...] = temp_stack_night
        # del mm










