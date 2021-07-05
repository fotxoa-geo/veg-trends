import numpy as np
from glob import glob
import os
from osgeo import gdal, osr
import matplotlib.pyplot as plt
import pymannkendall as mk
from scipy.stats import mannwhitneyu
from datetime import datetime
import pandas as pd
import matplotlib.gridspec as gridspec
from tqdm import tqdm


class trends:
    def __init__(self, base_directory: str, instrument: str):

        self.base_directory = base_directory
        self.instrument = instrument
        self.shapefile = os.path.join(base_directory, "data", "gis", "Jornada_GCS.shp")
        self.output_directory = os.path.join(base_directory, "output")
        self.daytime_temp = os.path.join(self.output_directory, 'jornada_day_temp' + ".tif")
        self.nighttime_temp = os.path.join(self.output_directory, 'jornada_night_temp' + ".tif")
        self.modis_dates = glob(os.path.join(base_directory, "data", "jornada", "MOD11A2", 'LST_*.tif'))
        self.output_figure_directory = os.path.join(base_directory, "output", 'figures')

    def climate_anomalies(self, year_start: int, year_end: int):
        # Read climate data
        ds_day = gdal.Open(self.daytime_temp, gdal.GA_ReadOnly)
        ds_night = gdal.Open(self.nighttime_temp, gdal.GA_ReadOnly)
        day_array = ds_day.ReadAsArray().transpose((1, 2, 0))
        night_array = ds_night.ReadAsArray().transpose((1, 2, 0))
        arrays = [day_array, night_array]

        # Data to write tif files
        originX, pixelWidth, b, originY, d, pixelHeight = ds_day.GetGeoTransform()
        driver = gdal.GetDriverByName('GTiff')

        # Temperature anomalies - by season (DJF, MAM, JJA, SON)
        length_analysis = (year_end - year_start) + 1
        for _array, array in enumerate(arrays):
            # Seasonal anomalies
            temp_anom_djf = np.zeros((day_array.shape[0], day_array.shape[1], length_analysis))
            temp_anom_mam = np.zeros((day_array.shape[0], day_array.shape[1], length_analysis))
            temp_anom_jja = np.zeros((day_array.shape[0], day_array.shape[1], length_analysis))
            temp_anom_son = np.zeros((day_array.shape[0], day_array.shape[1], length_analysis))

            anom_arrays = [temp_anom_djf, temp_anom_mam, temp_anom_jja, temp_anom_son]

            # Mann_kendall
            mk_original = np.zeros((day_array.shape[0], day_array.shape[1], 4))
            mk_modified = np.zeros((day_array.shape[0], day_array.shape[1], 4))

            # Data summary
            data_available_array = np.zeros((day_array.shape[0], day_array.shape[1]))

            # Average Climate
            temp_avg = np.zeros((day_array.shape[0], day_array.shape[1], 4))

            for _row, row in enumerate(array):
                for _col, col in enumerate(row):

                    # list for monthly average of temp
                    dates = []
                    temp = []

                    for _band, band in enumerate(col):
                        # Get date and sort into seasons
                        date = os.path.basename(self.modis_dates[_band]).split(".")[0].split("_")[1]
                        date_ym = datetime.strptime(date, "%Y%j").strftime("%Y%m")
                        temp.append(band)
                        dates.append(date_ym)

                    # Gap fill data
                    df = pd.DataFrame({'Date': dates, 'Temp': temp})
                    df = df.replace(0, np.nan)
                    df = df.dropna()

                    # If there is no data available, it is considered no data
                    if df.empty:
                        data_available_array[_row, _col] = -9999
                        temp_avg[_row, _col, :] = -9999
                        temp_anom_djf[_row, _col, :] = -9999
                        temp_anom_son[_row, _col, :] = -9999
                        temp_anom_jja[_row, _col, :] = -9999
                        temp_anom_mam[_row, _col, :] = -9999
                        mk_original[_row, _col, :] = -9999
                        mk_modified[_row, _col, :] = -9999
                        data_available_array[_row, _col] = -9999

                    # If data is not empty  Gapfill or drop?
                    else:
                        # Calculate data available as a % of 979 observations
                        data_available_array[_row, _col] = np.round(len(df)/979, 2)

                        # Resample data to monthly temperature averages, there should be no nan data at this point and
                        # monthly dataframes that we resample to should not be empty
                        df['Date'] = pd.to_datetime(df['Date'], format='%Y%m')
                        df_month = df.groupby(by='Date', as_index=False).mean()
                        df_month['year'] = df_month.Date.dt.year
                        df_month['month'] = df_month.Date.dt.month
                        year_min = df_month.year.min()
                        year_max = df_month.year.max()

                        # Season data
                        djf = []
                        mam = []
                        jja = []
                        son = []

                        season_dict = {'djf': (djf, (12, 1, 2)),
                                       'mam': (mam, (3, 4, 5)),
                                       'jja': (jja, (6, 7, 8)),
                                       'son': (son, (9, 10, 11))}

                        seasons_key = ['djf', 'mam', 'jja', 'son']

                        # Only use the 20 years of data, else pass
                        if year_max == 2021 and year_min == 2000:
                            for season in seasons_key:
                                for year in range(year_start, year_end + 1):
                                    if season == 'djf':
                                        met_winter = df_month[((df_month['month'] == 12) & (df_month['year'] == year-1)) |
                                                  ((df_month['month'] == 1) & (df_month['year'] == year)) |
                                                  ((df_month['month'] == 2) & (df_month['year'] == year))]
                                        winter_mean = met_winter['Temp'].mean()
                                        djf.append(winter_mean)

                                    else:
                                        met_season = df_month[((df_month['month'] == season_dict[season][1][0]) & (df_month['year'] == year)) |
                                                  ((df_month['month'] == season_dict[season][1][1]) & (df_month['year'] == year)) |
                                                  ((df_month['month'] == season_dict[season][1][2]) & (df_month['year'] == year))]
                                        met_mean = met_season['Temp'].mean()
                                        season_dict[season][0].append(met_mean)

                        else:
                            data_available_array[_row, _col] = -9999
                            temp_avg[_row, _col, :] = -9999
                            temp_anom_djf[_row, _col, :] = -9999
                            temp_anom_son[_row, _col, :] = -9999
                            temp_anom_jja[_row, _col, :] = -9999
                            temp_anom_mam[_row, _col, :] = -9999
                            mk_original[_row, _col, :] = -9999
                            mk_modified[_row, _col, :] = -9999

                        # Calculate the climate information for AOI
                        data_seasons = [djf, mam, jja, son]
                        for _data_season, (data_season, anom_array) in enumerate(zip(data_seasons, anom_arrays)):
                            data = np.array(data_season)
                            data_clean = data[~np.isnan(data)]

                            if len(data_season) == 0 or len(data_clean) ==0:
                                data_available_array[_row, _col] = -9999
                                temp_avg[_row, _col, :] = -9999
                                temp_anom_djf[_row, _col, :] = -9999
                                temp_anom_son[_row, _col, :] = -9999
                                temp_anom_jja[_row, _col, :] = -9999
                                temp_anom_mam[_row, _col, :] = -9999
                                mk_original[_row, _col, :] = -9999
                                mk_modified[_row, _col, :] = -9999

                            else:
                                # Average temperature maps
                                data = np.array(data_season)
                                data_clean = data[~np.isnan(data)]

                                mean_temp = np.mean(data_clean)
                                temp_avg[_row, _col, _data_season] = mean_temp

                                # Climate Anomalies for every year of analysis and mean
                                anomaly = data - mean_temp
                                anom_array[_row, _col, :] = anomaly

                                # Calculate Mann-Kendal Trends
                                trend, h, p, z, tau, s, var_s, slope, intercept = mk.original_test(data_clean, 0.05)
                                if trend == 'increasing':
                                    mk_original[_row, _col, _data_season] = 1
                                elif trend == 'decreasing':
                                    mk_original[_row, _col, _data_season] = -1
                                else:
                                    mk_original[_row, _col, _data_season] = 2

                                # Calculate modified mann_kendall trends
                                trend, h, p, z, tau, s, var_s, slope, intercept = mk.yue_wang_modification_test(data_clean, 0.05)
                                if trend == 'increasing':
                                    mk_modified[_row, _col, _data_season] = 1
                                elif trend == 'decreasing':
                                    mk_modified[_row, _col, _data_season] = -1
                                else:
                                    mk_modified[_row, _col, _data_season] = 2

            # Plot the data: only data with 4 subplots!
            data_to_plot = [temp_avg, mk_modified, mk_original]
            names = ['avg_temp', 'mk_mod', 'mk_original']

            if _array == 0:
                p_name = 'day'
            else:
                p_name = 'night '
            titles = ['Winter', 'Spring', 'Summer', 'Fall']

            for _data_plot, data_plot in tqdm(enumerate(data_to_plot), ncols=80, total=len(data_to_plot),
                                              desc='Loading Stats...'):
                data_plot[data_plot == 0] = np.nan
                data_plot[data_plot == -9999] = np.nan

                # Create figure
                fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(20, 20))
                axes = axes.flatten()

                for _ax, ax in enumerate(axes):
                    if _data_plot == 0:
                        vmin = np.nanmin(temp_avg)
                        vmax = np.nanmax(temp_avg)
                    else:
                        vmin = -1
                        vmax = 2

                    pcm = ax.imshow(data_plot[:, :, _ax], vmin=vmin, vmax=vmax)
                    ax.set_title(titles[_ax])

                fig.subplots_adjust(right=0.8)
                cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
                fig.colorbar(pcm, cax=cbar_ax)
                plt.savefig(os.path.join(self.output_figure_directory, p_name + '_' + names[_data_plot] +'.png'), dpi=300)
                plt.clf()

            # Plot anomaly data
            years = range(year_start, year_end + 1)
            for _anom, anom in tqdm(enumerate(anom_arrays), ncols=80, total=len(anom_arrays), desc='Loading Anomalies...'):
                anom[anom == 0] = np.nan
                anom[anom == -9999] = np.nan
                fig = plt.figure(figsize=(20, 20))
                ncols = 4
                nrows = 5

                gs = gridspec.GridSpec(ncols=ncols, nrows=nrows, wspace=0, hspace=0, figure=fig)

                for row in range(nrows):
                    for col in range(ncols):
                        ax = fig.add_subplot(gs[row, col])
                        idx = row + col * nrows
                        ax.set_title(str(years[idx]))

                        try:
                            pcm = ax.imshow(anom[:, :, idx], vmin=np.nanmin(anom), vmax=np.nanmax(anom),
                                            cmap='coolwarm')
                        except:
                            pass

                fig.subplots_adjust(right=0.8)
                cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
                fig.colorbar(pcm, cax=cbar_ax)
                plt.savefig(os.path.join(self.output_figure_directory, p_name + "_anomaly_" + seasons_key[_anom]
                                         + '.png'), dpi=300)
                plt.clf()

            # Create figure - data available
            fig, ax = plt.subplots(figsize=(20, 20))
            data_available_array[data_available_array == 0] = np.nan
            data_available_array[data_available_array == -9999] = np.nan

            pcm = ax.imshow(data_available_array[:, :], vmin=0, vmax=1, cmap='coolwarm')
            ax.set_title('% of data available (979 possible samples)')

            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
            fig.colorbar(pcm, cax=cbar_ax)
            plt.savefig(os.path.join(self.output_figure_directory, p_name + "_data_available" + '.png'), dpi=300)
            plt.clf()

    def mann_kendall(self):
        # Read data into arrays
        ds_day = gdal.Open(self.daytime_temp, gdal.GA_ReadOnly)
        ds_night = gdal.Open(self.nighttime_temp, gdal.GA_ReadOnly)
        day_array = ds_day.ReadAsArray().transpose((1, 2, 0))
        night_array = ds_night.ReadAsArray().transpose((1, 2, 0))
        arrays = [day_array, night_array]

        # Create output grids for Mann-Kendall results
        mk_day = np.zeros((day_array.shape[0], day_array.shape[1]))
        mk_night = np.zeros((day_array.shape[0], day_array.shape[1]))

        # Data to write tif files
        originX, pixelWidth, b, originY, d, pixelHeight = ds_day.GetGeoTransform()
        driver = gdal.GetDriverByName('GTiff')

        # Mann - Kendall Analysis
        for _array, array in enumerate(arrays):
            for _row, row in enumerate(array):
                for _col, col in enumerate(row):
                    data = array[_row, _col, :]
                    data = data[~np.isnan(data)]

                    # Mann - Kendall Test original
                    trend, h, p, z, tau, s, var_s, slope, intercept = mk.original_test(data, 0.05)
                    if trend == 'increasing':
                        if _array == 0:
                            mk_day[_row, _col] = 1
                        else:
                            mk_night[_row, _col] = 1
                    elif trend == 'decreasing':
                        print(_array, trend)
                        if _array == 0:
                            mk_day[_row, _col] = 2
                        else:
                            mk_night[_row, _col] = 2
                    else:
                        pass

        # Write to Geotif arrays
        # Save rasters
        day_raster = (os.path.join(self.output_directory, 'mk_day_temp' + ".tif"))
        outRaster = driver.Create(day_raster, day_array.shape[1], day_array.shape[0], 1, gdal.GDT_UInt16)
        outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
        outband = outRaster.GetRasterBand(1)
        outband.WriteArray(mk_day[:, :])
        outband.SetNoDataValue(0)

        # setteing srs from input tif file.
        prj = ds_day.GetProjection()
        outRasterSRS = osr.SpatialReference(wkt=prj)
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()

        night_raster = (os.path.join(self.output_directory, 'mk_night_temp' + ".tif"))
        outRaster = driver.Create(night_raster, day_array.shape[1], day_array.shape[0], 1, gdal.GDT_UInt16)
        outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
        outband = outRaster.GetRasterBand(1)
        outband.WriteArray(mk_night[:, :])
        outband.SetNoDataValue(0)

        # setteing srs from input tif file.
        prj = ds_day.GetProjection()
        outRasterSRS = osr.SpatialReference(wkt=prj)
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()
