import os

os.environ["USE_PYGEOS"] = "0"

import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
import time
from netCDF4 import Dataset
import matplotlib.pyplot as plt


def plot_forecast_timeseries(dataset, parameter, x, y):
    """
    Creates a graph with forecast timeseries (all ensemble members and median) of a specific xy point locations

    Args:
    - dataset (netCDF4.Dataset): a netcdf dataset class in with the netCDF4 package
    - parameter (str): one of the parameters provided by Deltares (qa, smdi, etdi)
    - x (float): x coordinate (EPSG:4326) 
    - y (float): y coordinate (EPSG:4326)

    Returns:
    - fig (matplotlib figure): figure object with timeseries.
    """
    timestamps = [
        pd.to_datetime(time.ctime(t * 60)) for t in dataset.variables["time"][:]
    ]
    lons = np.array(dataset.variables["x"])
    lats = np.array(dataset.variables["y"])

    fig, ax = plt.subplots(figsize=(16, 8))
    ax.set_xlabel("Date")
    ax.set_ylabel("index")
    ax.grid()

    for ensemble in range(51):
        values = parameter[
            :, ensemble, (np.abs(lons - x)).argmin(), (np.abs(lats - y)).argmin()
        ]
        ax.plot(timestamps, values, color="grey", alpha=0.25)

    medians = []
    for t in range(107):
        medians.append(
            np.median(
                parameter[
                    t, :, (np.abs(lons - x)).argmin(), (np.abs(lats - y)).argmin()
                ]
            )
        )
    ax.plot(timestamps, medians, color="firebrick", linewidth=3)

    return fig


def plot_forecast_spatial(
    dataset,
    parameter,
    timestamp,
    percentiles=[10, 50, 90],
    vmin=-2,
    vmax=2,
    colormap="coolwarm",
):
    """
    Creates a subplot with 3 spatial maps of a 3 percentile values of a certain dataset on a certain timestamp

    Args:
    - dataset (netCDF4.Dataset): a netcdf dataset class in with the netCDF4 package
    - parameter (str): one of the parameters provided by Deltares (qa, smdi, etdi)
    - timestamp (str): timestamp in 'yyyy-mm-dd' format
    - percentiles (list): list with the three percentile values to plot
    - vmin (float): minimum value to be used for colormap
    - vmax (float): maximum value to be used for colormap
    - colormap (str): colormap to be used for the maps

    Returns:
    - fig (matplotlib figure): figure object with spatial maps.
    """
    timestamps = [
        pd.to_datetime(time.ctime(t * 60)).strftime("%Y-%m-%d")
        for t in dataset.variables["time"][:]
    ]
    index = timestamps.index(timestamp)

    fig, axs = plt.subplots(1, 3, figsize=(15, 3))
    axs[0].imshow(
        np.flip(
            np.percentile(parameter[index, :, :, :], percentiles[0], axis=0), axis=0
        ),
        vmin=vmin,
        vmax=vmax,
        cmap=colormap,
    )
    axs[0].set_title("{} - p{}".format(timestamp, percentiles[0]))
    axs[0].grid()

    axs[1].imshow(
        np.flip(
            np.percentile(parameter[index, :, :, :], percentiles[1], axis=0), axis=0
        ),
        vmin=vmin,
        vmax=vmax,
        cmap=colormap,
    )
    axs[1].set_title("{} - p{}".format(timestamp, percentiles[1]))
    axs[1].grid()

    p90 = axs[2].imshow(
        np.flip(
            np.percentile(parameter[index, :, :, :], percentiles[2], axis=0), axis=0
        ),
        vmin=vmin,
        vmax=vmax,
        cmap=colormap,
    )
    axs[2].set_title("{} - p{}".format(timestamp, percentiles[2]))
    axs[2].grid()

    cbar = fig.colorbar(p90, ax=axs, orientation="vertical", fraction=0.1)
    cbar.set_label("Index values")
    return fig


def pixels_to_gdf(dataset):
    """
    Creates a geodataframe with points for each pixel center

    Args:
    - dataset (netCDF4.Dataset): a netcdf dataset class in with the netCDF4 package

    Returns:
    - gdf (geopandas.GeoDataFrame): geodataframe with points for each pixel center
    """
    lats = np.abs(dataset.variables["y"])
    lons = np.abs(dataset.variables["x"])

    points = []
    for x in lons:
        for y in lats:
            points.append(Point(x, y))

    gdf = gpd.GeoDataFrame(geometry=points)
    return gdf
