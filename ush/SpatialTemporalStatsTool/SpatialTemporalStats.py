import argparse
import os
from datetime import datetime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray
from shapely.geometry import Point, Polygon


class SpatialTemporalStats:
    def __init__(self):
        self.grid_gdf = None
        self.obs_gdf = None

    def generate_grid(self, resolution=1):
        self.resolution = resolution
        # Generate the latitude and longitude values using meshgrid
        grid_lons, grid_lats = np.meshgrid(
            np.arange(-180, 181, resolution), np.arange(-90, 91, resolution)
        )

        # Flatten the arrays to get coordinates
        grid_coords = np.vstack([grid_lons.flatten(), grid_lats.flatten()]).T

        # Create a GeoDataFrame from the coordinates
        self.grid_gdf = gpd.GeoDataFrame(
            geometry=[
                Polygon(
                    [
                        (lon, lat),
                        (lon + resolution, lat),
                        (lon + resolution, lat + resolution),
                        (lon, lat + resolution),
                    ]
                )
                for lon, lat in grid_coords
            ],
            crs="EPSG:4326",
        )  # CRS for WGS84
        self.grid_gdf["grid_id"] = np.arange(1, len(self.grid_gdf) + 1)

    def _extract_date_times(self, filenames):
        date_times = []
        for filename in filenames:
            # Split the filename by '.' to get the parts
            parts = filename.split(".")

            # Extract the last part which contains the date/time information
            date_time_part = parts[-2]

            # date/time format in filename is 'YYYYMMDDHH', can parse it accordingly
            year = int(date_time_part[:4])
            month = int(date_time_part[4:6])
            day = int(date_time_part[6:8])
            hour = int(date_time_part[8:10])

            # Construct the datetime object
            date_time = datetime(year, month, day, hour)

            date_times.append(date_time)

        return date_times

    def read_obs_values(
        self,
        obs_files_path,
        sensor,
        var_name,
        channel_no,
        start_date,
        end_date,
        filter_by_vars,
        QC_filter,
    ):
        self.sensor = sensor
        self.channel_no = channel_no
        # read all obs files
        all_files = os.listdir(obs_files_path)
        obs_files = [
            os.path.join(obs_files_path, file)
            for file in all_files
            if file.endswith(".nc4") and "diag_%s_ges" % sensor in file
        ]

        # get date time from file names
        files_date_times_df = pd.DataFrame()

        files_date_times = self._extract_date_times(obs_files)
        files_date_times_df["file_name"] = obs_files
        files_date_times_df["date_time"] = files_date_times
        files_date_times_df["date"] = pd.to_datetime(
            files_date_times_df["date_time"].dt.date
        )

        # read start date
        start_date = datetime.strptime(start_date, "%Y-%m-%d")
        end_date = datetime.strptime(end_date, "%Y-%m-%d")

        studied_cycle_files = files_date_times_df[
            (
                (files_date_times_df["date"] >= start_date)
                & ((files_date_times_df["date"] <= end_date))
            )
        ]["file_name"]

        studied_gdf_list = []
        for this_cycle_obs_file in studied_cycle_files:
            ds = xarray.open_dataset(this_cycle_obs_file)

            Combined_bool = ds["Channel_Index"].data == channel_no

            if QC_filter:
                QC_bool = ds["QC_Flag"].data == 0
                Combined_bool = Combined_bool * QC_bool

            # apply filters by variable
            for this_filter in filter_by_vars:
                filter_var_name, filter_operation, filter_value = this_filter
                if filter_operation == "lt":
                    this_filter_bool = ds[filter_var_name].data <= filter_value
                else:
                    this_filter_bool = ds[filter_var_name].data >= filter_value
                Combined_bool = (
                    Combined_bool * ~this_filter_bool
                )  # here we have to negate the above bool to make it right

            this_cycle_var_values = ds[var_name].data[Combined_bool]
            this_cycle_lat_values = ds["Latitude"].data[Combined_bool]
            this_cycle_long_values = ds["Longitude"].data[Combined_bool]
            this_cycle_long_values = np.where(
                this_cycle_long_values <= 180,
                this_cycle_long_values,
                this_cycle_long_values - 360,
            )

            geometry = [
                Point(xy) for xy in zip(this_cycle_long_values, this_cycle_lat_values)
            ]

            # Create a GeoDataFrame
            this_cycle_gdf = gpd.GeoDataFrame(geometry=geometry, crs="EPSG:4326")
            this_cycle_gdf["value"] = this_cycle_var_values

            studied_gdf_list.append(this_cycle_gdf)

        studied_gdf = pd.concat(studied_gdf_list)

        # Perform spatial join
        joined_gdf = gpd.sjoin(studied_gdf, self.grid_gdf, op="within", how="right")

        # Calculate average values of points in each polygon
        self.obs_gdf = self.grid_gdf.copy()
        self.obs_gdf[var_name + "_Average"] = joined_gdf.groupby("grid_id")[
            "value"
        ].mean()
        self.obs_gdf[var_name + "_RMS"] = joined_gdf.groupby("grid_id")["value"].apply(
            lambda x: np.sqrt((x**2).mean())
        )
        self.obs_gdf[var_name + "_Count"] = joined_gdf.groupby("grid_id")[
            "value"
        ].count()

        # convert count of zero to null. This will help also for plotting
        self.obs_gdf[var_name + "_Count"] = np.where(
            self.obs_gdf[var_name + "_Count"].values == 0,
            np.nan,
            self.obs_gdf[var_name + "_Count"].values,
        )

        return self.obs_gdf

    def plot_obs(self, selected_var_gdf, var_name, region, resolution, output_path):
        self.resolution = resolution
        var_names = [var_name + "_Average", var_name + "_Count", var_name + "_RMS"]

        for _, item in enumerate(var_names):
            plt.figure(figsize=(12, 8))
            ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
            # Add global map coastlines
            ax.add_feature(cfeature.GSHHSFeature(scale="auto"))
            filtered_gdf = selected_var_gdf.copy()

            if region == 1:
                # Plotting global region (no need for filtering)
                title = "Global"
                # filtered_gdf = selected_var_gdf

            elif region == 2:
                # Plotting polar region (+60 latitude and above)
                title = "Polar Region (+60 latitude and above)"
                filtered_gdf[item] = np.where(
                    filtered_gdf.geometry.apply(
                        lambda geom: self.is_polygon_in_polar_region(geom, 60)
                    ),
                    filtered_gdf[item],
                    np.nan,
                )

            elif region == 3:
                # Plotting northern mid-latitudes region (20 to 60 latitude)
                title = "Northern Mid-latitudes Region (20 to 60 latitude)"
                filtered_gdf[item] = np.where(
                    filtered_gdf.geometry.apply(
                        lambda geom: self.is_polygon_in_latitude_range(geom, 20, 60)
                    ),
                    filtered_gdf[item],
                    np.nan,
                )

            elif region == 4:
                # Plotting tropics region (-20 to 20 latitude)
                title = "Tropics Region (-20 to 20 latitude)"

                filtered_gdf[item] = np.where(
                    filtered_gdf.geometry.apply(
                        lambda geom: self.is_polygon_in_latitude_range(geom, -20, 20)
                    ),
                    filtered_gdf[item],
                    np.nan,
                )

            elif region == 5:
                # Plotting southern mid-latitudes region (-60 to -20 latitude)
                title = "Southern Mid-latitudes Region (-60 to -20 latitude)"
                filtered_gdf[item] = np.where(
                    filtered_gdf.geometry.apply(
                        lambda geom: self.is_polygon_in_latitude_range(geom, -60, -20)
                    ),
                    filtered_gdf[item],
                    np.nan,
                )

            elif region == 6:
                # Plotting southern polar region (less than -60 latitude)
                title = "Southern Polar Region (less than -60 latitude)"
                filtered_gdf[item] = np.where(
                    filtered_gdf.geometry.apply(lambda geom: geom.centroid.y < -60),
                    filtered_gdf[item],
                    np.nan,
                )
                # filtered_gdf = selected_var_gdf[
                #     selected_var_gdf.geometry.apply(lambda geom: geom.centroid.y < -60)
                # ]

            min_val, max_val, std_val, avg_val = (
                filtered_gdf[item].min(),
                filtered_gdf[item].max(),
                filtered_gdf[item].std(),
                filtered_gdf[item].mean(),
            )

            if item == "Obs_Minus_Forecast_adjusted_Average":
                max_val_cbar = 5.0 * std_val
                min_val_cbar = -5.0 * std_val
                cmap = "bwr"
            else:
                max_val_cbar = max_val
                min_val_cbar = min_val
                cmap = "jet"

            if item == "Obs_Minus_Forecast_adjusted_Count":
                cbar_label = "grid=%dx%d,   min=%.3lf,   max=%.3lf\n" % (
                    resolution,
                    resolution,
                    min_val,
                    max_val,
                )
            else:
                cbar_label = (
                    "grid=%dx%d,   min=%.3lf,   max=%.3lf,   bias=%.3lf,   std=%.3lf\n"
                    % (
                        resolution,
                        resolution,
                        min_val,
                        max_val,
                        avg_val,
                        std_val,
                    )
                )

            filtered_gdf.plot(
                ax=ax,
                cmap=cmap,
                vmin=min_val_cbar,
                vmax=max_val_cbar,
                column=item,
                legend=True,
                missing_kwds={"color": "lightgrey"},
                legend_kwds={
                    "orientation": "horizontal",
                    "shrink": 0.5,
                    "label": cbar_label,
                },
            )

            filtered_gdf.to_file(
                os.path.join(
                    output_path,
                    "%s_ch%d_%s_region_%d.gpkg"
                    % (self.sensor, self.channel_no, item, region),
                )
            )

            plt.title("%s\n%s ch:%d %s" % (title, self.sensor, self.channel_no, item))
            plt.savefig(
                os.path.join(
                    output_path,
                    "%s_ch%d_%s_region_%d.png"
                    % (self.sensor, self.channel_no, item, region),
                )
            )
            plt.close()

    def is_polygon_in_polar_region(self, polygon, latitude_threshold):
        """
        Check if a polygon is in the polar region based on a latitude threshold.
        """
        # Get the centroid of the polygon
        centroid = polygon.centroid

        # Extract the latitude of the centroid
        centroid_latitude = centroid.y

        # Check if the latitude is above the threshold
        return centroid_latitude >= latitude_threshold

    def is_polygon_in_latitude_range(self, polygon, min_latitude, max_latitude):
        """
        Check if a polygon is in the specified latitude range.
        """
        # Get the centroid of the polygon
        centroid = polygon.centroid

        # Extract the latitude of the centroid
        centroid_latitude = centroid.y

        # Check if the latitude is within the specified range
        return min_latitude <= centroid_latitude <= max_latitude

    def list_variable_names(self, file_path):
        ds = xarray.open_dataset(file_path)
        print(ds.info())

    def make_summary_plots(
        self,
        obs_files_path,
        sensor,
        var_name,
        start_date,
        end_date,
        QC_filter,
        output_path,
    ):
        self.sensor = sensor
        # read all obs files
        all_files = os.listdir(obs_files_path)
        obs_files = [
            os.path.join(obs_files_path, file)
            for file in all_files
            if file.endswith(".nc4") and "diag_%s_ges" % sensor in file
        ]

        # get date time from file names.
        # alternatively could get from attribute but that needs reading the entire nc4
        files_date_times_df = pd.DataFrame()

        files_date_times = self._extract_date_times(obs_files)
        files_date_times_df["file_name"] = obs_files
        files_date_times_df["date_time"] = files_date_times
        files_date_times_df["date"] = pd.to_datetime(
            files_date_times_df["date_time"].dt.date
        )

        # read start date
        start_date = datetime.strptime(start_date, "%Y-%m-%d")
        end_date = datetime.strptime(end_date, "%Y-%m-%d")

        studied_cycle_files = files_date_times_df[
            (
                (files_date_times_df["date"] >= start_date)
                & ((files_date_times_df["date"] <= end_date))
            )
        ]["file_name"]
        index = studied_cycle_files.index

        Summary_results = []

        # get unique channels from one of the files
        ds = xarray.open_dataset(studied_cycle_files[index[0]])
        unique_channels = np.unique(ds["Channel_Index"].data).tolist()
        print("Total Number of Channels ", len(unique_channels))
        Allchannels_data = {}
        for this_channel in unique_channels:
            Allchannels_data[this_channel] = np.empty(shape=(0,))
        for this_cycle_obs_file in studied_cycle_files:
            ds = xarray.open_dataset(this_cycle_obs_file)
            if QC_filter:
                QC_bool = ds["QC_Flag"].data == 0.0
            else:
                QC_bool = np.ones(
                    ds["QC_Flag"].data.shape, dtype=bool
                )  # this selects all obs as True
            for this_channel in unique_channels:
                channel_bool = ds["Channel_Index"].data == this_channel

                this_cycle_channel_var_values = ds[var_name].data[
                    channel_bool * QC_bool
                ]
                Allchannels_data[this_channel] = np.append(
                    Allchannels_data[this_channel], this_cycle_channel_var_values
                )

        for this_channel in unique_channels:
            this_channel_values = Allchannels_data[this_channel]
            squared_values = [x**2 for x in this_channel_values]
            mean_of_squares = sum(squared_values) / len(squared_values)
            rms_value = mean_of_squares**0.5
            Summary_results.append(
                [
                    this_channel,
                    np.size(this_channel_values),
                    np.std(this_channel_values),
                    np.mean(this_channel_values),
                    rms_value,
                ]
            )

        Summary_resultsDF = pd.DataFrame(
            Summary_results, columns=["channel", "count", "std", "mean", "rms"]
        )
        # Plotting
        plt.figure(figsize=(10, 6))
        plt.scatter(Summary_resultsDF["channel"], Summary_resultsDF["count"], s=50)
        plt.xlabel("Channel")
        plt.ylabel("Count")
        plt.title("%s %s" % ((self.sensor, var_name)))
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(
            os.path.join(
                output_path, "%s_%s_sumamryCounts.png" % (self.sensor, var_name)
            )
        )
        plt.close()

        # Plotting scatter plot for mean and std
        plt.figure(figsize=(15, 6))
        plt.scatter(
            Summary_resultsDF["channel"],
            Summary_resultsDF["mean"],
            s=50,
            c="green",
            label="Mean",
        )
        plt.scatter(
            Summary_resultsDF["channel"],
            Summary_resultsDF["std"],
            s=50,
            c="red",
            label="Std",
        )
        plt.scatter(
            Summary_resultsDF["channel"],
            Summary_resultsDF["rms"],
            s=50,
            label="Rms",
            facecolors="none",
            edgecolors="blue",
        )
        plt.xlabel("Channel")
        plt.ylabel("Statistics")
        plt.title("%s %s" % ((self.sensor, var_name)))
        plt.grid(True)
        plt.tight_layout()
        plt.legend()
        plt.savefig(
            os.path.join(output_path, "%s_%s_mean_std.png" % (self.sensor, var_name))
        )

        return Summary_resultsDF


def main(
    input_path,
    output_path,
    sensor,
    var_name,
    channel,
    grid_size,
    qc_flag,
    region,
    start_date,
    end_date,
    filter_by_vars,
):
    # Initialize SpatialTemporalStats object
    my_tool = SpatialTemporalStats()

    # Generate grid
    my_tool.generate_grid(grid_size)  # Call generate_grid method)
    print("grid created!")

    # Read observational values and perform analysis
    o_minus_f_gdf = my_tool.read_obs_values(
        input_path,
        sensor,
        var_name,
        channel,
        start_date,
        end_date,
        filter_by_vars,
        qc_flag,
    )

    print("read obs values!")

    # Plot observations
    print("creating plots...")

    my_tool.plot_obs(o_minus_f_gdf, var_name, region, grid_size, output_path)
    print("Time/Area stats plots created!")

    # Make summary plots
    print("Creating summary plots...")
    summary_results = my_tool.make_summary_plots(
        input_path, sensor, var_name, start_date, end_date, qc_flag, output_path
    )
    summary_results.to_csv(
        os.path.join(output_path, "%s_summary.csv" % sensor), index=False
    )
    print("Summary plots created!")


def parse_filter(s):
    try:
        var_name, comparison, threshold = s.split(",")
        if comparison not in ("lt", "gt"):
            raise ValueError("Comparison must be 'lt' or 'gt'")
        return (var_name, comparison, float(threshold))
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Filter must be in format 'var_name,comparison,threshold'"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Python Tool for Spatial and Temporal Analysis"
    )

    # Add arguments
    parser.add_argument(
        "-input",
        dest="input_path",
        help=r"REQUIRED: path to input config nc files",
        required=True,
        metavar="DIR",
        type=str,
    )
    parser.add_argument(
        "-output",
        dest="output_path",
        help=r"REQUIRED: path to output files",
        required=True,
        metavar="DIR",
        type=str,
    )

    parser.add_argument(
        "-sensor",
        dest="sensor",
        help=r"REQUIRED: sensor name",
        required=True,
        metavar="string",
        type=str,
    )
    parser.add_argument(
        "-var",
        dest="var_name",
        help=r"REQUIRED: variable name",
        required=True,
        metavar="string",
        type=str,
    )
    parser.add_argument(
        "-ch",
        dest="channel",
        help=r"REQUIRED the channel number",
        required=True,
        metavar="integer",
        type=int,
    )
    parser.add_argument(
        "-grid",
        dest="grid_size",
        help=r"optional: size of grid for plotting (choices: 0.5, 1, 2)",
        required=False,
        default=1,
        metavar="float",
        type=float,
    )
    parser.add_argument(
        "-no_qc_flag",
        dest="no_qc_flag",
        help=r"Optional: qc flag for filtering",
        action="store_true",
    )
    parser.add_argument(
        "-region",
        dest="region",
        help="REQUIRED: region for mapplot. 1: global, 2: polar region, 3: mid-latitudes region,"
        "4: tropics region, 5:southern mid-latitudes region, 6: southern polar region",
        required=False,
        default=0,
        metavar="integer",
        type=int,
    )
    parser.add_argument(
        "-sdate",
        dest="start_date",
        help=r"REQUIRED: start date of evaluation",
        required=False,
        default=0,
        metavar="string",
        type=str,
    )
    parser.add_argument(
        "-edate",
        dest="end_date",
        help=r"REQUIRED: end date of evaluation",
        required=False,
        default=0,
        metavar="string",
        type=str,
    )

    # New argument for filter criteria
    parser.add_argument(
        "-filter_by_vars",
        dest="filter_by_vars",
        help="Optional: Filtering criteria in format 'var_name,comparison,"
        "threshold'. Example: Land_Fraction,lt,0.9",
        nargs="+",
        type=parse_filter,
        default=[],
    )

    args = vars(parser.parse_args())

    input_path = args["input_path"]
    output_path = args["output_path"]
    sensor = args["sensor"]
    var_name = args["var_name"]
    channel = args["channel"]
    grid_size = args["grid_size"]
    region = args["region"]
    start_date = args["start_date"]
    end_date = args["end_date"]

    if args["no_qc_flag"]:
        qc_flag = False
    else:
        qc_flag = True

    # Accessing and printing the parsed filter criteria
    if args["filter_by_vars"]:
        for filter_criteria in args["filter_by_vars"]:
            print(
                f"Variable: {filter_criteria[0]},"
                f"Comparison: {filter_criteria[1]},"
                f"Threshold: {filter_criteria[2]}"
            )

    main(
        input_path,
        output_path,
        sensor,
        var_name,
        channel,
        grid_size,
        qc_flag,
        region,
        start_date,
        end_date,
        args["filter_by_vars"],
    )
