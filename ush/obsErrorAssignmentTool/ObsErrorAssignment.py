import argparse
import glob

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import gaussian_kde


def errorAssignment(nc4_files, channel, BinValue, qc_num):
    Bins = np.arange(0, 2.51, BinValue)
    BinsLength = len(Bins) - 1

    # create a dict to hold values for each bin
    QCData = {}
    AllData = {}
    All_o_f = np.array([])  # this keeps o-f raw data (no binning)
    All_clw = np.array([])  # this keeps clw raw data (no binning)

    for i in range(BinsLength):
        BinLow, BinHigh = str(round(Bins[i], 2)), str(round(Bins[i + 1], 2))
        QCData["%s_%s" % (BinLow, BinHigh)] = np.array([])
        AllData["%s_%s" % (BinLow, BinHigh)] = np.array([])

    for fileNo, filePath in enumerate(nc4_files):
        print(fileNo + 1)
        file = xr.open_dataset(filePath)
        channelNum = file["Channel_Index"].isin(channel)
        O_F = file["Obs_Minus_Forecast_adjusted"][channelNum].values
        QC_Flag = file["QC_Flag"][channelNum].values
        clw_mean = (
            file["clw_guess"][channelNum].values + file["clw_obs"][channelNum].values
        ) / 2
        # first remove O minus F outliers
        CheckBool = QC_Flag == qc_num
        O_F = O_F[CheckBool]
        clw_mean = clw_mean[CheckBool]
        # record O-F for use in PDF plot
        All_o_f = np.append(All_o_f, O_F)
        All_clw = np.append(All_clw, clw_mean)
        for i in range(BinsLength):
            BinLow, BinHigh = str(round(Bins[i], 2)), str(round(Bins[i + 1], 2))
            CheckBool = (clw_mean <= float(BinHigh)) & (clw_mean >= float(BinLow))
            ThisFile_O_F = O_F[CheckBool]
            QCData["%s_%s" % (BinLow, BinHigh)] = np.append(
                QCData["%s_%s" % (BinLow, BinHigh)], ThisFile_O_F
            )

        # collect data without QC
        BT = file["Observation"][channelNum].values
        O_F = file["Obs_Minus_Forecast_adjusted"][channelNum].values
        clw_mean = (
            file["clw_guess"][channelNum].values + file["clw_obs"][channelNum].values
        ) / 2

        CheckBool = (BT < 500) & (BT > 50)
        O_F = O_F[CheckBool]
        clw_mean = clw_mean[CheckBool]
        for i in range(BinsLength):
            BinLow, BinHigh = str(round(Bins[i], 2)), str(round(Bins[i + 1], 2))
            CheckBool = (clw_mean <= float(BinHigh)) & (clw_mean >= float(BinLow))
            ThisFile_O_F = O_F[CheckBool]
            AllData["%s_%s" % (BinLow, BinHigh)] = np.append(
                AllData["%s_%s" % (BinLow, BinHigh)], ThisFile_O_F
            )

    AllDataSummary = []
    for bins, bin_values in AllData.items():
        cle_mean = 0.5 * (float(bins.split("_")[0]) + float(bins.split("_")[1]))
        if len(bin_values) > 1:
            AllDataSummary.append(
                [
                    bins,
                    cle_mean,
                    np.mean(bin_values),
                    np.std(bin_values),
                    len(bin_values),
                ]
            )
    AllDataSummaryDF = pd.DataFrame(
        AllDataSummary,
        columns=[
            "Bin",
            "Mean cloud Amount",
            "FG Departure",
            "FG Departure std.Dev",
            "Count",
        ],
    )
    AllDataSummaryDF.dropna(inplace=True)

    QCDataSummary = []
    for bins, bin_values in QCData.items():
        cle_mean = 0.5 * (float(bins.split("_")[0]) + float(bins.split("_")[1]))
        if len(bin_values) > 1:
            QCDataSummary.append(
                [
                    bins,
                    cle_mean,
                    np.mean(bin_values),
                    np.std(bin_values),
                    len(bin_values),
                ]
            )
    QCDataSummaryDF = pd.DataFrame(
        QCDataSummary,
        columns=[
            "Bin",
            "Mean cloud Amount",
            "FG Departure",
            "FG Departure std.Dev",
            "Count",
        ],
    )
    QCDataSummaryDF.dropna(inplace=True)

    return QCDataSummaryDF, AllDataSummaryDF, All_o_f, All_clw


# find the second point
def findDesiredX(QCDataSummaryDF):
    QCDataSummaryDF["x"] = QCDataSummaryDF["Mean cloud Amount"]
    QCDataSummaryDF["y"] = QCDataSummaryDF["FG Departure std.Dev"]
    x = QCDataSummaryDF["x"].values
    y = QCDataSummaryDF["y"].values
    MaxOfY = QCDataSummaryDF["y"].max()

    # coordinates of the first point
    XOfFirst = x[np.argmin(x)]
    YOfFirst = y[np.argmin(x)]

    # find the point on west of max y
    XOfYMax = x[np.argmax(y)]
    CandidatePoints = QCDataSummaryDF[QCDataSummaryDF["x"] <= XOfYMax]
    CandidatePoints = CandidatePoints[CandidatePoints["x"] != XOfFirst]

    maxslope = 0
    for i, row in CandidatePoints.iterrows():
        tempslope = (row["y"] - YOfFirst) / (row["x"] - XOfFirst)
        if tempslope > maxslope:
            maxslope = tempslope

    clw_cld = (MaxOfY - YOfFirst) / maxslope
    return clw_cld + XOfFirst


def ErrorEstimation(QCDataSummaryDF):
    cLW_mean = QCDataSummaryDF["Mean cloud Amount"].values
    err = QCDataSummaryDF["FG Departure std.Dev"].values
    ErrCld = np.max(err)
    ErrClr = err[0]
    lastX = cLW_mean[-1]
    Clw_clr = cLW_mean[0]
    clw_cld = findDesiredX(QCDataSummaryDF)
    x = [0, Clw_clr, clw_cld, lastX]
    y = [ErrClr, ErrClr, ErrCld, ErrCld]
    errorParam = [ErrCld, ErrClr, Clw_clr, clw_cld]
    return x, y, errorParam


def read_previous_par(
    sensor_name, channel_num, path_cloudy_radiance_info, path_global_satinfo
):
    # first from cloudy_radiance file
    senseor_generic_name = sensor_name.split("_")[0]
    with open(path_cloudy_radiance_info) as TXTFile:
        for linenumer, line in enumerate(TXTFile):
            if line.startswith("obs_%s" % senseor_generic_name):
                for _ in range(4):
                    TXTFile.readline()
                while 1:
                    ch, cclr, ccld = TXTFile.readline().strip().split()
                    if int(ch) == channel_num:
                        break
    # now read global_satinfo
    satinfoDF = pd.read_csv(path_global_satinfo, delim_whitespace=True)
    satinfoDF = satinfoDF.rename(columns={"!sensor/instr/sat": "sensor_name"})
    satinfoDF = satinfoDF.set_index("sensor_name")
    satinfoDF = satinfoDF[satinfoDF["chan"] == channel_num]
    error = satinfoDF.loc[sensor_name, "error"]
    error_cld = satinfoDF.loc[sensor_name, "error_cld"]
    x1 = [float(cclr), float(ccld)]
    y1 = [float(error), float(error_cld)]
    return x1, y1


def plots(
    AllDataSummaryDF,
    QCDataSummaryDF,
    output_path,
    sensor,
    channel,
    x,
    y,
    x1,
    y1,
    All_o_f,
    AllErrors,
):
    fig = plt.figure(figsize=(8, 12))
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.scatter(
        QCDataSummaryDF["Mean cloud Amount"],
        QCDataSummaryDF["FG Departure std.Dev"],
        marker="o",
        color="dodgerblue",
        label="QC",
        s=4,
    )
    ax1.scatter(
        AllDataSummaryDF["Mean cloud Amount"],
        AllDataSummaryDF["FG Departure std.Dev"],
        marker="o",
        color="dimgrey",
        label="All",
        s=4,
    )
    ax1.scatter(x1, y1, color="red", label="original", marker="o", s=6)
    ax1.plot(x, y, color="green", label="New", linestyle="--")
    ax1.scatter(x[1:3], y[1:3], marker="X", color="green", label="New Parameters", s=15)

    ax1.set_xlabel("Mean Cloud Amount", fontsize=14)
    ax1.set_ylabel("FG Departure std.Dev", fontsize=14)
    ax1.legend(loc="upper left", markerscale=2, scatterpoints=1)
    inst, sat = sensor.split("_")
    ax1.set_title(
        "%s %s Channel %d" % (inst.upper(), sat.upper(), channel), fontsize=18
    )

    ax2 = plt.subplot(2, 2, 2)
    ax2.grid(True, which="both", color="lightgrey", linestyle="--")
    ax2.scatter(
        QCDataSummaryDF["Mean cloud Amount"],
        QCDataSummaryDF["Count"],
        marker="o",
        color="dodgerblue",
        label="QC",
        s=4,
    )
    ax2.scatter(
        AllDataSummaryDF["Mean cloud Amount"],
        AllDataSummaryDF["Count"],
        marker="o",
        color="dimgrey",
        label="All",
        s=4,
    )
    ax2.set_yscale("log")
    ax2.set_xlabel("Mean cloud Amount", fontsize=14)
    ax2.set_ylabel("Number of observations", fontsize=14)
    ax2.legend(loc="upper right", markerscale=2, scatterpoints=1)

    # Plot PDF of omf
    ax3 = fig.add_subplot(2, 2, 3)
    # Compute the PDF using kernel density estimation
    kde = gaussian_kde(All_o_f)
    x = np.linspace(All_o_f.min(), All_o_f.max(), 100)
    pdf = kde.evaluate(x)
    ax3.plot(x, pdf, label="Un-normalized", color="blue")
    # plot the pDF of normalized FG
    nomalizedFG = All_o_f / AllErrors
    kde_norm = gaussian_kde(nomalizedFG)
    x_norm = np.linspace(nomalizedFG.min(), nomalizedFG.max(), 100)
    pdf_norm = kde_norm.evaluate(x_norm)
    ax3.plot(x_norm, pdf_norm, label="normalized", color="red")
    ax3.legend(loc="upper right", markerscale=2, scatterpoints=1)
    ax3.set_xlabel("FG departure",fontsize=14)
    ax3.set_ylabel("PDF",fontsize=14)
    ax3.set_yscale("log")
    ax3.set_xlim(-10, 10)

    ax4 = fig.add_subplot(2, 2, 4)
    ax4.hist(AllErrors, bins=100, density=True, alpha=0.5, label="Un-normalized")
    ax4.set_xlabel("Errors",fontsize=14)
    ax4.legend(loc="upper right", markerscale=2, scatterpoints=1)

    plt.tight_layout()
    plt.savefig(output_path + "%s_ch%d.png" % (sensor, channel))
    plt.close()


def GetAllErrors(x, y, All_clw):
    Clw_clr, clw_cld = x[1], x[2]
    ErrClr, ErrCld = y[1], y[2]
    AllErrors = np.empty(len(All_clw))
    AllErrors = np.where(All_clw >= clw_cld, ErrCld, AllErrors)
    AllErrors = np.where(All_clw <= Clw_clr, ErrClr, AllErrors)

    # now assign the values in the middle
    slope = (ErrCld - ErrClr) / (clw_cld - Clw_clr)
    MiddleFormula = (All_clw - Clw_clr) * slope + ErrClr
    AllErrors = np.where(
        ((All_clw > Clw_clr) & (All_clw < clw_cld)), MiddleFormula, AllErrors
    )
    return AllErrors


def compute_obs_error_parameter(
    config_path,
    output_path,
    global_satinPath,
    cloudy_path,
    sensor,
    Channels,
    bin_size,
    qc_flag,
):
    print("+=================================================================+")
    print(
        "|                          compute_obs_error_parameter                 |"
    )
    print("+-----------------------------------------------------------------+")
    print("  ---(n) path to config nc files: " + str(config_path))
    print("  ---(o) path to output files: " + str(output_path))
    print("  ---(g) path to global_satinPath.txt: " + str(global_satinPath))
    print("  ---(l) path to cloudy_path.txt: " + str(cloudy_path))
    print("  ---(s) sensor name: " + str(sensor))
    print("  ---(c) channels: " + str(Channels))
    print("  ---(b) Bin Size: " + str(bin_size))
    print("  ---(q) QC Flag: " + str(qc_flag))
    # convert list of channels into a list of integers
    Channels = list(map(int, Channels.split(",")))
    AllErrorParam = []
    nc4_files = glob.glob(config_path + "diag_%s_ges*.nc4" % sensor)
    print(nc4_files)
    for i in Channels:
        print("working on channel %d" % i)
        QCDataSummaryDF, AllDataSummaryDF, All_o_f, All_clw = errorAssignment(
            nc4_files, i, bin_size, qc_flag
        )
        x, y, errorParam = ErrorEstimation(QCDataSummaryDF)
        AllErrors = GetAllErrors(x, y, All_clw)
        x1, y1 = read_previous_par(sensor, i, cloudy_path, global_satinPath)
        errorParam.append(i)
        AllErrorParam.append(errorParam)
        plots(
            AllDataSummaryDF,
            QCDataSummaryDF,
            output_path,
            sensor,
            i,
            x,
            y,
            x1,
            y1,
            All_o_f,
            AllErrors,
        )
    AllErrorParamDF = pd.DataFrame(
        AllErrorParam,
        columns=["Error_cld", "Error_clr", "clw_clr", "clw_cld", "channel"],
    )
    AllErrorParamDF.to_csv(
        output_path + "%s_obsError_parameters.csv" % sensor, index=False
    )


if __name__ == "__main__":
    # Create an argument parser
    parser = argparse.ArgumentParser(description="obs error parameters")

    # Add arguments
    parser.add_argument(
        "-n",
        dest="config_path",
        help=r"REQUIRED: path to config nc files",
        required=True,
        metavar="DIR",
        type=str,
    )
    parser.add_argument(
        "-o",
        dest="output_path",
        help=r"optional: path to output files",
        required=True,
        metavar="DIR",
        type=str,
    )
    parser.add_argument(
        "-g",
        dest="global_satinPath",
        help=r"REQUIRED: path to global_satinPath.txt",
        required=True,
        metavar="DIR",
        type=str,
    )
    parser.add_argument(
        "-l",
        dest="cloudy_path",
        help=r"REQUIRED: path to cloudy_path.txt",
        required=True,
        metavar="DIR",
        type=str,
    )
    parser.add_argument(
        "-s",
        dest="sensor",
        help=r"REQUIRED: sensor name",
        required=True,
        metavar="string",
        type=str,
    )
    parser.add_argument(
        "-c",
        dest="Channels",
        help=r'REQUIRED:list of channels.enclose the list in quotation,'
             r' each separated by a comma like "1,2,3,5"',
        required=True,
        metavar="string",
        type=str,
    )
    parser.add_argument(
        "-b",
        dest="bin_size",
        help=r"optional: size of bin for plotting",
        required=False,
        default=0.05,
        metavar="float",
        type=float,
    )
    parser.add_argument(
        "-q",
        dest="qc_flag",
        help=r"optional: qc flag for filtering",
        required=False,
        default=0,
        metavar="integer",
        type=int,
    )

    args = vars(parser.parse_args())

    config_path = args["config_path"]
    output_path = args["output_path"]
    global_satinPath = args["global_satinPath"]
    cloudy_path = args["cloudy_path"]
    sensor = args["sensor"]
    Channels = args["Channels"]
    bin_size = args["bin_size"]
    qc_flag = args["qc_flag"]

    compute_obs_error_parameter(
        config_path,
        output_path,
        global_satinPath,
        cloudy_path,
        sensor,
        Channels,
        bin_size,
        qc_flag,
    )
