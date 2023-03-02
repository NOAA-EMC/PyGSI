import numpy as _np
import os as _os
import emcpy.utils.dateutils as _dateutils
import pyGSI.filter_obs as _filter_obs
from netCDF4 import Dataset


def time_trace(
    datapath,
    date1,
    date2,
    expt_names,
    n_mem,
    delt,
    skip_enkf_hours=[],
    ob_types=["u"],
    codes_uv=[187],
    codes_tq=[287],
    hem=None,
    p_max=1050.0,
    p_min=100.0,
    lat_max=90.0,
    lat_min=0.0,
    lon_max=360.0,
    lon_min=0.0,
    error_max=40.0,
    error_min=0.000001,
):
    """
    Computes observation space stats into 3D arrays (n_ob_type, n_expt, 24) and plots the stats as a
    function of hour of day.

    Args:
      datapath        : (str) netCDF filename
      date1           : (str "YYYYMMDDHH") start date
      date2           : (str "YYYYMMDDHH") end date
      expt_names      : (list of str) experiment names
      n_mem           : (int) number of ensemble members
      delt            : (int) delta time in hours
      skip_enkf_hours : (list of int) hours (UTC) to skip for enkf
      ob_types        : (list of str) observation types (u,v,t,q,etc.)
      codes_uv        : (list of int) uv bufr codes used to filter obs
      codes_tq        : (list of int) tq bufr codes used to filter obs
      hem             : (str) name of hemisphere to plot stats (also None)
      p_max           : (float) maximum pressure (mb) for including observation in calculations
      p_min           : (float) minimum pressure (mb) for including observation in calculations
      lat_max         : (float) maximum latitude (deg N) for including observation in calculations
      lat_min         : (float) minimum latitude (deg N) for including observation in calculations
      lon_max         : (float) maximum latitude (deg E) for including observation in calculations
      lon_min         : (float) minimum latitude (deg E) for including observation in calculations
      error_max       : (float) maximum error standard deviation for including observation in calculations
      error_min       : (float) minimum error standard deviation for including observation in calculations

    Returns:
      dates         : (list str) list of date strings of the form YYYYMMDDHH based on date1 and date2
      bias          : (array float) mean of (forecast - observation)
      rms           : (array float) rms of (F-O)
      std_dev       : (array float) standard deviation of (F-O)
      spread        : (array float) ensemble spread (standard deviation)
      ob_error      : (array float) observation error standard deviation
      total_spread  : (array float) total spread (standard deviation)
      num_obs_total : (array float) total number of observations
      num_obs_assim : (array float) total number of observations assimilated
      cr            : (array float) consistency ratio (total spread/rmsi)**2
      ser           : (array float) spread error ratio (intraensemble std_dev/ rmse of ensemble mean fcst)

    References:
        consistency ratio (cr)
          a) https://journals.ametsoc.org/view/journals/atot/26/5/2008jtecha1156_1.xml (eq 3.4)
          b) https://journals.ametsoc.org/view/journals/mwre/150/8/MWR-D-21-0289.1.xml (eq 7)
        spread error ratio (ser)
          a) https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013GL057630 (eq 1)
        bufr codes
          a) https://emc.ncep.noaa.gov/mmb/data_processing/prepbufr.doc/table_2.htm
    """

    hours = 24
    dates = _dateutils.daterange(date1, date2, delt)
    n_expt = len(expt_names)
    n_ob_type = len(ob_types)
    bias = _np.zeros(shape=(n_ob_type, n_expt, 24))
    rms = _np.zeros(shape=(n_ob_type, n_expt, 24))
    std_dev = _np.zeros(shape=(n_ob_type, n_expt, 24))
    spread = _np.zeros(shape=(n_ob_type, n_expt, 24))
    total_spread = _np.zeros(shape=(n_ob_type, n_expt, 24))
    ob_error = _np.zeros(shape=(n_ob_type, n_expt, 24))
    num_obs_assim = _np.zeros(shape=(n_ob_type, n_expt, 24))
    num_obs_total = _np.zeros(shape=(n_ob_type, n_expt, 24))
    cr = _np.zeros(shape=(n_ob_type, n_expt, 24))
    ser = _np.zeros(shape=(n_ob_type, n_expt, 24))

    sum_innov = _np.zeros(shape=(n_ob_type, n_expt, 24))
    sum_innovsq = _np.zeros(shape=(n_ob_type, n_expt, 24))
    sum_fcst_ens_var = _np.zeros(shape=(n_ob_type, n_expt, 24))
    sum_ob_err_var = _np.zeros(shape=(n_ob_type, n_expt, 24))

    bbreak = False

    for date in dates:
        times = _dateutils.datetohrs(date)
        pdy = str(date[0:8])
        hour = int(str(date[8:10]))
        if any(item == hour for item in skip_enkf_hours):
            continue  # As of right now, the rrfs only runs the EnKF at 18-23Z so skip those off times
        for ob_type in ob_types:
            i_o = ob_types.index(ob_type)
            # https://emc.ncep.noaa.gov/mmb/data_processing/prepbufr.doc/table_2.htm
            if ob_type == "u" or ob_type == "v":
                codes = codes_uv
            elif ob_type == "t" or ob_type == "q":
                codes = codes_tq
            for expt_name in expt_names:
                print(f"{date} {expt_name} {ob_type}")
                i_e = expt_names.index(expt_name)

                # Initial read of all the ensemble diag files to get a common list of used observations.
                # Preemptively read analysis use flag of all ensemble diags. All ensembles will not
                # necessarily use the same observations (gross errors checks etc.). Change all use
                # values that are less than 1 to 0, then multiply each use numpy array by the previous
                # result. This will leave you with a list of all observations used by all ensemble members.
                mem = 1
                while mem <= n_mem:
                    memid = str(mem).zfill(4)
                    if ob_type == "u" or ob_type == "v":
                        diagfile = _os.path.join(datapath, f"{expt_name}/{date}/mem{memid}/diag_conv_uv_ges.{date}.nc4")
                    else:
                        diagfile = _os.path.join(datapath, f"{expt_name}/{date}/mem{memid}/diag_conv_{ob_type}_ges.{date}.nc4")
                    # Check for if there are any diag files for a particular cycle date.
                    # If there are none, move to the next cycle date.
                    exists = _os.path.exists(diagfile)
                    if not exists:
                        print(f"diag file for {expt_name} mem{mem:0>4} {date} doesn't exist. Skipping {date}.")
                        bbreak = True  # need to break here and one loop higher.
                        break

                    nc = Dataset(diagfile)
                    use = nc["Analysis_Use_Flag"][:]
                    use[use < 1] = 0
                    if mem == 1:
                        analysis_use = use
                    else:
                        analysis_use = analysis_use * use
                    mem = mem + 1
                    # End of initial read of all ensemble diag files.

                if bbreak:
                    bbreak = False  # reset break
                    break  # break and move to next cycle time.

                mem = 1
                while mem <= n_mem:
                    memid = str(mem).zfill(4)
                    if ob_type == "u" or ob_type == "v":
                        diagfile = _os.path.join(datapath, f"{expt_name}/{date}/mem{memid}/diag_conv_uv_ges.{date}.nc4")
                    else:
                        diagfile = _os.path.join(datapath, f"{expt_name}/{date}/mem{memid}/diag_conv_{ob_type}_ges.{date}.nc4")
                    nc = Dataset(diagfile)
                    if mem == 1:
                        use = analysis_use
                        code = nc["Observation_Type"][:]
                        lat = nc["Latitude"][:]
                        lon = nc["Longitude"][:]
                        pressure = nc["Pressure"][:]
                        errorinv = nc["Errinv_Final"][:]

                        if ob_type == "u":
                            ob = nc["u_Observation"][:]
                        elif ob_type == "v":
                            ob = nc["v_Observation"][:]
                        else:
                            ob = nc["Observation"][:]

                        if ob_type == "q":
                            ob = 1000.0 * ob  # convert from kg/kg to g/kg
                            errorinv = errorinv / 1000.0  # convert from kg/kg to g/kg

                        # consider where use flag==1 and bound by error/lat/lon/pres
                        used = _filter_obs.filter_obs(
                            code,
                            codes,
                            errorinv,
                            lat,
                            lon,
                            pressure,
                            use,
                            hem,
                            p_max=p_max,
                            p_min=p_min,
                            lat_max=lat_max,
                            lat_min=lat_min,
                            lon_max=lon_max,
                            lon_min=lon_min,
                            error_max=error_max,
                            error_min=error_min,
                        )

                        errorinv = errorinv[used]
                        error = 1.0 / errorinv
                        itot = len(ob)
                        ob = ob[used]
                        iasm = len(ob)
                        # end if mem==1

                    if ob_type == "u":
                        omf = nc["u_Obs_Minus_Forecast_adjusted"][:]
                    elif ob_type == "v":
                        omf = nc["v_Obs_Minus_Forecast_adjusted"][:]
                    else:
                        omf = nc["Obs_Minus_Forecast_adjusted"][:]
                    if ob_type == "q":
                        omf = 1000.0 * omf  # convert from kg/kg to g/kg

                    omf = omf[used]

                    if mem == 1:
                        fcst_ens = ob - omf
                        fcst_ens_var = (ob - omf) ** 2
                    else:
                        fcst_ens = fcst_ens + ob - omf
                        fcst_ens_var = fcst_ens_var + (ob - omf) ** 2
                    mem = mem + 1
                    # end while n_mem

                # if we broke out of the loop above before, because we were missing a mem.
                # this will break out of the loop here too.
                try:
                    error_var = error**2
                except NameError:
                    break

                fcst_ens_mean = fcst_ens / n_mem
                if n_mem > 1:
                    fcst_ens_var = (fcst_ens_var - n_mem * fcst_ens_mean**2) / (n_mem - 1)
                else:
                    fcst_ens_var = 0.0
                innov = ob - fcst_ens_mean
                num_obs_total[i_o, i_e, hour] = num_obs_total[i_o, i_e, hour] + itot
                num_obs_assim[i_o, i_e, hour] = num_obs_assim[i_o, i_e, hour] + iasm
                sum_innov[i_o, i_e, hour] = sum_innov[i_o, i_e, hour] + _np.sum(innov)
                sum_innovsq[i_o, i_e, hour] = sum_innovsq[i_o, i_e, hour] + _np.sum(innov**2)
                sum_fcst_ens_var[i_o, i_e, hour] = sum_fcst_ens_var[i_o, i_e, hour] + _np.sum(fcst_ens_var)
                sum_ob_err_var[i_o, i_e, hour] = sum_ob_err_var[i_o, i_e, hour] + _np.sum(error_var)

                # end n_times

                if num_obs_assim[i_o, i_e, hour] > 0:
                    mean_innov = sum_innov[i_o, i_e, hour] / num_obs_assim[i_o, i_e, hour]
                    bias[i_o, i_e, hour] = -1 * mean_innov
                    rms[i_o, i_e, hour] = _np.sqrt(sum_innovsq[i_o, i_e, hour] / num_obs_assim[i_o, i_e, hour])
                    mean_ob_err_var = sum_ob_err_var[i_o, i_e, hour] / num_obs_assim[i_o, i_e, hour]
                    ob_error[i_o, i_e, hour] = _np.sqrt(mean_ob_err_var)

                if num_obs_assim[i_o, i_e, hour] > 1:
                    innov_var = (sum_innovsq[i_o, i_e, hour] - num_obs_assim[i_o, i_e, hour] * mean_innov**2) / (num_obs_assim[i_o, i_e, hour] - 1.0)
                    std_dev[i_o, i_e, hour] = _np.sqrt(innov_var)
                    mean_fcst_var = sum_fcst_ens_var[i_o, i_e, hour] / num_obs_assim[i_o, i_e, hour]
                    spread[i_o, i_e, hour] = _np.sqrt(mean_fcst_var)
                    total_spread[i_o, i_e, hour] = _np.sqrt(mean_ob_err_var + mean_fcst_var)
                    cr[i_o, i_e, hour] = (total_spread[i_o, i_e, hour] / rms[i_o, i_e, hour]) ** 2
                    ser[i_o, i_e, hour] = spread[i_o, i_e, hour] / rms[i_o, i_e, hour]

                del errorinv
                del error
                del itot
                del ob
                del iasm
                del omf
                del code
                del lat
                del lon
                del pressure
                del use

            # end do n_expt
        # end do n_ob_type

    return dates, bias, rms, std_dev, spread, ob_error, total_spread, num_obs_total, num_obs_assim, cr, ser


def bincount(pindx, n_levs, variable, count):
    variable = _np.bincount(pindx, minlength=n_levs, weights=variable) / count[::-1]
    return variable

def profile(
    datapath,
    date1,
    date2,
    expt_names,
    n_mem,
    delt,
    skip_enkf_hours=[],
    ob_types=["u"],
    codes_uv=[187],
    codes_tq=[287],
    hem=None,
    p_max=1050.0,
    p_min=100.0,
    lat_max=90.0,
    lat_min=0.0,
    lon_max=360.0,
    lon_min=0.0,
    error_max=40.0,
    error_min=0.000001,
):
    """
    Computes observation space stats into 4D arrays (n_ob_type, n_expt, n_dates, n_levs) for plotting as a
    vertical profile.

    Args:
      datapath        : (str) netCDF filename
      date1           : (str "YYYYMMDDHH") start date
      date2           : (str "YYYYMMDDHH") end date
      expt_names      : (list of str) experiment names
      n_mem           : (int) number of ensemble members
      delt            : (int) delta time in hours
      skip_enkf_hours : (list of int) hours (UTC) to skip for enkf
      ob_types        : (list of str) observation types (u,v,t,q,etc.)
      codes_uv        : (list of int) uv bufr codes used to filter obs
      codes_tq        : (list of int) tq bufr codes used to filter obs
      hem             : (str) name of hemisphere to plot stats (also None)
      p_max           : (float) maximum pressure (mb) for including observation in calculations
      p_min           : (float) minimum pressure (mb) for including observation in calculations
      lat_max         : (float) maximum latitude (deg N) for including observation in calculations
      lat_min         : (float) minimum latitude (deg N) for including observation in calculations
      lon_max         : (float) maximum latitude (deg E) for including observation in calculations
      lon_min         : (float) minimum latitude (deg E) for including observation in calculations
      error_max       : (float) maximum error standard deviation for including observation in calculations
      error_min       : (float) minimum error standard deviation for including observation in calculations

    Returns:
      dates         : (list str) list of date strings of the form YYYYMMDDHH based on date1 and date2
      bias          : (array float) mean of (forecast - observation)
      rms           : (array float) rms of (F-O)
      std_dev       : (array float) standard deviation of (F-O)
      spread        : (array float) ensemble spread (standard deviation)
      ob_error      : (array float) observation error standard deviation
      total_spread  : (array float) total spread (standard deviation)
      num_obs_total : (array float) total number of observations
      num_obs_assim : (array float) total number of observations assimilated
      cr            : (array float) consistency ratio (total spread/rmsi)**2
      ser           : (array float) spread error ratio (intraensemble std_dev/ rmse of ensemble mean fcst)

    References:
        consistency ratio (cr)
          a) https://journals.ametsoc.org/view/journals/atot/26/5/2008jtecha1156_1.xml (eq 3.4)
          b) https://journals.ametsoc.org/view/journals/mwre/150/8/MWR-D-21-0289.1.xml (eq 7)
        spread error ratio (ser)
          a) https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013GL057630 (eq 1)
        bufr codes
          a) https://emc.ncep.noaa.gov/mmb/data_processing/prepbufr.doc/table_2.htm
    """

    dates = _dateutils.daterange(date1, date2, delt)
    n_dates = len(dates)
    n_expt = len(expt_names)
    n_ob_type = len(ob_types)
    n_skip = 0

    deltap = 50.0
    pbot = 975
    n_levs = 23
    levs = _np.zeros(n_levs, float)
    levs1 = _np.zeros(n_levs, float)
    levs2 = _np.zeros(n_levs, float)
    levs[0:18] = pbot - deltap * _np.arange(18)
    levs1[0:18] = levs[0:18] + 0.5 * deltap
    levs2[0:18] = levs[0:18] - 0.5 * deltap
    levs1[18] = levs2[17]
    levs2[18] = 70.0
    levs1[19] = 70.0
    levs2[19] = 50.0
    levs1[20] = 50.0
    levs2[20] = 30.0
    levs1[21] = 30.0
    levs2[21] = 10.0
    levs1[22] = 10.0
    levs2[22] = 0.0
    levs1[0] = 1200.0
    pbins = _np.zeros(n_levs + 1, float)
    pbins[0:n_levs] = levs1
    pbins[n_levs] = levs2[-1]
    for nlev in range(18, n_levs):
        levs[nlev] = 0.5 * (levs1[nlev] + levs2[nlev])

    bias = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    rms = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    std_dev = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    spread = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    total_spread = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    ob_error = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    num_obs_assim = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    num_obs_total = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    cr = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    ser = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))

    sum_innov = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    sum_innovsq = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    sum_fcst_ens_var = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    sum_fcst_ens = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    sum_ob_err_var = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    sum_ob_err = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))

    mean_innov = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    mean_ob_err_var = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    innov_var = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))
    mean_fcst_var = _np.zeros(shape=(n_ob_type, n_expt, n_dates, n_levs))

    iasms = _np.zeros(shape=(n_dates, n_levs))

    bbreak = False

    for date in dates:
        i_date = dates.index(date)
        times = _dateutils.datetohrs(date)
        pdy = str(date[0:8])
        hour = int(str(date[8:10]))
        if any(item == hour for item in skip_enkf_hours):
            continue  # As of right now, the rrfs only runs the EnKF at 18-23Z so skip those off times
        for ob_type in ob_types:
            i_o = ob_types.index(ob_type)
            # https://emc.ncep.noaa.gov/mmb/data_processing/prepbufr.doc/table_2.htm
            if ob_type == "u" or ob_type == "v":
                codes = codes_uv
            elif ob_type == "t" or ob_type == "q":
                codes = codes_tq
            for expt_name in expt_names:
                print(f"{date} {expt_name} {ob_type}")
                i_e = expt_names.index(expt_name)

                # Initial read of all the ensemble diag files to get a common list of used observations.
                # Preemptively read analysis use flag of all ensemble diags. All ensembles will not
                # necessarily use the same observations (gross errors checks etc.). Change all use
                # values that are less than 1 to 0, then multiply each use numpy array by the previous
                # result. This will leave you with a list of all observations used by all ensemble members.
                mem = 1
                while mem <= n_mem:
                    memid = str(mem).zfill(4)
                    if ob_type == "u" or ob_type == "v":
                        diagfile = _os.path.join(datapath, f"{expt_name}/{date}/mem{memid}/diag_conv_uv_ges.{date}.nc4")
                    else:
                        diagfile = _os.path.join(datapath, f"{expt_name}/{date}/mem{memid}/diag_conv_{ob_type}_ges.{date}.nc4")
                    # Check for if there are any diag files for a particular cycle date.
                    # If there are none, move to the next cycle date.
                    exists = _os.path.exists(diagfile)
                    if not exists:
                        print(f"diag file for {expt_name} mem{mem:0>4} {date} doesn't exist. Skipping {date}.")
                        n_skip = n_skip + 1
                        bbreak = True  # need to break here and one loop higher.
                        break

                    nc = Dataset(diagfile)
                    use = nc["Analysis_Use_Flag"][:]
                    use[use < 1] = 0
                    if mem == 1:
                        analysis_use = use
                    else:
                        analysis_use = analysis_use * use
                    mem = mem + 1
                    # End of initial read of all ensemble diag files.

                if bbreak:
                    bbreak = False  # reset break
                    break  # break and move to next cycle time.

                mem = 1
                while mem <= n_mem:
                    memid = str(mem).zfill(4)
                    if ob_type == "u" or ob_type == "v":
                        diagfile = _os.path.join(datapath, f"{expt_name}/{date}/mem{memid}/diag_conv_uv_ges.{date}.nc4")
                    else:
                        diagfile = _os.path.join(datapath, f"{expt_name}/{date}/mem{memid}/diag_conv_{ob_type}_ges.{date}.nc4")
                    nc = Dataset(diagfile)
                    if mem == 1:
                        use = analysis_use
                        code = nc["Observation_Type"][:]
                        lat = nc["Latitude"][:]
                        lon = nc["Longitude"][:]
                        pressure = nc["Pressure"][:]
                        errorinv = nc["Errinv_Final"][:]

                        if ob_type == "u":
                            ob = nc["u_Observation"][:]
                        elif ob_type == "v":
                            ob = nc["v_Observation"][:]
                        else:
                            ob = nc["Observation"][:]

                        if ob_type == "q":
                            ob = 1000.0 * ob  # convert from kg/kg to g/kg
                            errorinv = errorinv / 1000.0  # convert from kg/kg to g/kg

                        # consider where use flag==1 and bound by error/lat/lon/pres
                        used = _filter_obs.filter_obs(
                            code,
                            codes,
                            errorinv,
                            lat,
                            lon,
                            pressure,
                            use,
                            hem,
                            p_max=p_max,
                            p_min=p_min,
                            lat_max=lat_max,
                            lat_min=lat_min,
                            lon_max=lon_max,
                            lon_min=lon_min,
                            error_max=error_max,
                            error_min=error_min,
                        )

                        # Digitize itot (do this before "pressure = pressure[used]")
                        pindx = _np.digitize(pressure, pbins) - 1
                        itot, bin_edges = _np.histogram(pressure, pbins[::-1])

                        # Filter ob, pressure
                        ob = ob[used]
                        pressure = pressure[used]
                        errorinv = errorinv[used]
                        error = 1.0 / errorinv

                        # Digitize iasm (do this after "pressure = pressure[used]")
                        pindx = _np.digitize(pressure, pbins) - 1
                        iasm, bin_edges = _np.histogram(pressure, pbins[::-1])
                        # end if mem==1

                    if ob_type == "u":
                        omf = nc["u_Obs_Minus_Forecast_adjusted"][:]
                    elif ob_type == "v":
                        omf = nc["v_Obs_Minus_Forecast_adjusted"][:]
                    else:
                        omf = nc["Obs_Minus_Forecast_adjusted"][:]
                    if ob_type == "q":
                        omf = 1000.0 * omf  # convert from kg/kg to g/kg

                    omf = omf[used]

                    if mem == 1:
                        fcst_ens = ob - omf
                        fcst_ens_var = (ob - omf) ** 2
                    else:
                        fcst_ens = fcst_ens + ob - omf
                        fcst_ens_var = fcst_ens_var + (ob - omf) ** 2
                    mem = mem + 1
                    # end while n_mem

                # if we broke out of the loop above before, because we were missing a mem.
                # this will break out of the loop here too.
                try:
                    error_var = error**2
                except NameError:
                    break

                fcst_ens_mean = fcst_ens / n_mem
                if n_mem > 1:
                    fcst_ens_var = (fcst_ens_var - n_mem * fcst_ens_mean**2) / (n_mem - 1)
                else:
                    fcst_ens_var = 0.0
                innov = ob - fcst_ens_mean
                innovsq = innov**2
                num_obs_total[i_o, i_e, i_date] = num_obs_total[i_o, i_e, i_date] + itot[::-1]
                num_obs_assim[i_o, i_e, i_date] = num_obs_assim[i_o, i_e, i_date] + iasm[::-1]

                # Digitize vars
                innov = _np.bincount(pindx, minlength=n_levs, weights=innov)
                innovsq = _np.bincount(pindx, minlength=n_levs, weights=innovsq)
                error_var = _np.bincount(pindx, minlength=n_levs, weights=error_var)
                fcst_ens_var = _np.bincount(pindx, minlength=n_levs, weights=fcst_ens_var)

                sum_innov[i_o, i_e, i_date] = sum_innov[i_o, i_e, i_date] + innov
                sum_innovsq[i_o, i_e, i_date] = sum_innovsq[i_o, i_e, i_date] + innovsq
                sum_fcst_ens_var[i_o, i_e, i_date] = sum_fcst_ens_var[i_o, i_e, i_date] + fcst_ens_var
                sum_ob_err_var[i_o, i_e, i_date] = sum_ob_err_var[i_o, i_e, i_date] + error_var
                # end n_times

                if _np.sum(num_obs_assim[i_o, i_e, i_date]) > 0:
                    mean_innov = sum_innov[i_o, i_e, i_date] / num_obs_assim[i_o, i_e, i_date]
                    bias[i_o, i_e, i_date] = -1 * mean_innov
                    rms[i_o, i_e, i_date] = _np.sqrt(sum_innovsq[i_o, i_e, i_date] / num_obs_assim[i_o, i_e, i_date])
                    mean_ob_err_var = sum_ob_err_var[i_o, i_e, i_date] / num_obs_assim[i_o, i_e, i_date]
                    ob_error[i_o, i_e, i_date] = _np.sqrt(mean_ob_err_var)

                if _np.sum(num_obs_assim[i_o, i_e, i_date]) > 1:
                    innov_var = (sum_innovsq[i_o, i_e, i_date] - num_obs_assim[i_o, i_e, i_date] * mean_innov**2) / (num_obs_assim[i_o, i_e, i_date] - 1.0)
                    std_dev[i_o, i_e, i_date] = _np.sqrt(innov_var)
                    mean_fcst_var = sum_fcst_ens_var[i_o, i_e, i_date] / num_obs_assim[i_o, i_e, i_date]
                    spread[i_o, i_e, i_date] = _np.sqrt(mean_fcst_var)
                    total_spread[i_o, i_e, i_date] = _np.sqrt(mean_ob_err_var + mean_fcst_var)
                    cr[i_o, i_e, i_date] = (total_spread[i_o, i_e, i_date] / rms[i_o, i_e, i_date]) ** 2
                    ser[i_o, i_e, i_date] = spread[i_o, i_e, i_date] / rms[i_o, i_e, i_date]

                del errorinv
                del error
                del itot
                del ob
                del iasm
                del omf
                del code
                del lat
                del lon
                del pressure
                del use

            # end do n_expt
        # end do n_ob_type

    # Sum over all cycle times (axis=2)
    num_obs_total = _np.sum(num_obs_total, axis=2)
    num_obs_assim = _np.sum(num_obs_assim, axis=2)
    # Average over all cycle times (axis=2)
    n_dates = n_dates - n_skip/n_ob_type  # re-calc n_dates, we may have needed to skip some dates.
    bias = _np.sum(bias, axis=2) / n_dates
    rms = _np.sum(rms, axis=2) / n_dates
    ob_error = _np.sum(ob_error, axis=2) / n_dates
    std_dev = _np.sum(std_dev, axis=2) / n_dates
    spread = _np.sum(spread, axis=2) / n_dates
    total_spread = _np.sum(total_spread, axis=2) / n_dates
    cr = _np.sum(cr, axis=2) / n_dates
    ser = _np.sum(ser, axis=2) / n_dates

    return levs, levs2, levs1, dates, bias, rms, std_dev, spread, ob_error, total_spread, num_obs_total, num_obs_assim, cr, ser
