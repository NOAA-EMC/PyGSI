# coding: utf-8 -*-

'''
stats.py contains statistics utility functions
'''

__all__ = ['spectrum_stats']

import pandas as _pd


def spectrum_stats(df):

    '''
    Function to calculate major statistics for radiance data as a function of channel.
    Statisics are for observations that pass QC only.

    Args:
        df : (Pandas dataframe) A data frame containing raw radiance information as produced
             by the get_data method of Radiance.

    Returns:
        channel_stats : (Pandas dataframe) A data frame containing channel-by-channel statistics
    '''

    # Get channel numbers in dataframe
    sc = df.index.unique(level="Channel")

    # Initialise output dataframe
    channel_stats = _pd.DataFrame(index=sc,
                                  columns=["count", "omf_unadjusted_mean", "omf_adjusted_mean",
                                           "omf_unadjusted_stddev", "omf_adjusted_stddev"])

    # Sorting by channel increases efficiency and prevents a warning message
    dfs = df.sort_index()

    # Loop through channels and calculate stats
    for chan in sc:
        # Subset based on channel and QC indexes
        tmp = dfs.loc[(chan, 0.0)]
        # This next line makes sure that all rejected obs are avoided.
        # Using the inverse observation error is more reliable than the QC flag.
        tmp = tmp[(tmp['inverse_observation_error'] > 0.0)]
        channel_stats["count"].loc[chan] = len(tmp.axes[0])
        for var in ["omf_unadjusted", "omf_adjusted"]:
            channel_stats[var + '_mean'].loc[chan] = tmp[var].to_numpy().mean()
            channel_stats[var + '_stddev'].loc[chan] = tmp[var].to_numpy().std()

    return channel_stats
