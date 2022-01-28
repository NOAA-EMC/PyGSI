import numpy as np
import pandas as pd


def land_fraction(self, df):
    """
    Returns dataframe for only land data.
    """
    df = df[df['land_fraction'] > 0]

    return df


def water_fraction(self, df):
    """
    Returns dataframe for only water data.
    """
    df = df[df['water_fraction'] > 0]

    return df


def cloud_fraction(self, df):
    """
    Returns dataframe for only cloud data.
    """
    df = df[df['cloud_fraction'] > 0]

    return df


def vegetation_fraction(self, df):
    """
    Returns dataframe for only vegetation data.
    """
    df = df[df['vegetation_fraction'] > 0]

    return df


def ice_fraction(self, df):
    """
    Returns dataframe for only ice data.
    """
    df = df[df['ice_fraction'] > 0]

    return df


def snow_fraction(self, df):
    """
    Returns dataframe for only snow data.
    """
    df = df[df['snow_fraction'] > 0]

    return df
