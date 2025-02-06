#################################
## These functions are called in
## the kernel-generating Python
## scripts in this directory.
#################################
import numpy as np


def fill_t_q_profile(x, p):
    """
    Fill missing data in temperature or specific humidity profiles in zonal-mean climatology using the 
    value from the closest pressure level with available data.

    Parameters:
        x (np.array): Vector of temperature or specific humidity (with NaNs for missing data).
        p (np.array): Vector of pressure.

    Returns:
        np.array: Temperature or specific humidity profile with missing data filled in.
    """
    # Find indices of available and missing data
    data_ind = np.where(~np.isnan(x))[0]  # Levels with available data
    no_data_ind = np.where(np.isnan(x))[0]  # Levels with missing data
    
    # Pressure and data at levels with available data
    p_data = p[data_ind]
    x_data = x[data_ind]
    
    # Fill missing data
    for ind in no_data_ind:
        # Find the index of the closest pressure level with data
        closest_ind = np.argmin(np.abs(p_data - p[ind]))
        x[ind] = x_data[closest_ind]
    
    return x
# compute_cosSZA function
def compute_cosSZA(lat, lon, calday):
    """
    Calculate the cosine of the solar zenith angle (cosSZA).
    
    Args:
        lat (float or np.array): Latitude in degrees (positive for north, negative for south).
        lon (float or np.array): Longitude in degrees (positive for east, negative for west).
        calday (float or np.array): Calendar day of the year (1 to 365).
        
    Returns:
        np.array: Cosine of the solar zenith angle (cosSZA).
    """
    # Eccentricity factor
    eccf = 1 + 0.034 * np.cos(2 * np.pi * calday / 365)
    
    # Sine and cosine of latitude
    sin_lat = np.sin(np.pi * lat / 180)
    cos_lat = np.cos(np.pi * lat / 180)
    
    # Sun's latitude in radians
    sun_lat = (np.pi * 23.5 / 180) * np.cos(np.pi * (calday - 172) / 183)
    
    # Sun's longitude in radians
    sun_azim = 2 * np.pi * (calday + lon / 360)
    
    # Cosine of the solar zenith angle
    coszrs = -np.cos(sun_lat) * np.cos(sun_azim) * cos_lat + np.sin(sun_lat) * sin_lat
    
    return coszrs
