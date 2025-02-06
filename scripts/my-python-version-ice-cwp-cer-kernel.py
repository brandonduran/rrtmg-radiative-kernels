#!/usr/bin/env python
# coding: utf-8

# RRTM code expects arrays with (ncol, nlay)
# and with pressure decreasing from surface at element 0

#import pathlib
#print(pathlib.Path(__file__).parent.resolve())
import sys
from helper_function import fill_t_q_profile, compute_cosSZA
#import scripts.helper_functions as HF
sys.path.append("/home/bmduran/climlab-rrtmg/climlab-rrtmg")  # Update with the actual path to the climlab-rrtmg repository
from climlab_rrtmg import rrtmg_sw
import resource
import psutil
import os
import time
import numpy as np
import mat73

datadir = "../data/"

def get_memory_usage():
    process = psutil.Process(os.getpid())
    return f"Memory usage: {process.memory_info().rss / (1024 * 1024):.2f} MB"

#def get_memory_usage():
#    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss  # In KB
#    return f"Memory usage: {usage / 1024:.2f} MB"

# Define bins for ice water path (IWP) and ice-crystal effective radius (Reff)
IWP_bins_MODIS = np.array([0, 20, 50, 100, 200, 400, 1000, 20000])  # g m^-2
Reff_bins_MODIS = np.array([5, 10, 20, 30, 40, 50, 60, 90])  # um

# Calculate midpoints for IWP and Reff bins
IWP_mdpts_MODIS = (IWP_bins_MODIS[:-1] + IWP_bins_MODIS[1:]) / 2
Reff_mdpts_MODIS = (Reff_bins_MODIS[:-1] + Reff_bins_MODIS[1:]) / 2

# Load seasonal climatology of temperature and humidity from GCM simulations
#data = mat73.loadmat('/home/cawall/scripts/scripts_for_Brandon/data/CTL_seasonal_climatology_zonal_mean.mat')
data = mat73.loadmat(datadir+'CTL_seasonal_climatology_zonal_mean.mat')

lat = data['lat'].squeeze()  # Extract latitude
months = data['months'].squeeze()  # Extract months
p = data['p'].squeeze()  # Pressure levels

# Compute multi-model mean for q, t, and ts
q = np.mean(data['q'], axis=3)  # Shape reduction along the 4th dimension
t = np.mean(data['t'], axis=3)  # Shape reduction along the 4th dimension
ts = np.mean(data['ts'], axis=2)  # Shape reduction along the 3rd dimension

# Surface albedo values from Zelinka et al. 2012
albcs = np.array([0, 0.5, 1])

# Save data for loops (assigning to new variables)
t_save = t.copy()
q_save = q.copy()
ts_save = ts.copy()
lat_save = lat.copy()
months_save = months.copy()

def read_atm(which_atm):
    """
    Reads a standard atmospheric profile in the RRTM TAPE5 format.
    
    Parameters:
        which_atm (int): Specifies the atmospheric profile to use.
            1 - Tropical
            2 - Midlatitude Summer
            3 - Midlatitude Winter
            4 - Subarctic Summer (Not Supported)
            5 - Subarctic Winter
            6 - US Standard, 1976 (Not Supported)
    
    Returns:
        plev, play, Tlev, Tlay, Tsfc, ngas, gasvmr, semis, tag
    """
    # Directory path
    #dir_path = '/home/bmduran/RRTMG/'
    #dir_path = datadir
    
    # File mapping based on atmospheric type
    atmos_files = {
        1: ("runs_std_atm_lw/input_rrtm_TROP-clr", "TROP"),
        2: ("runs_std_atm_lw/input_rrtm_MLS-clr", "MLS"),
        3: ("runs_std_atm_lw/input_rrtm_MLW-clr", "MLW"),
        4: None,  # Subarctic Summer not supported
        5: ("runs_std_atm_lw/input_rrtm_SAW-clr", "SAW"),
        6: None   # US Standard not supported
    }

    if which_atm not in atmos_files or atmos_files[which_atm] is None:
        raise ValueError(f"Invalid or unsupported atmosphere type: {which_atm}")

    file_name, tag = atmos_files[which_atm]
    #file_path = dir_path + file_name
    file_path = datadir + file_name

    try:
        with open(file_path, 'r') as f:
            # Skip header until the "$" character is found
            while True:
                line = f.readline().strip()
                if line.startswith('$'):
                    atm_type = line[1:].strip()
                    print(f"Using trace gas sounding from {atm_type}")
                    break

            # Skip next line with run parameters
            f.readline()
            
            # Surface temperature
            Tsfc = float(f.readline().strip().split()[0])

            # Atmospheric layers and gas information
            a = list(map(int, f.readline().strip().split()[:3]))
            nlay, ngas = a[1], a[2]

            play = np.zeros(nlay)
            Tlay = np.zeros(nlay)
            plev = np.zeros(nlay + 1)
            Tlev = np.zeros(nlay + 1)
            gasvmr = np.zeros((ngas, nlay))

            # Read first layer
            line = list(map(float, f.readline().strip().split()[:9]))
            play[0], Tlay[0] = line[0], line[1]
            plev[0], Tlev[0] = line[4], line[5]
            plev[1], Tlev[1] = line[7], line[8]
            
            zlay = [line[6]]
            
            # Read gas volume mixing ratios
            line = list(map(float, f.readline().strip().split()))
            gasvmr[:, 0] = line[:ngas]

            num_molec = gasvmr[0, 0] > 1
            if num_molec:
                total_num = sum(line)
                gasvmr[:, 0] /= total_num

            # Read remaining layers
            for n in range(1, nlay):
                line = f.readline().strip()
                if len(list(line.split())) == 6:
                    a = list(map(float, line.split()))
                else:
                    a = list(map(float, line.split()[:3]))
                    a.extend([float(line.split()[3][:6]), float(line.split()[3][6:]), float(line.split()[4])]) #need to split index 3 , first 6 numbers, last 8 numbers
                play[n], Tlay[n] = a[0], a[1]
                plev[n + 1], Tlev[n + 1] = a[4], a[5]
                zlay.append(a[3])

                line = list(map(float, f.readline().strip().split()))
                gasvmr[:, n] = line[:ngas]
                if num_molec:
                    total_num = sum(line)
                    gasvmr[:, n] /= total_num

            semis = 1
            return plev, play, Tlev, Tlay, gasvmr

    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {file_path}")

print('before assigning p arrays')
#print(get_memory_usage())  # Call this at different points in your script
# --------- RRTMG parameters that don't change ---------

# Add one level for the surface
p = np.append(p, 1012.5)

# Switch from pressure-increasing indexing to pressure-decreasing indexing
p = np.flip(p)

# Layer midpoint pressure (hPa)
play = np.array([1.00625e+03, 9.62500e+02, 8.87500e+02, 8.25000e+02, 7.87500e+02,
        7.37500e+02, 6.90000e+02, 6.40000e+02, 5.80000e+02, 5.30000e+02,
        4.70000e+02, 4.20000e+02, 3.55000e+02, 3.05000e+02, 2.75000e+02,
        2.25000e+02, 1.90000e+02, 1.65000e+02, 1.25000e+02, 8.50000e+01,
        6.00000e+01, 4.00000e+01, 2.50000e+01, 1.50000e+01, 8.50000e+00,
        6.00000e+00, 4.00000e+00, 2.50000e+00, 1.50000e+00, 5.00000e-01])
play = play[np.newaxis, ...]

# Interface pressure (hPa)
plev = p
plev = plev[np.newaxis, ...]
#print(plev)
# Number of layers
nlay = len(play[0])

# Get GHG soundings
plev_std, play_std, tlev_std, tlay_std, gasvmr_std = read_atm(2)

# CO2, CH4, N2O mixing ratios from Zelinka et al. 2012
co2vmr = 330e-6 * np.ones(len(play[0]))
#print(co2vmr)
co2vmr = co2vmr[np.newaxis,...]
#print(co2vmr)
ch4vmr = 1.6e-6 * np.ones(len(play[0]))
ch4vmr = ch4vmr[np.newaxis,...]
n2ovmr = 0.28e-6 * np.ones(len(play[0]))
n2ovmr = n2ovmr[np.newaxis,...]

# Interpolation for other gas volume mixing ratios
o3vmr = np.interp(play, play_std, gasvmr_std[2, :])
covmr = np.interp(play, play_std, gasvmr_std[4, :])
o2vmr = np.interp(play, play_std, gasvmr_std[6, :])

# Molecular weights for gases
mwdry = 28.966
mwh2o = 18.016

# Other parameters
ncol = 1
icld = 1  # Cloud overlap assumption: 1 = random overlap
adjes = 1  # Use constant Earth-Sun distance (1 AU)
dyofyr = 1  # Dummy variable
scon = 1366  # Solar constant (W/m²) from Zelinka et al. 2012

# Cloud and aerosol flags
inflgsw = 2  # 1 or 2: input cloud physical properties; 0: input optical properties
iceflgsw = 3  # Ice flag: generalized effective particle size (Dge) from 1996
liqflgsw = 1  # 1: liquid clouds active, 0: inactive
iaer = 0  # Aerosols not active

# Start timer
start_time = time.time()

# RRTMG parameters
nbndsw = 14  # Number of shortwave bands
ngptsw = 112 # number of g-points
naerec = 6 # number of ECMWF aerosol types for ecaer input
nbndlw = 16
ngptlw = 140
#naerec = 5   # Number of aerosol records
isolvar = -1
bndsolvar = np.ones(nbndsw) # Solar variability scale factors for each shortwave band
indsolvar = np.ones(2) # Facular and sunspot amplitude scale factors (isolvar=1),
                                 # or Mg and SB indices (isolvar=2)
solcycfrac = 1. # Fraction of averaged solar cycle (0-1) at current time (isolvar=1)

#tauc, ssac, asmc, fsfc
# Initialize arrays with zeros
tauc = np.zeros((nbndsw, ncol, nlay))
ssac = np.zeros((nbndsw, ncol, nlay))
fsfc = np.zeros((nbndsw, ncol, nlay))
asmc = np.zeros((nbndsw, ncol, nlay))
tauaer = np.zeros((ncol, nlay, nbndsw))
ssaaer = np.zeros((ncol, nlay, nbndsw))
asmaer = np.zeros((ncol, nlay, nbndsw))
ecaer = np.zeros((ncol, nlay, naerec))

# Cloud-top pressure (CTP) assumption (this can be adjusted)
CTP = 250 # hPa

# Find the level that equals CTP
if CTP >= np.max(plev):
    CTP_ind = np.argmax(plev)  # Index of the maximum pressure level
else:
    CTP_ind = np.argmax(plev >= CTP) - 1  # Find level equal to or just above CTP

# Start timer
start_time = time.time()
print('before loop and before kernel initialized')

# Initialize array for kernel
SW_kernel = np.full(
    (
        len(albcs),
        len(lat_save),
        len(Reff_mdpts_MODIS),
        len(IWP_mdpts_MODIS),
        len(months_save),
    ),
    np.nan,
)

for albedo_ind in range(len(albcs)):
    for lat_ind in range(len(lat_save)):
        for month_ind in range(len(months_save)):
            # Extract values for the current latitude and month
            t = t_save[lat_ind, :, month_ind]
            q = q_save[lat_ind, :, month_ind]
            ts = ts_save[lat_ind, month_ind]
            lat = lat_save[lat_ind]
            months = int(months_save[month_ind])
            print(months)
            asdir = aldir = asdif = aldif = albcs[albedo_ind]

            # Calculate day of the year
            if months == 2:
                day_of_yr = (
                    np.datetime64(f"2000-{months:02d}-14")
                    - np.datetime64("2000-01-01")
                ).astype("timedelta64[D]").astype(int) + 1
            else:
                day_of_yr = (
                    np.datetime64(f"2000-{months:02d}-15")
                    - np.datetime64("2000-01-01")
                ).astype("timedelta64[D]").astype(int) + 1

            hrs = np.arange(0, 25)
            cosSZA_mean = []
            for n in range(len(hrs) - 1):
                dt = 0.01
                day_num = day_of_yr + np.arange(hrs[n], hrs[n + 1], dt) / 24
                cosSZA = compute_cosSZA(lat, 0, day_num)
                cosSZA[cosSZA < 0] = 0
                cosSZA_mean.append(np.nanmean(cosSZA))
            cosSZA_mean = np.array(cosSZA_mean)

            # Add one level for the surface
            t = np.append(t, ts)
            q = np.append(q, q[-1])

            # Switch from pressure-increasing to pressure-decreasing indexing
            q = np.flip(q)
            t = np.flip(t)

            # Fill missing values
            if np.any(np.isnan(t)):
                t = fill_t_q_profile(t, p)
            if np.any(np.isnan(q)):
                q = fill_t_q_profile(q, p)

            # RRTMG settings
            tlay = (t[:-1] + t[1:]) / 2  # Layer midpoint temperature (K)
            tlay = tlay[np.newaxis,...]
            tlev = t  # Interface temperature (K)
            tlev = tlev[np.newaxis,...]
            tsfc = ts  # Surface temperature (K)
            h2ovmr = (mwdry / mwh2o) * (q[:-1] + q[1:]) / 2  # Water vapor mole fraction
            h2ovmr = h2ovmr[np.newaxis,...]

            # Calculate SW CRE
            SWCRE = np.full((len(Reff_bins_MODIS), len(IWP_bins_MODIS)), np.nan)

            for n in range(len(Reff_bins_MODIS)):
                for m in range(len(IWP_bins_MODIS)):
                    Reff = Reff_bins_MODIS[n]  # um
                    IWP = IWP_bins_MODIS[m]  # g m-2

                    Dge = 8 / (3 * np.sqrt(3)) * Reff # generalized effective particle size (Fu 1996)

                    ciwp = np.zeros_like(play)
                    clwp = np.zeros_like(play)
                    reic = np.zeros_like(play)
                    relq = np.zeros_like(play)
                    cldfrac = np.zeros_like(play)
                    ciwp[CTP_ind] = IWP

                    if iceflgsw == 3:
                        reic[CTP_ind] = Dge
                    else:
                        reic[CTP_ind] = Reff

                    #cldfr[cicewp > 0] = 1
                    cldfrac[ciwp > 0] = 1

                    # Radiative transfer calculations over SZA
                    uflx_sw_TOA_save = []
                    uflxc_sw_TOA_save = []
                    permuteseed = 1 
                    irng = 1
                    for coszen in cosSZA_mean:
                        #  Call the Monte Carlo Independent Column Approximation (McICA, Pincus et al., JC, 2003)
                        (cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl,
                        ssacmcl, asmcmcl, fsfcmcl) = rrtmg_sw.climlab_mcica_subcol_sw(
                                        ncol, nlay, icld, permuteseed, irng, play,
                                        cldfrac, ciwp, clwp, reic, relq, tauc, ssac, asmc, fsfc)
                        
                        # a few different arguments seem to be required for the Python version
                        # isolvar (flag for solar variability method)
                        # indsolvar
                        # bndsolvar
                        # solcycfrac
                        
                        #  Call the RRTMG_SW driver to compute radiative fluxes
                        uflx_sw, _, _, uflxc_sw, _, _ = rrtmg_sw.climlab_rrtmg_sw(ncol,nlay,icld,
                            iaer, play, plev, tlay, tlev, tsfc,
                            h2ovmr,o3vmr,co2vmr,ch4vmr,n2ovmr,o2vmr,
                            asdir,asdif,aldir,aldif,
                            coszen,adjes,dyofyr,scon,isolvar,
                            inflgsw,iceflgsw,liqflgsw,cldfmcl,#also called cldfmcl?
                            #taucld,ssacld,asmcld,
                            taucmcl, ssacmcl, asmcmcl,fsfcmcl,
                            ciwpmcl, clwpmcl, reicmcl,relqmcl,
                            #fsfcld,cicewp,cliqwp,reice,reliq,
                            tauaer,ssaaer,asmaer,ecaer,
                            bndsolvar, indsolvar, solcycfrac
                        )
                        uflx_sw_TOA_save.append(uflx_sw[-1])
                        uflxc_sw_TOA_save.append(uflxc_sw[-1])
                    uflx_sw_TOA = np.nanmean(uflx_sw_TOA_save)
                    uflxc_sw_TOA = np.nanmean(uflxc_sw_TOA_save)

                    if IWP == 0: # a cloud doesn't exist
                        SWCRE[n, m] = 0
                    else:
                        SWCRE[n, m] = uflxc_sw_TOA - uflx_sw_TOA

            # Calculate SW kernel
            for n in range(len(Reff_mdpts_MODIS)):
                for m in range(len(IWP_mdpts_MODIS)):
                    SWCREa = SWCRE[n : n + 2, m : m + 2]
                    SW_kernel[albedo_ind, lat_ind, n, m, month_ind] = (
                        np.mean(SWCREa) / 100
                    )

            if lat_ind == 0 and albedo_ind == 0:
                print(month_ind)

        print(lat_ind)

# End timer
end_time = time.time()
print(f"Elapsed time: {end_time - start_time:.2f} seconds")

