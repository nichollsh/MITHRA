import glob
import os
import gc
import subprocess
import xarray as xr
import numpy as np
import scipy.interpolate as interp
from scipy.integrate import cumulative_trapezoid

import mithra.utils as utils

def scrape_npy(url:str):
    '''
    Scrape a stellar spectrum from the Centro de Astrobiología website.
    
    This spectrum is saved to an npy file in the data folder of this repository. 
    The npy file contains a 2D array. The first element is the wavelength array,
    and the second is the spectral flux array [erg s-1 cm-2 nm-1].

    Parameters
        url (str): url from which to download the file

    Returns
        None
    '''

    # download file 
    tmp = "/tmp/btsettl"
    utils.rmsafe(tmp)
    subprocess.run(["wget","-q","-O",tmp,url])

    # get name 
    with open(tmp,'r') as hdl:
        lines = hdl.readlines()
    teff = float(lines[1].split(" ")[3])
    logg = float(lines[2].split(" ")[3]) * 100.0
    out = os.path.join(utils.dirs["data"],"btsettl-cifist_%04d_%04d"%(teff,logg))
    utils.rmsafe(out)

    wl = []
    fl = []
    for line in lines[9:]:
        s = line.split()
        if len(s) == 2:
            wl.append(float(s[0]))
            fl.append(float(s[1]))
        else:
            break

    # skip lines, to reduce file size
    pitch = 2
        
    # save data
    wl = np.array(wl[::pitch])*0.1   # convert to nm
    fl = np.array(fl[::pitch])/0.1   # convert to erg s-1 cm-2 nm-1
    np.save(out, np.array([wl,fl]))

    del wl
    del fl
    del lines
    return 

def scrape_all():
    '''
    Iteratively download all stellar spectra from the BT-SETTL CIFIST grid from
    the Centro de Astrobiología website. Saved to the data folder.

    Parameters
        None

    Returns
        None
    '''

    for i in range(1,447,1):
        print("Downloading %03d..."%i)
        url = "http://svo2.cab.inta-csic.es/theory/newov2/ssap.php?model=bt-settl-cifist&fid=%d&format=ascii"%i
        scrape_npy(url)
    return 

def download_tar():
    '''
    Download BT-SETTL CIFIST grid from OSF page.

    This tar.gz file was created using the scrape functions above, and contains 
    all of the npy files containing wavelengths and fluxes. This tar.gz is
    unpackaged in the data folder of this repository.

    Parameters
        None

    Returns
        None
    '''

    # Setup
    out = "/tmp/btsettl.tar.gz"
    utils.rmsafe(out)
    cmd = ["osf","-p", "8r2sw", "fetch", "BTSETTL-CIFIST/npys.tar.gz",out]
    env = os.environ.copy()

    # Download
    print("Downloading tar.gz archive...")
    subprocess.run(cmd,env=env,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
    print("    Done")

    # Unzip
    print("Extracting archive...")
    print("    to: '%s'"%utils.dirs["data"])
    cmd = ["tar", "-xf", out, "-C", utils.dirs["data"]]
    subprocess.run(cmd)

    # Remove tar 
    cmd = ["rm",out]
    subprocess.run(cmd)
    print("    done")

    return 


def get_params_from_name(fpath:str):
    '''
    Get model gridpoint parameters (Teff, logg) from its file name.

    Parameters
        fpath (str): path to npy file

    Returns
        teff (float): effective temperature [K]
        logg (float): log surface gravity [log cm/s^2]
    '''

    name = fpath.split("/")[-1].split(".")[0]
    splt = name.split("_")
    teff = float(splt[-2])
    logg = float(splt[-1])/100.0
    return teff,logg

def list_files():
    '''
    List all BT-SETTL npy files in the data directory

    Parameters
        None

    Returns
        files (list): list of absolute file paths
    '''
    files = glob.glob(os.path.join(utils.dirs["data"], "btsettl-cifist*.npy")) 
    return list(files)

def integrate_spectrum(wl,fl):
    '''
    Integrate a spectrum over wavelength.

    Returns the cumulative integral at each point, starting from zero.

    Parameters
        wl (np.ndarray): wavelengths [nm]
        fl (np.ndarray): fluxes [erg s-1 cm-2 nm-1]

    Returns
        fi (np.ndarray): integrated fluxes [erg s-1 cm-2]
    '''

    fi = cumulative_trapezoid(fl,wl, initial=0.0)
    return fi

def create_interp(num_teff=0, num_logg=0, num_wave=40, teff_lims=(1.0, 2e5), logg_lims=(1.0, 50.0)):
    '''
    Create interpolated grid of teff, logg, and wavelength from npy files.

    Reads all npy files in the data folder, creating a single data structure of 
    flux versus wavelength, teff, and logg. This grid may have missing values.
    The data are interpolated within the same (or reduced) teff and logg grid, 
    so no extrapolation is performed but missing values are filled.

    Parameters
        num_teff (int): number of interpolated teff points (0 for same as input)
        num_logg (int): number of interpolated logg points (0 for same as input)
        num_wave (int): number of interpolatedwavelength points
        teff_lims (tuple): minimum and maximum teff values to use
        logg_lims (tuple): minimum and maximum logg values to use

    Returns
        vi (np.ndarray): array of fluxes in a 3D array (teff, logg, wavelength)
        xi (np.ndarray): array of teff [K] on same axes as above
        yi (np.ndarray): array of logg [log cm/s^2] on same axes as above
        zi (np.ndarray): array of wavelength [nm] on same axes as above
    '''

    # get wavelength data from first file 
    data = np.load(list_files()[0])
    min_wave, max_wave = np.amin(data[0]), np.amax(data[1])

    # limit wl range
    max_wave = min(max_wave, 1e5)  # 100 um
    min_wave = max(min_wave, 1.0)  # 1 nm
    target_wave = np.linspace(np.log10(min_wave+1), np.log10(max_wave-1), num_wave)
    len_ds = int(num_wave * 2)  # length of downsampled spectrum

    # flattened data from files
    flat_teff = []
    flat_logg = []
    flat_wave = []
    flat_flux = []
    print("Reading npy files...")
    for i,f in enumerate(list_files()):
        # get header
        teff,logg = get_params_from_name(f)

        # skip this point?
        if not(teff_lims[0] <= teff <= teff_lims[1]):
            continue 
        if not(logg_lims[0] <= logg <= logg_lims[1]):
            continue

        # load data
        data = np.load(f)
        w,f = np.log10(data[0]),data[1]

        # drop duplicate values and ensure sorted
        _,mask = np.unique(w, return_index=True)
        w = w[mask]
        f = f[mask]

        # interpolate (downsample) flux and wavelength arrays...
        # this is to ensure that all spectra cover the same wavelength range
        # and also to significantly improve the performance of the unstructured
        # interpolation step by reducing the number of points
        w_ds = np.linspace(w[0],w[-1], len_ds)
        f_ds = interp.pchip_interpolate(w,f,w_ds)

        # store
        flat_wave.extend(list(w_ds))
        flat_flux.extend(list(f_ds))
        flat_teff.extend(list(np.ones(len(w_ds))*teff))
        flat_logg.extend(list(np.ones(len(w_ds))*logg))

        del w_ds 
        del f_ds 
        del w 
        del f 
        del data 
        gc.collect()
    print("    done")

    # unique
    uniq_teff = np.unique(flat_teff)
    uniq_logg = np.unique(flat_logg)
    print("Source axes: (teff, logg, wave) = (%d, %d, %d)"%(len(uniq_teff), len(uniq_logg), num_wave))

    # output sizes
    if num_teff <= 1:
        num_teff = len(uniq_teff)
    if num_logg <= 1:
        num_logg = len(uniq_logg)

    # target grid to interpolate to
    xi = np.linspace(np.amin(uniq_teff), np.amax(uniq_teff), num_teff)
    yi = np.linspace(np.amin(uniq_logg), np.amax(uniq_logg), num_logg)
    zi = target_wave
    print("Interpolation target: (teff, logg, wave) = (%d, %d, %d)"%(len(xi), len(yi), len(zi)))

    print("Meshgrid",flush=True)
    xi,yi,zi = np.meshgrid(xi,yi,zi, indexing='ij')

    # interpolate
    print("Interpolating...")
    print("    please wait")
    vi = interp.griddata((flat_teff,flat_logg,flat_wave),flat_flux,(xi,yi,zi),method='linear')

    print("    done")
    return (vi, xi, yi, 10.0**zi)


def create_dataset(itp:tuple):
    '''
    Write interpolated grid to NetCDF file using Xarray.

    Parameters
       itp (tuple): tuple of four 3D output arrays from create_interp()

    Returns
        ds (xr.Dataset): dataset of interpolated data
        
    '''

    # get coords
    v,x,y,z = itp[0],itp[1],itp[2],itp[3]

    arr_teff = np.unique(x)
    arr_logg = np.unique(y)
    arr_wave = z[0,0,:]

    coords={
            "teff": arr_teff,
            "logg": arr_logg,
            "wave": arr_wave,
            }
    
    
    # get data
    data = {"flux":xr.DataArray(data=v, coords=coords, dims=["teff","logg","wave"])}
    
    # create dataset
    ds = xr.Dataset(data_vars=data, coords=coords)

    return ds

def write_dataset(ds:xr.Dataset, fp:str):
    '''
    Write interpolated grid to NetCDF file using Xarray.

    Parameters
       ds (xr.Dataset): dataset of interpolated data
       fp (str): path to nc file

    Returns
        None
        
    '''

    # save
    fp = os.path.abspath(fp)
    utils.rmsafe(fp)
    ds.to_netcdf(fp)
    return 

def read_dataset(fp:str)->xr.Dataset:
    '''
    Read interpolated grid from NetCDF file using Xarray

    Parameters
       fp (str): path to nc file

    Returns
        out (xr.Dataset): dataset of interpolated data
    '''

    print("Reading dataset...")
    if not os.path.exists(fp):
        raise FileNotFoundError(fp)
    with xr.open_dataset(fp) as ds:
        out = ds.copy(deep=True)
        print("    done")
    return out


def get_spec_from_dataset(ds:xr.Dataset, teff:float, logg:float):
    '''
    Get stellar spectrum from interpolated grid, searching for best point.

    Parameters
       ds (xr.Dataset): dataset containing interpolated spectra
       teff (float): target teff
       logg (float): target logg

    Returns
        wl (np.ndarray): wavelengths [nm]
        fl (np.ndarray): spectral fluxes [erg s-1 cm-2 nm-1]
        teff (float): optimal teff
        logg (logg): optimal logg
    '''

    close = ds.sel(teff=teff, logg=logg, method='nearest')
    wl   = np.array(close.coords["wave"].values)
    fl   = np.array(close["flux"].values)
    teff = float(close.coords["teff"].values)
    logg = float(close.coords["logg"].values)

    return wl, fl, teff, logg

def get_axes():
    '''
    Get coordinates of logg and teff axes from files in data folder.

    Parameters
        None

    Returns
        arr_teff (np.ndarray): teff points
        arr_logg (np.ndarray): logg points
    '''

    # get all files
    files = list_files()
    names = [f.split("/")[-1].split(".")[0] for f in files]

    # get all teff,logg
    arr_teff = []
    arr_logg = []
    for f in names:
        teff,logg = get_params_from_name(f)
        arr_teff.append(teff)
        arr_logg.append(logg)
    arr_teff = np.unique(arr_teff)
    arr_logg = np.unique(arr_logg)

    return arr_teff, arr_logg

def get_spec_from_npy(teff:float, logg:float, xmax:float=1.0e9):
    '''
    Read stellar spectrum from npy file, locating it according to teff and logg.

    Parameters
        teff (float): effective temperature [K]
        logg (float): log surface gravity [log cm/s^2]

    Returns
        wl (np.ndarray): wavelengths [nm]
        fl (np.ndarray): spectral fluxes [erg s-1 cm-2 nm-1]
    '''

    # get axes 
    arr_teff, arr_logg = get_axes()
    
    # get best indices
    grid_i = np.argmin(np.abs(arr_teff - teff))
    grid_j = np.argmin(np.abs(arr_logg - logg))
    grid_t = arr_teff[grid_i]
    grid_l = arr_logg[grid_j]

    # get data 
    best = os.path.join(utils.dirs["data"], "btsettl-cifist_%04d_%04d.npy"%(grid_t,grid_l*100.0))
    data = np.load(best)
    wl = data[0]
    fl = data[1]

    # truncate?
    if wl[-1] > xmax:
        imax = np.argmin(np.abs(wl - xmax))
        wl = wl[:imax]
        fl = fl[:imax]

    # return spectrum
    return wl,fl


def write_csv(fp:str, wl:float,fl:float):
    '''
    Write stellar spectrum to csv file.

    Parameters
        fp (string): path to output file
        wl (np.ndarray): wavelengths [nm]
        fl (np.ndarray): spectral fluxes [erg s-1 cm-2 nm-1]

    Returns
        None
    '''

    fp = os.path.abspath(fp)
    utils.rmsafe(fp)

    X = np.array([wl,fl]).T
    head = "Wavelength [nm] , Flux [erg s-1 cm-2 nm-1]"
    np.savetxt(fp, X, fmt="%.5e", delimiter=',',header=head)
    return 


def extend_planck(teff, wl, fl, xmax=1.0e9):
    dx = wl[-1]*0.1
    print(dx)
    xp = np.arange(wl[-1]+dx, xmax, dx)
    yp = np.zeros(np.shape(xp))
    for i,x in enumerate(xp):
        lam = x * 1.0e-9  # nm -> m

        # Calculate planck function value [W m-2 sr-1 m-1]
        # http://spiff.rit.edu/classes/phys317/lectures/planck.html
        yp[i] = 2.0 * utils.h_pl * utils.c_vac**2.0 / lam**5.0   *   1.0 / ( np.exp(utils.h_pl * utils.c_vac / (lam * utils.k_B * teff)) - 1.0)

        # Integrate solid angle (hemisphere), convert units
        yp[i] = yp[i] * np.pi * 1.0e-9 # [W m-2 nm-1]
        yp[i] = yp[i] * 1000.0 # [erg s-1 cm-2 nm-1]

    wl = np.concatenate((wl, xp))
    fl = np.concatenate((fl, yp))

    return wl,fl
