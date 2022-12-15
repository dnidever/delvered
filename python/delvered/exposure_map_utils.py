import numpy as np
from matplotlib.path import Path
from astropy.coordinates import SkyCoord
import astropy.units as u
from gala.coordinates import MagellanicStreamNidever08


def decam_array_map(coords1, coords2, exptimes, rmin=-33.0, rmax=30.0, res=0.25):
    '''
    This function creates a 2-dimensional numpy array that represent an exposure map of the given
    coordinates plus exposure times (these can be represented as teff) for a list of DECam exposures.
    Note that the returned array will be square. This method is not suitable for all-sky representations.

    Parameters:
    -----------
    coords1: array
        Contains the first coordinate set of all the exposures that contribute to exposure map (e.g. RA).
    coords2: array
        Contains the second coordinate set of all the exposures that contribute to exposure map (e.g. DEC).
    exptimes: array
        Array of the exposure times of the given coordinates. NaNs are supported.
    rmin: float, default=-33.0
        Minimum value for the map range in both given sky coordinates. The default is optimized to
        fit the DELVE-MC area (Drlica-Wagner et al. 2021) in Magellanic Stream coordinates.
    rmax: float, default=30.0
        Maximum value for the map range in both given sky coordinates. The default is optimized to
        fit the DELVE-MC area (Drlica-Wagner et al. 2021) in Magellanic Stream coordinates.
    res: float, default=0.25
        Resolution of the map in degrees.

    Returns:
    --------
    expmap: array
        2-dimensional array that represents the exposure map of the input.
    '''


    # DECam mask (deg)
    xx = [  -0.45599600,  -0.16197606,  -0.16945115,0.16443590,0.16692760,0.45845584,0.45596415,
            0.61044920,0.60795751,0.76742595,0.76493426,0.92440270,0.92191100, 1.0763961,
            1.0763961,0.92191100,0.91941931,0.76742595,0.76493426,0.61044920,0.61044920,
            0.45845584,0.45347246,  -0.45599600,  -0.45599600,  -0.61297275,  -0.61297275,  -0.76496611,
           -0.76745780,  -0.92194286,  -0.92443455,-1.0764279,-1.0814113,  -0.92194286,  -0.92194286,
           -0.76745780,  -0.76496611,  -0.61048106,  -0.61048106,  -0.45848770]
    yy = [ -0.97753922,  -0.97379390,  -0.80899989,  -0.80525457,  -0.97379390,  -0.97004858,  -0.80525457,
           -0.80525457,  -0.64046057,  -0.64420589,-0.47941188,  -0.47566657,  -0.14982387,  -0.14607856,
            0.15729222,0.16103754,0.49062555,0.49062555,0.64792892,0.65541956,0.81272292,
            0.82021356,0.98126225,0.98126225,0.81646824,0.81646824,0.64792892,0.65167424,
            0.49062555,0.49437087,0.16103754,0.15729222,  -0.14607856,  -0.14982387,  -0.47941188,
           -0.47941188,  -0.64046057,  -0.64795121,  -0.80525457,  -0.81274521]
    decam_vertices = np.array(tuple(zip(xx,yy)))

    final_exposures = np.array(tuple(zip(coords1, coords2)))
    pointings = np.tile(decam_vertices, (len(final_exposures),1,1))
    pointings_mag = np.tile(final_exposures, (1,1,40)).reshape(pointings.shape)
    final_coordinates = pointings + pointings_mag

    pointings_set = []
    for i, pointing in enumerate(final_coordinates):
        pointings_set.append(Path(pointing))

    bins = np.arange(rmin, rmax+res, res)
    points = (bins[1:]+bins[:-1])/2
    expmap = np.zeros((len(points), len(points)))
    xx,yy = np.meshgrid(points,points)
    for i, pointing in enumerate(pointings_set):
        flags = pointing.contains_points(np.hstack((xx.flatten()[:,np.newaxis],yy.flatten()[:,np.newaxis])))
        grid = np.zeros((len(points),len(points)),dtype='bool')
        grid[((xx.flatten()-rmin)/res).astype('int'),((yy.flatten()-rmin)/res).astype('int')] = flags
        if np.isnan(exptimes[i]): continue
        expmap += grid*exptimes[i]
        if i % 500 == 0: print(f'Added {str(i)} of {str(len(pointings))}')
        if i==len(pointings_set)-1: print(f'Added {str(i+1)} of {str(len(pointings))}')

    return expmap

def add_exposures(expmap, ra, dec, exptimes, teff=0.2):
    '''
    Adds the summed effective exposure time to a DELVE-MC exposure map.

    Parameters:
    -----------
    expmap: array
        Original exposure map.
    ra: array
        Contains the RA of the exposures to be added to exposure map, expressed in degrees.
    dec: array
        Contains the DEC of the exposures to be added to exposure map, expressed in degrees.
    exptimes: float or array
        Array of the exposure times of the given coordinates. NaNs are supported.
    teff: float or array
        Effective exposure value, called effective exposure time in DES vocabulary.

    Returns:
    --------
    expmap: array
        2-dimensional array with same shape as input and the added exposures.
    '''

    ms_coords = SkyCoord(ra*u.deg,dec*u.deg).transform_to(MagellanicStreamNidever08())

    extra_expmap = decam_array_map(ms_coords.L.deg, ms_coords.B.deg, exptimes*teff)

    return expmap+extra_expmap

def tot_2_eff(exptimes, skybright, transm, fwhm, filter='g', pix_scale=0.263):
    '''
    Calculates effective exposure time using atmosferic conditions, using Equation 3 from Neilsen et al. 2019.

    Parameters:
    -----------
    exptimes: float or array
        Exposure times from header.
    skybright: float or array
        Sky brightness for every exposure.
    transm: float or array
        Sky transmission (G-TRANSP for DECam headers).
    fwhm: float or array
        FWHM in the image, expressed in pixels.
    filter: str, default='g'
        Filter of the exposures, used to define the dark sky value at zenith. It only has values for DECam filters.
        It is based on DELVERED processing figures.
    pix_scale: float, default=0.263
        Pixel scale of the detector. Default value corresponds to a DECam CCD.

    Returns:
    --------
    teff: float or array
        Effective exposure time, i.e. tau * exptime.
    '''
    dark_sky = dict([('u',0.1), ('g',1), ('r',5), ('i',11), ('z',18), ('Y',16)])
    sky = skybright/exptimes # sky rate
    taus = transm**2/(fwhm*pix_scale/0.9)**2/(sky/dark_sky[filter])
    return exptimes*taus

if __name__=="__main__":
    print('This is the main program.')
