from astropy.coordinates import SkyCoord
from astropy.time import Time

import numpy as np
import pandas as pd
from math import pi
import matplotlib
import matplotlib.pyplot as plt
import healpy as hp
from ligo.skymap import postprocess
import ligo.skymap.plot.allsky # for astro degrees mollweide projection
from astropy.units import deg


from utility.arcasite import local_event, local_frame


def get_region_90(gwskymap: np.ndarray):
    #From the GW skymap, define the region containing 90% of the probability.
    #Returns the indices of these pixels.

    probs = gwskymap / np.sum(gwskymap)
    isort = np.argsort(probs)[::-1]

    cumprobs = np.cumsum(probs[isort])
    ibound = np.argwhere(cumprobs >= 0.9)[0][0]

    return isort[:(ibound+1)]


def extend_region_90(region90: np.ndarray, nside: int, extension_deg: float):
    """From the region containing 90% of the probability, extend this region by some amount.
    Return the indices of all the pixels in this extended region.
    """

    ra, dec = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    coords = SkyCoord(ra=ra*deg, dec=dec*deg)

    mindist_to_region = np.inf * np.ones(hp.nside2npix(nside))
    for pix in region90:
        dists = (coords.separation(coords[pix])).deg#angular_separation(coords, coords[pix]).deg
        mindist_to_region = np.amin([mindist_to_region, dists], axis=0)
    regionext = np.argwhere(mindist_to_region < extension_deg)

    return np.unique(np.append(region90, regionext))


def convert_region_to_local(region: np.ndarray, nside: int, start_time: Time, end_time: Time, detector: str):
    """From the defined region of interest (RoI), convert it to local coordinates for times in a time window of interest.
    Returns a skymap {p_i} with p_i = (how often the pixel i is in the RoI during the time window).
    """

    npix = hp.nside2npix(nside)

    # prepare equatorial skymap with 1 if inside region90, 0 outside
    skymap_region = np.zeros(npix)
    skymap_region[region] = 1

    # prepare empty local skymap and get the coordinates of each pixel
    skymap_region_local = np.zeros(npix)
    theta, phi = hp.pix2ang(nside, np.arange(npix))

    # probe uniformly (nbins) times between start and end time
    nbins = 100
    deltatime = (end_time - start_time) / (nbins - 1)
    for ibin in range(nbins):
        time = start_time + ibin * deltatime
        # get local coords using km3astro objects
        coords_local = local_event(azimuth=phi, zenith=theta, time=time, location=detector.lower())
        coords_eq = coords_local.transform_to("icrs")
        ipix_eq = hp.ang2pix(nside, coords_eq.ra.deg, coords_eq.dec.deg, lonlat=True)
        skymap_region_local += 1/nbins * skymap_region[ipix_eq]

    return skymap_region_local


def define_zenith_bands(skymap_region_local: np.ndarray, size_bands_deg: float):
    """Define the zenith bands that are covering the full local region probed by the RoI, with bands of fixed size.
    The first one starts at the lowest theta value.
    Returns list of list of pixels for each band, list of sizes of off bands, list of sizes of ON region in each band.
    """

    nside = hp.get_nside(skymap_region_local)
    theta, _ = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    theta = np.rad2deg(theta)

    # get the theta range to probe
    theta_region_local = theta[skymap_region_local > 0]
    thetamin = min(theta_region_local)
    thetamax = max(theta_region_local)

    # define the theta bands
    # /!\ This may be optimized as there may be a gap between min and max that it is currently still considered though it has no contribution from the RoI.
    thetabands = np.arange(thetamin, thetamax+size_bands_deg, size_bands_deg)
    thetabands = np.clip(thetabands, 0, 180)

    pixelbands = []
    omegaoff, omegaon = [], []
    for iband in range(len(thetabands)-1):
        band = (theta >= thetabands[iband]) & (theta < thetabands[iband+1])
        pixelbands.append(np.argwhere(band).flatten())
        omegaoff.append(hp.nside2pixarea(nside) * len(pixelbands[-1]))
        omegaon.append(hp.nside2pixarea(nside) * np.sum(skymap_region_local[band]))

    return pixelbands, np.array(omegaoff), np.array(omegaon)

def plot_skymap_healpix(skymap: np.ndarray, on_region: list, df: pd.DataFrame, evt_time: Time, outfig: str, id):
    """Plot the skymap overlayed with the ON region and the selected KM3NeT events."""

    plt.figure(figsize=(13, 7))
    ax = plt.axes([0.05, 0.05, 0.9, 0.85], projection="astro degrees mollweide")
    ax.set_title("Skymap "+id, fontsize=22, pad=25)
    ax.grid()
    #plt.title(ID+" skymap", fontsize=30, pad=20)
    ax.tick_params(labelsize=15)
    
    npix = len(skymap)
    nside = hp.npix2nside(npix)
    
    # Draw skymap
    im_skymap = ax.imshow_hpx(skymap, cmap=plt.get_cmap("Reds"))
    
    # Skymap contours
    cls_skymap = 100 * postprocess.find_greedy_credible_levels(skymap/np.sum(skymap))
    ax.contour_hpx(
        (cls_skymap, "ICRS"),
        colors=["red", "red"],
        linewidths=[1, 1],
        linestyles=[":", "-"],
        levels=(50, 90),
    )
    # Draw ON region
    mask_on = -1 * np.ones_like(skymap)
    mask_on[on_region] = 1
    ax.imshow_hpx(mask_on, cmap=matplotlib.colors.ListedColormap(["white", "blue"]), vmin=0.0, vmax=1.0, alpha=0.15)

    # Draw horizon (darkens the half of the sky that is above horizon)
    ra, dec = hp.pix2ang(nside, range(npix), lonlat=True)
    coords = SkyCoord(ra=ra*deg, dec=dec*deg, frame="icrs")
    orca_frame = local_frame(time=evt_time, location="arca")
    coords = coords.transform_to(orca_frame)
    horizon = -1 * np.ones_like(skymap)
    horizon[np.sin(coords.alt.rad) > 0] = 1
    ax.imshow_hpx(horizon, cmap=matplotlib.colors.ListedColormap(["white", "grey"]), vmin=0.0, vmax=1.0, alpha=0.05)

    # Draw KM3NeT events
    if df.empty:
        a = 5
    else:
        ra, dec = np.array(df['trackfit_ra'])*180/pi, np.array(df['trackfit_dec'])*180/pi
        coords = SkyCoord(ra=ra*deg, dec=dec*deg, frame='fk5')
        ax.scatter(coords.ra, coords.dec, marker="x", color="blue", s=15, transform=ax.get_transform('fk5'))
    
    # Draw custom legend
    handles = []
    handles.append(matplotlib.lines.Line2D([], [], color="red", linewidth=1, linestyle="-", label="50%/90% contour"))
    handles.append(matplotlib.patches.Patch(facecolor="blue", alpha=0.15, label="ON region"))
    handles.append(matplotlib.patches.Patch(facecolor="grey", alpha=0.05, label=r"Region above horizon at $t_0$"))
    handles.append(matplotlib.lines.Line2D([], [], color="blue", marker="x", linewidth=0, label=r"Events at t_${GW}$"))
    plt.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.06), ncols=4)
    
    ax.set_xlabel("Right ascension", fontsize=15, labelpad=13.5)
    ax.set_ylabel("Declination", fontsize=15)
    plt.savefig(outfig, dpi=300)

