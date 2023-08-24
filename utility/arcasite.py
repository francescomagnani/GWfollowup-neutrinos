from astropy.coordinates import (
    EarthLocation,
    SkyCoord,
    AltAz,
    Longitude,
    Latitude,
)
from astropy.time import Time
from astropy.units import deg, rad
import astropy.units as u

import utm
import numpy as np
from datetime import datetime

# implied below
HEMISPHERE = "north"
DATUM = "WGS84"

# ARCA Sicily site
# taken from det ID 42 (km3db CLI -> detx 49)
arca_northing = 4016800
arca_easting = 587600
arca_height = -3450  # m
arca_utm_zone_number = 33
arca_utm_zone_letter = "N"
arca_utm_zone = "{num}{let}".format(num=arca_utm_zone_number, let=arca_utm_zone_letter)
arca_latitude, arca_longitude = utm.to_latlon(
    arca_easting, arca_northing, arca_utm_zone_number, arca_utm_zone_letter
)

# ----------------------------------------------------------------

LOCATIONS = {
    "arca": EarthLocation.from_geodetic(
        lon=Longitude(arca_longitude * deg),
        lat=Latitude(arca_latitude * deg),
        height=arca_height,
    ),
}

def get_location(location):
    try:
        loc = LOCATIONS[location]
    except KeyError:
        raise KeyError("Invalid location, valid is 'arca'")
    return loc

def utm_zone(lat):
    """The UTM zone for a given latitude

    Parameters
    ----------
    lat : number
        Latitude in rad

    """
    return 1 + int((np.pi + lat) / (6 * np.pi / 180))


def longitude_of_central_meridian(utmzone):
    """The longitude of the central meridian for a given UTM zone.

    Parameters
    ----------
    utmzone : number
        The UTM zone.

    """
    zone_width = 6 * np.pi / 180
    return -np.pi + (utmzone - 1) * zone_width + zone_width / 2




def convergence_angle(lat, lon):
    """Calculate the converge angle on the UTM grid.

    Parameters
    ----------
    lon : number
        Longitude in rad
    lat : number
        Latitude in rad

    """
    latitude_deg = lat * u.deg

    if latitude_deg > 84 * u.deg or latitude_deg < -80 * u.deg:
        raise ValueError(
            "UTM coordinate system is only defined between -80deg S and 84deg N."
        )

    # detector position, longitude and latitude in rad
    # lambda  = longitude
    phi = lat

    # find UTM zone and central meridian

    # longitude of the central meridian of UTM zone in rad
    lambda0 = longitude_of_central_meridian(utm_zone(lon))
    omega = lon - lambda0

    # parameters of the Earth ellipsoid
    sma = 6378137  # semi-major axis in meters (WGS84)
    ecc = 0.0066943800  # eccentricity (WGS84)

    rho = sma * (1 - ecc) / pow(1 - ecc * np.sin(phi) ** 2, 3 / 2)
    nu = sma / np.sqrt(1 - ecc * np.sin(phi) ** 2)
    psi = nu / rho
    t = np.tan(phi)

    angle = (
        np.sin(phi) * omega
        - np.sin(phi) * omega**3 / 3 * pow(np.cos(phi), 2) * (2 * psi**2 - psi)
        - np.sin(phi)
        * omega**5
        / 15
        * pow(np.cos(phi), 4)
        * (
            psi**4 * (11 - 24 * t**2)
            - psi**3 * (11 - 36 * t**2)
            + 2 * psi**2 * (1 - 7 * t**2)
            + psi * t**2
        )
        - np.sin(phi)
        * omega**7
        / 315
        * pow(np.cos(phi), 6)
        * (17 - 26 * t**2 + 2 * t**4)
    )

    return angle

def np_to_astrotime(intime):
    """Convert numpy/pandas datetime64 to list[datetime]."""
    nptime = np.atleast_1d(intime)
    np_corr = (nptime - np.datetime64("1970-01-01T00:00:00")) / np.timedelta64(1, "s")
    return Time([datetime.utcfromtimestamp(t) for t in np_corr])

# --------------------------------------------------------------------------------------------------------------------------------

def local_frame(time, location):
    """Get the (horizontal) coordinate frame of your detector."""
    if not isinstance(time, Time):
        # if np.datetime64, convert to astro time
        time = np_to_astrotime(time)
    loc = get_location(location)
    frame = AltAz(obstime=time, location=loc)
    return frame


def local_event(azimuth, time, zenith, location, radian=True, **kwargs):
    """Create astropy events from detector coordinates."""
    zenith = np.atleast_1d(zenith).copy()
    azimuth = np.atleast_1d(azimuth).copy()
    if not radian:
        azimuth *= np.pi / 180
        zenith *= np.pi / 180
    altitude = zenith - np.pi / 2

    loc = get_location(location)
    # neutrino telescopes call the co-azimuth "azimuth"
    true_azimuth = (
        np.pi / 2 - azimuth + np.pi + convergence_angle(loc.lat.rad, loc.lon.rad)
    ) % (2 * np.pi)
    frame = local_frame(time, location=location)
    event = SkyCoord(alt=altitude * rad, az=true_azimuth * rad, frame=frame, **kwargs)
    return event