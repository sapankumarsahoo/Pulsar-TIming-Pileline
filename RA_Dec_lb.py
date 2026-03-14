from astropy.coordinates import SkyCoord
import astropy.units as u

def ra_dec_to_galactic(ra_hours, ra_minutes, dec):
    """
    Converts Right Ascension (RA) in hours and minutes and Declination (Dec) in degrees to
    Galactic Longitude (l) and Latitude (b).
    
    Parameters:
        ra_hours: float, Right Ascension hours.
        ra_minutes: float, Right Ascension minutes.
        dec: float, Declination in degrees.
        
    Returns:
        l: float, Galactic Longitude in degrees.
        b: float, Galactic Latitude in degrees.
    """
    # Create SkyCoord object with RA and Dec in hours, minutes, and degrees
    c = SkyCoord(ra=ra_hours*u.hour + ra_minutes*u.minute, dec=dec*u.degree, frame='icrs')

    # Convert to Galactic coordinates
    galactic = c.galactic

    # Extract Galactic Longitude and Latitude
    l = galactic.l.degree
    b = galactic.b.degree

    return l, b

# Example usage
ra_hours = 00
ra_minutes = 51
dec = -27
l, b = ra_dec_to_galactic(ra_hours, ra_minutes, dec)
print("Galactic Longitude (l):", l)
print("Galactic Latitude (b):", b)
