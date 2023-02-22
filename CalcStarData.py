# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 10:27:16 2023


CalcStar Data Object:
    
    This class takes data from a star database.
    
Structure: 
    
    There is one assigned method, which converts the star image coordinates to horizontal coordinates.


"""


from math import radians, sin, cos
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz
import utilities

# always pass numpy arrays, if list
# star_ra in decimal degrees, star_dec in decimal degrees -->database



class CalcStarData:
    def __init__(self, star_name, star_ra, star_dec, current_position, current_date, current_time):
        self.star_name = star_name
        self.star_ra = star_ra  #Note: convert decimal degrees to decimal hours
        self.star_dec = star_dec
        self.horizontal_coord = None  #Note: converted coordinates (image -> sky)
        self.true_coordinates = None  #Note: ??
        self.error_vector = None    # not used in here but used to store data
        #self.error_length = None    # not used in here but used to store data
        #self.star_size = None   # not used in here but used to store data
        self.current_position = current_position
        self.current_date = current_date
        self.current_time = current_time
        self.stars_dist_zenith = None
        self.height = 0
        self.rotate_point = utilities.rotate
        self.star_pos = np.array([])
        
    def equatorial_horizontal(self, time_obj, location_object):
        star = SkyCoord(ra=self.star_ra * u.degree, dec=self.star_dec * u.degree, frame='icrs')
         # read docs about astropy.coordinates
        star_altaz = star.transform_to(AltAz(obstime=time_obj, location=location_object))
        azimuth = star_altaz.az.degree; elevation = star_altaz.alt.degree
        horizontals = np.array([azimuth, elevation])
        self.star_pos.append(horizontals)
        return azimuth, elevation
      
    """
    METHOD: This method converts image coord. to horizontal coord. It takes inputs from the calibration, including arcs per px,
    north offset and image zenith. We should check where is this in use??
    
    """
    #horizontal 2 pixel
    def image_coordinates(self, arcseconds_per_px, north_offset, image_zenith):
        if arcseconds_per_px and north_offset and any(image_zenith) is not None:
            self.stars_dist_zenith = (((90 - self.star_pos[1]) * 3600) / arcseconds_per_px) # Note: px output. 3600 is the conversion between hour --> second
            angle = (-(360 - self.star_pos[0]) + 90) % 360  # Note: This is because of the different layouts of polar coordinate systems
            azimuth_vector = np.array([cos(radians(angle)), -sin(radians(angle))])
            zenith_star_vector_off = (self.stars_dist_zenith * azimuth_vector)
            zenith_star_pos = self.rotate_point(360 - north_offset, zenith_star_vector_off) + image_zenith #px
            self.horizontal_coord = round(zenith_star_pos[0]), round(zenith_star_pos[1])
            
            return round(zenith_star_pos[0]), round(zenith_star_pos[1])
        else:
            return None

print(dir(CalcStarData))