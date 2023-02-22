# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 11:54:44 2023

@author: Jessica

"""
import numpy as np
import math
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz
import cv2
#from astropy.coordinates import SkyCoord, EarthLocation, AltA


"""
This file lays out a variety of functions (utilities). 

Utilities file structure: 
    
        PART 1: Mathematical/Geometric
        PART 2: Conversions/Coordinate Transforms
        PART 3: Image/Frame Rescaling Commands
        
"""


    # \\\ PART 1: MATHEMATICAL/GEOMETRIC \\\ # 

#gets the distance between pos_1 and pos_2 in arb units/coord
def get_distance(pos_1, pos_2):
        x_1, y_1 = pos_1[:]
        x_2, y_2 = pos_2[:]
        distance = math.sqrt((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2)
        return distance

#gets the vector between pos_1 and pos_2 in arb units.coord
def get_vector(pos_1, pos_2):
        vector = np.subtract(pos_1, pos_2)
        return vector
    
#gets the foot point of a satellite trace in arb coords/coord
def get_foot_point(p, a, b):
    ap = p - a
    ab = b - a
    foot = a + np.dot(ap, ab) / np.dot(ab, ab) * ab
    return foot


#calcs unit_vector in arb units/coord
def unit_vector(vector):
    return vector / np.linalg.norm(vector)


#calculates angle between two vectors in degrees/coord
def angle_btw_vectors(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return math.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))


#gets angle corresponding to the slope between two coords in degrees 
"""def slope_angle(x, y, vector):
    if vector == None:
        vector = np.subtract(x, y) 
    elif (x - y).any() == None:
        vector= vector # this order of vector direction
    else:
        return None
    slope_angle = math.degrees(math.atan2(vector[0], - vector[1])) % 360                     #/TODO/
    return slope_angle"""


def slope_angle_points(x, y):
    vector = np.subtract(x, y) 
    slope_angle = math.degrees(math.atan2(vector[0], - vector[1])) % 360                     #/TODO/
    return slope_angle

def slope_angle_vector(vector):
    slope_angle = math.degrees(math.atan2(vector[0], - vector[1])) % 360                     #/TODO/
    return slope_angle

#gets the intersection points of two circles --> c0: centre (x0, y0), radius r0 // ---> c1: centre(x1, y1), radius r1
def get_circle_intersections(x0, y0, r0, x1, y1, r1):

    d = math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
    
    # discards circles which are non-intersecting
    if d > r0 + r1:
        return None
    # discards circles which lie inside one another
    if d < abs(r0 - r1):
        return None
    # discards coincident circles
    if d == 0 and r0 == r1:
        return None
    # gets intersection points of any other suitable circles
    else:
        a = (r0 ** 2 - r1 ** 2 + d ** 2) / (2 * d)
        h = math.sqrt(r0 ** 2 - a ** 2)
        x2 = x0 + a * (x1 - x0) / d
        y2 = y0 + a * (y1 - y0) / d
        x3 = x2 + h * (y1 - y0) / d
        y3 = y2 - h * (x1 - x0) / d
    
        x4 = x2 - h * (y1 - y0) / d
        y4 = y2 + h * (x1 - x0) / d
    
        return x3, y3, x4, y4
    

#rotates a vector counter clockwise by an angle theta (rad) and returns this
def rotate(theta, vector): 
    theta = np.radians(theta)
    cos, sin = np.cos(theta), np.sin(theta)
    rot_matrix = np.array([[cos, -sin], [sin, cos]])   
    return np.around(np.dot(rot_matrix, vector))

#defines a line direction by finding the slope of the line connecting two points
def line_direction(x1, y1, x2, y2):
    orientation = (y1 - y2) / (x1 - x2)
    if orientation >= 0:
        return math.degrees(math.atan(orientation))
    else:
        return math.degrees(math.atan(abs(orientation))) + 90


def objective(x, a, b, c, d):
    print(x, a, b, c, d)
    return (a * x) + (b * x ** 2) + (c * x ** 3) + d


    # \\\ PART 2: Conversions and Coordinate Transformations \\\ # 


"""

Conversion: No. of px to No. of arc.

[Used in: ]

"""
def get_arc_per_px(aligned_points, elevation_difference):
    # elevation_difference decimal degrees
    vector_length = np.linalg.norm(np.subtract(aligned_points[0], aligned_points[1]))
    arc_per_px = elevation_difference * 3600 / vector_length
    return arc_per_px


"""

This transfroms RaDec (equatorial) to Azmith, Elevation (horizontal aka. sky/viewing coordinates)

[Used in: ]

"""
"""def eq_2_horizontal(time_obj, location_obj):
    
     # tales in a time, and location object (from outside the class, because its faster)
     # and calculates the position of a star in the sky (horizontal coord)
     # at given time and date
       # This is very slow!!! use  Celest library or raw numpy
    star = SkyCoord(ra=star_ra * u.degree, dec=star_dec * u.degree, frame='icrs')
    # read docs about astropy.coordinates
    star_altaz = star.transform_to(AltAz(obstime=time_obj, location=location_obj))
    azimuth, elevation = star_altaz.az.degree, star_altaz.alt.degree
    #star_pos = np.array([azimuth, elevation])
    return azimuth, elevation"""


def eq_to_horizontal(ra, dec, time_obj, location_obj):
    star = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
    star_altaz = star.transform_to(AltAz(obstime=time_obj, location=location_obj))
    azimuth, altitude = star_altaz.az.degree, star_altaz.alt.degree
    return azimuth, altitude

"""

This transforms RaDec (equatorial) to Cartesian for use in analysing satellite trace & foot point

[Used in: ]

"""
    
def eq_2_cartesian(theta, phi):  # theta = right assertion, phi = declination; 3D cartesian 
     # this function takes normal declination -90 degrees
     x = math.cos(math.radians(theta)) * math.sin(math.radians(phi))
     y = math.sin(math.radians(theta)) * math.sin(math.radians(phi))
     z = math.cos(math.radians(phi))
     return x, y, z


"""

This transforms Azimuth, Elevation to pixels

[Used in: ]

"""
"""def horizontal_2_px(azimuth, elevation, north_offset, zenith_pos, arcsec):  # translate sky to image coordinates
     # This function uses the settings from the beginning
     radec_dist_zenith = (((90 - elevation) * 3600) / arcsec)
     angle = (-(360 - azimuth) + 90) % 360
     azimuth_vector = np.array([math.cos(math.radians(angle)), -math.sin(math.radians(angle))])
     zenith_star_vector_off = (radec_dist_zenith * azimuth_vector)
     zenith_star_pos = rotate(360 - north_offset, zenith_star_vector_off) + zenith_pos
     px_coord= round(zenith_star_pos[0]), round(zenith_star_pos[1])
     return px_coord"""
 
"""
This and the one in StarCalcData should be merged....all px-> RaDec

[Used in: ]

"""
def hor_2_px(star_pos, image_zenith, north_offset_ca, arcsec):
     stars_dist_zenith = (((90 - star_pos[1]) * 3600) / arcsec)
     angle = 90 + star_pos[0] % 360 #(-(360 - star_pos[0]) + 90) % 360
     azimuth_vector = np.array([math.cos(math.radians(angle)), -math.sin(math.radians(angle))])
     zenith_star_vector_off = (stars_dist_zenith * azimuth_vector)
     zenith_star_pos = rotate(360 - north_offset_ca, zenith_star_vector_off) + image_zenith
     px_coord = round(zenith_star_pos[0]), round(zenith_star_pos[1])
     return px_coord
    
"""
This transforms pixel coordinates to equatorial coordinates 

NOTE: still needs a lot of inputs...maybe that could be changed

[Used in: ]
"""
def px_2_eq(star_pos, time_obj, zenith_pos, north_offset, arcsec, location_obj): # convert an image coordinate to sky coordinates
        # only pass undistorted points !!!!
        # This funcion uses calibration settings from the beginning
        stars_direction_vector = np.subtract(star_pos, zenith_pos)  # order is very important
        translated_star_vector = rotate(north_offset, stars_direction_vector)

        stars_azimuth = 360 - slope_angle_vector(translated_star_vector)  # minus 360 for counter-clockwise results

        stars_dist = np.linalg.norm(np.subtract(star_pos, zenith_pos))
        stars_elevation = 90 - ((stars_dist * arcsec) / 3600)

        newAltAzcoordiantes = SkyCoord(alt=stars_elevation * u.degree, az=stars_azimuth * u.degree,
                                       obstime=time_obj,
                                       frame='altaz',
                                       location=location_obj).icrs

        return newAltAzcoordiantes.ra.degree, newAltAzcoordiantes.dec.degree

 


    # \\\ PART 3: IMAGE/FRAME RESCALING COMMANDS \\\ # 
    
def rescale_frame(frame, scale):
    width = int(frame.shape[1] * scale)
    height = int(frame.shape[0] * scale)
    dimension = (width, height)
    return cv2.resize(frame, dimension, interpolation=cv2.INTER_AREA)

