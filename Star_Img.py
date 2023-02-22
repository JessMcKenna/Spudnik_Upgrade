# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 15:21:29 2023

@author: Jessica
"""

import numpy as np
from utilities import rotate, slope_angle_vector

"""
StarImageObject: 
    
This class stores every star as a single object. As the counterpart to Star_Calc_oop, it specifically takes in 
image coordinates, corrects (undistorts) these using the distortion coefficients found during the calibration cycle, 
and lastly calculates sky coordinates (horizontal coordinates) from the input pixel coordinates. 


"""

class ImgStar:
    
    def __init__(self, image_coordinates, image_center):
        self.image_coordinates = image_coordinates
        self.image_center = image_center
        self.translated_coord = None   # converted coordinate (image -> sky)
        self.corrected_coordinates = None   # coordinate which has been undistorted
        self.detected_Flag = False  # used to note whether star has been successfully detected
        self.image_zenith = None   # values that are assigned once and stored in all objects
        self.north_offset = None
        self.arcseconds_per_px = None


    """
    This method corrects distorted image coordinates, using the distortion coefficients returned from the calibration.
    It returns the corrected px coord of a star_img_coord; ie. the actual pixel coordinates after undistorting it.
    """

    def correct_coordinates(self, distortion_coef, image_input=None):
        center = self.image_center
        if not isinstance(image_input, type(None)):
            star_vector = np.subtract(image_input, center)
        else:
            star_vector = np.subtract(self.image_coordinates, center) 

        #STEP 1: The distortion is described by a parameter which represents distortion as a quadratic fn of distance from img centre 
        dist = np.linalg.norm(star_vector) #This is the distance from image center (assumed to be (0,0)) to star (given in px coords).
        distortion = (distortion_coef[0] * dist) + (distortion_coef[1] * (dist ** 2)) + (distortion_coef[2] * (dist ** 3)) + (distortion_coef[3]) #These coefficients come from calib cycle
        distortion_param = distortion / dist  
        distortion_vector = star_vector - (star_vector * distortion_param) #The distortion vector points from img centre to undistorted point

        #STEP 2: The corrected point (px coords) is: px coords of image center plus distortion vector.
        corrected_point = center + distortion_vector 
        corrected_point = round(corrected_point[0]), round(corrected_point[1])  
        
        if isinstance(image_input, type(None)):
            self.corrected_coordinates = corrected_point
        return corrected_point   #From the image coords and the distortion parameters, we get the corrected point in px coords.
        
    
    
    
    ''' This method takes in image coordinates, undistorts the position and then converts to horizontal coordinates (azimuth, elevation).
    '''
    
    
    def horizontal_coord(self, distortion_coef, position=None):
        if not isinstance(position, type(None)):
            star_pos = self.correct_coordinates(distortion_coef, position)
        else:
            star_pos = self.correct_coordinates(distortion_coef)
        
        zenith_star_vector = np.subtract(star_pos, self.image_zenith)          # Note: here order is very important
        corr_zenith_star_vector = rotate(self.north_offset, zenith_star_vector) 
        #Step 1: Correct the vector pointing from zenith to the star, by rotating it clockwise through the 'north offset angle.
        
        zenith_star_dist = np.linalg.norm(corr_zenith_star_vector)
        star_azimuth = 360 - slope_angle_vector(corr_zenith_star_vector)  
        #Step 2: Measure the counterclockwise angle obtained from the previous step. This is star Azimuth.

        star_elevation = 90 - ((zenith_star_dist * self.arcseconds_per_px) / 3600)   
        #Step 3: 90 - (px distance * arcseconds/px) /(# arcs/deg) = (# degrees) elevation from the equator

        if isinstance(position, type(None)):
            self.translated_coord = np.array([star_azimuth, star_elevation])

        return np.array([star_azimuth, star_elevation]) #From the distortion coefficients and image coords we get the horizontal coords.
