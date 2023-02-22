# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 10:26:14 2023
@author: Jessica and Anna

Manual Star

This class defines an object with img_pos and eq_pos as its attributes.

It contains a single method to convert the equatorial positions (which were input manually as a string into the API) into something that
is workable - ie. does unit conversion and returns the viewing angles (aka. horizontal coordinates // alt/az // alt/elevation) of the star.

"""

import numpy as np
import utilities 

class ManualStar:

    def __init__(self, img_pos, eq_pos):
        self.img_pos = img_pos
        self.eq_pos = eq_pos   #wheres eq_pos list/array? time?
        
        
    def star_viewing_angles(self, time_obj, location_obj):
       ra = (self.eq_pos[0] + self.eq_pos[1] / 60 + self.eq_pos[2] / 3600) * 15 #right ascension in deg
       if self.eq_pos[3] < 0:
           dec = self.eq_pos[3] - (self.eq_pos[4] / 60) - (self.eq_pos[5] / 3600)
       else:
           dec = self.eq_pos[3] + (self.eq_pos[4] / 60) + (self.eq_pos[5] / 3600)
       
       azimuth, altitude = utilities.eq_to_horizontal(ra, dec, time_obj, location_obj)
       self.viewing_angles = np.array([azimuth, altitude])
       return azimuth, altitude
   