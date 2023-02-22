# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 17:50:05 2023

@author: Jessica
"""

import numpy as np
import cv2
import math
import utilities

'''
Zenith Approximation Object:
    
    This class allows calculation of the Zenith via a triangulation method. The triangulation method is based on the angles 
    between known stars and the approximate Zenith. This angle is referred to as theta. 
    
    The triangulation method works by first constructing all possible circles which contain a) any two stars, and b) the approx Zenith.
    The method uses geometric properties of inscribed triangles. 
    
    The intersections of the allowed circles are then found; thus, lines can be constructed between the intersection points of any two intersecting circles.
    Where all of these intersecting lines meet is the image Zenith.
    
Structure: 
    
    A circle class is assigned to the Zenith class. There are 4 main methods; 
    1) Construct Triangulation circles, 
    2) Find average Zenith intersection, 
    3) Get the North vector, 
    4) Visualise what's been done
    
'''

class Zenith:
    def __init__(self, approx_zenith_pos, star_angles, star_positions, index=None):
        self.approx_zenith_pos = approx_zenith_pos
        self.approx_zenith_pos_loop = approx_zenith_pos ###?
        
        self.star_angles = star_angles
        
        self.star_positions = star_positions
        self.star_px_coord = star_positions ###?
        
        self.distortion_star_positions = None
        self.zenith_error = None
        self.zenith_pos = None
        self.circle_objects = None
        self.average_zenith_intersection = None
        self.north_vector = None

         
    """
    METHOD: 
        
        This takes two star positions and the input Zenith approximation. It constructs a line between these stars, 
    bisects that line and then uses trig and the properties of inscribed triangles to find the centre of the circle, 
    as well as the radius of the circle that can be formed from the 3 given points. 
    It returns the center point and radius (arbitrary coords). 
    
    """
    
    def triangulation(self, p1, p2, theta):
        # Convert p1 and p2 to NumPy arrays
        p1 = np.array(p1); p2 = np.array(p2)
        # Compute the vector between p1 and p2
        p1p2_vector = p2 - p1
        # Compute a unit vector perpendicular to p1p2_vector
        p1p2_perp = np.array([p1p2_vector[1], -p1p2_vector[0]]) / np.linalg.norm(p1p2_vector)
        # Compute the midpoint of the line segment connecting p1 and p2
        p1p2_midpoint = p1p2_vector * 0.5 + p1
        # Compute the length of the line connecting the midpoint and the circle center
        mid_length = math.tan(math.radians(90 - theta)) * np.linalg.norm(p1p2_vector * 0.5)
        # Compute the radius of the constructed circle
        radius = round(np.linalg.norm(p1p2_vector * 0.5) / math.cos(math.radians(90 - theta)))
        # Compute the two candidate midpoints
        midpoints = [p1p2_perp * mid_length + p1p2_midpoint, -p1p2_perp * mid_length + p1p2_midpoint]
        # Round the midpoint coordinates to integers
        midpoints = np.round(midpoints).astype(int)
        return midpoints[0], midpoints[1], radius, p1p2_midpoint
           

    """
    METHOD: 
        
        This takes the constructed circles and finds their intersections using functions in the utilities file.
    We enumerate the angles between each star and the north direction; we also find the angles between stars (theta).
    If theta < 90 we construct a circle between between the 2 stars at an angle theta, which intersects the approx zenith.
    We append this to the list of circle objects.
    
    """
    def zenith_intersection(self):
       
        while True:     
            
            circle_objects = []
        
            for i, angle in enumerate(self.star_angles):
                for p, sub_angle in enumerate(self.star_angles):
                    if angle != sub_angle:
                        theta = abs(angle - sub_angle)
                        # Not possible to construct the circle with obtuse angles. Thus this will filter the too big triangulation-circles
                        if theta > 180:                
                            theta = 360 - theta 
                            
                        else:
                            circle_1, circle_2, radius, mid = self.triangulation(self.star_positions[i], self.star_positions[p], theta)

                            dist_1 = np.linalg.norm(np.subtract(circle_1, self.approx_zenith_pos_loop)) #norm(midpoint 1 - approx zenith pos)
                            dist_2 = np.linalg.norm(np.subtract(circle_2, self.approx_zenith_pos_loop)) #norm(midpoint 2 - approx zenith pos)

                            if dist_1 < dist_2:
                                circle_objects.append(self.Circle(circle_1, radius))
                            else:
                                circle_objects.append(self.Circle(circle_2, radius))
            #append the circle objects for which the distance from the centre to the approx zenith position is less 
            #circle_objects.append(self.triangulation(self.star_positions[i], self.star_positions[p], theta))
                                       
            self.circle_objects = circle_objects  
            
            #initiating empty lists to store the coords of the calculated zenith position and error
            zenith_x_pos = [] 
            zenith_y_pos = []
            
            zenith_dist_error = []
            
            #This enumerates the items in the circle_objects list and picks out the circles that could be intersecting. 
            #It finds their centres and then gets their intersecting points. 
            #It checks to see that the intersecting points are below a threshold distance from the Zenith approx
            
            
            for i, circle in enumerate(circle_objects): #circle is the circle obj
                 for j, sub_circle in enumerate(circle_objects):
                     if i != j:  #find intersection of circle and all other circles
                         pos_1 = circle.center
                         pos_2 = sub_circle.center
                         intersections = utilities.get_circle_intersections(pos_1[0], pos_1[1], round(circle.radius), pos_2[0], pos_2[1], round(sub_circle.radius))
                         
                         if intersections:
    
                             intersect_1 = intersections[0:2] #finds (x,y) of first intersection point ---> labelling [0:2] = [0,1] = (x1, y1)
                             intersect_2 = intersections[2:4] #finds (x,y) of second intersection point ---> labelling [2:4] = [2,3] = (x2, y2)
    
                             dist_1 = abs(np.linalg.norm(np.subtract(intersect_1, self.approx_zenith_pos_loop))) #distance from (x1,y1) to approx zenith
                             dist_2 = abs(np.linalg.norm(np.subtract(intersect_2, self.approx_zenith_pos_loop))) #distance from (x2,y2) to approx zenith
    
                             if dist_1 < dist_2 and dist_1 < 200: #arbitrary choice for distance 
                                 zenith_intersection = intersect_1
                                 zenith_dist_error.append(dist_1)
    
                                 zenith_x_pos.append(zenith_intersection[0])
                                 zenith_y_pos.append(zenith_intersection[1])
    
                             elif dist_1 > dist_2 and dist_2 < 200:
                                 zenith_intersection = intersect_2
                                 zenith_dist_error.append(dist_2)
    
                                 zenith_x_pos.append(zenith_intersection[0])
                                 zenith_y_pos.append(zenith_intersection[1])
                                 
                                 #This part of code could be made much nicer very easily!

            self.average_zenith_intersection = np.array([round(sum(zenith_x_pos) / len(zenith_x_pos)), round(sum(zenith_y_pos) / len(zenith_y_pos))])
            self.zenith_error = sum(zenith_dist_error) / len(zenith_dist_error)

             #if the distance between zenith intersection av and approx zenith is > 3 we set approx zenith equal to average and repeat
            if np.linalg.norm(np.subtract(self.average_zenith_intersection, self.approx_zenith_pos_loop)) > 3:
                 self.approx_zenith_pos_loop = self.average_zenith_intersection
             #or else we return av
            else:
                 return self.average_zenith_intersection 

        
    """
    METHOD: 
    
    This method returns the north vector. It does this by finding the maximum star angle; we then find the star position of the star indexed by this angle,
    We construct a vector between that position and the average zenith just calculated, and rotate anticlocwise by delta (260 - star angle) degrees to find the rotated position,
    
    We do the same with the min star angle, except that the rotation is by the min angle and in the clocwise direction.
    
    The north vector is just the average coords of these two rotated vectors.
    
    The star_angles are azimuth and the star_positions are [az, el] values (I think)
    
    """
    
    def get_north_vector(self):
        max_angle = max(self.star_angles) #max azimuth
        max_angle_delta = 360 - max_angle #
        max_pos = self.star_positions[self.star_angles.index(max_angle)]
        max_pos_vector = np.subtract(max_pos, self.average_zenith_intersection)
        max_pos_rotate = utilities.rotate(- max_angle_delta, max_pos_vector)  # clockwise rotation, since rotate is defined anticlockwise

        min_angle = min(self.star_angles)
        min_pos = self.star_positions[self.star_angles.index(min_angle)]
        min_pos_vector = np.subtract(min_pos, self.average_zenith_intersection) 
        min_pos_rotate = utilities.rotate(min_angle, min_pos_vector) #anticlockwise rotation 

        self.north_vector = np.array([round((max_pos_rotate[0] + min_pos_rotate[0]) / 2), round((max_pos_rotate[1] + min_pos_rotate[1]) / 2)])
        return self.north_vector


    """
    METHOD: 
        
    This method plots various items in order to visualise the approach to zenith calculation.
    
        a) It constructs circles around the zenith approx and the actual calculated zenith.
        b) It plots the intersecting circles for each circle objects.
        c) It also plots the north vector. 
    
    """
    def visualize(self, blank, show=None):

        for circle in self.circle_objects:
            center = circle.center
            radius = circle.radius
            cv2.circle(blank, center, radius, (255, 200, 0), thickness=2)

        intersect = self.average_zenith_intersection

        #make circles around the approximate and the actual (average) calculated zenith
        cv2.circle(blank, intersect, 20, (255, 255, 0), thickness=6)   #args: cv2.circle(image, center_coordinates, radius, color, thickness)
        cv2.circle(blank, self.approx_zenith_pos, 20, (0, 0, 255), thickness=6)   #args: cv2.circle(image, center_coordinates, radius, color, thickness)

        #draw a line to represent north vector
        if self.north_vector is not None:
            cv2.line(blank, intersect, intersect + self.north_vector, (255, 255, 255), thickness=4)

        if show:
            blank_rescale = utilities.rescale_frame(blank, 0.25)
            cv2.imshow("Alignment: ", blank_rescale)
            cv2.waitKey()
            return blank
        else:
            return blank
        

    class Circle:
        def __init__(self, center, radius):
            self.center = center
            self.radius = radius 