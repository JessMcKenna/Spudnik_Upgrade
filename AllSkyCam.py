# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 16:14:19 2023

@author: Jessica McKenna

"""

import glob
import cv2
import os
import numpy as np
import argparse
import matplotlib.pyplot as plt 



# /// FILENAMES AND DETAILS \\\ #

filename_0 = 'Combined_Frames.jpg'
filename_1 = 'Grayscaled.jpg'
filename_2 = 'Threshold.jpg'
filename_3 = 'CannyEdges.jpg'
filename_4 = 'Cropped_Frame.jpg'



# /// PART 1: COMBINING ALL IMAGES IN DIRECTORY \\\ #

directory = r'C:\Users\Viesis\Documents\Thesis\ImageProcessing\Working_Scripts_Processing_Allsky\SatPassing2'

# Import all image files with the .jpg extension
os.chdir(directory) 
print(os.listdir(directory)) 

files = glob.glob ("*.jpg")
image_data = []

for my_file in files:
    this_image = cv2.imread(my_file, 1)
    image_data.append(this_image)
 
# Calculate blended image
img = image_data[0]
for i in range(len(image_data)):
    if i == 0:
        pass
    else:
        alpha = 1.0/(i + 1)
        beta = 1.0 - alpha
        img = cv2.addWeighted(image_data[i], alpha, img, beta, 0.0)
        


#Save the resulting image to the processed image folder
path = r'C:/Users/Viesis/Documents/Thesis/ImageProcessing/Working_Scripts_Processing_Allsky/SatPassing2/ProcessedSatPassing'

cv2.imshow('Combined_Frames', img)
cv2.imwrite(os.path.join(path, filename_0),img)



# /// PART 2: PRE-PROCESSING OF THE COMBINED IMAGE \\\ #

##create mask 
h,w =img.shape[:2]
mask = np.zeros((h,w), np.uint8)

#grayscaling
img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY) 

#thresholding
ret, threshold = cv2.threshold(img_gray, 13, 255, cv2.THRESH_BINARY) #13 - 15 works for this parameter
cv2.imshow('threshold', threshold)
cv2.imwrite(os.path.join(path, filename_2), threshold)

#Extract only the circular region of the image and crop this

#edge detection
aperture_size = 5 
t_lower = 220
ratio = 2
kernel_size = 3

img_edge = cv2.imread(os.path.join(path, 'Threshold.jpg'))
canny_edges = cv2.Canny(img_edge, t_lower, t_lower*ratio, apertureSize=aperture_size)

cv2.imshow('edge detection', canny_edges)
cv2.imwrite(os.path.join(path, filename_3), canny_edges)


#circle transform
circles = cv2.HoughCircles(canny_edges, cv2.HOUGH_GRADIENT, 1, 10000, param1 = 100, param2 = 1, minRadius = 100, maxRadius = 250)

if circles is not None:
    circles = np.round(circles[0,:]).astype("int")
    for (x, y, r) in circles:
        cv2.circle(mask, (x,y), r, (255, 255, 255), thickness = -1)

        
cv2.imshow('circle detected', img_edge)

#Copy image using mask 
masked_data = cv2.bitwise_and(img_gray, img_gray, mask=mask)

# Apply Threshold
_,thresh = cv2.threshold(mask,1,255,cv2.THRESH_BINARY)

# Find Contour
contours = cv2.findContours(thresh,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)

#Code to close Window
cv2.imshow('Detected Edge',img_gray)
cv2.imshow('Cropped Image',masked_data)
cv2.imwrite(os.path.join(path, filename_4), masked_data)
cv2.waitKey(0)
cv2.destroyAllWindows()




# ///PART 3: Processing the cropped frame \\\ # 

# /// FILENAMES AND DETAILS \\\ #
filename_5 = 'Blurred_and_Sharpened.jpg'

# Applying the Canny Edge filter
cropped_img =  cv2.imread(os.path.join(path, 'Cropped_Frame.jpg')) #cv2.imread(os.path.join(path,'LD_Grayscaled.jpg'))  #take grayscale image, blur it

# Applying a Gaussian Blur
Gaussian_blurred = np.hstack([cv2.GaussianBlur(cropped_img,(3,3),0)])

cv2.imshow('Gaussian_Blurred', Gaussian_blurred)

# Sharpening with Laplacian Kernel
kernel1 = np.array([[-1, -1, -1], [-1, 8, -1], [-1, -1, -1]]) 
#Note; there are different kinds of these sharpening kernels and this one produced the clearest image

image_sharpened= cv2.filter2D(src=Gaussian_blurred, ddepth=-1, kernel=kernel1)

cv2.imshow('Gaussian_Blurred_Laplacian_Sharpened', image_sharpened)
cv2.imwrite(os.path.join(path, filename_5), image_sharpened)

cv2.waitKey()
cv2.destroyAllWindows()


#Satellite Track Detection

#Canny Edge 
canny_sat = cv2.imread(os.path.join(path, filename_5))

aperture_size = 5 #Aperture size  -- Should be odd number ( 3 - 7)High aperture picks up on background stars; but very low also picks them up
t_lower = 200
ratio = 2
kernel_size = 3

edges = cv2.Canny(canny_sat, t_lower, t_lower*ratio, apertureSize=aperture_size)

#Hough Lines
Rres = 1
Thetares = 1*np.pi/180
Threshold =10
minLineLength = None
maxLineGap = 50

#lines = cv2.HoughLinesP(edges, Rres, Thetares, Threshold, minLineLength, maxLineGap )
lines = cv2.HoughLinesP(edges, 1, np.pi / 180, 10, 15, 50, 2)

N = lines.shape[0]

for i in range(N):
    x1 = lines[i][0][0]
    y1 = lines[i][0][1]    
    x2 = lines[i][0][2]
    y2 = lines[i][0][3]    
    cv2.line(canny_sat,(x1,y1),(x2,y2),(255,0,0),2)

plt.figure()
plt.imshow(canny_sat)
plt.savefig(os.path.join(path, 'Line_Detection'))
plt.title('Linear Structures & Satellites')
plt.axis('off')

plt.show()

