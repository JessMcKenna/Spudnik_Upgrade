# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 22:52:44 2023

@author: Jessica McKenna

The purpose of this script is to automatically retrieve calibration data for nightsky images using Astrometry.net. 

REFERENCES:
    1) The script is based off this: https://stellartrip.net/wp-content/uploads/2020/12/3.1-astrometrie-image_en-1.html 
    2)                               https://astroquery.readthedocs.io/en/latest/astrometry_net/astrometry_net.html#handling-results


To make the plate solving more accurate, it is best to refine the parameters by indicating the approximate scale etc 

    
    
"""

import requests
import json
import os
import keyring
from astropy.io import fits
from PIL import Image
import matplotlib.image as mpimg
import numpy as np
from astroquery.astrometry_net import AstrometryNet
import matplotlib.pyplot as plt
from astropy.wcs import WCS

from matplotlib.colors import LogNorm

# Generate a LOGIN key in JSON (JSON-encoded string as a text field; request-json.)
apikey = 'bmhgeshkrzpspsrg'
#This sends the post request to the system 
R = requests.post('http://nova.astrometry.net/api/login', data={'request-json': json.dumps({"apikey": apikey})})
#Prints the session key for logging in to the system at Astrometry.net
print(R.text)
# Returns a session key which gives amount of time before session terminates

# Connect to the Astrometry.net server
AstrometryNet.api_key = 'bmhgeshkrzpspsrg'

rectified_images = r'C:\Users\Viesis\Documents\TheMotherShip\Calibrated_Images\Rectified_Data_Images\26thMarchUTCRectified\trial' # Where images will be stored 
  


plate_solved_images = []
# Loop through all files in the folder
for file_name in os.listdir(rectified_images):
    print(file_name)
    # Check if the file is an image file
    if not file_name.endswith('.jpg'):
        print(f'{file_name} not suitable for solving.')
        pass
    else:
        try_again = True
        submission_id = None
        # Get the full path to the file
        file_path = os.path.join(rectified_images, file_name)
        # Load the JPEG image
        image = Image.open(file_path)
        # Use the full path to the file when calling solve_from_image
        wcs_header = AstrometryNet.solve_from_image(file_path, submission_id=submission_id, solve_timeout=240, force_image_upload=True )
        
        if wcs_header:
            print('Header retrieved; WCS solution produced')
            sol_wcs = WCS(wcs_header)
    
            # Convert the JPEG image to a FITS file
            fits_file = 'my_image.fits'
            image_data = np.array(image)
            hdu = fits.PrimaryHDU(image_data)
            hdu.writeto(fits_file, overwrite=True)
    
            print(sol_wcs) # Print the header
            #Prepare plot
            img = mpimg.imread(file_path) #Image name_file
            fig = plt.figure(figsize=(20, 20))
    
            ax = plt.subplot(projection=sol_wcs) #Add wcs projection with wcs_gamcas file
            plt.imshow(img, origin='lower', cmap='gray', aspect='equal', norm = LogNorm())
    
            plt.title(f'Astrometry.net Plate Solution / {file_name}',fontsize=12, loc='center')
    
            ax.coords['ra'].set_axislabel('RA')
            ax.coords['dec'].set_axislabel('DEC')
    
            overlay = ax.get_coords_overlay('icrs')
            overlay.grid(color='green', ls='dotted', linewidth=3)
            
            # Add the following code to save the figure to the same directory
            plt.savefig(os.path.join(rectified_images, file_name.split('.')[0] + '_wcs.png'))

            plt.show()
            
            # Convert the solution to a FITS file stored locally in the folder
           # fits_file = file_name.split('.')[0] + '_wcs.fits'
            fits_file = os.path.join(rectified_images, file_name.split('.')[0] + '_wcs.fits')
            image_data = np.array(image)
            hdu = fits.PrimaryHDU(image_data, header=wcs_header)
            hdu.writeto(fits_file, overwrite=True)
            
            print(f'{fits_file} saved')

        else:
            print('Solve has failed; no WCS solution produced')
            # Code to execute when solve fails
            
 

"""

       
plate_solved_images = []
# Loop through all files in the folder
for file_name in os.listdir(rectified_images):  
    print(file_name)
    # Check if the file is an image file
    if not file_name.endswith('.jpg'):
        print(f'{file_name} not suitable for solving.')
        pass
    else:
        try_again = True
        submission_id = None
        
        image = Image.open(os.path.join(rectified_images, file_name))

        while try_again:
            try:
                if not submission_id:
                    wcs_header = AstrometryNet.solve_from_image(os.path.join(rectified_images, file_name), submission_id=submission_id, force_image_upload=True )
                else:
                    wcs_header = AstrometryNet.monitor_submission(submission_id, solve_timeout=120)
                
            except TimeoutError as e:
                submission_id = e.args[1]
            else:
                # got a result, so terminate
                try_again = False
        
    
    if wcs_header:
        print('Header retrieved; WCS solution produced')
        sol_wcs = WCS(wcs_header, relax = True)
        print(f'{file_name} has solution:',  sol_wcs) # Print the header
        plate_solved_images.append(file_name)
        
        #Prepare plot
        img = mpimg.imread(file_name) #Image name_file
        fig = plt.figure(figsize=(20, 20))

        ax = plt.subplot(projection=sol_wcs) #Add wcs projection with wcs_gamcas file
        plt.imshow(img, origin='lower', cmap='gray', aspect='equal', norm = LogNorm())

        plt.title(f'Astrometry.net Plate Solvingof {file_name}',fontsize=16, loc='center')

        ax.coords['ra'].set_axislabel('RA')
        ax.coords['dec'].set_axislabel('DEC')

        overlay = ax.get_coords_overlay('icrs')
        overlay.grid(color='red', ls='dotted')
        plt.savefig(os.path.splitext(file_name)[0] + '_PlateSolved.jpg')
        plt.show()
        
    else:
        print('Solve has failed; no WCS solution produced')
        # Code to execute when solve fails
        
        
 # Save WCS header to a new FITS file
 file_path = os.path.join(rectified_images, file_name)
 wcs_file_name = os.path.splitext(file_name)[0] + '_wcs.fits'
 original_file_path = os.path.join(file_path, wcs_file_name)
 

 hdu = fits.PrimaryHDU()
 hdu.header['SIMPLE'] = True  # add the SIMPLE keyword
 hdu.header = sol_wcs.to_header()
 hdu_list = fits.HDUList([hdu])
 hdu_list.writeto(wcs_file_name, overwrite=True)"""
