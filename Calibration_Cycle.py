# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 17:13:48 2023

@author: Jessica 
"""

''' This script is an API for the calibration of a camera. The calibration methods are briefly described in the report
    and in the file with the code updates. You have to manually identify stars in the image, this can be done by using
   https://nova.astrometry.net/  for example. '''
   

from tkinter import Tk, Button, Label, Entry, StringVar, mainloop
from os.path import exists
from numpy import arange
from utilities import rescale_frame
import numpy as np
import datetime
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation 
import pytz
import cv2
from matplotlib import pyplot
from scipy.optimize import curve_fit

from utilities import slope_angle_vector, objective, hor_2_px
from Star_Img import ImgStar
from ManualStar import ManualStar
from CalcStarDataRev import CalcStarData
from Zenith_Approximation import Zenith


""" 
Part 1: GENERAL
    Part 1 contains:
        
        A data storage section to initialise the variables to be entered into the API or calculated by the program.
        General information entry windows for the API and management of input data types (string splitting etc)
        Management of error notices, labels, etc.

        Two functions:
            a) Comd_general: This manages the general info input manually into the API
            b) Comd_star: This extracts the star coordinates input manuallly into the API (both px coords and eq coords)
 """
 
# /// DATA STORAGE /// #

# General Entry Data Storage (approximation)
file_path = None
position = None
current_date = None
current_time = None
arcsec = None
zenith_ap = None
north_offset_ap = None
height = None
rescale_setting = 0.9

# Alignment Data (calculation)
north_offset_ca = None
zenith_pos_ca = None

# Adjustment Data (accurate, adjustment)
zenith_pos = None
north_offset = None

# Distortion
distortion_parameters = None
distortion_flag = False

# Star Entry Data Storage
star_img_cord = []
star_eq_cord = []

# star objects
star_objects = []

window = Tk()
window.title('Alignment Cycle')
position_window = f'+{100}+{100}'
window_size = f'{800}x{800}'
window.resizable(0, 0)
window.geometry(position_window)
window.geometry(window_size)
window.configure(bg='#404040')

# Flags
next_sep_flag = False


# /// GENERAL FUNCTION /// # --> This manages the data input to the API

image_cord_entry = []
eq_cord_entry = []


def comd_general():
    global file_path
    global position
    global current_time
    global arcsec
    global current_date
    global zenith_ap
    global north_offset_ap
    global height
    global rescale_setting

    # only autofill
    global star_img_cord
    global star_eq_cord

    # // Autofill inputs //
    star_img_cord_temp = []
    star_eq_cord_temp = []

    autofill_path = general_entry[8].get()
    if autofill_path != "":
        if autofill_path and autofill_path[-3:] == "txt": #exists(autofill_path)
            setting_str = general_entry[9].get()
            if setting_str != "":
                flag = all(element.replace('.', "1").isdigit() for element in setting_str)
                setting_value = float(setting_str)
                if flag and setting_value < 1:
                    rescale_setting = setting_value
                else:
                    error_label.config(text="Setting Error")
                    error = True
            else:
                rescale_setting = 0.8

            error_label.config(text="")

            data = open(autofill_path, "r").readlines()
            data = [element.strip("\n").replace(" ", "") for element in data]
            data = list(filter(None, data))
            print(data[0])
            if data[0][-3:] == "jpg" or "JPG":
                 if data[0]:
                     file_path = data[0]
                    

            pos_data = data[1].split(",")
            if len(pos_data) == 2:
                position = (float(pos_data[0]), float(pos_data[1]))
                
            date_data = data[2].split(",")
            if len(date_data) == 3:
                current_date = (int(date_data[0]), int(date_data[1]), int(date_data[2]))
              
            time_data = data[3].split(",")
            if len(time_data) == 3:
                current_time = (int(time_data[0]), int(time_data[1]), int(time_data[2]))
               
            arcsec = float(data[4])
            
            zenith_ap_data = data[5].split(",")
            if len(zenith_ap_data) == 2:
                zenith_ap = (int(zenith_ap_data[0]), int(zenith_ap_data[1]))

            north_offset_ap = float(data[6])

            height = float(data[7])

            for index in range(8, len(data)):
                star = data[index].split("/")
                img_cord = star[0].split(",")
                img_cord = (int(img_cord[0]), int(img_cord[1]))
                star_img_cord_temp.append(img_cord)

                eq_cord = star[1].split(",")
                eq_cord = (
                int(eq_cord[0]), int(eq_cord[1]), int(eq_cord[2]), int(eq_cord[3]), int(eq_cord[4]), int(eq_cord[5]))
                star_eq_cord_temp.append(eq_cord)


            if len(star_eq_cord_temp) > 4:
                error_label_2.config(text="")
                star_img_cord = star_img_cord_temp
                star_eq_cord = star_eq_cord_temp
                return
            else:
                error_label_2.config(text="Autofill Enter more Stars")

        else:
         error_label.config(text="Autofill Path Error")

    
    # // Normal manual inputs/
    error = False
    unfill_error = False

    
    path = general_entry[0].get()
    if path != "":
        if path and path[-3:] == "jpg":
            file_path = path
            error_label.config(text="")
        else:
            error_label.config(text="Path Error")
            error = True
    else:
        unfill_error = True
        file_path = None
   

    position_str = general_entry[1].get()
    if position_str != "":
        pos = position_str.replace(" ", "")
        pos = pos.split(",")
        _ = [element.replace('.', "1") for element in pos]
        _ = [element.replace(',', "1") for element in _]
        flag = all(element.isdigit() for element in _)
        if flag and len(pos) == 2:
            position = float(pos[0]), float(pos[1])
        else:
            error_label.config(text="Position Error")
            error = True
    else:
        unfill_error = True
        position = None

    date_str = general_entry[2].get()
    if date_str != "":
        date = date_str.replace(" ", "")
        date = date.split(",")
        flag = all(element.isdigit() for element in date)
        if flag and len(date) == 3:
            current_date = int(date[0]), int(date[1]), int(date[2])
        else:
            error_label.config(text="Date Error")
            error = True
    else:
        unfill_error = True
        current_date = None

    time_str = general_entry[3].get()
    if time_str != "":
        time = time_str.replace(" ", "")
        time = time.split(",")
        flag = all(element.isdigit() for element in time)
        if flag and len(time) == 3:
            current_time = int(time[0]), int(time[1]), int(time[2])
        else:
            error_label.config(text="Time Error")
            error = True
    else:
        unfill_error = True
        current_date = None

    arc_str = general_entry[4].get()
    if arc_str != "":
        # if "," not in arc_str:
        flag = all(element.replace('.', "1").isdigit() for element in arc_str)
        if flag:
            arcsec = float(arc_str)
        else:
            error_label.config(text="Arcsec Error")
            error = True
    else:
        unfill_error = True
        arcsec = None

    zenith_str = general_entry[5].get()
    if zenith_str != "":
        zenith = zenith_str.replace(" ", "")
        zenith = zenith.split(",")
        flag = all(element.isdigit() for element in zenith)
        if flag and len(zenith) == 2:
            zenith_ap = int(zenith[0]), int(zenith[1])
        else:
            error_label.config(text="Zenith AP Error")
            error = True
    else:
        unfill_error = True
        zenith_ap = None

    offset_str = general_entry[6].get()
    if offset_str != "":
        # if "," not in offset_str:
        flag = all(element.replace('.', "1").isdigit() for element in offset_str)
        if flag:
            north_offset_ap = float(offset_str)
        else:
            error_label.config(text="offset Error")
            error = True
    else:
        unfill_error = True
        north_offset_ap = None

    height_str = general_entry[7].get()
    if height_str != "":
        # if "," not in height_str:
        flag = all(element.replace('.', "1").isdigit() for element in height_str)
        if flag:
            height = float(height_str)
        else:
            error_label.config(text="height Error")
            error = True
    else:
        unfill_error = True
        height = None

    setting_str = general_entry[9].get()
    if setting_str != "":
        flag = all(element.replace('.', "1").isdigit() for element in setting_str)
        if flag:
            setting_value = float(setting_str)
            if setting_value < 1:
                rescale_setting = setting_value
            else:
                error_label.config(text="Setting Error")
                error = True
        else:
            error_label.config(text="Setting Error")
            error = True
    else:
        rescale_setting = 0.8

    if not error:
        error_label.config(text="")

    if unfill_error:
        unfill_label.config(text="Unfilled Grid")
    else:
        unfill_label.config(text="")
 
        
# /// STAR FUNCTION /// # --> This gets the input star coordinates in the correct form
def comd_star():
    error = False
    unfill_error = False

    global star_img_cord
    global star_eq_cord

    img_cal_stars = []
    eq_cal_stars = []

    for index, cord in enumerate(image_cord_entry):
        eq_cord = None
        img_cord = None
        img_entry = cord.get()

        #Extract the input image coordinates & run error check
        if img_entry != "":
            values = img_entry.split(",")
            values = [elements.replace(" ", "") for elements in values]
            flag = all(element.isdigit() for element in values)
            if flag and len(values) == 2:
                img_cord = (int(values[0]), int(values[1]))
            else:
                error_label_2.config(text="Img Cord Error")
                error = True
        else:
            img_cord = None

        #extract the star index of the equatorial coordinates
        eq_entry = eq_cord_entry[index].get()
        
        if eq_entry != "":
            values = eq_entry.split(",")
            values = [elements.replace(" ", "") for elements in values]
            flag = all(element.isdigit() for element in values)
            if flag and len(values) == 6:
                eq_cord = (
                int(values[0]), int(values[1]), int(values[2]), int(values[3]), int(values[4]), int(values[5]))
            else:
                error_label_2.config(text="Eq Cord Error")
                error = True
        else:
            eq_cord = None

        if eq_cord and img_cord:
            img_cal_stars.append(img_cord) #the cal stars (px) are just the input coords
            eq_cal_stars.append(eq_cord) #the cal stars (eq) are also just the input coords

        if eq_cord and not img_cord or img_cord and not eq_cord:
            if not error:
                unfill_error = True

    if not error:
        if unfill_error:
            error_label_2.config(text="Uncomplete Grid")
        elif len(img_cal_stars) < 5:
            error_label_2.config(text="Enter more Stars")
            star_img_cord = img_cal_stars  # The star px coords are just the input ones
            star_eq_cord = eq_cal_stars    # The star eq coords are just the input ones
        else:
            error_label_2.config(text="")
            star_img_cord = img_cal_stars  # The star px coords are just the input ones
            star_eq_cord = eq_cal_stars    # The star eq coords are just the input ones
            
 
"""

PART 2: Zenith
    
    Part 2 this part does the alignment of the coordinate systems. 
    It also stores the position, date, times and performs 
    
    Defining function which: 
        a) stores the star objects and angles for the zenith alignment, 
        b) begins the zenith approximation
        c) creates the plots/visualisation 

"""


# /// ALIGNMENT FUNCTION /// # 
#--> This contains local position info required for coordinate conversion
#--> It reads in the image and creates a blank image of the same size
#--> Shows the input stars 

def comd_align():
    global star_objects
    global next_sep_flag
    global north_offset_ca
    global zenith_pos_ca
       
    star_objects_temp = []
       
    general_info = [file_path, position, current_date, current_time, arcsec, zenith_ap, north_offset_ap, height]
    
    print(general_info)
    
    if any(v is None for v in general_info):
        error_label_2.config(text="General Error")
        return
    elif len(star_img_cord) < 5:
        error_label_2.config(text="Enter more Stars")
        return
    else:
         error_label_2.config(text="")
    
    # /// POSITION AND TIME FOR LOCAL COORD SYSTEM \\\ # 
    
    lat, long = position
    location_obj = EarthLocation(lat=lat * u.deg, lon=long * u.deg, height=height * u.m) # utcoffset = +2 * u.hour
    year, month, day = current_date
    hour, minute, second = current_time
    datetime_obj = datetime.datetime(year, month, day, hour, minute, second)
    timezone = pytz.timezone("Europe/Oslo")
    aware = timezone.localize(datetime_obj)
    utc_offset = int(aware.utcoffset().total_seconds() / 3600)
    time_obj = Time(datetime_obj) - utc_offset * u.hour
 

     # /// STAR ANGLES FOR THE ZENITH ALIGNMENT \\\ #
    
    star_angle = []

    for index, star in enumerate(star_eq_cord):
        star_obj = ManualStar(star_img_cord[index], star) # creates an object --> img_coord, eq_coord          
        #creates instance of ManualStar object which has attributes img_pos, eq_pos and a method to get viewing angles
        angles = star_obj.star_viewing_angles(time_obj, location_obj)
        
        star_angle.append(angles[0]) #star_angle is a list of AZIMUTH angles
        
        star_objects_temp.append(star_obj)

    star_objects = star_objects_temp #Star Obj is a list of Azimuth of different stars. Star_ang is a list of AZ of stars
    
    #star_objects is now a list of all ManualStar object instances


    # /// ZENITH APPROXIMATION : STEP 1 \\\ #

    img = cv2.imread(file_path) #reads in the image file

    blank = np.zeros(img.shape, dtype='uint8')

    #This calls entirely on the Zenith_Approximation file.
    zenith_1 = Zenith(zenith_ap, star_angle, star_img_cord) #creates instance of the Zenith object which has as attributes an approx position, star_angles and star positions
    #arguments are zenith_ap, star_angle (azimuth?), star_img
    zenith_pos_ca = zenith_1.zenith_intersection()  # gets average zenith intersection
    north_vector = zenith_1.get_north_vector()      # returns the north vector
    north_offset_ca = 360 - slope_angle_vector(north_vector)

    # set the next_step_flag to activate the button for the adjustment cycle
    next_sep_flag = True

    # north offset is the counterclockwise angle between the vector [0, -1] (-1 because in images y=0 is on the top) and the north vector
    north_out_text.set(f'{north_offset_ca}')
    zenith_out_text.set(f'{zenith_pos_ca[0]}, {zenith_pos_ca[1]}')

    # plot the zenith intersections on the blank image
    blank = zenith_1.visualize(blank)  #plot the north vector, the calculated zenith, the zenith approx and the intersecting circles

    #for each item in the manually entered coords draw a circle
    for pos in star_img_cord:
        cv2.circle(blank, (int(pos[0]), int(pos[1])), int(2) * 10, (0, 255, 0), thickness=4) #BLUE #circle centred on the star img coord with a set radius
        cv2.line(blank, (int(pos[0]), int(pos[1])), zenith_pos_ca, (0, 255, 0), thickness=2) #BLUE #line connecting the calculated zenith with the star img coord

    adjust_next = Button(window, text="Adjustment", command=comd_destroy)
    adjust_next.configure(bg='#a6a6a6', fg='black')
    adjust_next.grid(row=17, column=4, columnspan=1, rowspan=1, sticky= 'wens', padx=5, pady=5)

    blank_resize = rescale_frame(blank, rescale_setting)
    cv2.imshow("INPUT_STARS_AND_ZENITH_INTERSECTION: ", blank_resize) ##Is this wrong??
    cv2.waitKey()
    #draw north vector on the image also



def comd_destroy():
    window.destroy()

# Star Entry Windows
label_ = Label(window, text="Image Coord")
label_.configure(bg='#404040', fg='white')
label_.grid(row=0, column=3, pady=5)
label_ = Label(window, text="Equatorial Coord")
label_.configure(bg='#404040', fg='white')
label_.grid(row=0, column=4, pady=5)


for label in range(1, 16):
    label_ = Label(window, text=f"Star {label} :")
    label_.configure(bg='#404040', fg='white')
    label_.grid(row=label, column=2, sticky='w', padx=0, pady=2)

    image_cord_label = Entry(window)
    image_cord_label.configure(bg='#999999', fg='black')
    image_cord_label.grid(row=label, column=3, padx=5)
    image_cord_entry.append(image_cord_label)

    eq_cord_label = Entry(window)
    eq_cord_label.configure(bg='#999999', fg='black')
    eq_cord_label.grid(row=label, column=4)
    eq_cord_entry.append(eq_cord_label)
 

# /// GENERAL INFORMATION ENTRY WINDOW \\\ # 

label_ = Label(window, text="General Information")
label_.configure(bg='#404040', fg='white')
label_.grid(row=0, column=1, pady=5)

label_names = ["File Path : ", "Position : ", "Date : ", "Time : ",
               "Arcsec_px : ", "Zenith AP : ", "North off AP : ", "Height : ", "Autofill file : ", "Image Rescale :"]
general_entry = []

for index, name in enumerate(label_names):
    label_ = Label(window, text=name)
    label_.configure(bg='#404040', fg='white')
    label_.grid(row=index + 1, column=0, padx=5)

    general_label = Entry(window)
    general_label.configure(bg='#999999', fg='black')
    general_label.grid(row=index + 1, column=1, padx=5)
    general_entry.append(general_label)

#OUTPUT WINDOWS
# Head
head_row = 15

output_label_2 = Label(window, text="Outputs")
output_label_2.configure(bg='#404040', fg='white')
output_label_2.grid(row=head_row, column=1, padx=5)

# Labels
label_ = Label(window, text="Zenith Pos : ")
label_.configure(bg='#404040', fg='white')
label_.grid(row=head_row + 1, column=0, padx=5)

zenith_out_text = StringVar()
zenith_out_label = Entry(window, textvariable=zenith_out_text)
zenith_out_label.configure(bg='#999999', fg='black')
zenith_out_label.grid(row=head_row + 1, column=1, padx=5)

label_ = Label(window, text="North Off : ")
label_.configure(bg='#404040', fg='white')
label_.grid(row=head_row + 2, column=0, padx=5)

north_out_text = StringVar()
north_out_label = Entry(window, textvariable=north_out_text)
north_out_label.configure(bg='#999999', fg='black')
north_out_label.grid(row=head_row + 2, column=1, padx=5)

# Error Labels (general)
error_label = Label(window, text="")
error_label.configure(bg='#404040', fg='white')
error_label.grid(row=len(label_names) + 1, column=1, padx=5)

unfill_label = Label(window, text="")
unfill_label.configure(bg='#404040', fg='white')
unfill_label.grid(row=len(label_names) + 2, column=1, padx=5)

# Star Error labels
error_label_2 = Label(window, text="")
error_label_2.configure(bg='#404040', fg='white')
error_label_2.grid(row=17, column=2, columnspan=2, padx=5)

# Command buttons
star_b = Button(window, text="Commit Stars", command=comd_star)
star_b.configure(bg='#a6a6a6', fg='black')
star_b.grid(row=16, column=3, columnspan=1, rowspan=1, sticky='wens', padx=5, pady=5)

general_b = Button(window, text="Commit Info", command=comd_general)
general_b.configure(bg='#a6a6a6', fg='black')
general_b.grid(row=len(label_names) + 3, column=1, columnspan=1, rowspan=1, sticky='wens', padx=5)

alignment_b = Button(window, text="Start Alignment", command=comd_align)
alignment_b.configure(bg='#a6a6a6', fg='black')
alignment_b.grid(row=16, column=4, columnspan=1, rowspan=1, sticky='wens', padx=5, pady=5)

mainloop()

if next_sep_flag:
    
    window = Tk()
    window.grid_propagate(False)
    window.title('Adjustment Cycle')
    position_window = f'+{800}+{800}'
    window_size = f'{800}x{800}'
    #window.resizable(0, 0)
    window.geometry(position_window)
    window.geometry(window_size)
    window.configure(bg='#404040')

    label_ = Label(window, text="Inputs")
    label_.configure(bg='#404040', fg='white')
    label_.grid(row=0, column=1, pady=5)

    label_ = Label(window, text="Adjust")
    label_.configure(bg='#404040', fg='white')
    label_.grid(row=0, column=2, pady=5)

    # Zenith X
    label_ = Label(window, text="Zenith X : ")
    label_.configure(bg='#404040', fg='white')
    label_.grid(row=1, column=0, pady=5)

    label_ = Label(window, text=f'{zenith_pos_ca[0]}')
    label_.configure(bg='#404040', fg='white')
    label_.grid(row=1, column=1, pady=5)

    entry_zenith_x_adjust = Entry(window)
    entry_zenith_x_adjust.configure(bg='#999999', fg='black')
    entry_zenith_x_adjust.grid(row=1, column=2, padx=5)

    # Zenith Y
    label_ = Label(window, text="Zenith Y : ")
    label_.configure(bg='#404040', fg='white')
    label_.grid(row=2, column=0, pady=5)

    label_ = Label(window, text=f'{zenith_pos_ca[1]}')
    label_.configure(bg='#404040', fg='white')
    label_.grid(row=2, column=1, pady=5)

    entry_zenith_y_adjust = Entry(window)
    entry_zenith_y_adjust.configure(bg='#999999', fg='black')
    entry_zenith_y_adjust.grid(row=2, column=2, padx=5)

    # North Offset
    label_ = Label(window, text="North Offset : ")
    label_.configure(bg='#404040', fg='white')
    label_.grid(row=3, column=0, pady=5)

    label_ = Label(window, text=f'{north_offset_ca}')
    label_.configure(bg='#404040', fg='white')
    label_.grid(row=3, column=1, pady=5)

    entry_offset_adjust = Entry(window)
    entry_offset_adjust.configure(bg='#999999', fg='black')
    entry_offset_adjust.grid(row=3, column=2, padx=5)

    # output values
    label_ = Label(window, text="Outputs")
    label_.configure(bg='#404040', fg='white')
    label_.grid(row=0, column=3, pady=5)

    # Labels
    zenith_out_x_text = StringVar()
    zenith_out_x_label = Entry(window, textvariable=zenith_out_x_text)
    zenith_out_x_label.configure(bg='#999999', fg='black')
    zenith_out_x_label.grid(row=1, column=3, padx=5)

    zenith_out_y_text = StringVar()
    zenith_out_y_label = Entry(window, textvariable=zenith_out_y_text)
    zenith_out_y_label.configure(bg='#999999', fg='black')
    zenith_out_y_label.grid(row=2, column=3, padx=5)

    north_out_text = StringVar()
    north_out_label = Entry(window, textvariable=north_out_text)
    north_out_label.configure(bg='#999999', fg='black')
    north_out_label.grid(row=3, column=3, padx=5)

    label_ = Label(window, text="Distortion Output: ")
    label_.configure(bg='#404040', fg='white')
    label_.grid(row=6, column=0, padx=5, pady=5)

    distortion_out_text = StringVar()
    distortion_out_label = Entry(window, textvariable=distortion_out_text, width=65)  # width=43 
    distortion_out_label.configure(bg='#999999', fg='black')
    distortion_out_label.grid(row=6, column=1, columnspan=3, padx=5)

    # Distortion error
    label_ = Label(window, text="Distortion Error: ")
    label_.configure(bg='#404040', fg='white')
    label_.grid(row=7, column=0, padx=5, pady=5)

    distortion_error_label = Label(window, text="")
    distortion_error_label.configure(bg='#404040', fg='white')
    distortion_error_label.grid(row=7, column=1, padx=5)

    # Angle error
    label_ = Label(window, text="Angle Error: ")
    label_.configure(bg='#404040', fg='white')
    label_.grid(row=7, column=2, padx=5, pady=5)

    angle_error_label = Label(window, text="")
    angle_error_label.configure(bg='#404040', fg='white')
    angle_error_label.grid(row=7, column=3, padx=5)


"""

PART 3: 
    
    ALIGNMENT
    
    a) The alignment function calculates the distortion & performs cubic fit (f'({a} * x) + ({b} * (x^2)) + ({c} * (x^3)) + ({d})') 
    to get parameters a, b, c, d. 
    b) It also calculates the errors associated with the fit
    

"""
# /// DISTORTION FUNCTION \\\ #
def comd_distortion():
    if distortion_flag:
        comd_alignment(distortion_adj=True)
        
        
# /// ALIGNMENT FUNCTION  \\\ # 
#--> process of aligning the input stars with the expected star positions for the current location and time
#--> The expected positions are based on precomputed star positions using astronomical algorithms.
def comd_alignment(distortion_adj=False):

    global distortion_parameters
    global distortion_flag
    global zenith_pos
    global north_offset

    zenith_x_adj = entry_zenith_x_adjust.get()
    if zenith_x_adj != "":
        zenith_x_adj = int(zenith_x_adj)
    else:
        zenith_x_adj = 0

    zenith_y_adj = entry_zenith_y_adjust.get()
    if zenith_y_adj != "":
        zenith_y_adj = int(zenith_y_adj)
    else:
        zenith_y_adj = 0

    offset_adj = entry_offset_adjust.get()
    if offset_adj != "":
        offset_adj = float(offset_adj)
    else:
        offset_adj = 0

    #calculated zenith position plus adjustment values (manual..?)
    zenith_pos_adj = zenith_pos_ca[0] + zenith_x_adj, zenith_pos_ca[1] + zenith_y_adj
    north_offset_adj = north_offset_ca + offset_adj

    # set the adjusted values globally for calculation
    zenith_pos = zenith_pos_adj
    north_offset = north_offset_adj

    # adjusted values output
    zenith_out_x_text.set(f'{zenith_pos_adj[0]}')
    zenith_out_y_text.set(f'{zenith_pos_adj[1]}')
    north_out_text.set(f'{north_offset_adj}')

    # red image for image-size
    img = cv2.imread(file_path)                  #reads in the image
    blank = np.zeros(img.shape, dtype='uint8')   #creates an empty image of same size

    # plot zenith
    cv2.circle(blank, zenith_pos_adj, 8, (255, 255, 255), thickness=6)
    image_center = int(img.shape[1] / 2), int(img.shape[0] / 2)  # image shape has different order

    px_dist = [0]
    dist_error = [0]
    angle_error = [0]

    for index, manual_star in enumerate(star_objects): #enumerated list of manual star objects
        # translate horizontal coordinates to image coordinates
        
        img_cord = manual_star.img_pos           #Manual Star 
        sky_cord = manual_star.viewing_angles    #defining altaz coords

        cv2.circle(blank, img_cord, 4, (255, 255, 255), thickness=6)  #circle around img coord

        calculated_px_position = hor_2_px(sky_cord, zenith_pos_adj, north_offset_adj, arcsec) #above definition is commented out, should be px_2_horizontal, but here says opposite, def already implemented in CalcStarData
        #????????????????

        #connect the position(horizontal) with img_coord of point, img zenith, image centre
        cv2.line(blank, calculated_px_position, img_cord, (255, 255, 0), thickness=2)
        cv2.line(blank, calculated_px_position, zenith_pos_adj, (150, 150, 150), thickness=2)
        cv2.line(blank, calculated_px_position, image_center, (255, 150, 0), thickness=2)

        # distortion error handling
        #DISTANCE 1
        dist = np.linalg.norm(np.subtract(image_center, img_cord)) #????
        #distance between img centre and img coord
        px_dist.append(dist) # figure out how far the star pixels are from the centre of the image
       
        #DISTANCE 2
        calc_pos_vec = np.subtract(image_center, calculated_px_position) #find a corresponding vector
        calc_px_dist = np.linalg.norm(calc_pos_vec)
        #distance between image centre and the calculated position (ie. calcuated from looking angles)
        
        error_value = abs(dist - calc_px_dist)
        dist_error.append(error_value)

        # angle error handling
        scalar = - abs(dist) / abs(calc_px_dist)
        error_pos = image_center + (calc_pos_vec * scalar)
        error_pos = int(error_pos[0]), int(error_pos[1])

        cv2.circle(blank, error_pos, 8, (255, 255, 255), thickness=2)              #WHITE

        error_value = np.linalg.norm(np.subtract(error_pos, img_cord))
        angle_error.append(error_value)
        
        print('Pixel distances:', px_dist) #this is a list of the distances between img center and the img coord of stars
        print('Dist Error', dist_error) #this is a list of the difference between the calculated pos and the actual pos

        # put the generated distortion function in here to check the results

        if distortion_adj is True and distortion_flag is True:
            distortion = (distortion_parameters[0] * dist) + (distortion_parameters[1] * (dist ** 2)) + (
                        distortion_parameters[2] * (dist ** 3)) + (distortion_parameters[3])
        else:
            distortion = 0

        # distortion error
        vector = np.subtract(img_cord, image_center)
        scalar = distortion / dist
        distortion_vector = vector - (vector * scalar)

        corrected_point = image_center + distortion_vector
        corrected_point = round(corrected_point[0]), round(corrected_point[1])

        cv2.circle(blank, corrected_point, 8, (0, 0, 255), thickness=2)            #RED - corrected point 
        cv2.circle(blank, calculated_px_position, 8, (0, 255, 0), thickness=2)     #GREEN

     # curve fit
    pyplot.scatter(px_dist, dist_error, s=10)            #scatter plot of the points generated
    coef, _ = curve_fit(objective, px_dist, dist_error)  #cubic fitting with arb distortion parameters
    a, b, c, d = coef
    distortion_out_text.set(f'{a}, {b}, {c}, {d}')

    # distortion parameters global variable
    distortion_parameters = coef
    distortion_flag = True

    # error from distortion (average variance in distortion function)
    dist_error_total = 0
    for y_index, x_pos in enumerate(px_dist):
        y_pos = objective(x_pos, a, b, c, d)
        dist_error_total += abs(y_pos - dist_error[y_index])

    average_dist_error = dist_error_total / len(px_dist)
    distortion_error_label.config(text=average_dist_error)
    # print(f'Distortion Function Error: {y_average_error}')

    # error from angle  (actual error caused by angle)
    average_angle_error = sum(angle_error) / (len(angle_error) - 1)  # -1 remove the 0
    angle_error_label.config(text=average_angle_error)

    x_line = arange(min(px_dist), max(px_dist), 1)
    y_line = objective(x_line, a, b, c, d)

    pyplot.plot(x_line, y_line, '--', color='red')
    pyplot.title(label="Distortion Function")
    pyplot.xlabel('Distance to image center (px)')
    pyplot.ylabel('Distortion (px)')

    blank_r = rescale_frame(blank, rescale_setting)
    cv2.imshow("Alignment", blank_r)
    pyplot.show()
    cv2.waitKey()
        

# error labels
error_label = Label(window, text="")
error_label.configure(bg='#404040', fg='white')
error_label.grid(row=4, column=1, padx=5)

# command buttons
align_check_button = Button(window, text="Alignment", command=comd_alignment)
align_check_button.configure(bg='#a6a6a6', fg='black')
align_check_button.grid(row=5, column=1, columnspan=1, rowspan=1, sticky='wens', padx=5, pady=5)



"""
PART 4: 
    
    CALIBRATION:
        
    a) The calibration test function searches for stars in the database
    b) Processes the image and detects stars this way 
    c) Matches the detected stars with with stars in the database.
    d) Plots the visualisation of this method

"""

# /// CALIBRATION FUNCTION \\\ # 
def comd_cal_test():

    global position
    global distortion_parameters
    global distortion_flag
    global zenith_pos
    global north_offset

    zenith_x_adj = entry_zenith_x_adjust.get()
    if zenith_x_adj != "":
        zenith_x_adj = int(zenith_x_adj)
    else:
        zenith_x_adj = 0

    zenith_y_adj = entry_zenith_y_adjust.get()
    if zenith_y_adj != "":
        zenith_y_adj = int(zenith_y_adj)
    else:
        zenith_y_adj = 0

    offset_adj = entry_offset_adjust.get()
    if offset_adj != "":
        offset_adj = float(offset_adj)
    else:
        offset_adj = 0

    zenith_pos_adj = zenith_pos_ca[0] + zenith_x_adj, zenith_pos_ca[1] + zenith_y_adj
    north_offset_adj = north_offset_ca + offset_adj

    # set the adjusted values globally for calculation
    zenith_pos = zenith_pos_adj
    north_offset = north_offset_adj

    # adjusted values output
    zenith_out_x_text.set(f'{zenith_pos_adj[0]}')
    zenith_out_y_text.set(f'{zenith_pos_adj[1]}')
    north_out_text.set(f'{north_offset_adj}')

    if not zenith_pos or not north_offset:
        return

    # open database
    ra = open(r'C:/Users/Viesis/Documents/Thesis/Neils/New/hygfull (2) (1).txt', 'r')
    dec = open(r'C:/Users/Viesis/Documents/Thesis/Neils/New/hygfull (1) (1).txt', 'r')
    mag = open(r'C:/Users/Viesis/Documents/Thesis/Neils/New/hygfull (2) (2).txt', 'r')

    magnitude_threshold = 6
    
    # Star Calc location and time
    lat, long = position
    location_object = EarthLocation(lat=lat * u.deg, lon=long * u.deg, height=height * u.m)

    utcoffset = +2 * u.hour    # why +2 ?
    year, month, day = current_date
    hour, minute, second = current_time
    time_obj = Time(datetime.datetime(year, month, day, hour, minute, second)) - utcoffset
    
    #print('Position is:', position)
   # print('Location is:', location_object)
   # print("Time is:", time_obj)
    
    # load the image
    img = cv2.imread(file_path)

    # /// STAR POS CALCULATION ///
    # open the database with right-assertion and declination
    star_calc_objects = []
    ra_contents = ra.readlines()
    dec_contents = dec.readlines()

    for i, line in enumerate(mag):
        mag_float = float(line)
        if mag_float < magnitude_threshold:
            star_name = None
            star_ra = float(ra_contents[i]) * 15  # why times 15?
            declination = float(dec_contents[i])
            #print('current date:', current_date)
            #print('current time:', current_time)
            star_obj = CalcStarData(star_name, star_ra, declination, position, current_date, current_time)
            star_obj.star_viewing_angles(time_obj, location_object)
            # translate angle to image coordinates
           # v2 = star_obj.image_coordinates(arcsec, north_offset, zenith_pos)
            star_calc_objects.append(star_obj)

    # /// STAR IMAGE DETECTION ///
    # filter, blur, edge-cascade
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    # ret, thresh1 = cv2.threshold(gray, 125, 255, cv2.THRESH_BINARY)   # try with and without thresholding
    blur = cv2.GaussianBlur(gray, (3, 3), 0)
    canny = cv2.Canny(blur, 125, 175)  # 125, 175

    # find contours
    contours, _ = cv2.findContours(canny, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

    # Circle detection on all detected contours, sort them by size
    detected_centers = []
    detected_radius = []
    
    for i, c in enumerate(contours):
        # don't create star objects in here to get rid of double stars
        contours_poly = cv2.approxPolyDP(c, 3, False)  # try true and false
        center, radi = cv2.minEnclosingCircle(contours_poly)
        if 1 < radi < 7:  # Filter for only bigger objects (stars)
            detected_centers.append(center)
            detected_radius.append(radi)

    # filters out some double detected stars (needs improvement)
    image_detected_stars = []
    image_detected_radius = []
    
    for i, center in enumerate(detected_centers):
        already_star = False
        if image_detected_stars:
            for star in image_detected_stars:
                dist = np.linalg.norm(np.subtract(center, star))
                if dist < 5:
                    already_star = True
                    break
            if not already_star:
                image_detected_stars.append(center)
                image_detected_radius.append(detected_radius[i])
        else:
            image_detected_stars.append(center)
            image_detected_radius.append(detected_radius[i])

    # initialise all image detected stars (detected_circles) as objects from Star class
    ImgStar.arcseconds_per_px = arcsec #all initiated as nonetype
    ImgStar.north_offset = north_offset
    ImgStar.image_zenith = zenith_pos

    star_img_objects = []
    
    image_center = (round(img.shape[1] * 0.5), round(img.shape[0] * 0.5))
    
    for star_center in image_detected_stars:
        img_star_object = ImgStar(star_center, image_center)
        img_star_object.horizontal_coord(distortion_parameters)
        #img_star_object.image_zenith= zenith_pos
        
        star_img_objects.append(img_star_object)
    
    # /// Match the Star Positions and indentify stars ///
    
    #matched_stars = []
    
    for i, img_star_obj in enumerate(star_img_objects):
        old_dist = 10 ** 10  # random high number
        img_star_obj.star_size = detected_radius[i]
        position_trans = img_star_obj.translated_coord
        
        for star in star_calc_objects:
            # given in degrees not px-pos
            if star.star_pos[0] - 1.2 <= position_trans[0] <= star.star_pos[0] + 1.2:  # Still need to handle stars close to north!!!!
                # az is way more accurate than elevation (deviation = 1)
                if star.star_pos[1] - 1.2 <= position_trans[1] <= star.star_pos[1] + 1.2:  # el is less accurate thus (deviation = 2)
                    img_star_obj.detected_Flag = True
                    dist = np.linalg.norm(np.subtract(star.star_pos, position_trans))
                    if dist < old_dist:
                        name = star.star_name
                    old_dist = dist
        if img_star_obj.detected_Flag:
            img_star_obj.star_name = name

    # /// Visualisation ///
    blank = np.zeros(img.shape, dtype='uint8')

    # Image center
    cv2.circle(blank, image_center, 14, (255, 255, 0), thickness=4)

    # Image detected stars
    for img_star in star_img_objects:
        pos = img_star.image_coordinates
        pos = round(pos[0]), round(pos[1])
        if img_star.detected_Flag:
            cv2.circle(blank, pos, round(img_star.star_size) * 2, (255, 255, 255), thickness=2)
            cv2.circle(blank, pos, 18, (0, 255, 0), thickness=1)
            cv2.line(blank, pos, image_center, (255, 255, 0), thickness=1)

        else:
            cv2.circle(blank, pos, round(img_star.star_size) * 2, (255, 255, 255), thickness=4)

    # Image corrected Stars
    for img_star in star_img_objects:
        pos = img_star.corrected_coordinates
        pos = round(pos[0]), round(pos[1])
        if img_star.detected_Flag:
            cv2.circle(blank, pos, 14, (255, 255, 0), thickness=4)

    # Calc stars
    for calc_star in star_calc_objects:
        pos = calc_star.horizontal_coord
        cv2.circle(blank, pos, 12, (0, 0, 255), thickness=2)

    # Plot the coordinate system
    for elevation in range(0, 18):
        if elevation % 2 == 0:
            cv2.circle(blank, zenith_pos, round(((elevation * 5) * 3600) / arcsec), (200, 200, 200),
                       thickness=2)
        else:
            cv2.circle(blank, zenith_pos, round(((elevation * 5) * 3600) / arcsec), (255, 255, 255),
                       thickness=1)
    cv2.circle(blank, zenith_pos, 20, (255, 255, 255), thickness=2)

    blank_r = rescale_frame(blank, rescale_setting)
    cv2.imshow("Star Identification ", blank_r)
    cv2.waitKey()

distortion_button = Button(window, text="Distortion", command=comd_distortion)
distortion_button.configure(bg='#a6a6a6', fg='black')
distortion_button.grid(row=5, column=2, columnspan=1, rowspan=1, sticky='wens', padx=5, pady=5)

distortion_button = Button(window, text="Calibration test", command=comd_cal_test)
distortion_button.configure(bg='#a6a6a6', fg='black')
distortion_button.grid(row=5, column=3, columnspan=1, rowspan=1, sticky='wens', padx=5, pady=5)

mainloop()
