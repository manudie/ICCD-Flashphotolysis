# -*- coding: utf-8 -*-
"""
This script was written to evaluate 1024x1024 px images (from the Andor iStar ICCD Camera).
It allows the user to generate a heatmap and/or a spectrum images from a specified area of lines in the exported sensor data.
It is essential, to export the data as an .asc (ASCII) file and append measurement information (tick the box in the export window (!)) in the Andor SOLIS Software. 
The exportet data than just needs to be in the same folder as this evaluation script. It is possible to place multiple files of raw .asc data into to folder, the script will iterate through every single one.
To decide wether to plot a heatmap, a spectrum or both just create an object from the class iccd_evaluation("String") and change "String" to "heatmap", "spectrum", or "both", e.g. plot1 = iccd_evaluation("both").
The spectrum is calculated as mean from rows in the image. By changing the values of self.mean_row_start and self.mean_row_end you can determine the span of rows, the  mean is calculated from.
By switching self.single_row to True, it is possible to generate the spectrum from a single row of pixels from the image, instead of mulitple rows. You can select the single row by changing the value of self.input_row.
Calling the function plot1.evaluate() then starts the evaluation. Note that only in case of "both" the area wherefrom the mean for the spectrum is created is shown in red on the heatmap. To avoid the red area, just plot the heatmap and spectrum individual. 
The script will create folders for the plotted images, aswell as for processed data and will sort the files into these folders after it finished evaluating them. 

"""

import os
import time
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

class iccd_evaluation():
    def __init__(self, mode):
        self.mode = mode
        self.file = ""
        self.single_row = False
        self.input_row = 255
        self.mean_row_start = 740
        self.mean_row_end = 760
        self.peak_height_min = 1500
        self.peak_distance = 20
        self.wavelenght_cal_1 = 404.6565
        self.wavelenght_cal_2 = 435.8335
        self.wavelenght_cal_3 = 546.0750

    def read_file(self):
        self.ascii_grid = pd.read_table(self.file, engine="python", skipfooter=29, index_col=1, usecols=[*range(1,1024)], header=None) 
        self.ascii_grid_transposed = self.ascii_grid.T
        self.image_filename = self.file[:len(self.file)-4]
        try:
            self.title = self.image_filename.rsplit("/",1)[1]
        except:
            try:
                self.title = self.image_filename.rsplit("\\",1)[1]
            except:
                self.title = self.image_filename
                
    def get_calibration(self):
        #peaks = find_peaks(self.spectrum_list, height=self.peak_height_min, distance=self.peak_distance)[0]
        #peak_cal_1 = peaks[0]
        #peak_cal_2 = peaks[1]
        #peak_cal_3 = peaks[2]        
        #b1 = (self.wavelenght_cal_2)/(((self.wavelenght_cal_1/peak_cal_1)-1)*(peak_cal_2+1))
        #m1 = (self.wavelenght_cal_1/peak_cal_1)-b
        #m = self.wavelenght_cal_2/(peak_cal_2-peak_cal_1)
        #b = self.wavelenght_cal_1-(m*peak_cal_1)
        #b = (self.wavelenght_cal_3-((self.wavelenght_cal_2*peak_cal_3)/peak_cal_2))/(1-(peak_cal_3/peak_cal_2))
        #m = (self.wavelenght_cal_2-b)/peak_cal_2
        #print(b)
        #print(m)
        m = 0.30423     #2-Punkt: 0.30453
        b = 386.33293   #2-Punkt: 386.19437
        index_list = np.arange(1,len(self.ascii_grid_transposed.columns)+1)
        self.wavelengths = list(map(lambda x: m*x+b, index_list))
    
    def plot_heatmap(self):
        plt.pcolor(self.ascii_grid_transposed)
        plt.colorbar(pad=0.01)
        plt.tick_params(
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False, 
            left=False,        # ticks along the top edge are off
            labelbottom=False, # labels along the bottom edge are off
            labelleft=False)
        plt.title(label=self.title, pad=-262, loc="left", color="white")
        if self.mode == "both":
            plt.axhspan(self.mean_row_start, self.mean_row_end, color="red", alpha=0.3, lw=0)
        plt.savefig(self.image_filename + "_heatmap", bbox_inches="tight")
        plt.clf()

    def plot_spectrum(self):
        if self.single_row == True:
            spectrum = self.ascii_grid_transposed.iloc[self.input_row]#[self.input_row:self.input_row+1]  # besser: to_numpy() 
        else: 
            spectrum = self.ascii_grid_transposed.iloc[self.mean_row_start:self.mean_row_end].mean()
        self.spectrum_list = spectrum.to_numpy()
        self.get_calibration()
        plt.plot(self.wavelengths, self.spectrum_list)
        plt.title(label=self.title, pad=-262, loc="left")
        plt.xlabel("wavelength λ / nm")
        plt.ylabel("counts")
        #plt.show()
        plt.savefig(self.image_filename + "_spectrum", bbox_inches="tight")
        plt.clf()

    def evaluate_mode(self):
        if self.mode == "spectrum":
            self.evaluate_spectrum()
        elif self.mode == "heatmap":
            self.evaluate_heatmap()
        elif self.mode == "both":
            self.evaluate_spectrum()
            self.evaluate_heatmap()

    def evaluate_heatmap(self):
        self.read_file()
        self.plot_heatmap()     

    def evaluate_spectrum(self):
        self.read_file()
        self.plot_spectrum() 
 
    def evaluate(self):
        directory = os.getcwd()
        if not os.path.isdir("processed_files"):
            os.mkdir("processed_files")
        if not os.path.isdir("iccd_heatmaps") and (self.mode == "heatmap" or self.mode == "both"):
            os.mkdir("iccd_heatmaps")
        if not os.path.isdir("spectra") and (self.mode == "spectrum" or self.mode == "both"):
            os.mkdir("spectra")
        time.sleep(1)
        for entry in os.scandir(directory):
            if entry.path.endswith(".asc") and entry.is_file():
                self.file = entry.path
                self.evaluate_mode()
                shutil.move(self.file, "processed_files")
                try:
                    if self.mode == "spectrum":
                        shutil.move(self.image_filename + "_spectrum.png", "spectra")
                    elif self.mode == "heatmap":
                        shutil.move(self.image_filename + "_heatmap.png", "iccd_heatmaps")
                    elif self.mode == "both":
                        shutil.move(self.image_filename + "_spectrum.png", "spectra")
                        shutil.move(self.image_filename + "_heatmap.png", "iccd_heatmaps")
                except:
                    #os.remove(self.image_filename + ".png")
                    print("Failed to safe data!")
        

if __name__ =='__main__':
    plot1 = iccd_evaluation("both")
    plot1.evaluate()
    