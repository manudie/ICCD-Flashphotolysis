# -*- coding: utf-8 -*-
"""
This script was written to evaluate 1024x1024 px images (from the Andor iStar ICCD Camera).
It allows the user to generate a heatmap and/or a spectrum images from a specified area of lines in the exported sensor data.
It is essential, to export the data as an .asc (ASCII) file and append acquisition information (tick the box in the export window (!)) in the Andor SOLIS Software. 
The exportet data than just needs to be in the same folder as this evaluation script. It is possible to place multiple files of raw .asc data into to folder, the script will iterate through every single one.
To decide wether to plot a heatmap, a spectrum or both just create an object from the class iccd_evaluation("String") and change "String" to "heatmap", "spectrum", or "both", e.g. plot1 = iccd_evaluation("both").
The spectrum is calculated as mean from rows in the image. By changing the values of self.mean_row_start and self.mean_row_end you can determine the span of rows, the  mean is calculated from.
By switching self.single_row to True, it is possible to generate the spectrum from a single row of pixels from the image, instead of mulitple rows. You can select the single row by changing the value of self.input_row.
Calling the function plot1.evaluate() then starts the evaluation. Note that only in case of "both" the area wherefrom the mean for the spectrum is created is shown in red on the heatmap. To avoid the red area, just plot the heatmap and spectrum individual. 
The script will create folders for the plotted images, aswell as for processed data and will sort the files into these folders after it finished evaluating them. 

"""

### import packages ###
import os
import json
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

### define class ###
class iccd_evaluation():
    
    def __init__(self, mode):
        ### constructor that declares class variables (self.x) ###
        self.mode = mode
        self.test_run = False            # If self.test_run is True, Images will not be safed, but jut displayed and also files will not be moved. Made for easier developing.
        self.file = ""
        self.readout_mode = ""          # Supported: "Full Resolution Image" ("FRI") or "Single Track" ("ST")
        self.single_row = False
        self.input_row = 255
        self.mean_row_start = 500
        self.mean_row_end = 900
        self.peak_height_min = 2000
        self.peak_distance = 20
        self.wavelenght_cal_1 = 404.6565
        self.wavelenght_cal_2 = 435.8335
        self.wavelenght_cal_3 = 546.0750
        self.directory = os.getcwd()

    
    def read_file(self):
        ### imports the ASCII (.asc) file exportet from the Andor SOLIS software and generates a pandas dataframe from it ###
        with open(self.file) as f:
            if "Single Track" in f.read():
                self.readout_mode = "ST"
                self.ascii_grid = pd.read_table(self.file, engine="python", skipfooter=29, index_col=0, header=None).dropna(axis=1, how="all")
            else:
                f.seek(0)
                self.readout_mode = ""
            if self.readout_mode != "ST":
                if "Full Resolution Image" in f.read():
                    self.readout_mode = "FRI"
                    self.ascii_grid = pd.read_table(self.file, engine="python", skipfooter=29, index_col=0, header=None).dropna(axis=1, how="all") #usecols=[*range(1,1024)],
                else: 
                    print('Unsupported readout mode! Choose "Full Resolution Image" or "Single Track"')
            f.close()
        #print(self.ascii_grid)
        self.ascii_grid_transposed = self.ascii_grid.T
        self.image_filename = self.file[:len(self.file)-4]
        try:
            self.title = self.image_filename.rsplit("/",1)[1]
        except:
            try:
                self.title = self.image_filename.rsplit("\\",1)[1]
            except:
                self.title = self.image_filename

    def calibrate(self):
        self.get_calibration_file()
        self.read_file()
        self.calculate_spectrum()
        peaks = find_peaks(self.spectrum_list, height=self.peak_height_min, distance=self.peak_distance)[0]
        print("Peaks found over " + str(self.peak_height_min) + " counts at columns: " + str(peaks))
        wavelenght_cal_list = [self.wavelenght_cal_1, self.wavelenght_cal_2, self.wavelenght_cal_3]
        m,b = np.polyfit(peaks,wavelenght_cal_list,1)
        calibration_dict = {
            "wavelenghts (y):": "columns (x):",
            self.wavelenght_cal_1: int(peaks[0]),
            self.wavelenght_cal_2: int(peaks[1]),
            self.wavelenght_cal_3: int(peaks[2]),
            "slope and intercept:": "values:",
            "m" : float(m),
            "b" : float(b)
        } 
        with open("calibration_file.json", "w", encoding="utf-8") as f:
            json.dump(calibration_dict, f, ensure_ascii=False, indent=4)
        

    def get_calibration(self):
        try:
            f = open("calibration_file.json")
            calibration_dict_in = json.load(f)
            m = calibration_dict_in.get("m")
            b = calibration_dict_in.get("b")
            f.close()
        except:
            print("No calibration file found! Spectrum will be generated with false values!")
            m = 0.3033110988290711   
            b = 386.25945492645576   
        index_list = np.arange(1,len(self.ascii_grid_transposed.columns)+1)
        self.wavelengths = list(map(lambda x: m*x+b, index_list))
    
    def plot_heatmap(self):
        ### plots the generated dataframe to a heatmap and saves it to the current folder ###
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
        if self.mode == "both" and self.readout_mode == "FRI":
            plt.axhspan(self.mean_row_start, self.mean_row_end, color="red", alpha=0.3, lw=0)
        if self.test_run == True: 
            plt.show()
        else:
            plt.savefig(self.image_filename + "_heatmap", bbox_inches="tight")
        plt.clf()

    def calculate_spectrum(self):
        if self.readout_mode == "ST":
            spectrum = self.ascii_grid_transposed.mean()
        elif self.readout_mode == "FRI":
            if self.single_row == True:
                spectrum = self.ascii_grid_transposed.iloc[self.input_row]
            else: 
                spectrum = self.ascii_grid_transposed.iloc[self.mean_row_start:self.mean_row_end].mean()
        self.spectrum_list = spectrum.to_numpy()

    def plot_spectrum(self):
        ### plots the generated dataframe to a spectrum with wavelengths on the x-axis from get_calibration() and saves it to the current folder ###
        self.calculate_spectrum()
        self.get_calibration()
        plt.plot(self.wavelengths, self.spectrum_list)
        plt.title(label=self.title, pad=-262, loc="left")
        plt.xlabel("wavelength λ / nm")
        plt.ylabel("counts")
        if self.test_run == True: 
            plt.show()
        else:
            plt.savefig(self.image_filename + "_spectrum", bbox_inches="tight")
        plt.clf()

    def evaluate(self):
        ### decides what to plot, given on the parameter ("heatmap", "spectrum", or "both") the class is called with ###
        if self.mode == "spectrum":
            self.plot_spectrum()
        elif self.mode == "heatmap":
            self.plot_heatmap()
        elif self.mode == "both":
            self.plot_spectrum()
            self.plot_heatmap()

    def iterate(self):
        self.make_dirs()
        for entry in os.scandir(self.directory):
            if entry.path.endswith(".asc") and entry.is_file():
                self.file = entry.path
                self.read_file()
                self.evaluate()
                self.move_files()
        
    def get_calibration_file(self):
        asc_file_count =len([f for f in os.listdir(self.directory) if f.endswith('.asc') and os.path.isfile(os.path.join(self.directory, f))])
        if asc_file_count == 1:
            file_list = [f for f in os.listdir(self.directory) if f.endswith(".asc") and os.path.isfile(os.path.join(self.directory, f))]
            self.file = "".join(file_list)
        elif asc_file_count > 1:
            print("Only put one calibration spectrum file (.asc) from calibration lamp in this folder to calibrate!")
        elif asc_file_count == 0:
            print("Put exactly one calibration spectrum file (.asc) from calibration lamp in this folder to calibrate!")

    def make_dirs(self):
        if not os.path.isdir("processed_files"):
            os.mkdir("processed_files")
        if not os.path.isdir("iccd_heatmaps") and (self.mode == "heatmap" or self.mode == "both"):
            os.mkdir("iccd_heatmaps")
        if not os.path.isdir("spectra") and (self.mode == "spectrum" or self.mode == "both"):
            os.mkdir("spectra")

    def move_files(self):
        ### checks if there are folders for the generated images and processed data and creates them, if not. 
        # Then scans and iterates through the current directory to evaluate and subsequently move the files to the folders one after another ###
        if self.test_run != True: 
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
    plot1 = iccd_evaluation("both")         # calling an object from the class iccd_evaluation("String") with the parameters "heatmap", "spectrum", or "both"
    plot1.iterate()                        # start evaluating process by calling the function evaluate()
    #plot1.calibrate()
   