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
from scipy import optimize
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
import matplotlib.ticker as ticker

### define class ###
class iccd_evaluation():
    
    def __init__(self, mode):
        ### constructor that declares class variables (self.x) ###
        self.mode = mode
        self.test_run = True            # If self.test_run is True, Images will not be safed, but displayed and also files will not be moved. Made for easier developing.
        self.show_cal_marks = True
        self.drop_first_measurement = True             # In kinetic series w/ single track, often the first line of data is false due to build up charge in the ccd. Setting self.drop_first_measurent to True drops this line of data. Note to acquire n+1 mesaurements! 
        self.file = ""
        self.readout_mode = ""          # Supported: "Full Resolution Image" ("FRI") or "Single Track" ("ST")
        self.single_row = False
        self.input_row = 255
        self.mean_row_start = 250
        self.mean_row_end = 550
        self.peak_height_min = 2500
        self.peak_distance = 20
        self.wavelenght_cal_1 = 404.6565
        self.wavelenght_cal_2 = 435.8335
        self.wavelenght_cal_3 = 546.0750
        self.directory = os.getcwd()
        self.AD_filename_list = ["_I_","_I0_","_D_"]
        self.AD_dict = {}
    
    def iterate(self):
        ### main function that scans and iterates through the current directory to evaluates each .asc file after another ###
        self.make_dirs()
        for entry in os.scandir(self.directory):
            if entry.path.endswith(".asc") and entry.is_file():
                self.file = entry.path
                self.read_file()
                if self.mode == "spectrum" or self.mode == "heatmap" or self.mode == "both":
                    self.evaluate()
                    self.move_files()
                elif self.mode == "DA" or self.mode == "A":
                    for el in self.AD_filename_list:
                        if el in self.file:
                            self.calculate_spectrum()
                            self.AD_dict.update({el.replace("_",""):self.spectrum_list})
                    self.move_files()        
        if self.mode == "DA" or self.mode == "A":
            self.calculate_diff_absorbance()
            self.plot_spectrum()
            self.move_files() 
            
    def evaluate(self):
        ### decides what to plot, given on the parameter ("heatmap", "spectrum", or "both") the class is called with ###
        if self.mode == "spectrum":
            self.calculate_spectrum()
            self.plot_spectrum()
        elif self.mode == "heatmap":
            self.plot_heatmap()
        elif self.mode == "both":
            self.calculate_spectrum()
            self.plot_spectrum()
            self.plot_heatmap()
    
    def read_file(self):
        ### imports the ASCII (.asc) file exportet from the Andor SOLIS software and generates a pandas dataframe from it ###
        with open(self.file) as f:
            if "Single Track" in f.read():
                self.readout_mode = "ST"
                self.ascii_grid = pd.read_table(self.file, engine="python", skipfooter=29, index_col=0, header=None).dropna(axis=1, how="all")
                self.ascii_grid = self.ascii_grid.drop([1,2,3,4],axis=0)
                if self.drop_first_measurement == True:
                    self.ascii_grid = self.ascii_grid.drop([1],axis=1)
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

    def make_dirs(self):
        ### checks if there are folders for the generated images and processed data and creates them, if not ###
        if not os.path.isdir("processed_files"):
            os.mkdir("processed_files")
        if not os.path.isdir("iccd_heatmaps") and (self.mode == "heatmap" or self.mode == "both"):
            os.mkdir("iccd_heatmaps")
        if not os.path.isdir("spectra") and (self.mode == "spectrum" or self.mode == "both" or self.mode == "DA" or self.mode == "A"):
            os.mkdir("spectra")

    def move_files(self):
        ### move the files to the folders one after another ###
        if self.test_run != True: 
            try:
                shutil.move(self.file, "processed_files")
            except:
                pass
            try:
                if self.mode == "spectrum" or self.mode == "DA" or self.mode == "A":
                    shutil.move(self.image_filename + "_spectrum.png", "spectra")
                elif self.mode == "heatmap":
                    shutil.move(self.image_filename + "_heatmap.png", "iccd_heatmaps")
                elif self.mode == "both":
                    shutil.move(self.image_filename + "_spectrum.png", "spectra")
                    shutil.move(self.image_filename + "_heatmap.png", "iccd_heatmaps")
            except:
                #os.remove(self.image_filename + ".png")
                #print("Failed to safe data!")
                pass
    
    def calculate_spectrum(self):
        ### calculates the mean of the data to a spectrum, based on the readout mode the data was aquired. 
        # Supported readout modes: "Full Resolution Image" ("FRI") or "Single Track" ("ST") ###
        if self.readout_mode == "ST":
            spectrum = self.ascii_grid_transposed.mean()
        elif self.readout_mode == "FRI":
            if self.single_row == True:
                spectrum = self.ascii_grid_transposed.iloc[self.input_row]
            else: 
                spectrum = self.ascii_grid_transposed.iloc[self.mean_row_start:self.mean_row_end].mean()
        self.spectrum_list = spectrum.to_numpy()
    
    def calculate_diff_absorbance(self):
        I = self.AD_dict.get("I")
        I0 = self.AD_dict.get("I0")
        D = self.AD_dict.get("D")
        self.diff_absorbance_list = -np.log10((I-D)/(I0-D))
    
    def plot_spectrum(self):
        ### plots the generated dataframe to a spectrum with wavelengths on the x-axis from get_calibration() and saves it to the current folder ###
        if self.mode == "DA" or self.mode == "A":
            spectrum = self.diff_absorbance_list
        else:
            spectrum = self.spectrum_list
        self.get_calibration()           
        plt.plot(self.wavelengths, spectrum, linewidth=1)
        plt.title(label=self.title, pad=-262, loc="left")
        plt.xlabel("wavelength λ / nm")
        if self.mode == "DA":
            plt.ylabel("difference absorbance ΔA")
        elif self.mode == "A":
            plt.ylabel("absorbance A")
        else:
            plt.ylabel("counts")
        if self.show_cal_marks == True:
            plt.vlines([self.wavelenght_cal_1, self.wavelenght_cal_2, self.wavelenght_cal_3], ymin=self.spectrum_list.min(), ymax=self.spectrum_list.max(), color="grey", linewidth=0.5)
        if self.test_run == True: 
            plt.show()
        else:
            plt.savefig(self.image_filename + "_spectrum", bbox_inches="tight")
        plt.clf()
    
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
                           
    def get_calibration_file(self):
        ### function that reads in only one calibration .asc file for calibrate() mode, instead of iterating through the directory in iterate() mode. ###
        asc_file_count =len([f for f in os.listdir(self.directory) if f.endswith('.asc') and os.path.isfile(os.path.join(self.directory, f))])
        if asc_file_count == 1:
            file_list = [f for f in os.listdir(self.directory) if f.endswith(".asc") and os.path.isfile(os.path.join(self.directory, f))]
            self.file = "".join(file_list)
        elif asc_file_count > 1:
            print("Only put one calibration spectrum file (.asc) from calibration lamp in this folder to calibrate!")
        elif asc_file_count == 0:
            print("Put exactly one calibration spectrum file (.asc) from calibration lamp in this folder to calibrate!")

    def get_calibration(self):
        ### Reads the calibration data from the .json calibration file or uses default values for the calibration, if there is no file. 
        # The function than maps the column numbers over the calibration funtion to calculate wavelenghts. 
        # Note, that you must use a valid and up to date calibration file to get correct wavelenght values! ###
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
        #print(self.wavelengths)
    
    def _1Lorentzian(self,x, amp, cen, wid):
        return amp*wid**2/((x-cen)**2+wid**2)

    def _3Lorentzian(self,x, amp1, cen1, wid1, amp2,cen2,wid2, amp3,cen3,wid3):
        return (amp1*wid1**2/((x-cen1)**2+wid1**2)) +\
                (amp2*wid2**2/((x-cen2)**2+wid2**2)) +\
                    (amp3*wid3**2/((x-cen3)**2+wid3**2))
    
    def calibrate(self):
        ### In calibrate() mode this function will search for the column numbers of three peaks in the spectrum of the calibration lamp. 
        # The peaks will be determined by their height, set by self.peak_height_min. 
        # It then uses a three-point linear fit to calculate values for the slope m and intercept b, which are stored together with the wavelenghts and depending column values in a .json file for later use. ###
        self.get_calibration_file()
        self.read_file()
        self.calculate_spectrum()
        peaks = find_peaks(self.spectrum_list, height=self.peak_height_min, distance=self.peak_distance)[0]
        print(peaks)
        if peaks.size > 3:
            peaks = find_peaks(self.spectrum_list, height=self.peak_height_min, distance=self.peak_distance, width=40)[0]
            print(peaks)
            
            amp1 = 4000
            cen1 = 80
            wid1 = 20

            amp2 = 13000
            cen2 = 150
            wid2 = 10

            amp3 = 50
            cen3 = 200
            wid3 = 10

            index_list = np.arange(1,len(self.ascii_grid_transposed.columns)+1) #own code
            print(self.ascii_grid_transposed)
            print(index_list)

            popt_3lorentz, pcov_3lorentz = optimize.curve_fit(self._3Lorentzian, index_list, self.spectrum_list, p0=[amp1, cen1, wid1, \
                                                                                    amp2, cen2, wid2, amp3, cen3, wid3])
            perr_3lorentz = np.sqrt(np.diag(pcov_3lorentz))

            pars_1 = popt_3lorentz[0:3]
            pars_2 = popt_3lorentz[3:6]
            pars_3 = popt_3lorentz[6:9]
            lorentz_peak_1 = self._1Lorentzian(index_list, *pars_1)
            lorentz_peak_2 = self._1Lorentzian(index_list, *pars_2)
            lorentz_peak_3 = self._1Lorentzian(index_list, *pars_3)     

            print(lorentz_peak_1)   
            print(lorentz_peak_2)
            print(lorentz_peak_3)
            
            fig = plt.figure(figsize=(4,3))
            gs = gridspec.GridSpec(1,1)
            ax1 = fig.add_subplot(gs[0])

            ax1.plot(index_list, self.spectrum_list)

            ax1.plot(index_list, lorentz_peak_1, "g")
            ax1.fill_between(index_list, lorentz_peak_1.min(), lorentz_peak_1, facecolor="green", alpha=0.2)
            ax1.plot(index_list, lorentz_peak_2, "r")
            ax1.fill_between(index_list, lorentz_peak_2.min(), lorentz_peak_2, facecolor="red", alpha=0.2)  
            ax1.plot(index_list, lorentz_peak_3, "y")           #c
            ax1.fill_between(index_list, lorentz_peak_3.min(), lorentz_peak_3, facecolor="yellow", alpha=0.2)           #cyan
            #ax1.set_xlim(-5,105)
            #ax1.set_ylim(-0.5,5)

            #ax1.set_xlabel("x_array",family="serif",  fontsize=12)
            #ax1.set_ylabel("y_array",family="serif",  fontsize=12)

            #ax1.xaxis.set_major_locator(ticker.MultipleLocator(50))
            #ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))

            ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax1.yaxis.set_minor_locator(AutoMinorLocator(2))

            ax1.tick_params(axis='both',which='major', direction="in", top="on", right="on", bottom="on", length=5, labelsize=8)
            ax1.tick_params(axis='both',which='minor', direction="in", top="on", right="on", bottom="on", length=3, labelsize=8)

            fig.tight_layout()

            #fig.show()
            fig.savefig("rawGaussian.png", format="png",dpi=1000)



        #peaks = np.array([6,162,741])
        print(peaks)
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


if __name__ =='__main__':
    plot1 = iccd_evaluation("spectrum")         # calling an object from the class iccd_evaluation("String") with the parameters "heatmap", "spectrum", or "both". 
                                            # You can also set the mode "DA" for calculating a difference absorbance spectrum.                
    plot1.calibrate()                       # Use this mode to peak search and generate a calibration file with new calibration values. There must only be one file with the data from the calibration lamp in the current folder for this mode!
    #plot1.iterate()                         # start evaluating process by calling the function evaluate()