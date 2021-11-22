# %%
# -*- coding: utf-8 -*-
"""
This script was written to evaluate 1024x1024 px images (from the Andor iStar ICCD Camera).
It allows the user to generate a heatmap and/or a spectrum images from a specified area of lines in the exported sensor data.
It is essential, to export the data as an .asc (ASCII) file and append acquisition information (tick the box in the export window (!)) in the Andor SOLIS Software. 
The exportet data than just needs to be in the same folder as this evaluation script. It is possible to place multiple files of raw .asc data into to folder, the script will iterate through every single one.
To decide wether to plot a heatmap, a spectrum or both just create an object from the class iccd_evaluation("String") and change "String" to "heatmap", "spectrum", or "both", e.g. plot1 = iccd_evaluation("both").
The spectrum is calculated as mean from rows in the image. By changing the values of self.mean_row_start and self.mean_row_end you can determine the span of rows, the  mean is calculated from.
By switching self.single_row to True, it is possible to generate the spectrum from a single row of pixels from the image, instead of mulitple rows. You can select the single row by changing the value of self.single_row.
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
from scipy.ndimage.filters import uniform_filter1d
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
import matplotlib.ticker as ticker
import matplotlib.pylab as pl

### define class ###
class iccd_evaluation():
    
    def __init__(self, mode):
        ##### Options to set by user #####
        self.test_run = False            # If self.test_run is True, Images will not be safed, but displayed and also files will not be moved. Made for easier developing.
        self.file_format_mode = "asc"           # "asc" for measurements from Andor iStar Camera or "txt" for measurements from OceanOptics Minispec which are already calibrated.
        self.drop_first_measurement = False             # In kinetic series w/ single track, often the first line of data is false due to build up charge in the ccd. Setting self.drop_first_measurent to True drops this line of data. Note to acquire n+1 mesaurements! 
        self.MA_filter = False
        self.layer_DA_spectra = True
        self.stack_DA_spectra = False
        self.legend_label_file = "legend_label_files/legend_labels_OD1.json"
        self.plot_title = False
        self.mean_row_start = 480          #250   
        self.mean_row_end = 880            #550
        self.single_row_mode = False
        self.single_row = 255
        self.calibration_file = "calibration_files/calibration_file.json"
        self.calibration_points = 2
        self.peak_height_min = 2500
        self.peak_distance = 20
        self.calibration_fit_peaks = True
        self.calibration_fit_mode = "gauss"           # "gauss", "lorentz" or "voigt"
        self.plot_calibration_fit = True         
        self.show_calibration_marks = False
        self.linewidth = 0.4 # Default: 0.6, MA: 1.0, txt: 0.6

        ##### Constructor that declares class variables (self.x) #####
        self.mode = mode
        self.file = ""
        self.readout_mode = ""          # Supported: "Full Resolution Image" ("FRI") or "Single Track" ("ST")
        self.wavelenght_cal_1 = 404.6565
        self.wavelenght_cal_2 = 435.8335
        self.wavelenght_cal_3 = 546.0750
        self.directory = os.getcwd()
        self.DA_I0D_list = ["_D_","_I0_","_I_"]
        self.DA_dict = {}
        self.DA_spectra_dict = {}
        self.calibration_mode = False
    
    def iterate(self):
        ### main function that scans and iterates through the current directory to evaluates each .asc file after another ###
        self.make_dirs()
        if self.mode == "spectrum" or self.mode == "heatmap" or self.mode == "both":
            for entry in os.scandir(self.directory):
                if entry.path.endswith("." + self.file_format_mode) and entry.is_file():
                    self.file = entry.path
                    self.read_file()
                    self.evaluate()
                    self.move_files()        
        elif self.mode == "DA" or self.mode == "A":
            DA_file_counter = 0
            filename_list = []
            for entry in os.scandir(self.directory):
                if entry.path.endswith("." + self.file_format_mode) and entry.is_file():
                    self.file = entry.path
                    for el in self.DA_I0D_list:
                        if el in self.file:
                            experiment = self.file.replace(el,"_")
                            for experiment_entry in os.scandir(self.directory):
                                if experiment_entry.path.endswith("." + self.file_format_mode) and experiment_entry.is_file():
                                    current_experiment = experiment_entry.path.replace(self.DA_I0D_list[DA_file_counter],"_")
                                    if current_experiment == experiment:
                                        DA_file_counter += 1
                                        filename_list.append(experiment_entry.path)
                                    if DA_file_counter == 3:
                                        break
                        if DA_file_counter == 3:
                            break
                if DA_file_counter == 3:
                    for el in filename_list:
                        self.file = el
                        self.read_file()
                        self.calculate_spectrum()
                        for el in self.DA_I0D_list:
                            if el in self.title:
                                self.DA_dict.update({el.replace("_",""):self.spectrum_list})
                        self.move_files() 
                    self.calculate_diff_absorbance()
                    if self.layer_DA_spectra == False and self.stack_DA_spectra == False:
                        self.plot_spectrum()
                    self.move_files() 
                    DA_file_counter = 0
                    filename_list = []
            if self.layer_DA_spectra == True or self.stack_DA_spectra == True:
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
            if self.file_format_mode == "asc":
                if "Single Track" in f.read():
                    self.readout_mode = "ST"
                    self.input_data_frame = pd.read_table(self.file, engine="python", skipfooter=29, index_col=0, header=None).dropna(axis=1, how="all")
                    self.input_data_frame = self.input_data_frame.drop([1,2,3,4],axis=0)
                    if self.drop_first_measurement == True:
                        self.input_data_frame = self.input_data_frame.drop([1],axis=1)
                else:
                    f.seek(0)
                    self.readout_mode = ""
                if self.readout_mode != "ST":
                    if "Full Resolution Image" in f.read():
                        self.readout_mode = "FRI"
                        self.input_data_frame = pd.read_table(self.file, engine="python", skipfooter=29, index_col=0, header=None).dropna(axis=1, how="all") #usecols=[*range(1,1024)],
                    else: 
                        print('Unsupported readout mode! Choose "Full Resolution Image" or "Single Track"')
            elif self.file_format_mode == "txt":
                self.input_data_frame = pd.read_csv(self.file, index_col=0, delimiter = "\t", decimal=",").dropna(axis=1, how="all")
                self.wavelengths = self.input_data_frame.index.values.tolist()
            f.close()
        #print(self.input_data_frame)
        self.input_data_frame_transposed = self.input_data_frame.T
        #print(self.input_data_frame_transposed)
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
        if self.file_format_mode == "asc":
            if self.readout_mode == "ST":
                spectrum = self.input_data_frame_transposed.mean()
            elif self.readout_mode == "FRI":
                if self.single_row == True:
                    spectrum = self.input_data_frame_transposed.iloc[self.single_row]
                else: 
                    spectrum = self.input_data_frame_transposed.iloc[self.mean_row_start:self.mean_row_end].mean()
        elif self.file_format_mode == "txt":
                spectrum = self.input_data_frame.mean(axis=1, numeric_only=True)
        self.spectrum_list = spectrum.to_numpy()
        if self.MA_filter == True:
            self.spectrum_list = uniform_filter1d(self.spectrum_list, size=5)
    
    def calculate_diff_absorbance(self):
        I = self.DA_dict.get("I")
        I0 = self.DA_dict.get("I0")
        D = self.DA_dict.get("D")
        self.diff_absorbance_list = -np.log10((I-D)/(I0-D))
        if self.MA_filter == True:
            self.diff_absorbance_list = uniform_filter1d(self.diff_absorbance_list, size=10)
        if self.layer_DA_spectra == True or self.stack_DA_spectra:
            self.DA_spectra_dict.update({self.title:self.diff_absorbance_list})
    
    def plot_spectrum(self):
        ### plots the generated dataframe to a spectrum with wavelengths on the x-axis from get_calibration() and saves it to the current folder ### 
        params = {
            'text.latex.preamble': ['\\usepackage{gensymb}'],
            'image.origin': 'lower',
            'image.interpolation': 'nearest',
            #'image.cmap': 'gray',
            'axes.grid': False,
            'savefig.dpi': 1000,  # to adjust notebook inline plot size
            'axes.labelsize': 11, 
            'axes.titlesize': 11,
            'font.size': 11, 
            'legend.fontsize': 6, 
            'xtick.labelsize': 11,
            'ytick.labelsize': 11,
            'text.usetex': True,
            #'figure.figsize': [3, 3],
            'font.family': 'serif',
        }
        plt.rcParams.update(params)        
        if self.file_format_mode == "asc":
            self.get_calibration() 
        y_label_drawn = False
        if self.layer_DA_spectra == False and self.stack_DA_spectra == False:
            fig = plt.figure(figsize=(4,3))
            gs = gridspec.GridSpec(1,1)
            ax1 = fig.add_subplot(gs[0]) 
            if self.mode == "DA" or self.mode == "A":
                spectrum = self.diff_absorbance_list
            else:
                spectrum = self.spectrum_list
            if self.file_format_mode == "asc":
                ax1.plot(self.wavelengths, spectrum, linewidth=self.linewidth, color="darkcyan")
            elif self.file_format_mode == "txt":
                trim_wl = []
                trim_spec = []
                for i, el in enumerate(self.wavelengths):
                    if el > 385 and el < 702:
                        trim_wl.append(el)
                        trim_spec.append(spectrum[i])
                ax1.plot(trim_wl, trim_spec, linewidth=self.linewidth, color="darkcyan")
            ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
            ax1.xaxis.set_major_locator(ticker.MultipleLocator(50))
            #ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))
        elif self.layer_DA_spectra == True or self.stack_DA_spectra == True:
            f = open(self.legend_label_file)
            label_dict = json.load(f)
            f.close()
            label_order_list = [*label_dict.keys()]
            self.DA_spectra_dict = dict(sorted(self.DA_spectra_dict.items(), key=lambda pair: label_order_list.index(pair[0])))
            colors = pl.cm.turbo(np.linspace(0,1,len(self.DA_spectra_dict)))
            two_colors = ["darkcyan","firebrick"]
            try:
                f = open(self.legend_label_file, encoding='utf8')
                label_dict = json.load(f)
                f.close()  
            except:
                print("No label file found!")             
            i = 0
            if self.layer_DA_spectra == True:
                fig = plt.figure(figsize=(4,3))
                gs = gridspec.GridSpec(1,1)
                ax1 = fig.add_subplot(gs[0]) 
                ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
                ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
                ax1.xaxis.set_major_locator(ticker.MultipleLocator(50))
                for key in self.DA_spectra_dict:
                    try:
                        spectrum = self.DA_spectra_dict[key]
                        if self.file_format_mode == "txt":
                            trim_wl = []
                            trim_spec = []
                            for j, el in enumerate(self.wavelengths):
                                if el > 385 and el < 702:
                                    trim_wl.append(el)
                                    trim_spec.append(spectrum[j])
                            ax1.plot(trim_wl, trim_spec, linewidth=self.linewidth, color=two_colors[i], label=label_dict[key])
                        elif self.file_format_mode == "asc":
                            if len(self.DA_spectra_dict) > 1:        
                                ax1.plot(self.wavelengths, spectrum, color=colors[i], linewidth=self.linewidth,  label=label_dict[key])
                            else:
                                ax1.plot(self.wavelengths, spectrum, color="darkcyan", linewidth=self.linewidth,  label=label_dict[key])
                    except:
                        ax1.plot(self.wavelengths, self.DA_spectra_dict[key], color=colors[i], linewidth=self.linewidth, label=key)
                    ax1.legend(loc="lower right", prop={'family':"serif", 'size':5.8})
                    #ax1.set_ylim([-0.035, 0.025])
                    i += 1       
            elif self.stack_DA_spectra == True:
                fig, axs = plt.subplots(len(self.DA_spectra_dict), sharex=True, sharey=True, figsize=(3,7))
                plt.clf()
                gs = gridspec.GridSpec(len(self.DA_spectra_dict), 1, hspace=0)
                for key in self.DA_spectra_dict:
                    axs[i] = fig.add_subplot(gs[i]) 
                    try:
                        axs[i].plot(self.wavelengths, self.DA_spectra_dict[key], color="darkcyan", linewidth=self.linewidth,  label=label_dict[key])
                    except:
                        axs[i].plot(self.wavelengths, self.DA_spectra_dict[key], color="darkcyan", linewidth=self.linewidth,  label=key)
                    axs[i].legend(loc="lower right", prop={'family':"serif", 'size':5.8})
                    axs[i].xaxis.set_minor_locator(AutoMinorLocator(2))
                    axs[i].yaxis.set_minor_locator(AutoMinorLocator(2))
                    axs[i].xaxis.set_major_locator(ticker.MultipleLocator(50))
                    axs[i].tick_params(axis='x',which='major', direction="in", top="on", right="on", bottom="on", length=5, labelsize=8)
                    axs[i].tick_params(axis='x',which='minor', direction="in", top="on", right="on", bottom="on", length=3, labelsize=8)
                    axs[i].tick_params(axis='y',which='major', direction="in", top="on", right="on", bottom="on", length=5, labelsize=6)
                    axs[i].tick_params(axis='y',which='minor', direction="in", top="on", right="on", bottom="on", length=3, labelsize=6)
                    if i == len(self.DA_spectra_dict)//2:
                        plt.ylabel("ΔA", family="serif", fontsize=8)
                        y_label_drawn = True
                    '''
                    try:
                        axs[i] = fig.add_subplot(gs[i]) 
                        axs[i].plot(self.wavelengths, self.DA_spectra_dict[key], color=colors[i], linewidth=self.linewidth,  label=label_dict[key])
                    except:
                        ax1.plot(self.wavelengths, self.DA_spectra_dict[key], color=colors[i], linewidth=self.linewidth, label=key)
                    '''
                    i += 1
        if self.plot_title == True:
            plt.title(label=self.title, pad=-262, loc="left", family="serif")#, fontsize=8)
        plt.xlabel("$\lambda$ / nm", family="serif")#, fontsize=8)
        if y_label_drawn == False:
            #plt.xticks(family="serif", fontsize=8)
            #plt.yticks(family="serif", fontsize=8)
            plt.tick_params(axis='both',which='major', direction="in", top="on", right="on", bottom="on", length=5)#, labelsize=8)
            plt.tick_params(axis='both',which='minor', direction="in", top="on", right="on", bottom="on", length=3)#, labelsize=8)
            if self.mode == "DA":
                plt.ylabel("$\Delta $A")#, family="serif", fontsize=8)
            elif self.mode == "A":
                plt.ylabel("$A$")#, family="serif", fontsize=8)
            else:
                plt.ylabel("counts", family="serif")#, fontsize=8)
        if self.calibration_mode == True:
            if self.calibration_fit_peaks == True and self.plot_calibration_fit == True and self.calibration_fit_mode == "gauss":
                ax1.plot(self.wavelengths, self.gauss_peak_1, "orange")
                ax1.fill_between(self.wavelengths, self.gauss_peak_1.min(), self.gauss_peak_1, facecolor="orange", alpha=0.2)
                ax1.plot(self.wavelengths, self.gauss_peak_2, "orange")
                ax1.fill_between(self.wavelengths, self.gauss_peak_2.min(), self.gauss_peak_2, facecolor="orange", alpha=0.2)  
            elif self.calibration_fit_peaks == True and self.plot_calibration_fit == True and self.calibration_fit_mode == "lorentz":
                '''
                ax1.plot(self.wavelengths, self.lorentz_peak_1, "indianred")
                ax1.fill_between(self.wavelengths, self.lorentz_peak_1.min(), self.lorentz_peak_1, facecolor="indianred", alpha=0.2)
                '''
                ax1.plot(self.wavelengths, self.lorentz_peak_2, "palegreen")
                ax1.fill_between(self.wavelengths, self.lorentz_peak_2.min(), self.lorentz_peak_2, facecolor="palegreen", alpha=0.2)  
                ax1.plot(self.wavelengths, self.lorentz_peak_3, "palegreen")           #c
                ax1.fill_between(self.wavelengths, self.lorentz_peak_3.min(), self.lorentz_peak_3, facecolor="palegreen", alpha=0.2)           #cyan
            elif self.calibration_fit_peaks == True and self.plot_calibration_fit == True and self.calibration_fit_mode == "voigt":    
                ax1.plot(self.wavelengths, self.voigt_peak_2, "blueviolet")
                ax1.fill_between(self.wavelengths, self.voigt_peak_2.min(), self.voigt_peak_2, facecolor="blueviolet", alpha=0.2)  
                ax1.plot(self.wavelengths, self.voigt_peak_3, "blueviolet")           #c
                ax1.fill_between(self.wavelengths, self.voigt_peak_3.min(), self.voigt_peak_3, facecolor="blueviolet", alpha=0.2)           #cyan
            ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
            ax1.xaxis.set_major_locator(ticker.MultipleLocator(50))
        if self.show_calibration_marks == True:
            if self.calibration_points == 2:
                plt.vlines([self.wavelenght_cal_2, self.wavelenght_cal_3], ymin=self.spectrum_list.min(), ymax=self.spectrum_list.max(), color="grey", linewidth=0.5)
            elif self.calibration_points == 3:
                plt.vlines([self.wavelenght_cal_1, self.wavelenght_cal_2, self.wavelenght_cal_3], ymin=self.spectrum_list.min(), ymax=self.spectrum_list.max(), color="grey", linewidth=0.5)
        fig.tight_layout()
        if self.test_run == True: 
            plt.show()
        else:
            fig.savefig(self.image_filename + "_spectrum", bbox_inches="tight", dpi=1000)
        plt.clf()

    def plot_heatmap(self):
        ### plots the generated dataframe to a heatmap and saves it to the current folder ###
        params = {
            'text.latex.preamble': ['\\usepackage{gensymb}'],
            'image.origin': 'lower',
            'image.interpolation': 'nearest',
            #'image.cmap': 'gray',
            'axes.grid': False,
            'savefig.dpi': 1000,  # to adjust notebook inline plot size
            'axes.labelsize': 9, 
            'axes.titlesize': 9,
            'font.size': 9, 
            'legend.fontsize': 6, 
            'xtick.labelsize': 9,
            'ytick.labelsize': 9,
            'text.usetex': True,
            #'figure.figsize': [3, 3],
            'font.family': 'serif',
        }
        plt.rcParams.update(params)
        plt.pcolor(self.input_data_frame_transposed)
        cb = plt.colorbar(pad=0.01)
        cb.set_label(label="counts")
        #cb.ax.tick_params(labelsize=10)
        plt.tick_params(
            #labelsize=10,
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False, 
            left=False,        # ticks along the top edge are off
            labelbottom=False, # labels along the bottom edge are off
            labelleft=False)
        if self.plot_title == True:
            plt.title(label=self.title, pad=-262, loc="left", color="white", family="serif", fontsize=8)
        if self.mode == "both" and self.readout_mode == "FRI":
            plt.axhspan(self.mean_row_start, self.mean_row_end, color="red", alpha=0.15, lw=0)
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
            f = open(self.calibration_file)
            calibration_dict_in = json.load(f)
            m = calibration_dict_in.get("m")
            b = calibration_dict_in.get("b")
            f.close()
        except:
            print("No calibration file found! Spectrum will be generated with false values!")
            m = 0.3033110988290711   
            b = 386.25945492645576   
        self.index_list = np.arange(1,len(self.input_data_frame_transposed.columns)+1)
        self.wavelengths = list(map(lambda x: m*x+b, self.index_list))
        #print(len(self.wavelengths))
        #print(self.wavelengths.index(532.4293391415051)+4)
    
    def _1gaussian(self, x, amp1,cen1,sigma1):
        return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))

    def _2gaussian(self, x, amp1,cen1,sigma1, amp2,cen2,sigma2):
        return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2))) + \
                amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen2)/sigma2)**2)))

    def _1Lorentzian(self,x, amp, cen, wid):
        return amp*wid**2/((x-cen)**2+wid**2)

    def _3Lorentzian(self, x, amp1, cen1, wid1, amp2,cen2,wid2, amp3,cen3,wid3):
        return (amp1*wid1**2/((x-cen1)**2+wid1**2)) +\
                (amp2*wid2**2/((x-cen2)**2+wid2**2)) +\
                    (amp3*wid3**2/((x-cen3)**2+wid3**2))
    
    def _1Voigt(self, x, ampG1, cenG1, sigmaG1, ampL1, cenL1, widL1):
        return (ampG1*(1/(sigmaG1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cenG1)**2)/((2*sigmaG1)**2)))) +\
              ((ampL1*widL1**2/((x-cenL1)**2+widL1**2)) )

    def _3Voigt(self, x, ampG1, cenG1, sigmaG1, ampL1, cenL1, widL1, ampG2, cenG2, sigmaG2, ampL2, cenL2, widL2, ampG3, cenG3, sigmaG3, ampL3, cenL3, widL3):
        return (ampG1*(1/(sigmaG1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cenG1)**2)/((2*sigmaG1)**2)))) + ((ampL1*widL1**2/((x-cenL1)**2+widL1**2)) ) +\
                (ampG2*(1/(sigmaG2*(np.sqrt(2*np.pi))))*(np.exp(-((x-cenG2)**2)/((2*sigmaG2)**2)))) + ((ampL2*widL2**2/((x-cenL2)**2+widL2**2)) ) +\
                    (ampG3*(1/(sigmaG3*(np.sqrt(2*np.pi))))*(np.exp(-((x-cenG3)**2)/((2*sigmaG3)**2)))) + ((ampL3*widL3**2/((x-cenL3)**2+widL3**2)) ) 

    def calibrate(self):
        ### In calibrate() mode this function will search for the column numbers of three peaks in the spectrum of the calibration lamp. 
        # The peaks will be determined by their height, set by self.peak_height_min. 
        # It then uses a three-point linear fit to calculate values for the slope m and intercept b, which are stored together with the wavelenghts and depending column values in a .json file for later use. ###
        self.calibration_mode = True
        self.get_calibration_file()
        self.read_file()
        self.calculate_spectrum()
        peaks = find_peaks(self.spectrum_list, height=self.peak_height_min, distance=self.peak_distance)[0]
        #print(peaks)
        if peaks.size > 3:
            peaks = find_peaks(self.spectrum_list, height=self.peak_height_min, distance=self.peak_distance, width=40)[0]
            #print(peaks)
            
            if self.calibration_fit_peaks == True:
                self.index_list = np.arange(1,len(self.input_data_frame_transposed.columns)+1) 
                if self.calibration_fit_mode == "gauss":
                    ########## gauss ##########
                    gauss_amp1 = 100
                    gauss_sigma1 = 10
                    gauss_cen1 = 190
                    gauss_amp2 = 7000
                    gauss_sigma2 = 5
                    gauss_cen2 = 520
                    popt_2gauss, pcov_2gauss = optimize.curve_fit(self._2gaussian, self.index_list, self.spectrum_list, p0=[gauss_amp1, gauss_cen1, gauss_sigma1, gauss_amp2, gauss_cen2, gauss_sigma2])
                    perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
                    gauss_pars_1 = popt_2gauss[0:3]
                    gauss_pars_2 = popt_2gauss[3:6]
                    self.gauss_peak_1 = self._1gaussian(self.index_list, *gauss_pars_1) + self.spectrum_list.min()
                    self.gauss_peak_2 = self._1gaussian(self.index_list, *gauss_pars_2) + self.spectrum_list.min()
                    fit_peak_list = [gauss_pars_1[1], gauss_pars_2[1]]
                elif self.calibration_fit_mode == "lorentz":
                    ########## lorentz ##########
                    lorentz_amp1 = 4000
                    lorentz_cen1 = 80
                    lorentz_wid1 = 20
                    lorentz_amp2 = 13000
                    lorentz_cen2 = 150
                    lorentz_wid2 = 10
                    lorentz_amp3 = 50
                    lorentz_cen3 = 200
                    lorentz_wid3 = 10
                    popt_3lorentz, pcov_3lorentz = optimize.curve_fit(self._3Lorentzian, self.index_list, self.spectrum_list, p0=[lorentz_amp1, lorentz_cen1, lorentz_wid1, lorentz_amp2, lorentz_cen2, lorentz_wid2, lorentz_amp3, lorentz_cen3, lorentz_wid3])
                    perr_3lorentz = np.sqrt(np.diag(pcov_3lorentz))
                    lorentz_pars_1 = popt_3lorentz[0:3]
                    lorentz_pars_2 = popt_3lorentz[3:6]
                    lorentz_pars_3 = popt_3lorentz[6:9]
                    self.lorentz_peak_1 = self._1Lorentzian(self.index_list, *lorentz_pars_1) + self.spectrum_list.min()
                    self.lorentz_peak_2 = self._1Lorentzian(self.index_list, *lorentz_pars_2) + self.spectrum_list.min()
                    self.lorentz_peak_3 = self._1Lorentzian(self.index_list, *lorentz_pars_3) + self.spectrum_list.min()     
                    fit_peak_list = [lorentz_pars_2[1], lorentz_pars_3[1]]
                elif self.calibration_fit_mode == "voigt":
                    ########## voigt ##########
                    voigt_ampG1 = 5000
                    voigt_cenG1 = 60
                    voigt_sigmaG1 = 5
                    voigt_ampL1 = 5000
                    voigt_cenL1 = 5
                    voigt_widL1 = 50
                    voigt_ampG2 = 100
                    voigt_cenG2 = 190
                    voigt_sigmaG2 = 10
                    voigt_ampL2 = 13000
                    voigt_cenL2 = 160
                    voigt_widL2 = 20
                    voigt_ampG3 = 7000
                    voigt_cenG3 = 520
                    voigt_sigmaG3 = 5
                    voigt_ampL3 = 7000
                    voigt_cenL3 = 520
                    voigt_widL3 = 20
                    popt_3voigt, pcov_1voigt = optimize.curve_fit(self._3Voigt, self.index_list, self.spectrum_list, p0=[voigt_ampG1, voigt_cenG1, voigt_sigmaG1, voigt_ampL1, voigt_cenL1, voigt_widL1, voigt_ampG2, voigt_cenG2, voigt_sigmaG2, voigt_ampL2, voigt_cenL2, voigt_widL2, voigt_ampG3, voigt_cenG3, voigt_sigmaG3, voigt_ampL3, voigt_cenL3, voigt_widL3])
                    perr_3voigt = np.sqrt(np.diag(pcov_1voigt))                    
                    voigt_pars_1 = popt_3voigt[0:6]
                    voigt_pars_2 = popt_3voigt[6:12]
                    voigt_pars_3 = popt_3voigt[12:18]                    
                    self.voigt_peak_1 = self._1Voigt(self.index_list, *voigt_pars_1) + self.spectrum_list.min()
                    self.voigt_peak_2 = self._1Voigt(self.index_list, *voigt_pars_2) + self.spectrum_list.min()
                    self.voigt_peak_3 = self._1Voigt(self.index_list, *voigt_pars_3) + self.spectrum_list.min()  
                    fit_peak_list = [voigt_pars_2[1], voigt_pars_3[1]]
        #print("Peaks found over " + str(self.peak_height_min) + " counts at columns: " + str(peaks))
        if self.calibration_fit_peaks == True:
            dict_peaks = fit_peak_list
        else:
            dict_peaks = peaks                
        wavelenght_cal_list_2 = [self.wavelenght_cal_2, self.wavelenght_cal_3]
        wavelenght_cal_list_3 = [self.wavelenght_cal_1, self.wavelenght_cal_2, self.wavelenght_cal_3]
        #m,b = np.polyfit(peaks,wavelenght_cal_list_3,1)
        m,b = np.polyfit(fit_peak_list,wavelenght_cal_list_2,1)
        if self.calibration_points == 2:
            calibration_dict = {
            "wavelenghts (y):": "columns (x):",
            self.wavelenght_cal_2: int(dict_peaks[0]),
            self.wavelenght_cal_3: int(dict_peaks[1]),
            "slope and intercept:": "values:",
            "m" : float(m),
            "b" : float(b)
        } 
        elif self.calibration_points == 3:
            calibration_dict = {
                "wavelenghts (y):": "columns (x):",
                self.wavelenght_cal_1: int(dict_peaks[0]),
                self.wavelenght_cal_2: int(dict_peaks[1]),
                self.wavelenght_cal_3: int(dict_peaks[2]),
                "slope and intercept:": "values:",
                "m" : float(m),
                "b" : float(b)
            }         
        with open(self.calibration_file, "w", encoding="utf-8") as f:
            json.dump(calibration_dict, f, ensure_ascii=False, indent=4)        

        self.plot_spectrum()

if __name__ =='__main__':
    plot1 = iccd_evaluation("DA")         # calling an object from the class iccd_evaluation("String") with the parameters "heatmap", "spectrum", or "both". 
                                            # You can also set the mode "DA" for calculating a difference absorbance spectrum.                
    #plot1.calibrate()                       # Use this mode to peak search and generate a calibration file with new calibration values. There must only be one file with the data from the calibration lamp in the current folder for this mode!
    plot1.iterate()                         # start evaluating process by calling the function evaluate()

# %%

