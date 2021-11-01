# -*- coding: utf-8 -*-

import os
import time
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class iccd_evaluation():
    def __init__(self, mode):
        self.mode = mode
        self.file = ""
        self.single_row = False
        self.input_row = 255
        self.mean_row_start = 300
        self.mean_row_end = 320

    def read_file(self):
        self.ascii_grid = pd.read_table(self.file, engine="python", skipfooter=29, index_col=1, usecols=[*range(1,1024)], header=None) 
        self.ascii_grid_transposed = self.ascii_grid.T
        self.image_filename = self.file[:len(self.file)-4]
        try:
            self.title = self.image_filename.rsplit("/",1)[1]
        except:
            self.title = self.image_filename
    
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
        plt.savefig(self.image_filename + "_heatmap", bbox_inches="tight")
        plt.clf()

    def plot_spectrum(self):
        if self.single_row == True:
            spectrum = self.ascii_grid_transposed.iloc[self.input_row]#[self.input_row:self.input_row+1]  # besser: to_numpy() 
        else: 
            spectrum = self.ascii_grid_transposed.iloc[self.mean_row_start:self.mean_row_end].mean()
        spectrum_list = spectrum.to_numpy()
        plt.plot(spectrum_list)
        plt.title(label=self.title, pad=-262, loc="left")
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
