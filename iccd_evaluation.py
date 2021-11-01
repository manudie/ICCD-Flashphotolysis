# -*- coding: utf-8 -*-

import os
import time
import shutil
import pandas as pd
import matplotlib.pyplot as plt

class iccd_heatmap():
    def __init__(self):
        self.file = ""#"MD08_550nm.asc"

    def read_file(self):
        self.ascii_grid = pd.read_table(self.file, engine="python", skipfooter=29, index_col=1, usecols=[*range(1,1024)], header=None) 
        #print(self.ascii_grid)
        self.ascii_grid_transposed = self.ascii_grid.T
        self.image_filename = self.file[:len(self.file)-4]
    
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
        plt.title(label=self.image_filename.rsplit("/",1)[1], pad=-262, loc="left", color="white")#, fontsize=16,  , pad=0.3)
        #plt.show()
        plt.savefig(self.image_filename, bbox_inches="tight")
        plt.clf()

    def evaluate(self):
        self.read_file()
        self.plot_heatmap()     
        
    def iteration(self):
        directory = os.getcwd()
        if not os.path.isdir("processed_files"):
            os.mkdir("processed_files")
        if not os.path.isdir("iccd_heatmaps"):
            os.mkdir("iccd_heatmaps")
        time.sleep(1)
        for entry in os.scandir(directory):
            if entry.path.endswith(".asc") and entry.is_file():
                self.file = entry.path
                self.evaluate()
                shutil.move(self.file, "processed_files")
                shutil.move(self.image_filename + ".png", "iccd_heatmaps")



if __name__ =='__main__':
    plot1 = iccd_heatmap()
    #plot1.evaluate()
    plot1.iteration()
