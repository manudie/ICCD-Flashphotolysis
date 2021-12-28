# ICCD-Flashphotolysis

This script was written to evaluate 1024x1024 px images (from the Andor iStar ICCD Camera).
It allows the user to generate a heatmap and/or a spectrum images from a specified area of lines in the exported sensor data.

It is essential, to export the data as an .asc (ASCII) file and append acquisition information (tick the box in the export window (!)) in the Andor SOLIS Software. 
The exportet data than just needs to be in the same folder as this evaluation script. It is possible to place multiple files of raw .asc data into to folder, the script will iterate through every single one.

To decide wether to plot a heatmap, a spectrum or both just create an object from the class iccd_evaluation("String") and change "String" to "heatmap", "spectrum", or "both", e.g. plot1 = iccd_evaluation("both").

![MD08_320nm_heatmap](https://user-images.githubusercontent.com/42771712/147589268-1f441812-98a7-4361-bc58-2b34a2f2a636.png)

The spectrum is calculated as mean from rows in the image. By changing the values of self.mean_row_start and self.mean_row_end you can determine the span of rows, the  mean is calculated from.

![MD08_320nm_spectrum](https://user-images.githubusercontent.com/42771712/147589277-a4788804-9ad4-4737-ac19-0be2c5c9e924.png)

By switching self.single_row to True, it is possible to generate the spectrum from a single row of pixels from the image, instead of mulitple rows. You can select the single row by changing the value of self.input_row.

Calling the function plot1.iterate() then starts the evaluation. Note that only in case of "both" the area wherefrom the mean for the spectrum is created is shown in red on the heatmap. To avoid the red area, just plot the heatmap and spectrum individual. 

The script will create folders for the plotted images, aswell as for processed data and will sort the files into these folders after it finished evaluating them. 

To calculate the difference absorbance create an object of the class with the string "DA" or "A". "A" labels the y-axis in case you want to display the absorbance. It is necessary that you feed the code with three measurement files containing "_I_", "_I0_" and "_D_" in their names, in both of these modes. 
![MD13_9_I_1_spectrum_MA](https://user-images.githubusercontent.com/42771712/147589292-7090d008-354b-4261-bf3c-f5dd4283d112.png)
In case you want to calcultate the difference absorbance "_I_" should be the laser excited 
measurement, "_I0_" the measurement without excitement and "_D_" the dark measurement. 
You can layer or stack these DA-Spectra by setting either self.layer_DA_spectra or self.stack_DA_spectra to true. If there also should be a legend plotted, you can add a .json file where you define a legend label for each measurement. You can then refere to this file in self.legend_label_file.
