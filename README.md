# ICCD-Blitzlichtphotolyse

This script was written to evaluate 1024x1024 px images (from the Andor iStar ICCD Camera).
It allows the user to generate a heatmap and/or a spectrum images from a specified area of lines in the exported sensor data.

It is essential, to export the data as an .asc (ASCII) file and append measurement information (tick the box in the export window (!)) in the Andor SOLIS Software. 
The exportet data than just needs to be in the same folder as this evaluation script. It is possible to place multiple files of raw .asc data into to folder, the script will iterate through every single one.

To decide wether to plot a heatmap, a spectrum or both just create an object from the class iccd_evaluation("String") and change "String" to "heatmap", "spectrum", or "both", e.g. plot1 = iccd_evaluation("both").

The spectrum is calculated as mean from rows in the image. By changing the values of self.mean_row_start and self.mean_row_end you can determine the span of rows, the  mean is calculated from.
By switching self.single_row to True, it is possible to generate the spectrum from a single row of pixels from the image, instead of mulitple rows. 

Calling the function plot1.evaluate() then starts the evaluation. Note that only in case of "both" the area wherefrom the mean for the spectrum is created is shown in red on the heatmap. To avoid the red area, just plot the heatmap and spectrum individual. 

The script will create folders for the plotted images, aswell as for processed data and will sort the files into these folders after it finished evaluating them. 
