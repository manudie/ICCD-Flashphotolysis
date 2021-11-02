# ICCD-Blitzlichtphotolyse

This script was written to evaluate 1024x1024 px images (from the Andor iStar ICCD Camera).
It allows the user to generate a heatmap and/or a spectrum images from a specified area of lines in the exported sensor data.

It is essential, to export the data as an .asc (ASCII) file and append acquisition information (tick the box in the export window (!)) in the Andor SOLIS Software. 
The exportet data than just needs to be in the same folder as this evaluation script. It is possible to place multiple files of raw .asc data into to folder, the script will iterate through every single one.

To decide wether to plot a heatmap, a spectrum or both just create an object from the class iccd_evaluation("String") and change "String" to "heatmap", "spectrum", or "both", e.g. plot1 = iccd_evaluation("both").

The spectrum is calculated as mean from rows in the image. By changing the values of self.mean_row_start and self.mean_row_end you can determine the span of rows, the  mean is calculated from.
By switching self.single_row to True, it is possible to generate the spectrum from a single row of pixels from the image, instead of mulitple rows. You can select the single row by changing the value of self.input_row.

Calling the function plot1.evaluate() then starts the evaluation. Note that only in case of "both" the area wherefrom the mean for the spectrum is created is shown in red on the heatmap. To avoid the red area, just plot the heatmap and spectrum individual. 

The script will create folders for the plotted images, aswell as for processed data and will sort the files into these folders after it finished evaluating them. 


![MD05_Kamera2Test_550nm_3ms_heatmap](https://user-images.githubusercontent.com/42771712/139681285-1a5b1f46-b010-4216-82d8-9e3995aa9dd0.png)

![MD05_Kamera2Test_550nm_3ms_spectrum](https://user-images.githubusercontent.com/42771712/139681283-c81d5039-95d3-4238-a660-e134d28a68f5.png)
