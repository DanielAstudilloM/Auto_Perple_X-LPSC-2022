Instructions

Make sure you have all requirements installed. All Perple_X files need to be contained in the bin folder.
To get the Perple_X software go to https://www.perplex.ethz.ch/

Data and samples should go in the compositions.txt file in that format. You can put as many samples as you want and it will automatically detect oxide compositions for SiO2, TiO2, Al2O3, FeO, MgO, MnO, Cr2O3, CaO, Na2O, K2O. Any additional data in the file will be ignored. Avoid using special characters for sample names and spaces will not be used for the results names.

To run the program run the run_script.py from the main folder, modify the Auto_build command in the script to use your files and pressure/temperature conditions. After calculations are finished, use the save_script.py to save your data to the results folder. 

Note that if you repeat the sample names after you already ran them once, your previous data will be deleted. Move your results before doing that.