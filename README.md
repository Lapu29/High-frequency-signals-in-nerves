This repository contains the Python code necessary for generating data and graphs for the paper "High-frequency signals: a comparison between the cable equation and telegrapher's 
equations in nerves". The .pyx file needs to be Cythonized using setup1.py. Running the Cython file will generate three .csv files, which are used by the script Figure3.py. All other
scripts are self-contained and will produce the figures in the paper, corresponding to the file name.
