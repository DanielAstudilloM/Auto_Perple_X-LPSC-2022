import os
import shutil
import subprocess
import pandas as pd
import numpy as np
import src.Autoperplex as AP


# Has Examples of how to run things


# names = AP.Auto_build("compositions.txt")
# AP.echo_file_mover(names)
# np_names = np.array(names)
# np.savetxt("last_source.txt",np_names,delimiter=",",fmt="%s")
# for n in names:
#     AP.build_command(n)

# AP.batcher(names)
# print(AP.bin_folder)

# AP.weramizer("A15C","C:/Users/cuent/Desktop/perplex2/tests/templates/",1800,10000)
# asd = AP.extractor("A15C")
# print(asd.head)
# AP.Auto_build("compositions.txt",[1000,2000],[1,30000],source_type="file")
# AP.frac("compositions.txt","A15C",[[1000,1100],[1100,1200]],[[10000,13000],[13000,16000]],[[1050,11000]])
# name = pd.read_csv("compositions.txt",delimiter="\t",index_col=0)
# name = name["A15C"]

asd = AP.extractor("A15C")
solid_fraction_vol = 100 - asd["melt(HGP)"].loc["vol %"]
print (solid_fraction_vol)
# asd["melt(HGP)"].iloc[4:13] = asd["melt(HGP)"].iloc[4:13]*(1-0.5) + name.iloc[0:9] * 0.5

# print(asd)