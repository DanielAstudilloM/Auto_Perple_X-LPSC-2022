import numpy as np
import pandas as pd
import multiprocessing as mp
import os
import shutil as sh
import subprocess


def Auto_build(source):
    """Automates the creation of building files to be processed through vertex (or other command)
    by using a template file and the 'build.exe' command
    
    Input: source must be a tab delimited file with all samples to be used and their composition
    un columns. First column must contain all element names and first row all sample names
    
    no return"""
    # Reading input file
    source_file = pd.read_csv(source,delimiter="\t",index_col=0)

    # Current template asks for SiO2, TiO2, Al2O3, Cr2O3, FeO, MnO, MgO, CaO, Na2O, K2O
    # Testing availability or mistakes.
    ind = source_file.index
    cols = source_file.columns
    if ("SiO2" in ind and "TiO2" in ind and "Al2O3" in ind and "Cr2O3" in ind and "FeO" in ind
    and "MgO" in ind and "CaO" in ind and "Na2O" in ind and "K2O" in ind):
        print("All components present, continuing")
    else:
        print("One or more components are missing or mispelled, check source file")
        return
    
    sio2 = source_file.loc["SiO2"]
    elements = source_file.loc["SiO2":"K2O"]
   
    # Importing template directories: fh= first half, sh = second half, elements go between
    template_1_fh = open("templates/template1_first_half_noMn.txt","r")
    fh = template_1_fh.read()
    template_1_sh = open("templates/template2_second_half.txt","r")
    sec_h = template_1_sh.read()
    template_1_sh.close()
    template_1_fh.close()

    # Starting loop for writing input files
    for i, val in enumerate(cols):
        val2 = val.replace(" ","")
        if os.path.exists(str(val2)+".txt"):
            os.remove(str(val2)+".txt")
        
        c_file = open(str(val2)+".txt","w")
        c_file.write(str(val2)+"\n")
        c_file.write(fh+"\n")
        for t in elements[val]:
            c_file.write(str(t)+"\n")
        c_file.write(sec_h)
        c_file.close()
    names = cols+".txt"
    names = names.tolist()
    names2 = []
    for n in names:
        n2 = n.replace(" ","")
        names2.append(n2)
    # print("starting pool")
    return(names2)

def build_command(file):
    """Runs a command on cmd that executes the perple_x build
    with automatic input from a constructed file
    
    input: file: must be a string with the name of the file"""
    
    if os.path.exists(file.replace(".txt",".dat")):
        os.remove(file.replace(".txt",".dat"))
    ring = "build <"+str(file)
    subprocess.run(ring,shell=True,stdout=subprocess.DEVNULL)

def batcher(source):

    c = os.cpu_count()-4
    for v in range(c+1):
        if os.path.exists("batch"+ str(v+1) + ".bat"):
            os.remove("batch"+ str(v+1) + ".bat")
        filex = open("batch"+ str(v+1) + ".bat","x")
        filex.close()
    counter = 1
    for h in source:
        counter = counter + 1
        if counter <= c+1 and counter > 1:
            filex = open("batch"+ str(counter) + ".bat","a")
            filex.write("echo "+ h.replace(".txt","") + " |vertex\n")
            filex.close()

        else:
            counter = 1
            filex = open("batch"+ str(counter) + ".bat","a")
            filex.write("echo "+ h.replace(".txt","") + " |vertex\n")
            filex.close()

names = Auto_build("compositions_noMn.txt")
np_names = np.array(names)
np.savetxt("last_source.txt",np_names,delimiter=",",fmt="%s")

for n in names:
    build_command(n)

batcher(names)