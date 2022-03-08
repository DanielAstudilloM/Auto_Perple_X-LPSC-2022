import numpy as np
import pandas as pd
import multiprocessing as mp
import os
import platform
import pathlib
import shutil as sh
import subprocess

dir_path = pathlib.Path(__file__).parent.parent.as_posix() 
bin_folder = dir_path + "/bin/"
src_folder = dir_path + "/src/"
results_folder = dir_path + "/results/"
op_system = platform.system()

def Auto_build(source,temp,pressure,source_type="file",mp="optimal"):

    """Automates the creation of building files to be processed through vertex (or other command)
    by using a template file and the 'build.exe' command
    
    Input: 
        source: tab delimited file with all samples to be used and their composition in columns.
                First column must contain all element names and first row all sample names.
        source_type: "file" as default, referring to the source as a file, or "frac" to use a pandas
                dataframe obtained through the extractor function from FOIS.
    
    Return:
        names2: List of file names produced, corresponding to the sample names (without spaces) as txt files"""

    # Reading input file
    if source_type == "file":
        source_file = pd.read_csv(source,delimiter="\t",index_col=0)
        ind = source_file.index
        cols = source_file.columns

    elif source_type == "frac":
        assert ("melt(HGP)" in source.columns)
        source_file = source["melt(HGP)"]
        cols = ["melt(HGP)"]
        ind = source_file.index
    # Current template asks for SiO2, TiO2, Al2O3, Cr2O3, FeO, MnO, MgO, CaO, Na2O, K2O
    # Testing availability or mistakes.
    
    if ("SiO2" in ind and "TiO2" in ind and "Al2O3" in ind and "Cr2O3" in ind and "FeO" in ind
     and "MnO" in ind and "MgO" in ind and "CaO" in ind and "Na2O" in ind and "K2O" in ind):
        print("All components present, continuing")
    else:
        print("One or more components are missing or mispelled, check source file")
        return
    
    sio2 = source_file.loc["SiO2"]
    elements = source_file.loc["SiO2":"K2O"]
   
    # Importing template directories: fh= first half, sh = second half, elements go between
    template_1_fh = open(dir_path+"/templates/template1_first_half.txt","r")
    fh = template_1_fh.read().replace("temp1",str(temp[0])).replace("temp2",str(temp[1])).replace("pressure1",str(pressure[0])).replace("pressure2",str(pressure[1]))
    template_1_sh = open(dir_path+"/templates/template2_second_half.txt","r")
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

        if source_type == "file":
            for t in elements[val]:
                c_file.write(str(t)+"\n")

        elif source_type == "frac":
            for t in elements:
                c_file.write(str(t)+"\n")
        
        c_file.write(sec_h)
        c_file.close()
    
    if source_type == "file":
        names = cols+".txt"
        names = names.tolist()

    elif source_type == "frac":
        names = cols[0]+".txt"
        names = [names]

    names2 = []
    for n in names:
        n2 = n.replace(" ","")
        names2.append(n2)

    echo_file_mover(names2, folder="bin")

    for l in names2:
        build_command(l)

    batcher (names2,multiprocessing=mp)

    return(names2)

def build_command(file):
    """Runs a command on cmd that executes the perple_x build
    with automatic input from a constructed file
    
    input: file: must be a string with the name of the file"""
    

    if os.path.exists(bin_folder + file.replace(".txt",".dat")):
        os.remove(bin_folder + file.replace(".txt",".dat"))

    ring = "cd " + bin_folder + "&& build <"+str(file)
    subprocess.run(ring,shell=True,stdout=subprocess.DEVNULL)

def batcher(source, multiprocessing = "optimal"):
    cores =os.cpu_count()
    if multiprocessing == "optimal":
        if cores == 1 or cores == 2:
            c = 0
        elif cores == 3 or cores == 4:
            c = 1
        elif cores == 5 or cores == 6:
            c = 2
        elif cores == 7:
            c = 3
        elif cores >= 8:
            c = cores - 4
    if multiprocessing == "intense":
        c = cores-1
    if multiprocessing == "single":
        c = 0
    
    if type(multiprocessing) == int:
        c = multiprocessing

    if "Windows" in op_system:
        for v in range(c+1):
            if os.path.exists(bin_folder+"batch"+ str(v+1) + ".bat"):
                os.remove(bin_folder+"batch"+ str(v+1) + ".bat")
            filex = open(bin_folder+"batch"+ str(v+1) + ".bat","x")
            filex.close()
        counter = 1
        for h in source:
            counter = counter + 1
            if counter <= c+1 and counter > 1:
                filex = open(bin_folder+"batch"+ str(counter) + ".bat","a")
                filex.write("echo "+ h.replace(".txt","") + " |vertex\n")
                filex.close()

            else:
                counter = 1
                filex = open(bin_folder+"batch"+ str(counter) + ".bat","a")
                filex.write("echo "+ h.replace(".txt","") + " |vertex\n")
                filex.close()
    else:
        for v in range(c+1):
            if os.path.exists(bin_folder+"batch"+ str(v+1) + ".sh"):
                os.remove(bin_folder+"batch"+ str(v+1) + ".sh")
            filex = open(bin_folder+"batch"+ str(v+1) + ".sh","x")
            filex.write("#!/bin/sh")
            filex.close()
        counter = 1
        for h in source:
            counter = counter + 1
            if counter <= c+1 and counter > 1:
                filex = open(bin_folder+"batch"+ str(counter) + ".sh","a")
                filex.write("echo "+ h.replace(".txt","") + " |vertex\n")
                filex.close()

            else:
                counter = 1
                filex = open(bin_folder+"batch"+ str(counter) + ".sh","a")
                filex.write("echo "+ h.replace(".txt","") + " |vertex\n")
                filex.close()
        
def echo_file_mover(source,folder="results"):
    """Moves the files generated by the Auto_build function to the results or bin folders
    
    Input:
        source: Return of the Auto_build function.
        folder: Either "results" or "bin" for the corresponding folder
    Return
        No return"""
    
    if folder=="results":
        copying_folder = results_folder
    elif folder == "bin":
        copying_folder = bin_folder
    else:
        print("location not available")

    for echo_file in source:
        # Removing file on copy location
        if os.path.exists(copying_folder+echo_file):
            os.remove(copying_folder+echo_file)

        # copying file
        sh.copy(echo_file,copying_folder)

        # removing file on main location
        os.remove(echo_file)

def build_file_mover(source,folder):
    """Moves the build file from folder to folder"""
    for name in source:
        sh.copyfile()


def weramizer (source, source_container, temperature, pressure):
    
    """Creates a parser file for WERAMI and uses it to get the data at defined conditions
    
    Inputs:
        source: name of the sample previously ran by VERTEX
        source_container: folder where results are stored
        temperature: temperature condition to get data from
        pressure: pressure condition to get data from
        
    Return:
        "source_name.wer": file where WERAMI output is redirected to.
        """
    # Moving files from results to bin
    sh.copy(source_container + str(source) + ".tof",bin_folder)
    sh.copy(source_container + str(source) + ".arf",bin_folder)
    sh.copy(source_container + str(source) + ".blk",bin_folder)
    sh.copy(source_container + str(source) + ".plt",bin_folder)
    sh.copy(source_container + str(source) + ".tim",bin_folder)
    sh.copy(source_container + str(source) + ".dat",bin_folder)
    sh.copy(source_container + str(source) + "_auto_refine.txt",bin_folder)


    
    filex = open ("weramizer.txt","w")
    filex.write(str(source)+ "\n1\n" + str(pressure) + " " + str(temperature) + "\n99 99\n0" )
    filex.close()

    sh.copy("weramizer.txt",bin_folder)
    print("cd " + bin_folder)
    ring = "cd " + bin_folder + " && werami <weramizer.txt >>" + str(source) + ".wer"
    subprocess.call(ring, shell=True)

    sh.copy(bin_folder + str(source) + ".wer", results_folder)

    # Delete files from bin folder
    os.remove(bin_folder + str(source) + ".tof")
    os.remove(bin_folder + str(source) + ".arf")
    os.remove(bin_folder + str(source) + ".blk")
    os.remove(bin_folder + str(source) + ".plt")
    os.remove(bin_folder + str(source) + ".tim")
    os.remove(bin_folder + str(source) + "_auto_refine.txt")
    os.remove(bin_folder + str(source) + ".wer")
    os.remove(bin_folder + str(source) + ".dat")
    os.remove(bin_folder + "weramizer.txt")
    os.remove(bin_folder + str(source) + "_1.txt")
    os.remove(bin_folder + str(source) + "_seismic_data.txt")

def mover (source, destination,frac=False):
    if not os.path.exists(destination):
        os.makedirs(destination)
    
    if frac == False:
        sh.copy(bin_folder + str(source) + ".tof",destination)
        sh.copy(bin_folder + str(source) + ".arf",destination)
        sh.copy(bin_folder + str(source) + ".blk",destination)
        sh.copy(bin_folder + str(source) + ".plt",destination)
        sh.copy(bin_folder + str(source) + ".tim",destination)
        sh.copy(bin_folder + str(source) + ".dat",destination)
        sh.copy(bin_folder + str(source) + ".txt",destination)
        sh.copy(bin_folder + str(source) + "_auto_refine.txt",destination)
        sh.copy(bin_folder + str(source) + "_seismic_data.txt",destination)

        os.remove(bin_folder + str(source) + ".tof")
        os.remove(bin_folder + str(source) + ".arf")
        os.remove(bin_folder + str(source) + ".blk")
        os.remove(bin_folder + str(source) + ".plt")
        os.remove(bin_folder + str(source) + ".tim")
        os.remove(bin_folder + str(source) + ".txt")
        os.remove(bin_folder + str(source) + "_auto_refine.txt")
        os.remove(bin_folder + str(source) + ".dat")
        os.remove(bin_folder + str(source) + "_seismic_data.txt")

    if frac == True:
        sh.copy(bin_folder + "melt(HGP)" + ".tof",destination + str(source) + ".tof")
        sh.copy(bin_folder + "melt(HGP)" + ".arf",destination + str(source) + ".arf")
        sh.copy(bin_folder + "melt(HGP)" + ".blk",destination + str(source) + ".blk")
        sh.copy(bin_folder + "melt(HGP)" + ".plt",destination + str(source) + ".plt")
        sh.copy(bin_folder + "melt(HGP)" + ".tim",destination + str(source) + ".tim")
        sh.copy(bin_folder + "melt(HGP)" + ".dat",destination + str(source) + ".dat")
        sh.copy(bin_folder + "melt(HGP)" + ".txt",destination + str(source) + ".txt")
        sh.copy(bin_folder + "melt(HGP)" + "_auto_refine.txt",destination + str(source) + "_auto_refine.txt")
        sh.copy(bin_folder + "melt(HGP)" + "_seismic_data.txt",destination + str(source) + "_seismic_data.txt")

        os.remove(bin_folder + "melt(HGP)" + ".tof")
        os.remove(bin_folder + "melt(HGP)" + ".arf")
        os.remove(bin_folder + "melt(HGP)" + ".blk")
        os.remove(bin_folder + "melt(HGP)" + ".plt")
        os.remove(bin_folder + "melt(HGP)" + ".tim")
        os.remove(bin_folder + "melt(HGP)" + ".txt")
        os.remove(bin_folder + "melt(HGP)" + "_auto_refine.txt")
        os.remove(bin_folder + "melt(HGP)" + ".dat")
        os.remove(bin_folder + "melt(HGP)" + "_seismic_data.txt")
def fortran_output_is_stupid(source):
    
    """Transforms output from WERAMI and replaces spaces for single tabs for easy readability
    
    Input:
        source: path of WERAMI output file
        
    Return:
        transformed_output: string containing all the transformed text of the source file.
        lines: list with all lines from transformed_output."""
    
    # Opening and reading file
    filex = open(source + ".wer","r")
    transformed_output = filex.read()
    filex.close()

    #Transforming spaces into tabs
    c = 20
    while c > 1:
        transformed_output = transformed_output.replace(" "*c, "\t")
        c=c-1
    
    # Separating lines
    lines = transformed_output.splitlines()
    return (transformed_output,lines)


def extractor(source, source_container="results"):
    
    """Extracts information from phases in the reformatted output of WERAMI
    
    Input:
        source: name of the WERAMI output file
        source_container: folder containing the .wer file produced by weramizer
                          defaults to "bin" meaning the bin folder
                          "results" sets it to results folder
                          Complete paths may be set as well. Example "C:/user/tomato/desktop/perp_res/"
    Returns:
        hola.csv: file with extracted data in results folder.
        new_comp: Pandas Dataframe with information to use with Auto_build for recycling of data
        """

    if source_container == "results":
        source_folder = results_folder
    elif source_container == "bin":
        source_folder = bin_folder
    else:
        source_folder = source_container
    
    transformed_output, lines = fortran_output_is_stupid(source_folder+source)
    
    #Finder
    for index,line in enumerate(lines):
        if "Phase Compositions" in line:
            line_start = index
        if "Phase speciation" in line:
            line_end = index
            break
    
    file_name = source.replace(".wer",".csv")
    filex = open(source_folder+file_name,"w")
    filex.write(lines[line_start+1].replace("\twt","Phase\twt")+"\n")
    for x in range (line_end-line_start-3):
        filex.write(lines[line_start+(x+2)].replace(" ","")+"\n")
    filex.close()
    new_comp = pd.read_csv(results_folder+file_name,delimiter = "\t", index_col=0)
    return (new_comp.T)

def custom_fractionation(source,sample,temps,pressures,probe,mix=0):
    """Initiates a fractionating sequence
    
    Input
        source: csv file containing composition/compositions.

        sample: name of the sample in source file to be run.

        temps: array (example: [[1,2], [3,4]]) containing temperature steps to be used.

        pressures: array (example: [[1,2], [3,4]]) containing pressure steps to be used.

        probe: array (example: [[1,2], [3,4]]) that contain temperature and pressure to used
             with weramizer and get the composition to be used in all steps
        
        mix: default = 0, states a mixture proportion from 0 to 1 of the remnant melt with
             the original melt for the composition of the following step.
             
    Return
        No return (yet)"""

    # Assertions
    assert(len(temps) == len(pressures)),\
            "temps and pressures should be the same length"
    assert (len(probe) == len(temps)-1),\
            "probe list should be one item less than pressures and temps"
    assert mix >= 0 and mix <=1, "mix proportion not between 0-100"

    name = pd.read_csv(source,delimiter="\t",index_col=0)
    print(name.head())
    name = name[sample]
    name.to_csv("frac.tmp",sep="\t")

    names = Auto_build("frac.tmp",temps[0],pressures[0],mp="single")

    os.remove("frac.tmp")
    if "Windows" in op_system:
        subprocess.call("cd "+str(bin_folder)+" && call batch1.bat",shell=True)
    else:
        subprocess.call("cd "+str(bin_folder)+" && source batch1.sh",shell=True)

    destination = results_folder+str(sample)+"_frac/"+"T"+str(temps[0])+"_P"+str(pressures[0])+"/"
    mover(sample,destination)
    wer_result = weramizer(sample,destination,probe[0][0],probe[0][1])
    new_input = extractor(sample)
    assert "melt(HGP)" in source.columns, "probe data does not contain melt"
    new_input["melt(HGP)"].iloc[4:13] = new_input["melt(HGP)"].iloc[4:13]*(1-mix) + name.iloc[0:9] * mix
    # print(new_input)

    # loop
    for i in range(len(temps)-1):
        print(i)
        # Auto_Build in frac mode using melt(HGP) from extracted source mixed with initial composition
        loop = Auto_build(new_input,temps[i+1],pressures[i+1],source_type="frac",mp="single")
        if "Windows" in op_system:
            subprocess.call("cd "+str(bin_folder)+" && call batch1.bat",shell=True)
        else:
            subprocess.call("cd "+str(bin_folder)+" && source batch1.sh",shell=True)
            
        destination = results_folder+str(sample)+"_frac/"+"T"+str(temps[i+1])+"_P"+str(pressures[i+1])+"/"
        mover(sample,destination,frac=True)
        if not i == len(temps)-1:
            wer_result = weramizer(sample,destination,probe[i+1][0],probe[i+1][1])
            new_input = extractor(sample)
            assert "melt(HGP)" in source.columns, "probe data does not contain melt"
            new_input["melt(HGP)"].iloc[4:13] = new_input["melt(HGP)"].iloc[4:13]*(1-mix) + name.iloc[0:9] * mix

def LMO_type_fractionation(source,sample,T_limits,P_limits,T_step,P_step,mix="equilibrated",trapped_melt=1,
                            density="moon",gravity="moon",radius="moon",**kwargs):
    """Initiates a fractionation sequence following the rules for magma ocean crystallization:
        A shell of determined thickness at the initial conditions is equilibrated at the conditions
        of its upper edge, solid fraction is removed with an amount of trapped melt, the remaining melt
        is equilibrated with the remaining magma ocean and used for the following shell that starts at
        the top of the solid fraction depth.
        
        Input:
            source: csv file containing composition/compositions.

            sample: name of the sample in source file to be run.
            
            T_limits: tuple, (T of initial depth, T at final depth)
            
            P_limits: tuple, (P at initial depth, P at final depth)
            
            T_step: temperature difference for each modelled shell, choosing "gradient" establishes
                a continous gradient between T_limits and calculates T from the depth of each shell
            
            P_step: pressure difference at the limits of the shell, i.e. pressure "thickness" of each shell
            
            mix: default "equilibrated", mixes the remnant melt from any run with the initial composition
                or the reequilibrated magma ocean. "equilibrated" states the mixture is proportional to 
                the volume of the remnant melt and the remaining magma ocean (e.g. 50% remnant melt of a ring 
                that was 1% of the total ocean mixed with the remaining 99% of the initial composition). 
                "constant" states that the sourrounding ocean does not requilibrate and the mixture always 
                mixes with the initial composition at a certain rate (requires mix_rate variable in kwargs) 
                
            trapped_melt: number from 0 to 100, defaults to 1, percent of melt trapped within each solid fraction of shell.
            
        Return:
            no return """
        
        
    assert "equilibrated" in mix or "constant" in mix, \
        "Only equilibrated or constant values are supported for mix"
    
    if "gradient" not in T_step:
        assert T_limits[1]-T_limits[0]<=T_step, "T_step is larger than limit differences"
    
    assert P_limits[1]-P_limits[0]<=P_step, "P_step is larger than limit differences"

    if "constant" in mix:
        assert "mix_rate" in kwargs, "missing mix_rate argument"

    # Isolate selected composition for Autobuild
    name = pd.read_csv(source,delimiter="\t",index_col=0)
    print(name.head())
    name = name[sample]

    # Initial variables
    initial_P = P_limits[0]
    initial_T = T_limits[0]
    final_T = T_limits[1]
    final_P = P_limits[1]
    first_shell_P = initial_P - P_step

    if T_step == "gradient":
        gradient = (initial_T-final_T)/(initial_P-final_P)
        first_shell_T = gradient * (first_shell_P)
    else:
        first_shell_T = initial_T - T_step
    
    if density == "moon":
        rho = 3500
    elif type(density) == float or type(density) == int and density > 0:
        rho = density
    else:
        raise (TypeError("only moon density is supported right now"))
    
    if gravity == "moon":
        g = 1.7
    elif type(gravity) == float or type(gravity) == int and gravity > 0:
        g = gravity
    else:
        raise (TypeError("only moon gravity is supported right now"))
    
    if radius == "moon":
        R = 1737e3
    elif type(radius) == float or type(radius) == int and radius > 0:
        R = radius
    else:
        raise (TypeError("only moon radius is supported right now"))
    
    initial_shell_volume = 4/3*np.pi *((radius-(first_shell_P/(rho*g)))**3 - (radius-(initial_P/(rho*g)))**3)
    # Import and Auto_build isolated composition
    name.to_csv("frac.tmp",sep="\t")
    names = Auto_build("frac.tmp",[initial_T,first_shell_T],[initial_P,first_shell_P],mp="single")

    # Remove temporal file of isolated composition
    os.remove("frac.tmp")

    # In-thread vertex run of the built file
    if "Windows" in op_system:
        subprocess.call("cd "+str(bin_folder)+" && call batch1.bat",shell=True)
    else:
        subprocess.call("cd "+str(bin_folder)+" && source batch1.sh",shell=True)

    # Creting destination folder and moving results
    destination = results_folder+str(sample)+"_LMO/"+"T"+str(initial_T)+"_P"+str(initial_P)+"/"
    mover(sample,destination)

    # Extracting data using werami
    wer_result = weramizer(sample,destination,initial_T+T_step,initial_P+P_step)
    new_input = extractor(sample)

    # Finding solid fraction
    assert "melt(HGP)" in source.columns, "probe data does not contain melt"
    melt_fraction_vol = new_input["melt(HGP)"].loc["vol %"]
    solid_fraction_vol = 100 - melt_fraction_vol

    # Establishing new shell conditions
    new_shell_P = initial_P - (solid_fraction_vol * P_step / 100)
    if T_step == "gradient":
        new_shell_T = gradient * new_shell_P
    else:
        new_shell_T = initial_T + T_step
    
    # Mixing: proportion of remnant melt and ocean
    remaining_ocean_volume = 4/3*np.pi *((radius-(final_P/(rho*g)))**3 - (radius-(first_shell_P/(rho*g)))**3)
    remaining_melt_volume = initial_shell_volume * melt_fraction_vol / 100
    total_remaining_vol = remaining_melt_volume + remaining_ocean_volume
    remaining_ocean_percent = remaining_ocean_volume / total_remaining_vol
    remaining_melt_percent = remaining_melt_volume / total_remaining_vol
    new_input["melt(HGP)"].iloc[4:13] = new_input["melt(HGP)"].iloc[4:13]*remaining_melt_percent + name.iloc[0:9] * remaining_ocean_percent
   
    # Loop variables
    new_shell_base_P = new_shell_P
    new_shell_top_P = new_shell_P - P_step
    new_shell_base_T = new_shell_T
    new_shell_top_T = new_shell_top_P * gradient
    shell = 2
    ocean_composition = new_input["melt(HGP)"].iloc[4:13]
    # Loop
    while new_shell_top_P <= final_P:

        # Auto_Build in frac mode using new mix
        loop = Auto_build(new_input,[new_shell_base_T,new_shell_top_T],[new_shell_base_P,new_shell_top_P],source_type="frac",mp="single")
        if "Windows" in op_system:
            subprocess.call("cd "+str(bin_folder)+" && call batch1.bat",shell=True)
        else:
            subprocess.call("cd "+str(bin_folder)+" && source batch1.sh",shell=True)
            
        destination = results_folder+str(sample)+"_LMO/ring" + str(shell) + "/"
        mover(sample,destination,frac=True)
        wer_result = weramizer(sample,destination,new_shell_top_T,new_shell_top_P)
        new_input = extractor(sample)

        # Finding solid fraction
        assert "melt(HGP)" in source.columns, "probe data does not contain melt"
        melt_fraction_vol = new_input["melt(HGP)"].loc["vol %"]
        solid_fraction_vol = 100 - melt_fraction_vol

        # Establishing new shell conditions
        new_shell_P = initial_P - (solid_fraction_vol * P_step / 100)
        if T_step == "gradient":
            new_shell_T = gradient * new_shell_P
        else:
            new_shell_T = initial_T + T_step
        
        # Mixing: proportion of remnant melt and ocean
        remaining_ocean_volume = 4/3*np.pi *((radius-(final_P/(rho*g)))**3 - (radius-(new_shell_top_P/(rho*g)))**3)
        remaining_melt_volume = 4/3*np.pi *((radius-(new_shell_top_P/(rho*g)))**3 -
                                                    (radius-(new_shell_base_P/(rho*g)))**3) * melt_fraction_vol / 100
        total_remaining_vol = remaining_melt_volume + remaining_ocean_volume
        remaining_ocean_percent = remaining_ocean_volume / total_remaining_vol
        remaining_melt_percent = remaining_melt_volume / total_remaining_vol
        new_input["melt(HGP)"].iloc[4:13] = new_input["melt(HGP)"].iloc[4:13]*remaining_melt_percent + ocean_composition * remaining_ocean_percent
        
        # New loop variables
        new_shell_base_P = new_shell_P
        new_shell_top_P = new_shell_P - P_step
        new_shell_base_T = new_shell_T
        new_shell_top_T = new_shell_top_P * gradient
        shell = shell + 1
        ocean_composition = new_input["melt(HGP)"].iloc[4:13]
