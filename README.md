# rmsd_plot
Plot the MD simulation
import matplotlib.pyplot as plt
import numpy as np
import csv
import glob
import fileinput
import pandas as pd
import matplotlib.cm as cm
#loading file
exclude = []
files_rmsd = []
for x in range (1,21):
    if x in exclude:
        continue
    files_rmsd.append("/work9/haont/copy-update/500/glu/trial"+str(x)+"/trial"+str(x)+"rmsd.xvg")
filelist_rmsd = np.array(files_rmsd)

files_t1r2 = []
for x in range (1,21):
    if x in exclude:
        continue
    files_t1r2.append("/work9/haont/copy-update/500/glu/trial"+str(x)+"/trial"+str(x)+"t1r2.xvg")
filelist_t1r2 = np.array(files_t1r2)

files_t1r3 = []
for x in range (1,21):
    if x in exclude:
        continue
    files_t1r3.append("/work9/haont/copy-update/500/glu/trial"+str(x)+"/trial"+str(x)+"t1r3.xvg")
filelist_t1r3 = np.array(files_t1r3)
#create a list of GLUA files
files_glua = []
for x in range (1,21):
    if x in exclude:
        continue
    files_glua.append("/work9/haont/copy-update/500/glu/trial"+str(x)+"/trial"+str(x)+"rmsd_GLUA.xvg")
filelist_glua = np.array(files_glua)

#create a list of GLUB files
files_glub = []
for x in range (1,21):
    if x in exclude:
        continue
    files_glub.append("/work9/haont/copy-update/500/glu/trial"+str(x)+"/trial"+str(x)+"rmsd_GLUB.xvg")
filelist_glub = np.array(files_glub)

#create a list of dis2018 files
files_dis2018 = []
for x in range (1,21):
    if x in exclude:
        continue
    files_dis2018.append("/work9/haont/copy-update/500/glu/trial"+str(x)+"/trial"+str(x)+"distave2018.xvg")
filelist_dis2018 = np.array(files_dis2018)

#create a list of dis2019 files
files_dis2019 = []
for x in range (1,21):
    if x in exclude:
        continue
    files_dis2019.append("/work9/haont/copy-update/500/glu/trial"+str(x)+"/trial"+str(x)+"distave2019.xvg")
filelist_dis2019 = np.array(files_dis2019)

#create a list of dis2118 files
files_dis2118 = []
for x in range (1,21):
    if x in exclude:
        continue
    files_dis2118.append("/work9/haont/copy-update/500/glu/trial"+str(x)+"/trial"+str(x)+"distave2118.xvg")
filelist_dis2118 = np.array(files_dis2118)

#create a list of dis2119 files
files_dis2119 = []
for x in range (1,21):
    if x in exclude:
        continue
    files_dis2119.append("/work9/haont/copy-update/500/glu/trial"+str(x)+"/trial"+str(x)+"distave2119.xvg")
filelist_dis2119 = np.array(files_dis2119)

# viridis to assign different color to each line plot
# get_cmap to obtain color map
color_map = cm.get_cmap('viridis')
#obtain the color by passing the value of i/100 as an argument

#subplots function create the 2D array of subplots with 10 rows 10 columns
#flaten the subplot by ravel function
fig, axs = plt.subplots(nrows=5,ncols=4, figsize=(40, 20))
axs = axs.ravel()


for i in range(0,len(filelist_rmsd)):
    row = i // 2
    col = i % 2
    df1 = pd.read_csv(filelist_rmsd[i], skiprows=18, delim_whitespace=True, header=None)
    df2 = pd.read_csv(filelist_t1r2[i], skiprows=18, delim_whitespace=True, header=None)
    df3 = pd.read_csv(filelist_t1r3[i], skiprows=18, delim_whitespace=True, header=None)
    color1 = color_map(i/(len(filelist_rmsd)))
    color2 = color_map(i/(len(filelist_t1r2)))
    color3 = color_map(i/(len(filelist_t1r3)))
    x1 = df1.iloc[:, 0]
    y1 = df1.iloc[:, 1]
    x2 = df2.iloc[:, 0]
    y2 = df2.iloc[:, 1]
    x3 = df3.iloc[:, 0]
    y3 = df3.iloc[:, 1]
    axs[i].plot(x1, y1, label='protein')
    axs[i].plot(x2, y2, label='t1r2')
    axs[i].plot(x3, y3, label='t1r3')
    axs[i].set_ylim(0, 1.5)
    axs[i].set_xlabel('Time(ns')
    axs[i].set_ylabel("RMSD (nm)")
    axs[i].legend()
    axs[i].set_title(f'trial {i+1}')
    axs[i].tick_params(axis='y', labelsize=24)
    axs[i].tick_params(axis='x', labelsize=24)

plt.suptitle('RMSD protein')
plt.tight_layout()
plt.savefig("glu_rmsd_all.png")


# viridis to assign different color to each line plot
# get_cmap to obtain color map
color_map = cm.get_cmap('viridis')
#obtain the color by passing the value of i/100 as an argument

#subplots function create the 2D array of subplots with 10 rows 10 columns
#flaten the subplot by ravel function
fig, axs = plt.subplots(5,4, figsize=(40, 20))
axs = axs.ravel()


for i in range(0,len(filelist_rmsd)):
    df = pd.read_csv(filelist_rmsd[i], skiprows=18, delim_whitespace=True, header=None)
    color = color_map(i/(len(filelist_dis2119)))
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]

    axs[i].plot(x, y, color = color, label='Trial {}'.format(i+1))
    axs[i].set_ylim(0, 1.5)
    axs[i].set_xlabel('Time(ns')
    axs[i].set_ylabel("RMSD (nm)")
    axs[i].tick_params(axis='y', labelsize=24)
    axs[i].tick_params(axis='x', labelsize=24)
    axs[i].legend()

plt.suptitle('RMSD protein')
plt.tight_layout()
plt.savefig("glu_rmsd.png")
# viridis to assign different color to each line plot
# get_cmap to obtain color map
color_map = cm.get_cmap('viridis')
#obtain the color by passing the value of i/100 as an argument

#subplots function create the 2D array of subplots with 10 rows 10 columns
#flaten the subplot by ravel function
fig, axs = plt.subplots(5,4, figsize=(40, 20))
axs = axs.ravel()


for i in range(0,len(filelist_rmsd)):
    df = pd.read_csv(filelist_rmsd[i], skiprows=18, delim_whitespace=True, header=None)
    color = color_map(i/(len(filelist_dis2119)))
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]

    axs[i].plot(x, y, color = color, label='Trial {}'.format(i+1))
    axs[i].set_ylim(0, 1.5)
    axs[i].set_xlabel('Time(ns')
    axs[i].set_ylabel("RMSD (nm)")
    axs[i].tick_params(axis='y', labelsize=24)
    axs[i].tick_params(axis='x', labelsize=24)
    axs[i].legend()

plt.suptitle('RMSD protein')
plt.tight_layout()
plt.savefig("glu_rmsd.png")

# viridis to assign different color to each line plot
# get_cmap to obtain color map
color_map = cm.get_cmap('viridis')
#obtain the color by passing the value of i/100 as an argument

#subplots function create the 2D array of subplots with 10 rows 10 columns
#flaten the subplot by ravel function
fig, axs = plt.subplots(5,4, figsize=(40, 20))
axs = axs.ravel()


for i in range(0,len(filelist_t1r2)):
    df = pd.read_csv(filelist_t1r2[i], skiprows=18, delim_whitespace=True, header=None)
    color = color_map(i/(len(filelist_dis2119)))
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]

    axs[i].plot(x, y, color = color, label='Trial {}'.format(i+1))
    axs[i].set_ylim(0, 1.5)
    axs[i].set_xlabel('Time(ns')
    axs[i].set_ylabel("RMSD  T1R2LBD(nm)")
    axs[i].tick_params(axis='y', labelsize=24)
    axs[i].tick_params(axis='x', labelsize=24)
    axs[i].legend()

plt.suptitle('RMSD protein')
plt.tight_layout()
plt.savefig("glu_t1r2.png")


# viridis to assign different color to each line plot
# get_cmap to obtain color map
color_map = cm.get_cmap('viridis')
#obtain the color by passing the value of i/100 as an argument

#subplots function create the 2D array of subplots with 10 rows 10 columns
#flaten the subplot by ravel function
fig, axs = plt.subplots(5,4, figsize=(40, 20))
axs = axs.ravel()


for i in range(0,len(filelist_rmsd)):
    df = pd.read_csv(filelist_rmsd[i], skiprows=18, delim_whitespace=True, header=None)
    color = color_map(i/(len(filelist_dis2119)))
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]

    axs[i].plot(x, y, color = color, label='Trial {}'.format(i+1))
    axs[i].set_ylim(0, 1.5)
    axs[i].set_xlabel('Time(ns')
    axs[i].set_ylabel("RMSD (nm)")
    axs[i].tick_params(axis='y', labelsize=24)
    axs[i].tick_params(axis='x', labelsize=24)
    axs[i].legend()

plt.suptitle('RMSD T1TR3LBD')
plt.tight_layout()
plt.savefig("glu_t1r3.png")


# plot GLUA

#subplots function create the 2D array of subplots with 10 rows 10 columns
#flaten the subplot by ravel function
fig, axs = plt.subplots(5,4, figsize=(40, 20))
axs = axs.ravel()


for i in range(0,len(filelist_glua)):
    df = pd.read_csv(filelist_glua[i], skiprows=18, delim_whitespace=True, header=None)
    color = color_map(i/(len(filelist_dis2119)))
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]

    axs[i].plot(x, y, color = color, label='Trial {}'.format(i+1))
    axs[i].set_ylim(0, 12)
    axs[i].set_xlabel('Time(ns')
    axs[i].set_ylabel("RMSD (nm)")
    axs[i].tick_params(axis='y', labelsize=24)
    axs[i].tick_params(axis='x', labelsize=24)
    axs[i].legend()

plt.suptitle("RMSD GLUA")
plt.tight_layout()
plt.savefig("glu_glua.png")

# plot GLUB
#subplots function create the 2D array of subplots with 10 rows 10 columns
#flaten the subplot by ravel function
fig, axs = plt.subplots(5,4, figsize=(40, 20))
axs = axs.ravel()


for i in range(0,len(filelist_glub)):
    df = pd.read_csv(filelist_glub[i], skiprows=18, delim_whitespace=True, header=None)
    color = color_map(i/(len(filelist_dis2119)))
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]

    axs[i].plot(x, y, color = color, label='Trial {}'.format(i+1))
    axs[i].set_ylim(0, 12)
    axs[i].set_xlabel('Time(ns')
    axs[i].set_ylabel("RMSD (nm)")
    axs[i].tick_params(axis='y', labelsize=24)
    axs[i].tick_params(axis='x', labelsize=24)
    axs[i].legend()

plt.suptitle("RMSD GLUB")
plt.tight_layout()
plt.savefig("glu_glub.png")

# plot dis2018
#subplots function create the 2D array of subplots with 10 rows 10 columns
#flaten the subplot by ravel function
fig, axs = plt.subplots(5,4, figsize=(40, 20))
axs = axs.ravel()


for i in range(0,len(filelist_dis2018)):
    df = pd.read_csv(filelist_dis2018[i], skiprows=24, delim_whitespace=True, header=None)
    color = color_map(i/(len(filelist_dis2119)))
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]

    axs[i].plot(x, y, color = color, label='Trial {}'.format(i+1))
    axs[i].set_ylim(0, 12)
    axs[i].set_xlabel('Time(ns')
    axs[i].set_ylabel("Distance (nm)")     
    axs[i].tick_params(axis='y', labelsize=24)     
    axs[i].tick_params(axis='x', labelsize=24)
    axs[i].legend()

plt.suptitle("COM distance GLUA - T1R2")
plt.tight_layout()
plt.savefig("glu_dis2018.png")

# plot dis2019
#subplots function create the 2D array of subplots with 10 rows 10 columns
#flaten the subplot by ravel function
fig, axs = plt.subplots(5,4, figsize=(40, 20))
axs = axs.ravel()


for i in range(0,len(filelist_dis2019)):
    df = pd.read_csv(filelist_dis2019[i], skiprows=24, delim_whitespace=True, header=None)
    color = color_map(i/(len(filelist_dis2119)))
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]

    axs[i].plot(x, y, color = color, label='Trial {}'.format(i+1))
    axs[i].set_ylim(0, 12)
    axs[i].set_xlabel('Time(ns)')
    axs[i].set_ylabel("Distance (nm)")
    axs[i].tick_params(axis='y', labelsize=24)
    axs[i].tick_params(axis='x', labelsize=24)
    axs[i].legend()

    
plt.suptitle("COM distance GLUA - T1R3")
plt.tight_layout()
plt.savefig("glu_dis2019.png")

# plot dis2118
#subplots function create the 2D array of subplots with 10 rows 10 columns
#flaten the subplot by ravel function
fig, axs = plt.subplots(5,4, figsize=(40, 20))
axs = axs.ravel()


for i in range(0,len(filelist_dis2118)):
    df = pd.read_csv(filelist_dis2118[i], skiprows=24, delim_whitespace=True, header=None)
    color = color_map(i/(len(filelist_dis2119)))
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]

    axs[i].plot(x, y, color = color, label='Trial {}'.format(i+1))
    axs[i].set_ylim(0, 12)
    axs[i].set_xlabel('Time(ns)')
    axs[i].set_ylabel("Distance (nm)")     
    axs[i].tick_params(axis='y', labelsize=24)     
    axs[i].tick_params(axis='x', labelsize=24)
    axs[i].legend()
    

plt.suptitle("COM distance GLUB - T1R2")
plt.tight_layout()
plt.savefig("glu_dis2118.png")


# plot dis2119
#subplots function create the 2D array of subplots with 10 rows 10 columns
#flaten the subplot by ravel function
fig, axs = plt.subplots(5,4, figsize=(40, 20))
axs = axs.ravel()


for i in range(0,len(filelist_dis2119)):
    df = pd.read_csv(filelist_dis2119[i], skiprows=24, delim_whitespace=True, header=None)
    color = color_map(i/(len(filelist_dis2119)))
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]

    axs[i].plot(x, y, color = color, label='Trial {}'.format(i+1))
    axs[i].set_ylim(0, 12)
    axs[i].set_xlabel('Time(ns)')
    axs[i].set_ylabel("Distance (nm)")
    axs[i].tick_params(axis='y', labelsize=24)
    axs[i].tick_params(axis='x', labelsize=24)
    axs[i].legend()

plt.suptitle("COM distance GLUB - T1R3")
plt.tight_layout()
plt.savefig("glu_dis2119.png")



