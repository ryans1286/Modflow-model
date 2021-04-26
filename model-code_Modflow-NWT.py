# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 13:37:36 2021

@author: rstrickland
"""

'''
YOU MUST DOWNLOAD THE MODFLOW-NWT EXECUTABLE TO RUN THIS FILE
'''

import numpy as np
import flopy as fp
import matplotlib.pyplot as plt

#MODFLOW-NWT executable path
path = ""

#NAME THE MODEL
modelName = "simulation-name"

#INPUT MODEL DIMENSIONS

Lx = 1000. #distances in meters
Ly = 1000.

ztop = 100. #maximum height in meters
zbot = 0.

nrow = 100 
ncol = 100
nlay = 1

recharge = 0.001 #define the recharge rate in units of m/year

#define the hydraulic conductivity in the horizontal and vertical direction
horiz_cond = 1.
vert_cond = 1.

drow = Lx / nrow #cell dimensions in meters
dcol = Ly / ncol

#fill arrays with the top and bottom of each cell
botm = np.ones(shape = (nrow, ncol)) * zbot
top = np.ones(shape = (nrow, ncol)) * ztop


#CREATE SURFACE TOPOGRAPHY
for i in range(nrow):
    #this loop creates an inverted V-shape moraine
    if i <= 20:
        top[i, :] = top[i, :] - 0.5 * (20-i)
    else:
        top[i, :] = top[i, :] - 1 * (i-20)


#CREATE SUBSURFACE TOPOGRAPHY

#this surface is impermeable and defines the bottom of the aquifer
botm[1:20, :] = 70
botm[20:60, :] = 50
botm[0,:] = 60


#DEFINE BOUNDARY CONDITIONS

ibound = np.ones((nlay, nrow, ncol), dtype = np.int32) #create a boundary condition array
ibound[:, 0, :] = -1 #make first row a constant head boundary, defining the lake condition

strt = np.ones((nlay, nrow, ncol), dtype = np.float32) #create starting head array
strt[:, :, :] = top #set the starting head to the elevation of the cell

#DEFINE SURFACE DRAINAGE

#Allows groundwater to discharge on the surface
drainList = []
for i in range(1, nrow):
    #The tops of all cells are drains, excluding the first row. 
    for j in range(ncol):
        rowList = [0, i, j, top[i][j], 10000] #set the conductance (10000) to a very high value
        drainList.append(rowList)
    
drnPeriod = {0 : drainList} #make it a dictionary, should apply to all periods

#INSTANTIATE ALL THE MODULES NEEDED TO RUN THE MODEL

mf = fp.modflow.Modflow(modelName, exe_name = path, version = "mfnwt") #create the model and define the executable
drn = fp.modflow.ModflowDrn(mf, stress_period_data = drnPeriod) #surface drainage
bas = fp.modflow.ModflowBas(mf, ibound = ibound, strt = strt) #cell heads

#discretizes the model
dis = fp.modflow.ModflowDis(mf, nlay, nrow, ncol, delr = drow, delc = dcol, top = top, botm = botm, nstp = 1, itmuni = 5, steady = True)

spd = {(0, 0): ["print budget", "save head", "save budget"]} #stores stress period data
upw = fp.modflow.ModflowUpw(mf, hk = horiz_cond, vka = vert_cond) #define conductivity
rch = fp.modflow.ModflowRch(mf, nrchop = 1, rech = recharge) #define recharge
nwt = fp.modflow.ModflowNwt(mf) #define the MODFLOW solver
oc = fp.modflow.ModflowOc(mf, stress_period_data = spd, compact = True) #defines output storage

mf.write_input() #write the input files to run the code


#RUN THE MODEL

success, buff = mf.run_model()
if not success:
    raise Exception("MODFLOW did not terminate normally.")
    
