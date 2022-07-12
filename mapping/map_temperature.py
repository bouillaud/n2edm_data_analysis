import numpy as np
import os.path
import matplotlib.pyplot as plt
from map_tools import *

# plot parameters
plt.rcParams["font.size"]=10
plt.rcParams["lines.markersize"]=3
plt.rcParams["figure.figsize"]=[12,8]
plt.rcParams["figure.dpi"]=200

class Temperature:
    
    def __init__(self, run, cyc, calib=True):
        self.path = "../maps_bin/{}_{}_000_temperature.EDMdat".format(str(run).zfill(6), str(cyc).zfill(6)) 
        if os.path.exists(self.path):
            self.getData()
            self.file_exists=True
        else:
            self.file_exists=False
            return None
        if calib==True:
            self.calibrate()
            self.calibrated=True
        else:
            self.calibrated=False
        self.tempGrad = (self.tempTop - self.tempBot) / 157 # in °C/cm
        self.avgTop, self.avgBot = np.mean(self.tempTop), np.mean(self.tempBot)
        self.avgBack, self.avgBot2 = np.mean(self.tempBack), np.mean(self.tempBot2)
        self.stdTop, self.stdBot = np.std(self.tempTop), np.std(self.tempBot)
        self.stdBack, self.stdBot2 = np.std(self.tempBack), np.std(self.tempBot2)
        self.avgGrad = np.mean(self.tempGrad)
        self.stdGrad = np.std(self.tempGrad)
        ### new format
        self.avg = np.mean([self.tempBot, self.tempTop, self.tempBot2, self.tempBack], axis=1)
        self.std = np.std([self.tempBot, self.tempTop, self.tempBot2, self.tempBack], axis=1)
    
    def getData(self, dtformat='v2', ncol=50):
            # define data type according to the header. there is an extra column ('other') not specified in the header.
            if dtformat=='v2':
                dtype = np.dtype([("time", np.uint64)] 
                                 + [("temp {}".format(i), np.double) for i in range(1, ncol-1)] 
                                 + [("digital inputs", np.uint64)] 
                                 + [("other", np.uint64)])

            with open(self.path, "rb") as f:
                data = np.fromfile(f, dtype=dtype)
                
            self.t = 1e-9*np.array(data["time"] - data["time"][0])
            self.tempBot = data["temp 33"]
            self.tempTop = data["temp 34"]
            self.tempBot2 = data["temp 35"]
            self.tempBack = data["temp 36"]
            self.data = data
            
    def calibrate(self, offsets=[-0.019, -0.005, 0.001]):
        self.tempBot = self.tempBot ## calibration reference
        self.tempTop = self.tempTop + offsets[0]
        self.tempBot2 = self.tempBot2 + offsets[1]
        self.tempBack = self.tempBack + offsets[2]
            
class TemperatureRun:
    
    def __init__(self, run, cycrange, calib=True):
        cycles = {}
        tempsTop = []
        tempsBot = []
        tempsBack = []
        tempsBot2 = []
        tempsGrad = []
        ### new format
        avgs = []
        stds = []
        for cyc in cycrange:
            cycle = Temperature(run, cyc, calib=calib)
            if cycle.file_exists:
                cycles[cyc] = cycle
                tempsBot.append(cycle.avgBot)
                tempsTop.append(cycle.avgTop)
                tempsBot2.append(cycle.avgBot2)
                tempsBack.append(cycle.avgBack)
                tempsGrad.append(cycle.avgGrad)
                ### new format
                avgs.append(cycle.avg)
                stds.append(cycle.std)
        self.cycles = cycles
        self.tempsTop = np.array(tempsTop)
        self.tempsBot = np.array(tempsBot)
        self.tempsBack = np.array(tempsBack)
        self.tempsBot2 = np.array(tempsBot2)
        self.tempsGrad = np.array(tempsGrad)
        ### new format
        self.cycAvgs = np.transpose(avgs)
        self.cycStds = np.transpose(stds)
        self.avgs = np.mean(avgs, axis=0)
        self.stds = np.mean(stds, axis=0)
        
class TemperatureRunSet:
    
    def __init__(self, runrange, cycrange, calib=True):
        runs = {}
        tempsTop = []
        tempsBot = []
        tempsBack = []
        tempsBot2 = []
        tempsGrad = []
        ### new format
        avgs = []
        stds = []
        for r in runrange:
            run = TemperatureRun(r, cycrange, calib=calib)
            if len(run.tempsTop)>0:
                runs[r] = run
                tempsBot.append(np.mean(run.tempsBot))
                tempsTop.append(np.mean(run.tempsTop))
                tempsBot2.append(np.mean(run.tempsBot2))
                tempsBack.append(np.mean(run.tempsBack))
                tempsGrad.append(np.mean(run.tempsGrad))
                ### new format
                avgs.append(run.avgs)
                stds.append(run.stds)
        self.runs = runs
        self.tempsTop = np.array(tempsTop)
        self.tempsBot = np.array(tempsBot)
        self.tempsBack = np.array(tempsBack)
        self.tempsBot2 = np.array(tempsBot2)
        self.tempsGrad = np.array(tempsGrad)
        ### new format
        self.runAvgs = np.transpose(avgs)
        self.runStds = np.transpose(stds)
            