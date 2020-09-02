import matplotlib.pyplot as plt
import matplotlib.animation as animation
from time import time
from random import randrange
import warnings
warnings.filterwarnings('ignore')
import os
import sys
import datetime
from math import ceil
import csv

MicrobeContents = {} #dictionary contains names of all dictionaries with bacterial contents

#data is structures so that an item is the key, and the values are nested lists.
#Values are [[bool, bool, dictionary ID][row,column or None][conc over time or time]]. 
# [0] bool is for RT monitoring, and [1] bool is for output monitoring
data = {"Time" : [[False, True],[None]]}
settings = {"RowColumnDisplay" : [1,1], "Display" : True, "Output" : False}

class DataHandler:

    OutputFile = None #the file information is writen to
    OutputCSV = None

    def IdentifyDictionary(self):
        #this function identifies which metabolites/properties belong to which dictionary in the MicrobeContents
        #which should reduce processing time so other functions do not need to check each time to pull data.
        [[data[item][0].append(dictionary) for dictionary in MicrobeContents.keys() if item in MicrobeContents[dictionary]] for item in data.keys()]

    def CreateOutputCSV(self):
        DataHandler.OutputCSV = open(str("./VirtualMicrobeOutputs/VirtualMicrobe_Output_"+ datetime.datetime.now().strftime("%d_%b_%Y_%H_%M")), 'w', newline="")
        DataHandler.OutputFile = csv.writer(DataHandler.OutputCSV)
        DataHandler.OutputFile.writerow(item for item in data.keys() if data[item][0][1])

    def OutputDataToCSV(self):
        DataHandler.OutputFile.writerow(data[item][2][-1:][0] for item in data.keys() if data[item][0][1])

    def CloseOutputFile(self):
        DataHandler.OutputCSV.close()

    def AxesGenerator(self): #this function generates the row, columns for each axis in the plot
        for row in range(3):
            for col in range(2):
                yield [row, col]

    def ConfigSettings(self, row):

        GenAxes = DataHandler.AxesGenerator(self)

        if row.startswith("@MonitorRT:") and len(row[11:].strip()) > 0:
            try:
                for item in row[11:].split(","):
                    data[item.strip()] = [[True, False], next(GenAxes)]
            except:
                data[row[11:].strip()] = [[True, False], next(GenAxes)]
            finally:
                pass

            if len(data.keys()) > 7: #ensures user has not provided too many items to monitor, reason it is seven is to include "Time"
                print("Error: Too many items provided to be monitored. Maximum is six.")
                sys.exit()

            else: #this determines how many rows and columns are required for the RT figure
                if len(data.keys()) > 2: #if there is more than one item to monitor, then specifies the columns and rows
                    settings["RowColumnDisplay"] = [ceil((len(data.keys())-1)/2), 2]
                elif len(data.keys()) == 1:
                    settings["Display"] = False #will not display RT figures if none are given
                else:
                    settings["RowColumnDisplay"] = [1,1]

                
        elif row.startswith("@Monitor") and len(row[9:].strip()) > 0:

            try:
                for item in row[9:].split(","):
                    if item.strip() not in data:
                        data[item.strip()] = [[False, True],[None]]
                    else:
                        data[item.strip()][0][1] = True

            except:
                if row[9:].strip() not in data:
                    data[row[9:].strip()] = [[False, True],[None]]
                else:
                    data[row[9:].strip()][0][1] = True

            DataHandler.CreateOutputCSV(self)
            settings["Output"] = True

        elif row.startswith("@Settings") and len(row[10:].strip()) > 0:
            try:
                for item in row[10:].split(","):
                    name, prop = item.split(" ")
                    settings[name.strip()] = int(prop)
                
            except:
                name, prop = row[11:].split(" ")
                settings[name.strip()] = int(prop.strip())
        
    def TakeMonitoredData(self): #this function takes the data from the metabolite/property pool and adds them to data to be displayed or outputted

        for MonitoredItem in data.keys():
            if MonitoredItem == "Time":
                try:
                    data["Time"][2].append(round(time() - DataPlot.TimeStarted))
                except:
                    data["Time"].append([round(time() - DataPlot.TimeStarted)])
            else:
                try:
                    #this adds the metabolite conc from the pool and adds to the data dictionary
                    data[MonitoredItem][2].append(MicrobeContents[data[MonitoredItem][0][2]][MonitoredItem])
                except: 
                    data[MonitoredItem].append([MicrobeContents[data[MonitoredItem][0][2]][MonitoredItem]])

    def GetConfigs(self):
        if os.path.exists("./Virtual Microbe Configuration File.txt"):
            Config = open("./Virtual Microbe Configuration File.txt")

        else:
            print("""Error: Could not detect configuration file.
            Please ensure you have the configuration file in current directory.
            (Requires text file: Virtual Microbe Configuration File.txt)""")
            sys.exit()

        PlaceHolder = ""

        for row in Config:
            if row.startswith("/"):
                PlaceHolder = row[1:].strip()
                MicrobeContents[PlaceHolder] = { }
            elif row.startswith("@"):
                DataHandler.ConfigSettings(self, row)
            elif row.startswith('"'):
                item, value = row.split(",")
                MicrobeContents[PlaceHolder][item[1:-1]] = int(value.strip())

        DataHandler.IdentifyDictionary(self)

        Config.close()

#Randomly alters the metabolite concentrations by 1
def mixup():
    for dictionary in MicrobeContents.keys():
        for key, value in MicrobeContents[dictionary].items():
            if 5 < randrange(0, 10):
                MicrobeContents[dictionary][key] = value + 1
            else:
                MicrobeContents[dictionary][key] = value - 1

class DataPlot:

    TimeStarted = 0
    ani = None #animation
    fig = None #figure
    axes = None

    def CreatePlot(self):
        DataPlot.fig, DataPlot.axes = plt.subplots(nrows = settings["RowColumnDisplay"][0], ncols = settings["RowColumnDisplay"][1], constrained_layout=True)

    def InitiatePlot(self):
        DataPlot.TimeStarted = time()
        DataPlot.ani = animation.FuncAnimation(DataPlot.fig, DataPlot.PlotMetabolites, interval=5000, fargs=(self))
        if settings["Display"]:
            plt.show()

    def PlotMetabolites(self):

        mixup()
        DataHandler.TakeMonitoredData(self=None)
    
        for metabolite in data.keys():
            if data[metabolite][0][0]: #checks if metabolite is meant to be monitored in RT

                #the data[metabolite][1][0] and [1][1] is the row and column of the subplot
                DataPlot.axes[data[metabolite][1][0], data[metabolite][1][1]].clear()

                #Figure formatting
                DataPlot.axes[data[metabolite][1][0], data[metabolite][1][1]].plot(data["Time"][2][-10:], data[metabolite][2][-10:])
                DataPlot.axes[data[metabolite][1][0], data[metabolite][1][1]].set_title(metabolite, fontsize=12)
                DataPlot.axes[data[metabolite][1][0], data[metabolite][1][1]].set_xlabel("Time (seconds)", fontsize=10)
                DataPlot.axes[data[metabolite][1][0], data[metabolite][1][1]].set_ylabel("Conc (units)", fontsize=10)

        if settings["Output"]:
            DataHandler.OutputDataToCSV(self=None)

DataHandler.GetConfigs(self=None)

if settings["Display"] or settings["Output"]:
    DataPlot.CreatePlot(self=None)
    DataPlot.InitiatePlot(self=None)
    
   
