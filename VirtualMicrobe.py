import matplotlib.pyplot as plt
import matplotlib.animation as animation
from time import time
from random import randrange
from os import path
import datetime
from math import ceil
import csv

"""Virtual microbe with a fucntional, minimalistic metabolism and cellular functions.
This is a concept project :)"""

MicrobeContents = {}
MicrobeFacets = {"Toxicity" : 100, "Objective" : "Replicate", "TerminalElectronAcceptor" : ["Oxygen"]}

#Data is structured so that an item is the key, and the values are nested lists.
#Values are [[bool, bool][row,column or None][conc over time or time]]. 
# [0] bool is for RT monitoring, and [1] bool is for output monitoring
Data = {"Time" : [[False, True],[None]]}

Settings = {"RowColumnDisplay" : [1,1], "Display" : True, "Output" : False}

StoiMat = {'Substrate': ['LPS', 'PG', 'Glucose', 'NAD', 'NADH', 'ADP', 'ATP', 'nitrate', 'amino acids', 'RNA', 'precursor', 'e-', 'Hydrogen', 'Oxygen', 'nitrite', 'PMF', 'DNA', 'RNAstrand', 'peptides', 'CellWallIntegrity', 'nucleotides', 'Genome'], 
    'DNARepair': {'ADP': 4, 'ATP': -4, 'precursor': -4, 'Genome': 1}, 
    'WallSynthesis': {'LPS': -1, 'PG': -1, 'ADP': 2, 'ATP': -2, 'CellWallIntegrity': 1},
 'WallDecay': {'LPS': 1, 'PG': 1, 'CellWallIntegrity': -1},
  'AAsynthesis': {'ADP': 1, 'ATP': -1, 'amino acids': 1, 'precursor': -1},
   'RnaSynthesis': {'ADP': 1, 'ATP': -1, 'RNA': 1, 'precursor': -1},
    'PrecursorSynthesis': {'Glucose': -1, 'ADP': 1, 'ATP': -1, 'precursor': 1}, 
    'WallComponents': {'LPS': 1, 'PG': 1, 'ADP': 1, 'ATP': -1, 'precursor': -2},
 'ATPhydrolysis': {'ADP': 1, 'ATP': -1},
  'Gluconeogenesis': {'Glucose': 1, 'ADP': 1, 'ATP': -1, 'precursor': -1},
   'ADPsynthesis': {'ADP': 1, 'precursor': -1},
    'NucleotideSynthesis': {'ADP': 1, 'ATP': -1, 'precursor': -1, 'nucleotides': 1},
     'ProteinSynthesis': {'ADP': 20, 'ATP': -20, 'amino acids': -20, 'RNAstrand': -1, 'peptides': 20},
      'RNAsynthesis': {'ADP': 20, 'ATP': -20, 'RNA': -20, 'RNAstrand': 1},
       'DNAsynthesis': {'ADP': 1, 'ATP': -1, 'DNA': 1, 'peptides': -1, 'nucleotides': -1},
        'Replicate': {'LPS': -10, 'PG': -10, 'ADP': 10, 'ATP': -10, 'DNA': -10},
         'ElectronTransportChain': {'NAD': 1, 'NADH': -1, 'PMF': 1},
          'Gly_TCA': {'Glucose': -1, 'NAD': -12, 'NADH': 12, 'Hydrogen': -12},
           'OxidativePhosphorylation': {'ADP': -1, 'ATP': 1, 'Hydrogen': 2, 'PMF': -1}} #Stoichiometry Matrix


def React(Reaction): #this function performs the conversion of metabolites for the reactions
    print(Reaction)
    for substrate, flux in StoiMat[Reaction].items(): 
        MicrobeContents[substrate] += flux

class Transport:
     
    def In_Out(self, molecule):
        MicrobeContents["ex_" + molecule] += 1
        MicrobeContents[molecule] -= 1

    def Out_In(self, molecule):
        MicrobeContents[molecule] += 1
        MicrobeContents["ex_" + molecule] -= 1
        
class DecisionsUnit:

    def __init__(self):

        self.Objectives = []
        self.Stressers = {
                        "Hydrogen" : ">= 100",
                        "Genome" : "< 100",
                        "e-" : "> 0",
                        "CellWallIntegrity" : "< 100"}
        self.Reaction_queue = []

    def __call__(self):
        self.Set_objectives()
        self.Analyse_metabolome()

    def CellStress(self): #This will need modification later

        for metabolite, value in MicrobeContents.items():
            if value > int(MicrobeFacets["Toxicity"]):
                self.Objectives.append(metabolite)
        
        for stresser, condition in self.Stressers.items():
            if eval("MicrobeContents[stresser]" + condition):
                self.Objectives.append(stresser)
            else:
                self.Objectives.append(stresser) #tester
                            
    def Set_objectives(self): #checks physiological state and sets objectives accordingly
        #this could be placed elsewhere
        self.CellStress()
        self.Objectives.append(MicrobeFacets["Objective"])

    def Identify_actions(self, objectif):

        metabolic_requirement = {}
        #this checks if sufficient metabolites are present or if any will build up to toxicity
        for metabolite, flux in objectif.items():
            if MicrobeContents[metabolite] + flux <= 0:
                metabolic_requirement[metabolite] = MicrobeContents[metabolite] + flux
            elif MicrobeContents[metabolite] + flux >= int(MicrobeFacets["Toxicity"]):
                metabolic_requirement[metabolite] = (MicrobeContents[metabolite] + flux) - int(MicrobeFacets["Toxicity"])
        return metabolic_requirement #returns a dictionary of metabolites and their actions required...
    
    def Identify_reactions(self, metabolite, deficit=-10): #identifies which reactions produce/consume a metabolite
        Possibilities = []

        for reaction in StoiMat.keys():
            if reaction != "Substrate" and metabolite in StoiMat[reaction]:
                if deficit > 0 and StoiMat[reaction][metabolite] < 0:
                    Possibilities.append(reaction) 
                elif deficit < 0 and StoiMat[reaction][metabolite] > 0:
                    Possibilities.append(reaction)
                
        return Possibilities #returns a list of reactions

    def Analyse_metabolome(self):
        print(self.Objectives)

        for objective in self.Objectives:
            if objective in StoiMat.keys():
                print(self.Identify_actions(StoiMat[objective]))
                for metabolite, deficit in self.Identify_actions(StoiMat[objective]):
                    print(self.Identify_reactions(metabolite, deficit))
            else:
                #print(self.Identify_reactions(objective))
                for reaction in self.Identify_reactions(objective):
                    print(self.Identify_actions(StoiMat[reaction]))

class CoreFunctions:
    
    def __str__(self):
        """This is a core function class for the programme. 
        It retrieves essential information and checks what reactions occur and helps to ensures the system works"""

    def __call__(self):
        print("Initiating Virtual Microbe's core functions... \n")
        #self.RetrieveSM()
        # self.GetEquations()
        print("Stoichiometry matrix uploaded and checked\n")
        print("Reading configuration file...")
        self.GetConfigs()
        print("Done\n")
        print("Checking configurations...")
        # self.Check_metabolites_present()
        print("Done\n")
        print("Virtual Microbe's core functions have launched successfully\n")

    def RetrieveSM(self): #SM = stoichiometry matrix
        if path.exists("./VirtualMicrobeStoichiometryMatrix.csv"):
            csv_file = open("./VirtualMicrobeStoichiometryMatrix.csv", encoding="utf-8")
            SM = csv.reader(csv_file, delimiter=',')
        else:
            raise SystemExit("""Error: Could not detect Stoichiometry Matrix file.
            Please ensure file is in same directory and named correctly.
            (Requires CSV file: VirtualMicrobeStoichiometryMatrix.csv)""")
        
        #extracts data from csv file
        for row in SM:
            if row[1] == "" or row[0] == "~":
                pass
            elif row[0].replace("\ufeff", "") == "":
                StoiMat["Substrate"] = [flux for flux in row[1:] if flux != ""]
            else:
                StoiMat[row[0]] = {}
                [StoiMat[row[0]].update({StoiMat["Substrate"][int(order)] : int(row[int(order)+1])}) for order in range(len(row[1:])) if row[int(order)+1] != "0" and row[int(order)+1] != ""]
        csv_file.close()
                    
    def GetEquations(self):
        #prints the stoichiometry matrix in a human readable form as balanced equations
        print("Reactions from the stoichiometry matrix:")

        def Balance(reaction):
            balance = 0
            for flux in StoiMat[reaction].values():
                balance += flux
            return balance

        for reaction in StoiMat.keys():

            equation = "->"
            
            if reaction == "Substrate":
                pass
            
            else:
                for substrate, flux in StoiMat[reaction].items():
                    if flux < 0:
                        equation = str(flux).replace("-", "") + " " + substrate + " " + equation if equation.startswith("->") else str(flux).replace("-", "") + " " + substrate + " + " + equation
                        
                    elif flux > 0:
                        equation = equation + " " + str(flux) + " " + substrate if equation.endswith("->") else equation + " + " + str(flux) + " " + substrate
    

                print(reaction + " is ", "BALANCED\n" if Balance(reaction) == 0 else "UNBALANCED\n", equation + "\n") #checks whether equations are balanced

    def GetMetabolites(self):
        print("List of metabolites: " + str([item for item in StoiMat["Substrate"]])[1:-1].replace("'", ""))

    def AxesGenerator(self): #this function generates the row, columns for each axis in the plot
        for row in range(3):
            for col in range(2):
                yield [row, col]

    def GetConfigs(self):
        if path.exists("./Virtual Microbe Configuration File.txt"):
            Config = open("./Virtual Microbe Configuration File.txt")

        else:
            raise SystemExit("""Error: Could not detect configuration file.
            Please ensure you have the configuration file in current directory.
            (Requires text file: Virtual Microbe Configuration File.txt)""")

        for row in Config:

            if row.startswith("@"):
                self.ConfigSettings(row)

            elif row.startswith("!"):
                self.Facets_Physiology(row)

            elif row.startswith('"'):
                item, value = row.split(",")
                MicrobeContents[item[1:-1]] = int(value.strip())

        Config.close()
 
    def ConfigSettings(self, row):

        GenAxes = self.AxesGenerator()

        if row.startswith("@MonitorRT:"):

            if len(row[11:].strip()) == 0: #if nothing has been given...
                Settings["Display"] = False
            
            elif len(row[11:].split(",")) > 6:
                raise Exception("Error: Too many items provided to be monitored. Maximum is six.")

            else:
                try:
                    for item in row[11:].split(","):
                        Data[item.strip()] = [[True, False], next(GenAxes)]
                except:
                    Data[row[11:].strip()] = [[True, False], next(GenAxes)]

            #this determines how many rows and columns are required for the RT figure
            #if len(Data.keys()) > 2: #if there is more than one item to monitor, then specifies the columns and rows
            Settings["RowColumnDisplay"] = [ceil((len(Data.keys())-1)/2), 2]


        elif row.startswith("@Monitor"):

            if len(row[9:].strip()) == 0:
                Settings["Output"] = False
            else:
                try:
                    for item in row[9:].split(","):
                        if item.strip() not in Data:
                            Data[item.strip()] = [[False, True],[None]]
                        else:
                            Data[item.strip()][0][1] = True

                except:
                    if row[9:].strip() not in Data:
                        Data[row[9:].strip()] = [[False, True],[None]]
                    else:
                        Data[row[9:].strip()][0][1] = True

                DataHandler.CreateOutputCSV(self)
                Settings["Output"] = True

        elif row.startswith("@Settings") and len(row[10:].strip()) > 0:
            try:
                for item in row[10:].split(","):
                    name, prop = item.split(" ")
                    Settings[name.strip()] = int(prop)
                
            except:
                name, prop = row[11:].split(" ")
                Settings[name.strip()] = int(prop.strip())
   
    def Facets_Physiology(self, row):
        facet, value = row.split(":")
        MicrobeFacets[facet[1:]] = value.strip()

    def Check_metabolites_present(self):
        for metabolite in MicrobeContents.keys():
            if metabolite not in StoiMat["Substrate"] and not metabolite.startswith("ex_"):
                raise Exception(metabolite + " not found in stoichiometry matrix file. Please add to matrix.")
        
        for metabolite in StoiMat["Substrate"]:
            if metabolite not in MicrobeContents:
                raise Exception(metabolite + " not found in configuration file. Please add to configs and specify concentration.")

    def mixup(self):
        for key, value in MicrobeContents.items():
            if 5 < randrange(0, 10):
                MicrobeContents[key] = value + 1
            else:
                MicrobeContents[key] = value - 1

class DataHandler:

    OutputFile = None #the file information is writen to
    OutputCSV = None
    TimeStarted = 0

    def CreateOutputCSV(self):
        self.OutputCSV = open(str("./VirtualMicrobeOutputs/VirtualMicrobe_Output_"+ datetime.datetime.now().strftime("%d_%b_%Y_%H_%M")), 'w', newline="")
        self.OutputFile = csv.writer(self.OutputCSV)
        self.OutputFile.writerow(item for item in Data.keys() if Data[item][0][1])

    def OutputDataToCSV(self):
        self.OutputFile.writerow(Data[item][2][-1:][0] for item in Data.keys() if Data[item][0][1])

    def CloseOutputFile(self):
        self.OutputCSV.close()
    
    def TakeMonitoredData(self): #this function takes the data from the metabolite/property pool and adds them to data to be displayed or outputted

        for MonitoredItem in Data.keys():
            if MonitoredItem == "Time":
                try:
                    Data["Time"][2].append(round(time() - DataHandler.TimeStarted))
                except:
                    Data["Time"].append([round(time() - DataHandler.TimeStarted)])
            else:
                try:
                    #this adds the metabolite conc from the pool and adds to the data dictionary
                    Data[MonitoredItem][2].append(MicrobeContents[MonitoredItem])
                except: 
                    Data[MonitoredItem].append([MicrobeContents[MonitoredItem]])

class DataPlot:

    def __init__(self):
        self.ani = None #animation
        self.fig = None #figure
        self.axs = None
    
    def __call__(self):
        self.CreatePlot()
        self.InitiatePlot()

    def CreatePlot(self):
        self.fig, self.axs = plt.subplots(nrows = Settings["RowColumnDisplay"][0], ncols = Settings["RowColumnDisplay"][1], constrained_layout=True)

    def InitiatePlot(self):
        DataHandler.TimeStarted = time()
        self.ani = animation.FuncAnimation(self.fig, self.PlotMetabolites, interval=5000)
        if Settings["Display"]:
            plt.show()

    def PlotMetabolites(self, frames):

        CoreFunctions.mixup(self)
        
        if Settings["Display"] or Settings["Output"]:
            DataHandler.TakeMonitoredData(self)
    
        for metabolite in Data.keys():
            if Data[metabolite][0][0]: #checks if metabolite is meant to be monitored in RT

                #the data[metabolite][1][0] and [1][1] is the row and column of the subplot
                self.axs[Data[metabolite][1][0], Data[metabolite][1][1]].clear()

                #Figure formatting
                self.axs[Data[metabolite][1][0], Data[metabolite][1][1]].plot(Data["Time"][2][-10:], Data[metabolite][2][-10:])
                self.axs[Data[metabolite][1][0], Data[metabolite][1][1]].set_title(metabolite, fontsize=12)
                self.axs[Data[metabolite][1][0], Data[metabolite][1][1]].set_xlabel("Time (seconds)", fontsize=10)
                self.axs[Data[metabolite][1][0], Data[metabolite][1][1]].set_ylabel("Conc (units)", fontsize=10)

        if Settings["Output"]:
            DataHandler.OutputDataToCSV(self)

if __name__ == "__main__":
    CoreFunc = CoreFunctions()
    CoreFunc()
    print(MicrobeContents)
    Decisions = DecisionsUnit()
    Decisions()


    if Settings["Display"] or Settings["Output"]:
        print("Loading display...\n")
        DP = DataPlot()
        DP()
        
