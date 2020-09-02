import csv

#This script takes a CSV file containing the stoichiometry matrix of a metabolic network and imports it into the virtual microbe

StoiMat = {}#Stoichiometry Matrix

csv_file = open(".\VirtualMicrobeStoichiometryMatrix.csv", encoding="utf-8")

SM = csv.reader(csv_file, delimiter=',')

#extracts data from csv file
for row in SM:
    if row[1] == "" or row[0] == "~":
        pass
    else:
        StoiMat["Substrate" if row[0].replace("\ufeff", "") == "" else row[0]] = [sub for sub in row[1:] if sub != ""] 
csv_file.close()
#print(StoiMat)

#checks all fluxes are present in the matrix
for reaction, flux in zip(StoiMat.keys(), StoiMat.values()):
    if len(flux) != len(StoiMat["Substrate"]):
        print(reaction + " is missing a flux value, please check your csv file. Ensure zero flux is denoted with 0.")

#prints the stoichiometry matrix in a human readable form as balanced equations
for reaction, FluxList in zip(StoiMat.keys(), StoiMat.values()):
    equation = " -> "
    Balance = 0
    if reaction == "Substrate" or reaction == "Replicate":
        pass
    else:
        for flux in range(len(FluxList)):
            if FluxList[flux].startswith("-"):
                equation = " " + FluxList[flux].replace("-", "") + StoiMat["Substrate"][flux] + equation
                Balance += int(FluxList[flux])
            elif FluxList[flux] != "0":
                equation = equation + " " + FluxList[flux] + StoiMat["Substrate"][flux]
                Balance += int(FluxList[flux])
        print("BALANCED" if Balance == 0 else "UNBALANCED", " " + reaction + ":" + equation) #checks whether equations are balanced


#block to determine whether a metabolite builds up or decomposes (is used or produced by only one reaction)
for element in range(len(StoiMat["Substrate"])):
    elem_balance = 0
    for reac in list(StoiMat)[1:]:
        elem_balance += int(StoiMat[reac][element])
    print(StoiMat["Substrate"][element], " is not balanced by " if elem_balance != 0 else "is balanced ", + elem_balance if elem_balance != 0 else "")
    #print(StoiMat["Substrate"][element] + " is not balanced by " + str(elem_balance) if elem_balance != 0 else "")