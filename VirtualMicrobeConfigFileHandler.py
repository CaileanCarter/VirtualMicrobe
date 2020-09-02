import os
import sys

Metabolites = {}
Extracellular_substrates = {}
Bacterial_properties = {}

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
    elif row.startswith('"'):
        item, value = row.split(",")
        try:
            exec(PlaceHolder + "[item[1:len(item)-1]] = int(value.strip())")
        except:
            exec(PlaceHolder + "= dict(item = int(value))")

Config.close()