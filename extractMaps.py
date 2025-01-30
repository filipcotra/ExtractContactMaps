import sys
import re
import pandas as pd

# Input should be one file containing contact information.
fileName = sys.argv[1]
contactFile = open(fileName, "r")
# Finding the PDB code for the file.
pdbCode = re.search('(?<=Results\/).*(?=.pdb)', fileName).group(0)


# Setting a constant to define the threshold for what area
# is considered a contact.
CONTACT_THRESHOLD = 25.0
# Setting a constant to define the pandas dataframe columns,
# as they will be identical for each.
COL_NAMES = ["PAR", "CONTACT", "AREA"]

# Creating dictionaries to track contact maps, which will be
# stored in pandas dataframes. Key values will be strings
# representing the type of contact map: BB2BB, SC2SC, and SC2BB.
contactMaps = {"bb2bb": pd.DataFrame(columns = COL_NAMES), "sc2sc": pd.DataFrame(columns = COL_NAMES),
               "sc2bb": pd.DataFrame(columns = COL_NAMES)}
# Tracking the sequence.
sequence = []

# Making a function to write to file.
def writeDf(df, dftype):
    # Building a header to write.
    header = f"#CMVIEW GRAPH FILE ver: 1.0\n#SEQUENCE: {sequence}\n#PDB: {pdbCode}\n" +\
             f"#PDB CHAIN CODE: CUSTOM\n#CHAIN: CUSTOM\n#MODEL: 1\n#CT: {dftype}\n#CUTOFF: 8.0"
    # Writing to files.
    outputName = f"{pdbCode}.{dftype}.CM.txt"
    outputFile = open(outputName, "w")
    outputFile.write(f"{header}\n")
    outputFile.write(df.to_csv(sep = '\t', header = False, index = False))
    outputFile.close()

# Looping through the file to build a contact map.
for line in contactFile:
    # Skipping non-informative lines.
    if "SAS" in line or "#" in line or line == "\n" or line == "":
        continue
    # Collecting relevant information about each contact.
    lineList = line.split()
    try:
        # Information about the source particle.
        parName = lineList[1]  # Name of the particle (V0, for example)
        parAA = parName[0]  # One-letter AA code for the particle
        parType = "bb" if parName[-1] == "0" else "sc"  # 0 indicates backbone, 1 indicates sidechain
        resNum = lineList[2]  # Residue number of the particle (Like 1)
        # Adding to sequence if unique.
        parRes = f"{parAA}{resNum}"
        if parRes not in sequence:
            sequence.append(parRes)
        # Information about the contact particle.
        contactName = lineList[5]  # Name of the contact particle (L0, for example)
        contactAA = contactName[0]  # One-letter AA code for the contact particle
        contactType = "bb" if contactName[-1] == "0" else "sc"  # 0 indicates backbone, 1 indicates sidechain
        contactResNum = lineList[6]  # Residue of the contact particle (Like 2)
        # General contact information.
        contactArea = lineList[8]  # Area of the contact (>25 indicates a contact)
        sectorInfo = lineList[-1]  # Sector of the contact particle relative to the current particle
    except:
        continue
    # Getting the sector information.
    if "=" in sectorInfo:
        sectorNum = int(re.search('(?<=\=).+', lineList[-1]).group(0).rstrip().lstrip())
    else:
        sectorNum = int(sectorInfo.rstrip().lstrip())
    # Tracking any s = 0 values. If found, just move on.
    if sectorNum == 0:
        continue
    # Tracking contact information.
    if float(contactArea) < CONTACT_THRESHOLD:
        continue
    if parType == "bb":
        if contactType == "bb":
            newRow = pd.DataFrame(data = [[resNum, contactResNum, contactArea]], columns = COL_NAMES)
            contactMaps["bb2bb"] = pd.concat([contactMaps["bb2bb"], newRow], ignore_index=True)
    else:
        if contactType == "bb":
            newRow = pd.DataFrame(data = [[resNum, contactResNum, contactArea]], columns = COL_NAMES)
            contactMaps["sc2bb"] = pd.concat([contactMaps["sc2bb"], newRow], ignore_index=True)
        else:
            newRow = pd.DataFrame(data = [[resNum, contactResNum, contactArea]], columns = COL_NAMES)
            contactMaps["sc2sc"] = pd.concat([contactMaps["sc2sc"], newRow], ignore_index=True)
contactFile.close()

# Editing sequence to remove particle type indicators.
i = 0
while i < len(sequence):
    sequence[i] = sequence[i][0]
    i += 1
sequence = "".join(sequence)
# Writing to files.
writeDf(contactMaps["bb2bb"], "bb2bb")
writeDf(contactMaps["sc2sc"], "sc2sc")
writeDf(contactMaps["sc2bb"], "sc2bb")
