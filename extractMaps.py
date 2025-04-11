import sys;
import re;
import json;

# Input should be one file containing contact information.
fileName = sys.argv[1];
contactFiles = open(fileName, "r");
# Setting a constant to define the threshold for what area
# is considered a contact.
CONTACT_THRESHOLD = 20.0

# Making a function to write to file.
def writeCM(BB2BB_map, SC2SC_map, SC2BB_map, sequence):
    # Building a header to write.
    header = f"# {sequence}";
    # Writing to files.
    outputName = f"{pdbCode}.CM.txt";
    outputFile = open(outputName, "w");
    outputFile.write(f"{header}\n");
    outputFile.write(f"{json.dumps(BB2BB_map)}\n{json.dumps(SC2SC_map)}\n{json.dumps(SC2BB_map)}");
    outputFile.close();

for fileName in contactFiles:
    fileName = fileName.rstrip();
    # Finding the PDB code for the file.
    pdbCode = re.search('(?<=Results\/).*(?=.pdb)', fileName).group(0);
    sequence = [];
    BB2BB_map = {};
    SC2SC_map = {};
    SC2BB_map = {};
    contactFile = open(fileName, "r");
    # Looping through the file to build a contact map.
    for line in contactFile:
        # Skipping non-informative lines.
        if "SAS" in line or "#" in line or line == "\n" or line == "":
            continue;
        # Collecting relevant information about each contact.
        lineList = line.split();
        try:
            # Information about the source particle.
            parName = lineList[1];  # Name of the particle (V0, for example)
            parAA = parName[0];  # One-letter AA code for the particle
            parType = "bb" if parName[-1] == "0" else "sc";  # 0 indicates backbone, 1 indicates sidechain
            resNum = lineList[2];  # Residue number of the particle (Like 1)
            # Information about the contact particle.
            contactName = lineList[5];  # Name of the contact particle (L0, for example)
            contactAA = contactName[0];  # One-letter AA code for the contact particle
            contactType = "bb" if contactName[-1] == "0" else "sc";  # 0 indicates backbone, 1 indicates sidechain
            contactResNum = lineList[6];  # Residue of the contact particle (Like 2)
            # General contact information.
            contactArea = lineList[8];  # Area of the contact (>25 indicates a contact)
            sectorInfo = lineList[-1];  # Sector of the contact particle relative to the current particle
        except:
            continue;
        # Getting the sector information.
        if "=" in sectorInfo:
            sectorNum = int(re.search('(?<=\=).+', lineList[-1]).group(0).rstrip().lstrip());
        else:
            sectorNum = int(sectorInfo.rstrip().lstrip());
        # Tracking any s = 0 values. If found, just move on.
        if sectorNum == 0:
            continue;
        # Tracking contact information.
        if float(contactArea) < CONTACT_THRESHOLD:
            continue;
        if parType == "bb":
            if contactType == "bb":
                if resNum not in BB2BB_map.keys():
                    BB2BB_map[resNum] = [];
                BB2BB_map[resNum].append(contactResNum);
        else:
            if contactType == "bb":
                if resNum not in SC2BB_map.keys():
                    SC2BB_map[resNum] = [];
                SC2BB_map[resNum].append(contactResNum);
            else:
                if resNum not in SC2SC_map.keys():
                    SC2SC_map[resNum] = [];
                SC2SC_map[resNum].append(contactResNum);
        # Adding to sequence if unique.
        parRes = f"{parAA}{resNum}"
        if parRes not in sequence:
            sequence.append(parRes);
    sequenceStr = "".join([AA[0] for AA in sequence]);
    writeCM(BB2BB_map, SC2SC_map, SC2BB_map, sequenceStr);
    contactFile.close();
contactFiles.close();