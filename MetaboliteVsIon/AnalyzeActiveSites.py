import os
import re
import pandas as pd

CWD = os.getcwd()
PARENT_DIR = os.path.dirname(CWD)

PATH_TO_DATAFILES = os.path.join(PARENT_DIR, "datafiles")
PATH_TO_UNIPROT_ENTRIES = os.path.join(PATH_TO_DATAFILES, "uniprot_entries")
PATH_TO_EVAL_FILES = os.path.join(PATH_TO_DATAFILES, "eval_files")
PATH_TO_MSA_FILE = os.path.join(PATH_TO_DATAFILES, "muscle_data", "alignment.aln")

def convert_uniprot_data_to_position(entry):
    value_str = entry[2]
    # Check if the value string contains periods
    if re.search(r'\.+', value_str):
        # Split the string using regular expression to handle any number of periods
        start_str, end_str = re.split(r'\.+', value_str)

        # Convert the start and end strings to integers
        start = int(start_str)
        end = int(end_str)

        # Generate the range of values
        return list(range(start, end + 1))
    else:
        # Convert the string to an integer
        return int(value_str)


def read_uniprot_file_to_analyze_active_sites(directory, filename):
    data = []
    act_str = 'ACT_SITE'
    bind_str = 'BINDING'
    for file in os.listdir(directory):
        if file == filename:
            filepath = os.path.join(directory, filename)
            with open(filepath, 'r') as inF:
                lines = inF.readlines()
                for i in range(len(lines) - 1):  # Iterate through lines except the last one
                    line = lines[i].strip()
                    next_line = lines[i + 1].strip()

                    if line.startswith('AC'):
                        uniProt_ID = line.split()[1].replace(';', '')

                    if line.startswith('FT'):
                        entry = line.split()               
                        # Process the current FT line (active site or binding site)
                        if entry[1] == act_str:
                            position = convert_uniprot_data_to_position(entry)
                            data.append({
                                "UniProt_ID": uniProt_ID,
                                "Type": act_str,
                                "Position": position, #if there are two numbers like [450, 455] then the active sites goes from 450 to 455
                                "Description": re.search(r'(?<==")([^"]+)', ' '.join(next_line.split()[1:])).group(1)
                            })                         
                        
                        if entry[1] == bind_str:
                            position = convert_uniprot_data_to_position(entry)
                            data.append({
                                "UniProt_ID": uniProt_ID,
                                "Type": bind_str,
                                "Position": position,
                                "Description": re.search(r'="([^"]+)"', next_line.split()[1]).group(1)
                            })
    activeSitesDF = pd.DataFrame(data)
    return activeSitesDF

def FindActiveSitesInMSA(activeSitesDF, MSAFile):
    print(activeSitesDF)
    with open(MSAFile, 'r') as inF:
        lines = inF.readlines()
        for i in range(len(lines)):
            line = lines[i]
            print(line)

                    # new idea, find the binding site and then check all the values above and below it and see if they are the same. 
                    # or just go to the bottom line and see if that is a star which is much faster
                    #add a conservation score in the DF database wich will have . : * or - Use this to see if any
                    #binding sites don't have a * and then see what it is.

if __name__ == "__main__":
    file_names = [file for file in os.listdir(PATH_TO_UNIPROT_ENTRIES) if os.path.isfile(os.path.join(PATH_TO_UNIPROT_ENTRIES, file))]
    # print(file_names)
    dataFrames = []
    for file in file_names:
        dataFrames.append(read_uniprot_file_to_analyze_active_sites(PATH_TO_UNIPROT_ENTRIES, file))
    final_dataframe = pd.concat(dataFrames, ignore_index=True)


    FindActiveSitesInMSA(read_uniprot_file_to_analyze_active_sites(PATH_TO_UNIPROT_ENTRIES, "A0A0H4VJ04_9SPHN.txt"), PATH_TO_MSA_FILE)

    # print(final_dataframe)
