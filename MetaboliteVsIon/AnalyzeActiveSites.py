import os
import re
import pandas as pd

from cluster_analysis_3d import CalculateResidueDistanceWithDataframeInput

CWD = os.getcwd()

PATH_TO_DATAFILES = os.path.join(CWD, "datafiles")
PATH_TO_UNIPROT_ENTRIES = os.path.join(PATH_TO_DATAFILES, "uniprot_entries")
PATH_TO_EVAL_FILES = os.path.join(PATH_TO_DATAFILES, "eval_files")
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
    df = pd.DataFrame(data)
    print(df)

    return df

if __name__ == "__main__":
    file_names = [file for file in os.listdir(PATH_TO_UNIPROT_ENTRIES) if os.path.isfile(os.path.join(PATH_TO_UNIPROT_ENTRIES, file))]
    print(file_names)
    for file in file_names:
        important_positions = read_uniprot_file_to_analyze_active_sites(PATH_TO_UNIPROT_ENTRIES, file)
