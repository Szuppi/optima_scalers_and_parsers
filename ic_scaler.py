import pandas as pd
import jinja2
import time
from ic_scaler_parser import find_sample_ICs, generateBounds_proteins, generateICs, compileDataRow
import os


def scaleData(source_data, variables, ics):
    dataDf = source_data.copy()
    for v in variables:
        if 'STD' in v:
            if v[0:-3] in ics.keys():
                dataDf[v] = dataDf[v]*ics[v[0:-3]]
                # Scaling standard devs
        if v in ics.keys():
            dataDf[v] = dataDf[v]*ics[v]
            # Scaling relative measurement values
    return dataDf


def compileDataTable(ics, source_data, inp):
    variables = source_data.columns
    # hacking Stress here
    ics['RAP'] = inp
    dataDf = scaleData(source_data, variables, ics)
    dataPoints = []
    for i, row in dataDf.iterrows():
        dataPoints.append(compileDataRow(variables, row.values))
    return dataPoints


def generateOutput(ics, variables, dataPoints):
    file_loader = jinja2.FileSystemLoader('.')
    env = jinja2.Environment(loader=file_loader)
    template = env.get_template('data_w_std.xml')
    # megszorozza a számolt hibával, a maximum értékét a mérésnek
    output = template.render(ics=ics, variables=variables,
                             dataPoints=dataPoints)
    return output


def generateFileName(file_index, directory, name, maxdigit=4):
    padded_number = str(file_index).zfill(maxdigit)
    file_name = 'Holczer2019_'+'rap'+'_'+name+'_'+padded_number+'.xml'
    path = os.path.join(directory, file_name)
    return path


# Define the function to generate a file with given content
def generate_file(file_index, directory, bounds, source_data):
    variables = source_data.columns
    species = [v.upper() for v in variables if v.upper() != 'TIME' and 'STD' not in v.upper()]
    # species vector (list) will contain only "real" variables from the sample

    random_ics = generateICs(bounds)
    species_ics = find_sample_ICs(species, random_ics)
    name = ['100nm', '2000nm'] # RAP cc.?
    inp = 100*10**(-12)

    dataPoints = compileDataTable(species_ics, source_data, inp)

    vars_to_xml = species
    actual_ics = species_ics.copy()
    for s in species:
        if s in species_ics.keys():
            actual_ics[s] = actual_ics[s]*source_data[s][0] # This is simply dataPoints time = 0 row values

    # actual_ics = sample prot amounts at t=0, vagy legyen benne REF és CASPASEA is?
    # vars_to_xml rap.csv file-bol a valtozok, ez kell legyen?
    output = generateOutput(actual_ics, vars_to_xml, dataPoints)
    filename = generateFileName(file_index, directory, name[0])

    # elmenti generált ic-ket df-be
    with open(filename, 'w') as f:
        f.write(output)
    return filename


# Directory to save files
output_directory = 'holczer2019/rap'
# Create the directory if it does not exist
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Scaling min_max to get mol/cm^3 and generating intervals for the proteins
df = pd.read_csv('min_max_ranges.csv')
df_scaled = df.select_dtypes(include='number') * 1e-12
bounds_proteins = generateBounds_proteins(df_scaled) 

# data file
data = pd.read_csv('rap.csv')

# Changes every col name to upper case
data.columns = [col.upper() for col in data.columns]

# print("originalData: \n", data, "\n")

start = time.time()
for i in range(1, 11):
    file_index = i
    generate_file(file_index, output_directory, bounds_proteins, data)
print("job finished in:", time.time()-start)