import numpy as np


# Finds the ic-s in random_ics that correspond to a protein in the sample
# If the sample prot. is not in random_ics, then ic = 0
def find_sample_ICs(species, random_ics):
    random_sample_ics = dict()
    for s in species:
        if s in random_ics.keys():
            random_sample_ics[s] = random_ics[s]
        else:
            random_sample_ics[s] = 0
    # Ezek nem tudom, hogy ide kellenek-e
    random_sample_ics["REF"] = 1.0
    random_sample_ics["CASPASEA"] = 0
    return random_sample_ics


# Generates an initial cond. for all proteins in the network (proteins from min_max)
def generateICs(bounds):
    random_ics = dict()
    for s in bounds:
        if bounds[s][1] == 0:
            random_ics[s] = 0.0
        else:
            lb, ub = bounds[s]
            random_ics[s] = np.random.uniform(lb, ub)          
    # Ezek nem tudom, hogy ide kellenek-e
    random_ics["REF"] = 1.0
    random_ics["CASPASEA"] = 0
    return random_ics


def generateBounds_proteins(df):
    bounds = dict()
    for col in df.iloc[:, 4:]:
        lb, ub = df[col][:2]
        if lb == '' or lb < 0:
            lb = 0
        if ub == '' or ub < 0:
            ub = 0 
        col = col.replace('x_', '')
        bounds[col.upper()] = [lb, ub]
    return bounds


def compileDataRow(variables, dataPoints):
    meas = ""
    for v in variables:
        meas = meas+"<%s>" % v + "{:.4e}" + "</%s>" % v
    start = "<dataPoint>"
    close = "</dataPoint>"
    row = start+meas.format(*dataPoints)+close
    return row
