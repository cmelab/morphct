# ---=== Directory and File Structure ===---
inputDir = '/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/inputCGMorphs'
outputDir = '/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/outputFiles'


# ---=== Input Morphology Details ===---
morphology = 'p1-L15-f0.0-P0.1-T1.5-e0.5.xml'
inputSigma = 3.0
overwriteCurrentData = False

# ---=== Execution Modules ===---

executeFinegraining = False
executeMolecularDynamics = False
executeExtractMolecules = False
executeObtainChromophores = False
executeZINDO = False
executeCalculateTransferIntegrals = False
executeCalculateMobility = True

# ---=== Fine Graining Parameters ===---

CGToTemplateDirs = {\
'A':'/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/templates',\
'B':'/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/templates',\
'C':'/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/templates',\
}
CGToTemplateFiles = {\
'A':'mid3HT.xml',\
'B':'mid3HT.xml',\
'C':'mid3HT.xml',\
}
CGToTemplateAAIDs = {\
'A':[0, 1, 2, 3, 4, 24],\
'B':[5, 6, 7, 18, 19, 20, 21, 22, 23],\
'C':[8, 9, 10, 11, 12, 13, 14, 15, 16, 17],\
}
CGToTemplateBonds = {\
'bondB':['C2-C3', 2, 5],\
'bondC':['C5-C6', 7, 8],\
}
rigidBodySites = {\
'A':[0, 1, 2, 3, 4],\
}
additionalConstraints = [\
['C1-C10', 3, 25],\
['C1-C10-C9', 3, 25, 26],\
['C1-C10-S1', 3, 25, 29],\
['S1-C1-C10', 4, 3, 25],\
['C2-C1-C10', 2, 3, 25],\
['C1-C10-C9-C2', 3, 25, 26, 27],\
['C1-C10-S1-C1', 3, 25, 29, 28],\
['S1-C1-C10-S1', 4, 3, 25, 29],\
['C2-C1-C10-S1', 2, 3, 25, 29],\
]
moleculeTerminatingUnits = [\
['H1',0,0,0],\
]
moleculeTerminatingBonds = [\
]
moleculeTerminatingConnections = [\
['C10-H1', 0],\
['C1-H1', 3],\
]

# ---=== Forcefield Parameters ===---

pairRCut = 10.0
pairDPDGammaVal = 0.0
# --== Lennard-Jones Pair ==--
ljCoeffs = [\
['C1-C1', 0.070, 3.550],\
['C1-C2', 0.070, 3.550],\
['C1-C3', 0.068, 3.525],\
['C1-C4', 0.068, 3.525],\
['C1-C5', 0.068, 3.525],\
['C1-C6', 0.068, 3.525],\
['C1-C7', 0.068, 3.525],\
['C1-C8', 0.068, 3.525],\
['C1-C9', 0.070, 3.550],\
['C1-C10', 0.070, 3.550],\
['C1-H1', 0.046, 2.979],\
['C1-S1', 0.132, 3.550],\
['C2-C2', 0.070, 3.550],\
['C2-C3', 0.068, 3.525],\
['C2-C4', 0.068, 3.525],\
['C2-C5', 0.068, 3.525],\
['C2-C6', 0.068, 3.525],\
['C2-C7', 0.068, 3.525],\
['C2-C8', 0.068, 3.525],\
['C2-C9', 0.070, 3.550],\
['C2-C10', 0.070, 3.550],\
['C2-H1', 0.046, 2.979],\
['C2-S1', 0.132, 3.550],\
['C3-C3', 0.066, 3.500],\
['C3-C4', 0.066, 3.500],\
['C3-C5', 0.066, 3.500],\
['C3-C6', 0.066, 3.500],\
['C3-C7', 0.066, 3.500],\
['C3-C8', 0.066, 3.500],\
['C3-C9', 0.068, 3.525],\
['C3-C10', 0.068, 3.525],\
['C3-H1', 0.044, 2.958],\
['C3-S1', 0.128, 3.525],\
['C4-C4', 0.066, 3.500],\
['C4-C5', 0.066, 3.500],\
['C4-C6', 0.066, 3.500],\
['C4-C7', 0.066, 3.500],\
['C4-C8', 0.066, 3.500],\
['C4-C9', 0.068, 3.525],\
['C4-C10', 0.068, 3.525],\
['C4-H1', 0.044, 2.958],\
['C4-S1', 0.128, 3.525],\
['C5-C5', 0.066, 3.500],\
['C5-C6', 0.066, 3.500],\
['C5-C7', 0.066, 3.500],\
['C5-C8', 0.066, 3.500],\
['C5-C9', 0.068, 3.525],\
['C5-C10', 0.068, 3.525],\
['C5-H1', 0.044, 2.958],\
['C5-S1', 0.128, 3.525],\
['C6-C6', 0.066, 3.500],\
['C6-C7', 0.066, 3.500],\
['C6-C8', 0.066, 3.500],\
['C6-C9', 0.068, 3.525],\
['C6-C10', 0.068, 3.525],\
['C6-H1', 0.044, 2.958],\
['C6-S1', 0.128, 3.525],\
['C7-C7', 0.066, 3.500],\
['C7-C8', 0.066, 3.500],\
['C7-C9', 0.068, 3.525],\
['C7-C10', 0.068, 3.525],\
['C7-H1', 0.044, 2.958],\
['C7-S1', 0.128, 3.525],\
['C8-C8', 0.066, 3.500],\
['C8-C9', 0.068, 3.525],\
['C8-C10', 0.068, 3.525],\
['C8-H1', 0.044, 2.958],\
['C8-S1', 0.128, 3.525],\
['C9-C9', 0.070, 3.550],\
['C9-C10', 0.070, 3.550],\
['C9-H1', 0.046, 2.931],\
['C9-S1', 0.132, 3.550],\
['C10-C10', 0.070, 3.550],\
['C10-H1', 0.046, 2.979],\
['C10-S1', 0.132, 3.550],\
['H1-H1', 0.030, 2.500],\
['H1-S1', 0.087, 2.979],\
['S1-S1', 0.250, 3.550],\
]
# --== Disipative Particle Dynamics Pair ==--
dpdCoeffs = [\
['C1-C1', 100.0, 2.96],\
['C1-C2', 100.0, 2.96],\
['C1-C3', 100.0, 2.96],\
['C1-C4', 100.0, 2.96],\
['C1-C5', 100.0, 2.96],\
['C1-C6', 100.0, 2.96],\
['C1-C7', 100.0, 2.96],\
['C1-C8', 100.0, 2.96],\
['C1-C9', 100.0, 2.96],\
['C1-C10', 100.0, 2.96],\
['C1-H1', 100.0, 2.96],\
['C1-S1', 100.0, 2.96],\
['C2-C2', 100.0, 2.96],\
['C2-C3', 100.0, 2.96],\
['C2-C4', 100.0, 2.96],\
['C2-C5', 100.0, 2.96],\
['C2-C6', 100.0, 2.96],\
['C2-C7', 100.0, 2.96],\
['C2-C8', 100.0, 2.96],\
['C2-C9', 100.0, 2.96],\
['C2-C10', 100.0, 2.96],\
['C2-H1', 100.0, 2.96],\
['C2-S1', 100.0, 2.96],\
['C3-C3', 100.0, 2.96],\
['C3-C4', 100.0, 2.96],\
['C3-C5', 100.0, 2.96],\
['C3-C6', 100.0, 2.96],\
['C3-C7', 100.0, 2.96],\
['C3-C8', 100.0, 2.96],\
['C3-C9', 100.0, 2.96],\
['C3-C10', 100.0, 2.96],\
['C3-H1', 100.0, 2.96],\
['C3-S1', 100.0, 2.96],\
['C4-C4', 100.0, 2.96],\
['C4-C5', 100.0, 2.96],\
['C4-C6', 100.0, 2.96],\
['C4-C7', 100.0, 2.96],\
['C4-C8', 100.0, 2.96],\
['C4-C9', 100.0, 2.96],\
['C4-C10', 100.0, 2.96],\
['C4-H1', 100.0, 2.96],\
['C4-S1', 100.0, 2.96],\
['C5-C5', 100.0, 2.96],\
['C5-C6', 100.0, 2.96],\
['C5-C7', 100.0, 2.96],\
['C5-C8', 100.0, 2.96],\
['C5-C9', 100.0, 2.96],\
['C5-C10', 100.0, 2.96],\
['C5-H1', 100.0, 2.96],\
['C5-S1', 100.0, 2.96],\
['C6-C6', 100.0, 2.96],\
['C6-C7', 100.0, 2.96],\
['C6-C8', 100.0, 2.96],\
['C6-C9', 100.0, 2.96],\
['C6-C10', 100.0, 2.96],\
['C6-H1', 100.0, 2.96],\
['C6-S1', 100.0, 2.96],\
['C7-C7', 100.0, 2.96],\
['C7-C8', 100.0, 2.96],\
['C7-C9', 100.0, 2.96],\
['C7-C10', 100.0, 2.96],\
['C7-H1', 100.0, 2.96],\
['C7-S1', 100.0, 2.96],\
['C8-C8', 100.0, 2.96],\
['C8-C9', 100.0, 2.96],\
['C8-C10', 100.0, 2.96],\
['C8-H1', 100.0, 2.96],\
['C8-S1', 100.0, 2.96],\
['C9-C9', 100.0, 2.96],\
['C9-C10', 100.0, 2.96],\
['C9-H1', 100.0, 2.96],\
['C9-S1', 100.0, 2.96],\
['C10-C10', 100.0, 2.96],\
['C10-H1', 100.0, 2.96],\
['C10-S1', 100.0, 2.96],\
['H1-H1', 100.0, 2.96],\
['H1-S1', 100.0, 2.96],\
['S1-S1', 100.0, 2.96],\
]
# --== Bond ==--
bondCoeffs = [\
['C1-C2', 1028.54, 1.37368],\
['C1-S1', 582.50, 1.73373],\
['C10-C9', 1028.54, 1.37368],\
['C10-S1', 582.50, 1.73373],\
['C2-C3', 599.64, 1.50884],\
['C2-C9', 906.2, 1.43277],\
['C3-C4', 536.00, 1.54158],\
['C3-H1', 655.09, 1.09827],\
['C4-C5', 536.00, 1.54158],\
['C4-H1', 680.00, 1.09527],\
['C5-C6', 536.00, 1.54158],\
['C5-H1', 680.00, 1.09527],\
['C6-C7', 536.00, 1.54158],\
['C6-H1', 680.00, 1.09527],\
['C7-C8', 536.00, 1.54158],\
['C7-H1', 680.00, 1.09527],\
['C8-H1', 680.00, 1.09527],\
['C9-H1', 741.26, 1.0822],\
['C1-C10', 784.58, 1.45],\
['C1-H1', 741.26, 1.0822],\
['C10-H1', 741.26, 1.0822],\
]
# --== Angle ==--
angleCoeffs = [\
['C1-C2-C3', 166.32, 2.17388],\
['C2-C1-S1', 86.36, 1.92496],\
['C2-C3-C4', 120.14, 2.01481],\
['C2-C3-H1', 74.06, 1.90571],\
['C2-C9-H1', 35.263, 2.15897],\
['C3-C4-C5', 58.35, 1.96699],\
['C3-C4-H1', 37.5, 1.93208],\
['C4-C3-H1', 37.5, 1.93208],\
['C4-C5-C6', 58.35, 1.96699],\
['C4-C5-H1', 37.5, 1.93208],\
['C5-C4-H1', 37.5, 1.93208],\
['C5-C6-C7', 58.35, 1.96699],\
['C5-C6-H1', 37.5, 1.93208],\
['C6-C5-H1', 37.5, 1.93208],\
['C6-C7-C8', 58.35, 1.96699],\
['C6-C7-H1', 37.5, 1.93208],\
['C7-C6-H1', 37.5, 1.93208],\
['C7-C8-H1', 37.5, 1.93208],\
['C8-C7-H1', 37.5, 1.93208],\
['C9-C2-C1', 39.582, 1.97784],\
['C9-C2-C3', 166.545, 2.15335],\
['C9-C10-S1', 86.36, 1.92496],\
['C10-C9-C2', 39.582, 1.97784],\
['C10-C9-H1', 35.263, 2.14639],\
['C10-S1-C1', 86.36, 1.61921],\
['H1-C3-H1', 33.0, 1.88146],\
['H1-C4-H1', 33.0, 1.88146],\
['H1-C5-H1', 33.0, 1.88146],\
['H1-C6-H1', 33.0, 1.88146],\
['H1-C7-H1', 33.0, 1.88146],\
['H1-C8-H1', 33.0, 1.88146],\
['C1-C10-C9', 54.694, 2.27137],\
['C1-C10-S1', 41.74, 2.08687],\
['S1-C1-C10', 41.74, 2.08687],\
['C2-C1-C10', 54.694, 2.27137],\
]
# --== Dihedral ==--
dihedralCoeffs = [\
['C9-C10-S1-C1', 126.32, -109.81, -19.738, -25.303, 28.53],\
['C10-C9-C2-C1', 126.32, -109.81, -19.738, -25.303, 28.53],\
['C10-C9-C2-C3', 117.65, 238.26, 205.96, 112.81, 27.467],\
['C10-S1-C1-C2', 126.32, -109.81, -19.738, -25.303, 28.53],\
['C2-C3-C4-C5', 2.4469, -6.3946, 10.747, 30.695, 11.139],\
['C2-C9-C10-S1', 126.32, -109.81, -19.738, -25.303, 28.53],\
['C3-C4-C5-C6', 1.9475, -3.7121, 1.388, 8.6305, 1.6008],\
['C9-C2-C3-C4', 0.3175, 1.127, 14.143, -22.297, 6.7188],\
['C4-C5-C6-C7', 1.8922, -3.49094, 1.4665, 7.1418, 0.2859],\
['C5-C6-C7-C8', 1.9788, -3.8476, 1.1614, 7.419, 0.4146],\
['C1-C10-C9-C2', 75.595, 116.0, 42.679, -1.528, -3.8137],\
['C1-C10-S1-C1', 158.7, 418.34, 521.33, 376.73, 115.12],\
['S1-C1-C10-S1', 2.9533, 0.1571, -4.2326, 0.3979, 1.8855],\
['C2-C1-C10-S1', 2.9533, -0.1571, -4.2326, -0.3979, 1.8855],\
]
# --== Improper ==--
improperCoeffs = [\
]

# ---=== Molecular Dynamics Phase Parameters ===---
numberOfPhases = 8
temperatures = [1.0]
taus = [1.0]
pairTypes = ['none', 'dpd', 'lj', 'lj', 'lj', 'lj', 'lj', 'lj']
bondTypes = ['harmonic']
angleTypes = ['harmonic']
dihedralTypes = ['table']
integrationTargets = ['all']
timesteps = [1E-3, 1E-3, 1E-9, 5E-9, 1E-8, 1E-7, 1E-6, 1E-5]
durations = [1E5, 1E4, 1E3, 1E3, 1E3, 1E4, 1E5, 1E5]
terminationConditions = ['KEmin', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt']
groupAnchorings = ['all', 'all', 'all', 'all', 'all', 'all', 'all', 'none']

# ---=== Chromophore Parameters ===---

CGSiteSpecies = {\
'A':'Donor',\
'B':'None',\
'C':'None',\
}
maximumHopDistance = 10.0
removeORCAInputs = False
removeORCAOutputs = False
chromophoreLength = 3

# ---=== Chromophore Energy Scaling Parameters ===---

literatureHOMO = -5.0
literatureLUMO = None
targetDoSSTDHOMO = 0.1
targetDoSSTDLUMO = None


# ---=== Kinetic Monte Carlo Parameters ===---

systemTemperature = 290
numberOfCarriersPerSimulationTime = 10000
hopLimit = 0
simulationTimes = [1.00e-11, 2.15e-11, 4.64e-11, 1.00e-10, 2.15e-10, 4.64e-11, 1.00e-9]
recordCarrierHistory = True
reorganisationEnergy = 0.3063
combineKMCResults = True

# ---=== Begin run ===---

parameterFile = __file__

if __name__ == "__main__":
    import runMorphCT
    import sys

    sys.path.append('./code')

    import helperFunctions

    procIDs = helperFunctions.getCPUCores()
    parameterNames = [i for i in dir() if (not i.startswith('__')) and (i not in ['runMorphCT', 'helperFunctions', 'sys'])]
    parameters = {}
    for name in parameterNames:
        parameters[name] = locals()[name]
    runMorphCT.simulation(**parameters) # Execute MorphCT using these simulation parameters
