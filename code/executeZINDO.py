import os
import sys
import numpy as np
import helperFunctions
import subprocess as sp
import multiprocessing as mp
import cPickle as pickle


def createInputFiles(chromophoreList, AAMorphologyDict, parameterDict):
    # Singles first
    for chromophore in chromophoreList:
        # Include the molecule terminating units on the required atoms of the chromophore
        terminatingGroupPositions = terminateMonomers(chromophore, parameterDict, AAMorphologyDict)
        writeOrcaInp(AAMorphologyDict, chromophore.AAIDs, [chromophore.image] * len(chromophore.AAIDs), terminatingGroupPositions, [chromophore.image] * len(terminatingGroupPositions), parameterDict['outputDir'] + '/' + parameterDict['morphology'][:-4] + chromophore.orcaInput)
    print ""
    # Determine how many pairs there are first:
    numberOfPairs = 0
    for chromo in chromophoreList:
        for neighbour in chromo.neighbours:
            if int(neighbour[0]) > chromo.ID:
                numberOfPairs += 1
    numberOfPairs = np.sum([len(chromo.neighbours) for chromo in chromophoreList])
    print "There are", numberOfPairs/2, "total neighbour pairs to consider."  # /2 because the forwards and backwards hops are identical
    # Then consider each chromophore against every other chromophore
    for chromophore1 in chromophoreList:
        neighboursID = [neighbour[0] for neighbour in chromophore1.neighbours]
        neighboursImage = [neighbour[1] for neighbour in chromophore1.neighbours]
        for chromophore2 in chromophoreList:
            # Skip if chromophore2 is not one of chromophore1's neighbours
            # Also skip if chromophore2's ID is < chromophore1's ID to prevent duplicates
            if (chromophore2.ID not in neighboursID) or (chromophore2.ID < chromophore1.ID):
                continue
            # Find the correct relative image for the neighbour chromophore
            chromophore2RelativeImage = neighboursImage[neighboursID.index(chromophore2.ID)]
            chromophore2Transformation = list(np.array(chromophore1.image) - np.array(chromophore2.image) + np.array(chromophore2RelativeImage))
            # Find the dimer AAIDs and relative images for each atom
            AAIDs = chromophore1.AAIDs + chromophore2.AAIDs
            images = [[0, 0, 0] for i in range(len(chromophore1.AAIDs))] + [chromophore2Transformation for i in range(len(chromophore2.AAIDs))]
            # Now add the terminating groups to both chromophores
            terminatingGroupPositions1 = terminateMonomers(chromophore1, parameterDict, AAMorphologyDict)
            terminatingGroupPositions2 = terminateMonomers(chromophore2, parameterDict, AAMorphologyDict)
            # We don't want to add the terminating hydrogens for adjacent monomers, so remove the ones that are within a particular distance
            terminatingGroupPositions1, terminatingGroupPositions2 = removeAdjacentTerminators(terminatingGroupPositions1, terminatingGroupPositions2)
            terminatingGroupImages1 = [[0, 0, 0] for i in range(len(terminatingGroupPositions1))]
            terminatingGroupImages2 = [chromophore2Transformation for i in range(len(terminatingGroupPositions2))]
            # Update the ORCA input name
            inputName = chromophore1.orcaInput.replace('.inp', '-%04d.inp' % (chromophore2.ID)).replace('single', 'pair')
            # Write the dimer input file
            writeOrcaInp(AAMorphologyDict, AAIDs, images, terminatingGroupPositions1 + terminatingGroupPositions2, terminatingGroupImages1 + terminatingGroupImages2, parameterDict['outputDir'] + '/' + parameterDict['morphology'][:-4] + inputName)
    print ""


def removeAdjacentTerminators(group1, group2):
    popList = [[], []]
    for index1, terminatingHydrogen1 in enumerate(group1):
        for index2, terminatingHydrogen2 in enumerate(group2):
            separation = np.linalg.norm(terminatingHydrogen2 - terminatingHydrogen1)
            if separation < 1.2:
                popList[0].append(index1)
                popList[1].append(index2)
    for groupNo, group in enumerate(popList):
        for index in sorted(group, reverse=True):
            [group1, group2][groupNo].pop(index)
    return group1, group2


def writeOrcaInp(AAMorphologyDict, AAIDs, images, terminatingGroupPosns, terminatingGroupImages, inputName):
    linesToWrite = []
    allAtomTypes = []
    allPositions = []
    # Format the atom positions ready for ORCA
    for index, atomID in enumerate(AAIDs):
        # Cut the integer bit off the atomType
        allAtomTypes.append(''.join([i for i in AAMorphologyDict['type'][atomID] if not i.isdigit()]))
        # Add in the correct periodic images to the position
        allPositions.append(AAMorphologyDict['unwrapped_position'][atomID] + np.array([(images[index][i] * [AAMorphologyDict['lx'], AAMorphologyDict['ly'], AAMorphologyDict['lz']][i]) for i in range(3)]))
    # Now add in the terminating Hydrogens
    for index, position in enumerate(terminatingGroupPosns):
        # Cut the integer bit off the atomType
        allAtomTypes.append('H')
        # Add in the correct periodic images to the position
        allPositions.append(position + np.array([(terminatingGroupImages[index][i] * [AAMorphologyDict['lx'], AAMorphologyDict['ly'], AAMorphologyDict['lz']][i]) for i in range(3)]))
    # Now geometrically centralize all of the atoms that are to be included in this input file to make it easier on ORCA
    centralPosition = np.array([np.average(np.array(allPositions)[:, 0]), np.average(np.array(allPositions)[:, 1]), np.average(np.array(allPositions)[:, 2])])
    # Create the lines to be written in the input file
    for index, position in enumerate(allPositions):
        linesToWrite.append(" %s  %.5f  %.5f  %.5f\n" % (allAtomTypes[index], position[0] - centralPosition[0], position[1] - centralPosition[1], position[2] - centralPosition[2]))
    # Load the ORCA input template
    with open(os.getcwd() + '/templates/template.inp', 'r') as templateFile:
        inpFileLines = templateFile.readlines()
    # Insert the linesToWrite
    inpFileLines[-1:-1] = linesToWrite
    # Write the ORCA input file
    with open(inputName, 'w+') as orcaFile:
        orcaFile.writelines(inpFileLines)
    print "\rOrca Input File written as", inputName[helperFunctions.findIndex(inputName, '/')[-1] + 1:],


def terminateMonomers(chromophore, parameterDict, AAMorphologyDict):
    # Get the connections to the terminating groups from the parXX.py
    terminatingBonds = [bond for bond in parameterDict['moleculeTerminatingConnections']]
    # Remove any termination connections that already exist (i.e. terminating unit at the end of the molecule)
    popList = []
    for bondNo, bond in enumerate(terminatingBonds):
        if bond[0] in np.array(np.array(chromophore.bonds)[:, 0]):
            popList.append(bondNo)
    for index in sorted(popList, reverse=True):
        terminatingBonds.pop(index)
    # Because we didn't reorder the AAID list at any point, the integer in terminatingBond[1] should correspond
    # to the correct AAID in chromo.AAIDs
    AAIDsToAttachTo = [chromophore.AAIDs[index] for index in map(int, list(np.array(terminatingBonds)[:, 1]))]
    # Now work out the positions of any bonded atoms for each of these terminating atoms to work out where we
    # should put the hydrogen
    newHydrogenPositions = []
    for terminatingAtomID in AAIDsToAttachTo:
        # To do this, find the relative positions of the bonded atoms...
        thisAtomPosition = AAMorphologyDict['unwrapped_position'][terminatingAtomID]
        averagePositionOfBondedAtoms = np.array([0.0, 0.0, 0.0])
        for bond in chromophore.bonds:
            if (bond[1] == terminatingAtomID):
                averagePositionOfBondedAtoms += np.array(AAMorphologyDict['unwrapped_position'][bond[2]]) - np.array(thisAtomPosition)
            elif (bond[2] == terminatingAtomID):
                averagePositionOfBondedAtoms += np.array(AAMorphologyDict['unwrapped_position'][bond[1]]) - np.array(thisAtomPosition)
        # ... then reverse that vector and whack on a hydrogen 1.06 A away in that direction
        newHydrogenPositions.append(thisAtomPosition + (- averagePositionOfBondedAtoms / np.linalg.norm(averagePositionOfBondedAtoms) * 1.06))
    # Return terminatingGroups (positions of those hydrogens to be added to the ORCA input)
    return newHydrogenPositions


def getORCAJobs(inputDir, procIDs):
    # First delete any previous log files as we're about to start again with the ZINDO/S calculations
    try:
        os.unlink(inputDir.replace('/inputORCA', '/*.log'))
    except OSError:
        pass
    # Obtain a list of files to run
    singleORCAFileList = os.listdir(inputDir + '/single')
    pairORCAFileList = os.listdir(inputDir + '/pair')
    ORCAFilesToRun = []
    for fileName in singleORCAFileList:
        if fileName[-4:] == '.inp':
            ORCAFilesToRun.append(inputDir + '/single/' + fileName)
    for fileName in pairORCAFileList:
        if fileName[-4:] == '.inp':
            ORCAFilesToRun.append(inputDir + '/pair/' + fileName)
    ORCAFilesToRun.sort()
    # Do not run any jobs that have already have an output file (and so have at least started to run if not finished)
    popList = []
    for jobNo, job in enumerate(ORCAFilesToRun):
        try:
            with open(job.replace('inputORCA', 'outputORCA').replace('.inp', '.out'), 'r'):
                popList.append(jobNo)
        except IOError:
            pass
    popList.sort(reverse=True)
    for popIndex in popList:
        ORCAFilesToRun.pop(popIndex)
    if len(ORCAFilesToRun) == 0:
        return []
    # Create a jobslist for each procID
    jobsList = [ORCAFilesToRun[i:i + (int(np.ceil(len(ORCAFilesToRun) / len(procIDs)))) + 1] for i in xrange(0, len(ORCAFilesToRun), int(np.ceil(len(ORCAFilesToRun) / float(len(procIDs)))))]
    return jobsList


def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList):
    createInputFiles(chromophoreList, AAMorphologyDict, parameterDict)
    inputDir = parameterDict['outputDir'] + '/' + parameterDict['morphology'][:-4] + '/chromophores/inputORCA'
    procIDs = parameterDict['procIDs']
    jobsList = getORCAJobs(inputDir, procIDs)
    numberOfInputs = sum([len(ORCAFilesToRun) for ORCAFilesToRun in jobsList])
    print "Found", numberOfInputs, "ORCA files to run."
    if numberOfInputs > 0:
        # Create pickle file containing the jobs sorted by ProcID to be picked up by singleCoreRunORCA.py
        pickleName = inputDir.replace('inputORCA', 'ORCAJobs.pickle')
        with open(pickleName, 'w+') as pickleFile:
            pickle.dump(jobsList, pickleFile)
        print "ORCA jobs list written to", pickleName
        if len(jobsList) <= len(procIDs):
            procIDs = procIDs[:len(jobsList)]
        runningJobs = []
        # Open the required processes to execute the ORCA jobs
        for CPURank, jobs in enumerate(jobsList):
            print 'python ' + os.getcwd() + '/code/singleCoreRunORCA.py ' + parameterDict['outputDir'] + '/' + parameterDict['morphology'][:-4] + ' ' + str(CPURank) + ' &'
            runningJobs.append(sp.Popen(['python', str(os.getcwd()) + '/code/singleCoreRunORCA.py', parameterDict['outputDir'] + '/' + parameterDict['morphology'][:-4], str(CPURank)]))
        # Wait for all jobs to complete
        [p.wait() for p in runningJobs]
        # Delete the job pickle
        os.system('rm ' + pickleName)
    return AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList


if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print "Please specify the pickle file to load to continue the pipeline from this point."
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList = helperFunctions.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList)