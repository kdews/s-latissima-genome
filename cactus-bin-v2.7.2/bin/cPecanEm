#!/usr/bin/env python3

"""Script to train pair-HMMs of cPecan.
"""

import math
import os
import random
import argparse
from toil.common import Toil
from toil.job import Job
from sonLib.bioio import setLoggingFromOptions, logger, system, popenCatch, cigarRead, cigarWrite, nameValue, prettyXml, fastaRead
import numpy
import xml.etree.cElementTree as ET
from itertools import product
from functools import reduce

SYMBOL_NUMBER=4

class Hmm:
    def __init__(self, modelType):
        self.modelType=modelType
        self.stateNumber = { "fiveState":5, "fiveStateAsymmetric":5, "threeState":3, "threeStateAsymmetric":3}[modelType]
        self.transitions = [0.0] * self.stateNumber**2
        self.emissions = [0.0] * (SYMBOL_NUMBER**2 * self.stateNumber)
        self.likelihood = 0.0
        self.runningLikelihoods = []

    def _modelTypeInt(self):
        return { "fiveState":0, "fiveStateAsymmetric":1, "threeState":2, "threeStateAsymmetric":3}[self.modelType]

    def write(self, file):
        f = open(file, 'w')
        f.write(("%s " % self._modelTypeInt()) + " ".join(map(str, self.transitions)) + (" %s\n" % self.likelihood))
        f.write(" ".join(map(str, self.emissions)) + "\n")
        f.close()

    def addExpectationsFile(self, file):
        fH = open(file, 'r')
        l = list(map(float, fH.readline().split()))
        assert len(l) == len(self.transitions)+2
        assert int(l[0]) == self._modelTypeInt()
        self.likelihood += l[-1]
        self.transitions = [sum(x) for x in zip(self.transitions, l[1:-1])]
        assert len(self.transitions) == self.stateNumber**2
        l = list(map(float, fH.readline().split()))
        assert len(l) == len(self.emissions)
        self.emissions = [sum(x) for x in zip(self.emissions, l)]
        assert len(self.emissions) == self.stateNumber * SYMBOL_NUMBER**2
        self.runningLikelihoods = list(map(float, fH.readline().split())) #This allows us to keep track of running likelihoods
        fH.close()
        return self

    @staticmethod
    def loadHmm(file):
        fH = open(file, 'r')
        l = fH.readline().split()
        assert len(l) > 0
        fH.close()
        return Hmm({ 0:"fiveState", 1:"fiveStateAsymmetric", 2:"threeState", 3:"threeStateAsymmetric"}[int(l[0])]).addExpectationsFile(file)

    def normalise(self):
        """Normalises the EM probs.
        """
        for fromState in range(self.stateNumber):
            i = self.stateNumber * fromState
            j = sum(self.transitions[i:i+self.stateNumber])
            for toState in range(self.stateNumber):
                self.transitions[i + toState] = self.transitions[i + toState] / j
        for state in range(self.stateNumber):
            i = state * SYMBOL_NUMBER * SYMBOL_NUMBER
            j = sum(self.emissions[i:i+SYMBOL_NUMBER * SYMBOL_NUMBER])
            for emission in range(SYMBOL_NUMBER * SYMBOL_NUMBER):
                self.emissions[i + emission] = self.emissions[i + emission] / j

    def randomise(self):
        """Randomise the values in the HMM to small values.
        """
        self.transitions = [random.random() for x in range(self.stateNumber*self.stateNumber)]
        self.emissions = [random.random() for x in range(self.stateNumber*SYMBOL_NUMBER*SYMBOL_NUMBER)]
        self.normalise()

    def equalise(self):
        """Initialise the hmm with all equal probabilities.
        """
        self.transitions = [1.0/self.stateNumber] * (self.stateNumber**2)
        self.emissions = [1.0/(SYMBOL_NUMBER*SYMBOL_NUMBER)] * (self.stateNumber * SYMBOL_NUMBER**2)

    def setEmissionsToJukesCantor(self, divergence):
        i = (0.25 + 0.75*math.exp(-4.0*divergence/3.0))/4.0
        j = (0.25 - 0.25*math.exp(-4.0*divergence/3.0))/4.0
        for state in range(self.stateNumber):
            for x in range(SYMBOL_NUMBER):
                for y in range(SYMBOL_NUMBER):
                    self.emissions[(state * SYMBOL_NUMBER**2) + x * SYMBOL_NUMBER + y] = i if x == y else j

    def tieEmissions(self):
        """Sets the emissions to reflect overall divergence, but not to distinguish between different base identity
        """
        for state in range(self.stateNumber):
            a = self.emissions[state*SYMBOL_NUMBER**2:(state+1)*SYMBOL_NUMBER**2]
            identityExpectation = sum([float(a[i]) if (i % SYMBOL_NUMBER) == (i / SYMBOL_NUMBER) else 0.0 for i in range(SYMBOL_NUMBER**2)])
            a = [identityExpectation/SYMBOL_NUMBER if (i % SYMBOL_NUMBER) == (i / SYMBOL_NUMBER) else (1.0 - identityExpectation)/(SYMBOL_NUMBER**2 - SYMBOL_NUMBER) for i in range(SYMBOL_NUMBER**2)]
            assert sum(a) + 0.001 > 1.0 and sum(a) - 0.001 < 1.0
            self.emissions[state*SYMBOL_NUMBER**2:(state+1)*SYMBOL_NUMBER**2] = a
        assert len(self.emissions) == self.stateNumber * SYMBOL_NUMBER**2

def expectationMaximisation(job, sequences, alignments, outputModel, options):
    #Iteratively run cPecanRealign to get expectations and load model.
    if options.inputModel != None: #Read in the model
        job.log("Loading the model from the input file %s" % options.inputModel)
        hmm = Hmm.loadHmm(options.inputModel)
        job.log("Loaded the model, has type %s" % hmm.modelType)
        hmm.normalise()
    else:
        job.log("Making model of type %s" % options.modelType)
        hmm = Hmm(options.modelType)
        if options.randomStart: #Make random parameters
            job.log("Using random starting parameters")
            hmm.randomise()
        else:
            hmm.equalise()
    if options.setJukesCantorStartingEmissions != None:
        hmm.setEmissionsToJukesCantor(float(options.setJukesCantorStartingEmissions))

    #Write out the first version of the output model
    hmm.write(outputModel)

    #Make a set of split alignment files
    alignmentsLength = 0
    splitAlignmentFiles = []
    fH = None
    for cigar in cigarRead(alignments):
        if fH == None:
            splitAlignmentFiles.append(os.path.join(job.getGlobalTempDir(), "alignments_%s.cigar" % len(splitAlignmentFiles)))
            fH = open(splitAlignmentFiles[-1], 'w')
        alignmentsLength += (abs(cigar.start1 - cigar.end1) + abs(cigar.start2 - cigar.end2))/2.0
        cigarWrite(fH, cigar)
        if alignmentsLength > options.maxAlignmentLengthPerJob:
            fH.close()
            fH = None
            splitAlignmentFiles[-1] = (splitAlignmentFiles[-1], alignmentsLength)
            alignmentsLength = 0
    if fH != None:
        fH.close()
        splitAlignmentFiles[-1] = (splitAlignmentFiles[-1], alignmentsLength)

    #Sample the alignment files so that we do EM on no more than options.maxAlignmentLengthToSample bases
    random.shuffle(splitAlignmentFiles) #This ensures we don't just take the first N alignments in the provided alignments file
    sampledSplitAlignmentFiles = []
    totalSampledAlignmentLength = 0.0
    for alignmentsFile, alignmentsLength in splitAlignmentFiles:
        totalSampledAlignmentLength += alignmentsLength
        sampledSplitAlignmentFiles.append(alignmentsFile)
        if totalSampledAlignmentLength >= options.maxAlignmentLengthToSample:
            break
    job.log("We sampled: %s bases of alignment length and %s alignment files, of a possible %s base and %s files" % \
                       (totalSampledAlignmentLength, len(sampledSplitAlignmentFiles), sum([x[1] for x in splitAlignmentFiles]), len(splitAlignmentFiles)))
    splitAlignmentFiles = sampledSplitAlignmentFiles

    #Files to store expectations in
    expectationsFiles = [os.path.join(job.getGlobalTempDir(), "expectation_%i.txt" % i) for i in range(len(splitAlignmentFiles))]
    assert len(splitAlignmentFiles) == len(expectationsFiles)

    job.addFollowOnJobFn(expectationMaximisation2, sequences, splitAlignmentFiles, outputModel, expectationsFiles, 0, [], options)

def expectationMaximisation2(job, sequences, splitAlignments, modelsFile, expectationsFiles, iteration, runningLikelihoods, options):
    if iteration < options.iterations:
        # FIXME RefactoringTool: Line 168: You should use a for loop here
        list(map(lambda x : job.addChildTargetFn(calculateExpectations,
                    args=(sequences, x[0], None if (options.useDefaultModelAsStart and iteration == 0) else modelsFile, x[1], options)),
            list(zip(splitAlignments, expectationsFiles))))
        job.addFollowOnJobFn(calculateMaximisation, sequences, splitAlignments, modelsFile, expectationsFiles, iteration, runningLikelihoods, options)
    else:
        #Write out the likelihoods to the bottom of the file
        fH = open(modelsFile, 'a')
        fH.write("\t".join(map(str, runningLikelihoods)) + "\n")
        fH.close()

def calculateExpectations(job, sequences, alignments, modelsFile, expectationsFile, options):
    #Run cPecanRealign
    system("cat %s | cPecanRealign --logLevel DEBUG %s %s --outputExpectations=%s %s" % (alignments, sequences, nameValue("loadHmm", modelsFile, str), expectationsFile, options.optionsToRealign))

def calculateMaximisation(job, sequences, splitAlignments, modelsFile, expectationsFiles, iteration, runningLikelihoods, options):
    #Load and merge the models files
    if len(expectationsFiles) > 0:
        hmm = Hmm.loadHmm(expectationsFiles[0])
        for expectationsFile in expectationsFiles[1:]:
            hmm.addExpectationsFile(expectationsFile)
        hmm.normalise()
        job.log("On %i iteration got likelihood: %s for model-type: %s, model-file %s" % (iteration, hmm.likelihood, hmm.modelType, modelsFile))
        runningLikelihoods.append(hmm.likelihood)
        job.log("On %i iteration got transitions: %s for model-type: %s, model-file %s" % (iteration, " ".join(map(str, hmm.transitions)), hmm.modelType, modelsFile))
        #If not train emissions then load up the old emissions and replace
        if options.trainEmissions:
            if options.tieEmissions:
                hmm.tieEmissions()
            job.log("On %i iteration got emissions: %s for model-type: %s, model-file %s" % (iteration, " ".join(map(str, hmm.emissions)), hmm.modelType, modelsFile))
        else:
            hmm.emissions = Hmm.loadHmm(modelsFile).emissions
            job.log("On %i using the original emissions" % iteration)

        #Write out
        hmm.write(modelsFile)

    #Generate a new set of alignments, if necessary
    if options.updateTheBand:
        # FIXME: RefactoringTool: Line 206: You should use a for loop here
        list(map(lambda alignments : job.addChildTargetFn(calculateAlignments, args=(sequences, alignments, modelsFile, options)), splitAlignments))

    #Start the next iteration
    job.addFollowOnJobFn(expectationMaximisation2, sequences, splitAlignments, modelsFile, expectationsFiles, iteration+1, runningLikelihoods, options)
    #Call em2

def calculateAlignments(job, sequences, alignments, modelsFile, options):
    temporaryAlignmentFile=os.path.join(job.getLocalTempDir(), "realign.cigar")
    system("cat %s | cPecanRealign --logLevel DEBUG %s --loadHmm=%s %s > %s" % (alignments, sequences, modelsFile, options.optionsToRealign, temporaryAlignmentFile))
    system("mv -f %s %s" % (temporaryAlignmentFile, alignments))

def expectationMaximisationTrials(job, sequences, alignments, outputModel, options):
    if (options.inputModel is not None) or (not options.randomStart): #No multiple iterations
        job.addFollowOnJobFn(expectationMaximisation, sequences, alignments, outputModel, options)
    else:
        job.log("Running %s random restart trials to find best hmm" % options.trials)
        trialModels = [os.path.join(job.getGlobalTempDir(), "trial_%s.hmm" % i) for i in range(options.trials)]
        # FIXME: RefactoringTool: Line 223: You should use a for loop here
        list(map(lambda trialModel : job.addChildTargetFn(expectationMaximisation, sequences, alignments, trialModel, options), trialModels))
        job.addFollowOnJobFn(expectationMaximisationTrials2, sequences, trialModels, outputModel, options)

def expectationMaximisationTrials2(job, sequences, trialModels, outputModel, options):
    trialHmms = [Hmm.loadHmm(x) for x in trialModels]
    if options.outputTrialHmms: #Write out the different trial hmms
        for i in range(options.trials):
            trialHmms[i].write(outputModel + ("_%i" % i))
    #Pick the trial hmm with highest likelihood
    hmm = max(trialHmms, key=lambda x : x.likelihood)
    hmm.write(outputModel)
    if options.outputXMLModelFile != None:
        open(options.outputXMLModelFile, 'w').write(prettyXml(hmmsXML(trialHmms)))
    if options.blastScoringMatrixFile != None:
        matchProbs, gapOpen, gapExtend = makeBlastScoringMatrix(hmm, [x[1] for x in reduce(lambda x, y : list(x) + list(y), list(map(fastaRead, sequences.split())))])
        fH = open(options.blastScoringMatrixFile, 'w')
        writeLastzScoringMatrix(fH, matchProbs, gapOpen, gapExtend)
        fH.close()
    job.log("Summary of trials:" + prettyXml(hmmsXML(trialHmms)))
    job.log("Hmm with highest likelihood: %s" % hmm.likelihood)

def hmmsXML(hmms):
    """Converts a bunch of hmms into a stupid little XML file that can be used to plot relevant statistics in a sensible way
    """
    #Get distributions on likelihood convergence and report.
    if len(hmms) == 0:
        raise RuntimeError("No hmms to summarise")

    #Do some checks that they are all the same time of hmm
    stateNumber = hmms[0].stateNumber
    modelType = hmms[0].modelType
    for hmm in hmms[1:]:
        if hmm.modelType != modelType:
            raise RuntimeError("Hmms not all of the same type")
        if hmm.stateNumber != stateNumber:
            raise RuntimeError("Hmms do not all have the same number of states")

    parent = ET.Element("hmms", { "modelType":str(modelType), "stateNumber":str(stateNumber)})

    #For each hmm
    for hmm in hmms:
        child = ET.SubElement(parent, "hmm")
        child.attrib["likelihood"] = str(hmm.likelihood)
        child.attrib["runningLikelihoods"] = "\t".join(map(str, hmm.runningLikelihoods))
        child.attrib["transitions"] = "\t".join(map(str, hmm.transitions))
        child.attrib["emissions"] = "\t".join(map(str, hmm.emissions))

    #Plot aggregate distributions

    ##Get the distribution on likelihoods.
    likelihoods = [x.likelihood for x in hmms]
    parent.attrib["maxLikelihood"] = str(max(likelihoods))
    parent.attrib["likelihoods"] = "\t".join(map(str, likelihoods))
    parent.attrib["likelihoodAvg"] = str(numpy.average(likelihoods))
    parent.attrib["likelihoodStdDev"] = str(numpy.std(likelihoods))

    def statFn(values, node):
        node.attrib["max"] = str(max(values))
        node.attrib["avg"] = str(numpy.average(values))
        node.attrib["std"] = str(numpy.std(values))
        node.attrib["min"] = str(min(values))
        node.attrib["distribution"] = "\t".join(map(str, values))

    #For each transitions report ML estimate, and distribution of parameters + variance.
    for fromState in range(stateNumber):
        for toState in range(stateNumber):
            statFn([x.transitions[fromState*stateNumber + toState] for x in hmms],
                   ET.SubElement(parent, "transition", {"from":str(fromState), "to":str(toState)}))

    #For each emission report ML estimate, and distribution of parameters + variance.
    for state in range(stateNumber):
        for x in range(4):
            for y in range(4):
                statFn([z.emissions[state * 16 + x * 4 + y] for z in hmms],
                       ET.SubElement(parent, "emission", {"state":str(state), "x":"ACGT"[x], "y":"ACGT"[y]}))

    return parent

def makeBlastScoringMatrix(hmm, sequences):
    """Converts an hmm into a lastz style scoring matrix
    """
    #convert to a three state hmm
    hmm2 = Hmm("threeState")
    hmm2.transitions = hmm.transitions[:3] + hmm.transitions[hmm.stateNumber*1:hmm.stateNumber*1+3] + hmm.transitions[hmm.stateNumber*2:hmm.stateNumber*2+3]
    hmm2.emissions = hmm.emissions[:3 * SYMBOL_NUMBER**2]
    hmm2.normalise()
    hmm = hmm2

    #Get gap distribution, assuming we include reverse complement sequences then it's fraction of GCs
    gcFraction = sum([sum([1.0 if y in 'GC' else 0.0 for y in x]) for x in sequences]) / sum(map(len, sequences))
    logger.debug("Got the GC fraction in the sequences for making the scoring matrix: %s" % gcFraction)
    baseProb = lambda x : gcFraction/2.0 if x in (1,2) else (1.0 - gcFraction)/2.0

    #Calculate match matrix
    logger.debug("Original match probs: %s" % " ".join(map(str, hmm.emissions[:SYMBOL_NUMBER**2])))
    matchProbs = [ hmm.emissions[x * SYMBOL_NUMBER + y] / (baseProb(x) * baseProb(y)) for x, y in product(list(range(SYMBOL_NUMBER)), list(range(SYMBOL_NUMBER))) ]
    logger.debug("Blast emission match probs: %s" % " ".join(map(str, matchProbs)))
    matchContinue = hmm.transitions[0]
    #The 6.94 is the 1/100th the sum of the lastz scoring matrix
    nProb = math.sqrt(math.exp((6.94+sum([math.log(x * matchContinue) for x in matchProbs]))/len(matchProbs)))
    logger.debug("N prob is: %s" % nProb) #Note it may go above 1!
    weight=100
    matchProbs = [weight*math.log((x * matchContinue) / nProb**2) for x in matchProbs]
    logger.debug("Blast match probs, %s: %s" % (sum(matchProbs)/4.0, " ".join(map(str, matchProbs))))

    #Calculate gap open
    gapOpen = weight*math.log((0.5 * (hmm.transitions[1]/nProb + hmm.transitions[2]/nProb)) * \
    ((hmm.transitions[hmm.stateNumber*1 + 0] + hmm.transitions[hmm.stateNumber*2 + 0])/(2*nProb**2)) * \
    ((nProb**2)/matchContinue))
    logger.debug("Gap open: %s" % gapOpen)

    #Calculate gap extend
    gapContinue = weight*math.log(0.5 * (hmm.transitions[hmm.stateNumber*1 + 1]/nProb + hmm.transitions[hmm.stateNumber*2 + 2]/nProb))
    logger.debug("Gap continue: %s" % gapContinue)

    return matchProbs, gapOpen, gapContinue

def writeLastzScoringMatrix(fileHandle, matchProbs, gapOpen, gapExtend):
    """# This matches the default scoring set for BLASTZ

    bad_score          = X:-1000  # used for sub['X'][*] and sub[*]['X']
    fill_score         = -100     # used when sub[*][*] is not defined
    gap_open_penalty   =  400
    gap_extend_penalty =   30

         A     C     G     T
    A   91  -114   -31  -123
    C -114   100  -125   -31
    G  -31  -125   100  -114
    T -123   -31  -114    91
    """
    fileHandle.write("gap_open_penalty = %s\n" % int(round(-gapOpen)))
    fileHandle.write("gap_extend_penalty = %s\n" % int(round(-gapExtend)))
    bases = "ACGT"
    fileHandle.write("\t\t" + "\t".join(bases) + "\n")
    for x in range(4):
        fileHandle.write("\t%s\t%s\n" % (bases[x], "\t".join([str(int(round(x))) for x in matchProbs[x*SYMBOL_NUMBER:((x+1)*SYMBOL_NUMBER)]])))

class Options:
    """Dictionary representing options, can be used for running pipeline within another job.
    """
    def __init__(self):
        self.modelType="fiveState"
        self.inputModel=None
        self.iterations=10
        self.trials=3
        self.outputTrialHmms = False
        self.randomStart=False
        self.optionsToRealign="--diagonalExpansion=10 --splitMatrixBiggerThanThis=3000"
        self.updateTheBand=False
        self.maxAlignmentLengthPerJob=1000000
        self.maxAlignmentLengthToSample=50000000
        self.useDefaultModelAsStart = False
        self.setJukesCantorStartingEmissions=None
        self.tieEmissions = False
        self.trainEmissions=False
        self.outputXMLModelFile = None
        self.blastScoringMatrixFile = None

def addExpectationMaximisationOptions(parser, options):
    group = parser.add_argument_group( "Expectation Maximisation Options", "These are options are used in doing expectation maximisation on the reads.")
    group.add_argument("--inputModel", default=options.inputModel, help="Input model")
    group.add_argument("--outputModel", default="hmm.txt", help="File to write the model in")
    group.add_argument("--outputXMLModelFile", default=options.outputXMLModelFile, help="File to write XML representation of model in - useful for stats")
    group.add_argument("--modelType", default=options.modelType, help="Specify the model type, currently either fiveState, threeState, threeStateAsymmetric")
    group.add_argument("--iterations", default=options.iterations, help="Number of iterations of EM", type=int)
    group.add_argument("--trials", default=options.trials, help="Number of independent EM trials. The model with the highest likelihood will be reported. Will only work if randomStart=True", type=int)
    group.add_argument("--outputTrialHmms", default=options.outputTrialHmms, help="Writes out the final trained hmm for each trial, as outputModel + _i", action="store_true")
    group.add_argument("--randomStart", default=options.randomStart, help="Iterate start model with small random values, else all values are equal", action="store_true")
    group.add_argument("--optionsToRealign", default=options.optionsToRealign, help="Further options to pass to cPecanRalign when computing expectation values, should be passed as a quoted string")
    group.add_argument("--updateTheBand", default=options.updateTheBand, help="After each iteration of EM update the set of alignments by realigning them, so allowing stochastic updating of the constraints. This does not alter the input alignments file", action="store_true")
    group.add_argument("--maxAlignmentLengthPerJob", default=options.maxAlignmentLengthPerJob, help="Maximum total alignment length of alignments to include in one job during EM..", type=int)
    group.add_argument("--maxAlignmentLengthToSample", default=options.maxAlignmentLengthToSample, help="Maximum total alignment length of alignments to include in doing EM. Alignments are randomly sampled without replacement to achieve maximum. "
                       "The alignment length of an alignment is the avg. of the lengths of the query and target substrings covered by the alignment.", type=int)
    group.add_argument("--useDefaultModelAsStart", default=options.useDefaultModelAsStart, help="Use the default BAR hmm model as the starting point", action="store_true")
    group.add_argument("--setJukesCantorStartingEmissions", default=options.setJukesCantorStartingEmissions, help="[double] Set the starting hmm emissions by jukes cantor expectation, using given subs/site estimate", type=float)
    group.add_argument("--trainEmissions", default=options.trainEmissions, help="Train the emissions as well as the transitions.", action="store_true")
    group.add_argument("--tieEmissions", default=options.tieEmissions, help="Normalise all emissions to reflect overall level of diversity, but be tied to not reflect differences between different bases, other than identity/difference.", action="store_true")
    group.add_argument("--blastScoringMatrixFile", default=options.blastScoringMatrixFile, help="Calculate a BLAST scoring matrix from the HMM, output in Lastz/Blastz format")

def main():
    #Parse the inputs args/options
    parser = argparse.ArgumentParser(description="train pair-HMMs of cPecan")
    parser.add_argument("--sequences", dest="sequences", help="Quoted list of fasta files containing sequences")
    parser.add_argument("--alignments", dest="alignments", help="Cigar file ")
    addExpectationMaximisationOptions(parser, Options())
    Job.Runner.addToilOptions(parser)
    options = parser.parse_args()

    #Log the inputs
    logger.info("Got '%s' sequences, '%s' alignments file, '%s' output model and '%s' iterations of training" % (options.sequences, options.alignments, options.outputModel, options.iterations))

    with Toil(options) as toil:
        if options.restart:
            toil.restart()
        else:
            toil.start(Job.wrapJobFn(expectationMaximisationTrials, options.sequences, options.alignments, options.outputModel, options))


if __name__ == '__main__':
    from cPecan.cPecanEm import *
    main()
