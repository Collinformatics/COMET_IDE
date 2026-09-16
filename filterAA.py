# PURPOSE: This code will load in your extracted substrates for processing

# IMPORTANT: Process all of your data with extractSubstrates before using this script


from functions import getFileNames, NGS
import os
import sys



# ===================================== User Inputs ======================================
# Input 1: Select Dataset
inEnzymeName = 'Mpro2'
inPathFolder = os.path.join('Enzymes', inEnzymeName)
inSaveFigures = True
inSetFigureTimer = False

# Input 2: Computational Parameters
inFixResidues = True
inFixedResidue = ['Q']
inFixedPosition = [5]
inExcludeResidues = False
inExcludedResidue = [['F','Y']]
inExcludedPosition = [5]
inMinimumSubstrateCount = 1
inShowSampleSize = True
inCodonSequence = 'NNS' # Baseline probs of degenerate codons (can be N, S, or K)
inUseCodonProb = False # Use AA prob from inCodonSequence to calculate enrichment
inAvgInitialProb = False
inDropResidue = ['R9'] # To drop 9th to last AA: ['R9'], For nothing: []

# Input 3: Making Figures
inBlockFigures = False
inPlotAADistribution = False
inPlotEntropy = True
inPlotEnrichmentMap = True
inPlotEnrichmentMapScaled = False
inPlotLogo = True
inPlotWeblogo = True
inPlotMotifEnrichment = True
inPlotWordCloud = True
inPlotBarGraphs = True
inPlotPCA = False
inPlotCounts = False
inPlotPositionalProbDist = False # For understanding shannon entropy
if inBlockFigures:
    inPlotAADistribution = False
    inPlotEntropy = False
    inPlotEnrichmentMap = False
    inPlotEnrichmentMapScaled = False
    inPlotLogo = False
    inPlotWeblogo = False
    inPlotMotifEnrichment = False
    inPlotWordCloud = False # Word cloud
    inPlotMotifEnrichment = False
    inPlotBarGraphs = False
    inPlotPCA = False
    inPlotCounts = False

# Input 4: Inspecting The data
inPrintNumber = 10
inFindSequences = False
inFindSeq = ['AVLQS', 'VILQS','VILQT','VILQS','VILHS','VIMQS','VPLQS','NILQS']
inFindSeq = ['CC', 'CI', 'CV', 'VV', 'II', 'IC', 'IV', 'VV', 'VC', 'VI']

# Input 5: CSV
inSaveCSV = False # Save substrates in a csv file
inMinSubsCSV = 10000 # Minimum counts for saved substrates
inSubLengthCSV = False # If: False, use full seq, If: 6, use 6 AA seq
inUseBgSubs = True
inExcludeSeq = 'LQ'
inMaxBgSubstrateCount = 20
inModulo = 2500 # Increase to select fewer Bg substrates
inScaleModulo = False # Continually increase modulus value to keep fewer low count bg subs

# Input 6: Plot Heatmap
inShowEnrichmentScores = True
inShowEnrichmentAsSquares = False
inYLabelEnrichmentMap = 2 # 0 for full Residue name, 1 for 3-letter code, 2 for 1 letter

# Input 7: Plot Sequence Motif
inBigLettersOnTop = False
inLimitYAxis = False

# Input 8: Word Cloud
inLimitWords = True
inTotalWords = 50

# Input 9: Bar Graphs
inNSequences = 50
inPlotAllSequences = False

# Input 10: PCA
inPCAMotif = False
inNumberOfPCs = 2
inTotalSubsPCA = 10000
inIncludeSubCountsESM = False
inExtractPopulations = False
inPlotEntropyPCAPopulations = False
inAdjustZeroCounts = False # Prevent counts of 0 in PCA EM & Motif

# Input 11: Delete Datapoint
inDeleteSubstrates = [''] # Delete individual substrates
inMaxSubstrateCounts = False #3*10**5 # Set a maximum substrate count



# =================================== Setup Parameters ===================================
# Colors:
white = '\033[38;2;255;255;255m'
greyDark = '\033[38;2;144;144;144m'
purple = '\033[38;2;189;22;255m'
magenta = '\033[38;2;255;0;128m'
pink = '\033[38;2;255;0;242m'
cyan = '\033[38;2;22;255;212m'
green = '\033[38;2;5;232;49m'
greenLight = '\033[38;2;204;255;188m'
greenDark = '\033[38;2;30;121;13m'
yellow = '\033[38;2;255;217;24m'
orange = '\033[38;2;247;151;31m'
red = '\033[91m'
resetColor = '\033[0m'

# Load: Dataset labels
enzymeName, filesInitial, filesFinal, labelAAPos = getFileNames(enzyme=inEnzymeName)

inPlotMotifEnrichmentNBars = True



# =================================== Initialize Class ===================================
ngs = NGS(
    enzyme=inEnzymeName, enzymeName=enzymeName, substrateLength=len(labelAAPos),
    filterSubs=inFixResidues, fixedAA=inFixedResidue, fixedPosition=inFixedPosition,
    excludeAAs=inExcludeResidues, excludeAA=inExcludedResidue,
    excludePosition=inExcludedPosition, minCounts=inMinimumSubstrateCount,
    minEntropy=None, figEMSquares=inShowEnrichmentAsSquares,
    xAxisLabels=labelAAPos, printNumber=inPrintNumber, showNValues=inShowSampleSize,
    bigAAonTop=inBigLettersOnTop, findMotif=False, folderPath=inPathFolder,
    filesInit=filesInitial, filesFinal=filesFinal, plotPosS=inPlotEntropy,
    plotFigEM=inPlotEnrichmentMap, plotFigEMScaled=inPlotEnrichmentMapScaled,
    plotFigLogo=inPlotLogo, plotFigWebLogo=inPlotWeblogo,
    plotFigMotifEnrich=inPlotMotifEnrichment, plotFigWords=inPlotWordCloud,
    wordLimit=inLimitWords, wordsTotal=inTotalWords, plotFigBars=inPlotBarGraphs,
    NSubBars=inNSequences, plotFigPCA=inPlotPCA, numPCs=inNumberOfPCs,
    NSubsPCA=inTotalSubsPCA, plotSuffixTree=False,
    saveFigures=inSaveFigures, setFigureTimer=inSetFigureTimer
)



# ====================================== Load data =======================================
# Get dataset tag
ngs.getDatasetTag(useCodonProb=inUseCodonProb, codon=inCodonSequence)

# Load: Counts
countsInitial, countsInitialTotal = ngs.loadCounts(
    filter=False, fileType='Initial Sort', dropColumn=inDropResidue
)

# Load: Substrates
if inFindSequences or inUseBgSubs:
    substratesInitial, totalSubsInitial = ngs.loadUnfilteredSubs(loadInitial=True,
                                                                 dropColumn=inDropResidue)

# Calculate: Initial sort probabilities
rfInitial = ngs.calculateRF(counts=countsInitial, N=countsInitialTotal,
                              fileType='Initial Sort', calcAvg=inAvgInitialProb)

substratesInitial = None
loadUnfilteredSubs = False
filePathFixedCountsFinal, filePathFixedSubsFinal = None, None
substratesFinal, countsFinal, countsFinalTotal = None, None, None
if inFixResidues or inExcludeResidues:
    filePathFixedSubsFinal, filePathFixedCountsFinal = (
        ngs.getFilePath(datasetTag=ngs.datasetTagSave, sortType='FinalSort'))

    # Verify that the file exists
    if (os.path.exists(filePathFixedSubsFinal) and
            os.path.exists(filePathFixedCountsFinal)):

        # Load: Counts
        countsFinal, countsFinalTotal = ngs.loadCounts(
            filter=True, fileType='Final Sort',
            datasetTag=ngs.datasetTagSave, dropColumn=inDropResidue
        )

        # Load: Substrates
        substratesFinal, totalSubsFinal = ngs.loadSubstratesFiltered()

        if countsFinalTotal != totalSubsFinal:
            print(f'{orange}ERROR: '
                  f'The total number of Loaded Counts ({cyan}{countsFinalTotal:,}'
                  f'{orange}) =/= number of Total Substrates '
                  f'({cyan}{totalSubsFinal:,}{orange})\n')
            sys.exit()
    else:
        loadUnfilteredSubs = True

        # Load: Substrates
        substratesFinal, totalSubsFinal = ngs.loadUnfilteredSubs(loadFinal=True,
                                                                dropColumn=inDropResidue)
else:
    substratesFinal, totalSubsFinal = ngs.loadUnfilteredSubs(loadFinal=True,
                                                             dropColumn=inDropResidue)

    # Load: Counts
    countsFinal, countsFinalTotal = ngs.loadCounts(fileType='Final Sort', filter=False,
                                                   dropColumn=inDropResidue)


def limCounts(substrates, maxCounts=5*10**5):
    print('================================= Limit Counts '
          '==================================')
    print(f'Max Counts: {red}{maxCounts:,}{resetColor}')
    subs = {}
    for substrate, counts in substrates.items():
        if counts <= maxCounts:
            subs[substrate] = counts

    print('Substrates:')
    for i, (substrate, counts) in enumerate(subs.items()):
        print(f'    {pink}{substrate}{resetColor}, Counts: {red}{counts:,}{resetColor}')
        if i == inPrintNumber:
            break
    totalCounts = sum(subs.values())
    print(f'\nTotal Substrates: {red}{totalCounts:,}{resetColor}\n'
          f'Unique Sequences: {red}{len(subs.values()):,}{resetColor}\n\n')

    return subs, totalCounts

if inMaxSubstrateCounts:
    substratesFinal, countsFinalTotal = limCounts(substrates=substratesFinal,
                                                  maxCounts=inMaxSubstrateCounts)



# ================================== Evaluate The Data ===================================
saveSubs = False
if inFixResidues and loadUnfilteredSubs:
    saveSubs = True
    substratesFinal, countsFinalTotal = ngs.fixResidue(
        substrates=substratesFinal, fixedString=ngs.datasetTag,
        printRankedSubs=True, sortType='Final Sort')
elif inExcludeResidues and loadUnfilteredSubs:
    saveSubs = True
    substratesFinal, countsFinalTotal = ngs.exclResidue(
        substrates=substratesFinal, fixedString=ngs.datasetTag,
        printRankedSubs=True, sortType='Final Sort')


# Delete datapoint
if not isinstance(inDeleteSubstrates, list):
    inDeleteSubstrates = list(inDeleteSubstrates)
if inDeleteSubstrates[0]:
    print('=============================== Delete Substrates '
          '===============================')
    print(f'Delete Substrates:')
    for substrate in inDeleteSubstrates:
        print(f'    {greenLight}{substrate}{resetColor}')
        substratesFinal.pop(substrate, None)

    print('\nSubs:')
    for i, (substrate, count) in enumerate(substratesFinal.items()):
        print(f'    {pink}{substrate}{resetColor}, Counts: {red}{count:,}{resetColor}')
        if i > 5:
            break
    print('\n')


if countsFinal is None:
    # Count fixed substrates
    countsFinal, countsFinalTotal = ngs.countResidues(substrates=substratesFinal,
                                                      datasetType='Final Sort')

# Save the data
if saveSubs and not inMaxSubstrateCounts:
    ngs.saveData(substrates=substratesFinal, counts=countsFinal)

# # Filter counts matrix
# if inDropResidue:
#     countsFinal = ngs.dropColumnsFromMatrix(countMatrix=countsFinal,
#                                             datasetType='Final Sort',
#                                             dropColumn=inDropResidue)

# Display current sample size
ngs.recordSampleSize(NInitial=countsInitialTotal, NFinal=countsFinalTotal,
                     NFinalUnique=len(substratesFinal.keys()))

# Calculate: RF
rfFinal = ngs.calculateRF(counts=countsFinal, N=countsFinalTotal,
                            fileType='Final Sort')

if inPlotAADistribution:
    # Plot: AA probabilities in initial and final sorts
    ngs.plotLibraryAADist(rfInitial=rfInitial, rfFinal=rfFinal,
                            codonType=inCodonSequence, datasetTag=ngs.datasetTag)

    # Evaluate: Degenerate codon probabilities
    probCodon = ngs.calculateProbCodon(codonSeq=inCodonSequence)

    # Plot: Codon probabilities
    ngs.plotLibraryAADist(rfInitial=rfFinal, rfFinal=probCodon,
                            codonType=inCodonSequence, datasetTag=inCodonSequence,
                            skipInitial=True)

# Calculate: Positional entropy
entropy = ngs.calculateEntropy(rf=rfFinal)

# Calculate: Enrichment scores
if inUseCodonProb:
    # Evaluate: Degenerate codon probabilities
    probCodon = ngs.calculateProbCodon(codonSeq=inCodonSequence)
    enrichmentScores = ngs.calculateEnrichment(rfInitial=probCodon,
                                               rfFinal=rfFinal)
else:
    enrichmentScores = ngs.calculateEnrichment(rfInitial=rfInitial,
                                               rfFinal=rfFinal)

# Create csv
if inSaveCSV:
    if inUseBgSubs:
        ngs.saveSubstrateCSV(
            seqs=substratesFinal, initialRF=rfInitial, finalRF=rfFinal,
            minCounts=inMinSubsCSV, seqsBg=substratesInitial, excludeAA=inExcludeSeq,
            maxCountsBg=inMaxBgSubstrateCount, mod=inModulo, modScale=inScaleModulo,
            chopSeq=inSubLengthCSV
        )
    else:
        ngs.saveSubstrateCSV(
            seqs=substratesFinal, initialRF=rfInitial,
            finalRF=rfFinal, minCounts=inMinSubsCSV, chopSeq=inSubLengthCSV
        )

# Plot: PCA
if inPlotPCA:
    if inPCAMotif:
        substrates = ngs.getMotif(substrates=substratesFinal)
        datasetTag = ngs.datasetTagMotif
        saveTag = ngs.datasetTagSave
        labelPos = ngs.xAxisLabelsMotif
        subsPCA = substrates

    else:
        datasetTag = ngs.datasetTag
        saveTag = datasetTag
        labelPos = labelAAPos
        subsPCA = substratesFinal

    # Generate ESM embeddings
    tokensESM, subsESM, subCountsESM = ngs.ESM(
        substrates=subsPCA, subLabel=labelPos, useSubCounts=inIncludeSubCountsESM
    )

    # Cluster substrates
    subPopulations = ngs.plotPCA(
        substrates=subsPCA, data=tokensESM, indices=subsESM, N=subCountsESM
    )

    # Plot: Substrate clusters
    if subPopulations is not None:
        clusterCount = len(subPopulations)
        for index, subCluster in enumerate(subPopulations):
            # Plot data
            ngs.plotSubstratePopulations(
                substrates=subCluster, probInit=rfInitial, clusterIndex=index,
                numClusters=clusterCount, datasetTag=datasetTag, saveTag=saveTag)

# Plot: Word cloud
if inPlotWordCloud:
    ngs.plotWordCloud(substrates=substratesFinal)

# Plot: Bar graphs
if inPlotBarGraphs:
    ngs.plotBarGraph(substrates=substratesFinal, dataType='Counts', plotAllSubs=False)
    if inPlotAllSequences:
        ngs.plotBarGraph(substrates=substratesFinal, dataType='Counts', plotAllSubs=True)
    ngs.plotBarGraph(substrates=substratesFinal, dataType='RF', plotAllSubs=False)
    if inPlotAllSequences:
        ngs.plotBarGraph(substrates=substratesFinal, dataType='RF', plotAllSubs=True)

# Find sequences
if inFindSequences:
    ngs.findSequence(substrates=substratesFinal,
                     sequence=inFindSeq,
                     sortType='Final Sort')

# Plot counts
if inPlotCounts:
    # Plot the data
    ngs.plotMatrix(data=countsFinal, figLabel=ngs.datasetTag,
                   totalCounts=countsFinalTotal)

# Plot AA distributions
if inPlotPositionalProbDist:
    ngs.plotPositionalProbDist(probability=rfFinal, entropyScores=entropy,
                               sortType='Final Sort', datasetTag=ngs.datasetTag)
