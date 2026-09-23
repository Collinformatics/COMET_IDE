# PURPOSE: This script contains the functions that you will need to process your NGS data
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
import csv
import esm
import gzip
from itertools import combinations, product
import json
import logomaker
import math
import matplotlib
from matplotlib.lines import Line2D
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap, LinearSegmentedColormap, Normalize
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import RectangleSelector
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import os
import pandas as pd
import random
import seaborn as sns
import sys
import threading
import warnings
from wordcloud import WordCloud



# ================================== Setup Residue List ==================================
defaultResidues = (('Alanine', 'Ala', 'A'), ('Arginine', 'Arg', 'R'),
                   ('Asparagine', 'Asn', 'N'), ('Aspartic Acid', 'Asp', 'D'),
                   ('Cysteine', 'Cys', 'C'),  ('Glutamic Acid', 'Glu', 'E'),
                   ('Glutamine', 'Gln', 'Q'),('Glycine', 'Gly', 'G'),
                   ('Histidine', 'His ', 'H'),('Isoleucine', 'Ile', 'I'),
                   ('Leucine', 'Leu', 'L'), ('Lysine', 'Lys', 'K'),
                   ('Methionine', 'Met', 'M'), ('Phenylalanine', 'Phe', 'F'),
                   ('Proline', 'Pro', 'P'), ('Serine', 'Ser', 'S'),
                   ('Threonine', 'Thr', 'T'),('Tryptophan', 'Typ', 'W'),
                   ('Tyrosine', 'Tyr', 'Y'), ('Valine', 'Val', 'V'))



# ===================================== Set Options ======================================
# Colors: Console
greyDark = '\033[38;2;144;144;144m'
white = '\033[38;2;255;255;255m'
purple = '\033[38;2;189;22;255m'
magenta = '\033[38;2;255;0;128m'
pink = '\033[38;2;255;0;242m'
cyan = '\033[38;2;22;255;212m'
blue = '\033[38;5;51m'
green = '\033[38;2;5;232;49m'
greenLight = '\033[38;2;204;255;188m'
greenLightB = '\033[38;2;204;255;188m'
greenDark = '\033[38;2;30;121;13m'
yellow = '\033[38;2;255;217;24m'
orange = '\033[38;2;247;151;31m'
red = '\033[91m'
resetColor = '\033[0m'



# =================================== Define Functions ===================================
def getFileNames(enzyme):
    if enzyme.lower() == 'eln' or enzyme.lower() == 'hne':
        enzyme = 'Human Neutrophil Elastase'
        inFileNamesInitialSort = ['ELN-I_S1_L001', 'ELN-I_S1_L002']
        inFileNamesFinalSort = ['ELN-R4_S2_L001', 'ELN-R4_S2_L002']
        inAAPositions = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8']
    elif enzyme.lower() == 'ide':
        enzyme = 'IDE'
        inFileNamesInitialSort = ['IDE-I_S3_L001']
        inFileNamesFinalSort = ['IDE-F_S5_L001']
        inAAPositions = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8']
    elif enzyme.lower() == 'ide prev':
        enzyme = 'IDE Prev'
        inFileNamesInitialSort = ['IDE-S1_L001', 'IDE-S1_L002',
                                  'IDE-S1_L003', 'IDE-S1_L004']
        inFileNamesFinalSort = ['IDE-S2_L001', 'IDE-S2_L002',
                                'IDE-S2_L003', 'IDE-S2_L004']
        inAAPositions = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8']
    elif enzyme.lower() == 'mpro' or enzyme.lower() == 'mpro1':
        enzyme = f'SARS-CoV M{'ᵖʳᵒ'}'
        # inFileNamesInitialSort = ['Mpro1-I_S1_L001', 'Mpro1-I_S1_L002',
        #                           'Mpro1-I_S1_L003', 'Mpro1-I_S1_L004']
        # inFileNamesFinalSort = ['Mpro1-R4_S3_L001', 'Mpro1-R4_S3_L002',
        #                         'Mpro1-R4_S3_L003', 'Mpro1-R4_S3_L004']
        inFileNamesInitialSort = ['Mpro-I_S1_L001']
        inFileNamesFinalSort = ['Mpro-R4_S3_L001']
        inAAPositions = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8']
    elif enzyme.lower() == 'mpro2':
        enzyme = f'SARS-CoV-2 M{'ᵖʳᵒ'}'
        # inFileNamesInitialSort = ['Mpro2-I_S1_L001', 'Mpro2-I_S1_L002',
        #                           'Mpro2-I_S1_L003', 'Mpro2-I_S1_L004']
        # inFileNamesFinalSort = ['Mpro2-R4_S3_L001', 'Mpro2-R4_S3_L002',
        #                         'Mpro2-R4_S3_L003', 'Mpro2-R4_S3_L004']
        inFileNamesInitialSort = ['Mpro2-I_S1_L001']
        inFileNamesFinalSort = ['Mpro2-R4_S3_L001']
        inAAPositions = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8']
    elif enzyme.lower() == 'mpro2-lq':
        enzyme = f'SARS-CoV-2 M{'ᵖʳᵒ'} LQ-NNS'
        inFileNamesInitialSort = ['Mpro2-LQ-Initial_S3_L001', 'Mpro2-LQ-Initial_S3_L002',
                                  'Mpro2-LQ-Initial_S3_L003', 'Mpro2-LQ-Initial_S3_L004']
        inFileNamesFinalSort = ['Mpro2-LQ-R5_S4_L001', 'Mpro2-LQ-R5_S4_L002',
                                'Mpro2-LQ-R5_S4_L003', 'Mpro2-LQ-R5_S4_L004']
        inAAPositions = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8']
    elif enzyme.lower() == 'mmp7':
        enzyme = 'Matrix Metalloproteinase-7'
        inFileNamesInitialSort = ['MMP7-I_S3_L001', 'MMP7-I_S3_L002']
        inFileNamesFinalSort = ['MMP7-R4_S4_L001', 'MMP7-R4_S4_L002']
        inAAPositions = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8']
    elif enzyme.lower() == 'fyn':
        enzyme = 'Fyn'
        inFileNamesInitialSort = ['Fyn-I_S6_L001', 'Fyn-I_S6_L002']
        inFileNamesFinalSort = ['Fyn-F_S1_L001', 'Fyn-F_S1_L002']
        inAAPositions = ['-4', '-3', '-2', '-1', '0', '1', '2', '3', '4']
    elif enzyme.lower() == 'src':
        enzyme = 'Src'
        inFileNamesInitialSort = ['Src-I_S4_L001', 'Src-I_S4_L002']
        inFileNamesFinalSort = ['Src-F_S2_L001', 'Src-F_S2_L002']
        inAAPositions = ['-4', '-3', '-2', '-1', '0', '1', '2', '3', '4']
    elif enzyme.lower() == 'den':
        enzyme = 'Dengue Virus - NS2B/NS3'
        inFileNamesInitialSort = ['DEN-I_S1_L001']
        inFileNamesFinalSort = ['DEN-F_S2_L001']
        inAAPositions = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9']
    elif enzyme.lower() == 'veevprev' or enzyme.lower() == 'veprev':
        enzyme = 'Venezuelan Equine Encephalitis Virus - nsP2'
        inFileNamesInitialSort = ['VEEV-I_S1']
        inFileNamesFinalSort = ['VEEV-R4_S2']
        inAAPositions = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10']
    elif enzyme.lower() == 'veev' or enzyme.lower() == 've':
        enzyme = 'Venezuelan Equine Encephalitis Virus - nsP2'
        inFileNamesInitialSort = ['VEEV-I_S3']
        inFileNamesFinalSort = ['VEEV-F_S4']
        inAAPositions = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9']
    elif enzyme.lower() == 'wnv':
        enzyme = 'West Nile Virus - NS2B/NS3'
        inFileNamesInitialSort = ['WNV-I_S3_L001']
        inFileNamesFinalSort = ['WNV-F_S4_L001']
        inAAPositions = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9']
    elif enzyme.lower() == 'zk':
        enzyme = 'Zika Virus - NS2B/NS3'
        inFileNamesInitialSort = ['ZK-I_S1_L001']
        inFileNamesFinalSort = ['ZK-F_S2_L001']
        inAAPositions = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9']
    else:
        print(f'{orange}ERROR: There are no file names for {cyan}{enzyme}{orange}\n'
              f'       Add information to the "{cyan}filePaths{resetColor}" function in '
              f'{cyan}{os.path.basename(__file__)}{resetColor}\n')
        sys.exit(1)

    return enzyme, inFileNamesInitialSort, inFileNamesFinalSort, inAAPositions


def pressKey(event):
    if event.key == 'escape':
        plt.close()
    elif event.key == 'e':
        sys.exit()
    elif event.key == 'r':
        python = sys.executable # Doesn't seem to work on windows?
        os.execl(python, python, *sys.argv)


def includeCommas(x):
    return f'{x:,}'



class NGS:
    def __init__(self, enzyme, enzymeName, substrateLength, filterSubs, fixedAA,
                 fixedPosition, excludeAAs, excludeAA, excludePosition, minCounts,
                 minEntropy, figEMSquares, xAxisLabels, printNumber, showNValues,
                 bigAAonTop, findMotif, folderPath, filesInit, filesFinal, plotPosS,
                 plotFigEM, plotFigEMScaled, plotFigLogo, plotFigWebLogo, plotFigWords,
                 wordLimit, wordsTotal, plotFigBars, NSubBars, plotFigPCA, numPCs,
                 NSubsPCA, plotSuffixTree, saveFigures, setFigureTimer,
                 translateDNA=False, minPhred=20, xAxisLabelsMotif=None,
                 motifFilter=False, releasedCounts=False, plotFigMotifEnrich=False):
        if not isinstance(fixedAA, list):
            fixedAA = [fixedAA]
        if not isinstance(fixedPosition, list):
            fixedPosition = [fixedPosition]
        if not isinstance(excludeAA, list):
            excludeAA = [excludeAA]
        if not isinstance(excludePosition, list):
            excludePosition = [excludePosition]

        # Parameters: Figures
        self.labelSizeTitle = 20 # Set fontsize
        self.labelSizeAxis = 18 # Set fontsize
        self.labelSizeTicks = 16 # Set fontsize
        self.labelSizeEM = 11 # Set fontsize
        self.lineThickness = 1.5
        self.tickLength = 4
        if figEMSquares:
            self.figSizeEM = (5, 8)  # (width, height)
        else:
            self.figSizeEM = (9.5, 8)
        self.figSize = (9.5, 8)
        self.figSizeMini = (self.figSize[0], 6)
        self.entropy, self.entropyMax = None, None
        self.eMap, self.eMapScaled = None, None
        self.eMapReleased, self.eMapReleasedScaled = None, None
        self.heights, self.weblogo = None, None
        self.plotFigEntropy = plotPosS
        self.plotFigEM = plotFigEM
        self.plotFigEMScaled = plotFigEMScaled
        self.plotFigLogo = plotFigLogo
        self.plotFigWebLogo = plotFigWebLogo
        self.plotFigMotifEnrich = plotFigMotifEnrich
        self.plotFigWords = plotFigWords
        self.wordsLimit = wordLimit
        self.wordsTotal = wordsTotal
        self.plotFigBars = plotFigBars
        self.NSubBars = NSubBars
        self.plotFigPCA = plotFigPCA
        self.numPCs = numPCs
        self.NSubsPCA = NSubsPCA
        self.plotSuffixTree = plotSuffixTree
        self.datasetTag = None
        self.datasetTagSave = None
        self.datasetTagMotif = None
        self.title = ''
        self.titleCombined = ''
        self.titleMotif = ''
        self.titleReleased = ''
        self.titleWeblogo = ''
        self.titleWords = ''
        self.substrateLength = substrateLength
        self.figEMSquares = figEMSquares
        self.residueLabelType = 2  # 0 = full AA name, 1 = 3-letter code, 2 = 1 letter
        self.residues = defaultResidues
        self.letters = [residue[2] for residue in self.residues]
        self.colorsAA = NGS.residueColors()
        self.showSampleSize = showNValues
        self.nSubsInitial = 0
        self.nSubsFinal = 0
        self.bigAAonTop = bigAAonTop
        self.findMotif = findMotif
        self.motifFilter = motifFilter
        self.initialize = True  # filterMotif.py: Set to False after NGS.calculateEntropy()
        self.motifTag = ''
        self.subFrame = None
        self.motifIndex = None
        self.motifIndexExtracted = []
        self.saveFigures = saveFigures
        self.setFigureTimer = setFigureTimer
        self.figureTimerDuration = 0.5
        self.saveFigureIteration = 0
        self.figureResolution = 300

        # Parameters: Dataset
        self.enzyme = enzyme
        self.enzymeName = enzymeName
        self.filterSubs = filterSubs
        self.fixedAA = fixedAA
        self.fixedPos = fixedPosition
        self.excludeAAs = excludeAAs
        self.excludeAA = excludeAA
        self.excludePos = excludePosition
        self.minSubCount = minCounts
        self.minEntropy = minEntropy
        self.xAxisLabels = xAxisLabels
        self.xAxisLabelsMotif = xAxisLabelsMotif
        self.motifLen = None
        self.motifTags = None
        self.printNumber = printNumber
        self.selectedSubstrates = []
        self.selectedDatapoints = []
        self.rectangles = []
        self.initialize = True
        self.maxValue = 0
        self.releasedCounts = releasedCounts

        # Parameters: DNA Processing
        self.translateDNA = translateDNA # Only set as True when processing DNA seqs
        self.minQS = minPhred # Minium Phred quality score for extracted substrates
        self.fileSize = []
        self.countExtractedSubs = []
        self.percentUnusableDNASeqs = []

        # Parameters: Files
        self.filesInit = filesInit
        self.filesFinal = filesFinal
        self.pathFolder = folderPath
        self.pathData = os.path.join(self.pathFolder, 'Data')
        self.pathSaveFigs = os.path.join(self.pathFolder, 'Figures')
        self.pathFilteredSubs = None
        self.pathFilteredCounts = None

        # Parameters: Misc
        self.roundVal = 3
        np.set_printoptions(suppress=True) # Prevent data from printing in sci notation
        np.seterr(divide='ignore')

        # Verify directory paths
        if not os.path.exists(self.pathFolder):
            os.makedirs(self.pathFolder, exist_ok=True)
            # print(f'{orange}ERROR: Folder not found\n'
            #       f'Check input: "{cyan}inPathFolder{orange}"\n'
            #       f'     inPathFolder = {self.pathFolder}\n')
            # sys.exit(1)
        if self.pathData is not None:
            if not os.path.exists(self.pathData):
                os.makedirs(self.pathData, exist_ok=True)
        if self.pathSaveFigs is not None:
            if not os.path.exists(self.pathSaveFigs):
                os.makedirs(self.pathSaveFigs, exist_ok=True)

        # Set options
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 1000)
        pd.options.display.float_format = '{:,.4f}'.format



    @staticmethod
    def alert(soundPath):
        # This function can be used to play .mp3 files
        # Used it to let you know when a process has been completed
        from playsound import playsound

        if os.path.exists(soundPath):
            threading.Thread(target=playsound, args=(soundPath,)).start()
        else:
            print(f'{orange}ERROR: The alerts sound was not found at\n'
                  f'     {soundPath}{resetColor}\n\n')



    @staticmethod
    def residueColors():
        color = ['darkgreen', 'firebrick', 'deepskyblue', 'pink', 'navy', 'black', 'gold']
                 # Aliphatic, Acidic, Basic, Hydroxyl, Amide, Aromatic, Sulfur

        return {
            'A': color[0],
            'R': color[2],
            'N': color[4],
            'D': color[1],
            'C': color[6],
            'E': color[1],
            'Q': color[4],
            'G': color[0],
            'H': color[2],
            'I': color[0],
            'L': color[0],
            'K': color[2],
            'M': color[6],
            'F': color[5],
            'P': color[0],
            'S': color[3],
            'T': color[3],
            'W': color[5],
            'Y': color[5],
            'V': color[0]
        }



    @staticmethod
    def createCustomColorMap(colorType):
        colorType = colorType.lower()
        if colorType == 'counts':
            useGreen = True
            if useGreen:
                # Green
                colors = ['#FFFFFF','#ABFF9B','#39FF14','#2E9418','#2E9418',
                          '#005000']
            else:
                # Orange
                colors = ['white','white','#FF76FA','#FF50F9','#FF00F2',
                          '#CA00DF','#BD16FF']
        elif colorType == 'standard deviation':
            colors = ['white','white','#FF76FA','#FF50F9','#FF00F2','#CA00DF','#BD16FF']
        elif colorType == 'word cloud':
            # ,'#F2A900','#2E8B57','black'
            colors = ['#BF5700','#BF5700','#F8971F','#FAA338',
                      '#00C01E','#1D680D','#003000','black']
            # colors = ['#008631','#39E75F','#BF5700','#F8971F','black']
        elif colorType == 'em':
            colors = ['navy','royalblue','dodgerblue','lightskyblue','white','white',
                      'lightcoral','red','firebrick','darkred']
        else:
            print(f'{orange}ERROR: Cannot create colormap. '
                  f'Unrecognized colorType parameter: {cyan}{colorType}{resetColor}\n')
            sys.exit(1)

        # Create colormap
        if len(colors) == 1:
            colorList = [(0, colors[0]), (1, colors[0])]
        else:
            colorList = [(i / (len(colors) - 1), color) for i, color in enumerate(colors)]
        return LinearSegmentedColormap.from_list('custom_colormap', colorList)


    @staticmethod
    def roundup(val, upperLim=True):
        # Round decimal up to units of 5
        vNew = val * 10
        vDiv = vNew / 5.0
        if upperLim:
            val = np.ceil(vDiv) * 5
        else:
            val = np.floor(vDiv) * 5
        val /= 10
        # print(f'Value: {val}, {vNew}, {vDiv}')
        # print(f'* New Value: {val}\n')
        return val


    @staticmethod
    def zScore(dist, AA):
        x = dist.loc[AA]
        dist = np.asarray(dist)
        mu = dist.mean()
        sigma = dist.std(ddof=0) # population std
        return (x - mu) / sigma



    @staticmethod
    def scoreSubstrate(seq, matrix, sqroot=False):
        score = 0
        for index in range(len(seq)):
            # Evaluate substrate
            AA = seq[index]
            pos = matrix.columns[index]
            value = matrix.loc[AA, pos]
            if score == 0:
                score = value
            else:
                score *= value
        if sqroot:
            score = np.sqrt(score)
        return score



    @staticmethod
    def rankScores(values):
        if not isinstance(values, list):
            values = list(values)
        ranked = pd.Series(values).rank(ascending=False, method='min')
        return [int(x) for x in ranked]




    def loadAndTranslate(self, filePath, fileName, fileType, fixedSubs,
                         startSeq, endSeq, printQS, forwardRead):
        if forwardRead is None:
            print(f'{orange}ERROR: The file {cyan}{fileName}{orange} does not contain '
                  f'a forward read (R1) or reverse read (R2) label.\n\n'
                  f'Please update the file name, or expand the conditional statements '
                  f'in extractSubs.py that call the function NGS.loadAndTranslate() '
                  f'so that the file can be processed.')
            sys.exit(1)


        # Define file location
        if not fileName.endswith(fileType):
            fileName += f'.{fileType}'
        fileLocation = os.path.join(filePath, f'{fileName}')

        if fileType == 'fasta':
            printQS = False

        # Determine the read direction
        subSequence = None
        # Extract the substrates
        subSequence = self.translate(path=fileLocation, fileName=fileName,
                                     fileType=fileType, fixData=fixedSubs,
                                     startSeq=startSeq, endSeq=endSeq,
                                     printQS=printQS, forwardRead=forwardRead)

        return subSequence



    def translate(self, path, fileName, fileType, fixData, startSeq, endSeq, printQS,
                  forwardRead):
        if forwardRead:
            read = 'R1'
            print('============================ Translate: Forward Read '
                  '============================')
        else:
            read = 'R2'
            print('============================ Translate: Reverse Read '
                  '============================')
        subSequence = {}
        totalSeqsDNA = 0
        printedSeqs = 0


        # Evaluate the file path
        gZipped = False
        if not os.path.isfile(path):
            pathZipped = path + '.gz'
            if os.path.isfile(pathZipped):
                gZipped = True
                path = pathZipped
            else:
                print(f'{orange}ERROR: File location does not lead to a file\n'
                      f'     {path}\n'
                      f'     {pathZipped}\n')
                sys.exit(1)
        print(f'File Location:\n'
              f'     {greenDark}{path}{resetColor}\n\n')

        # Define fixed substrate tag
        if fixData:
            self.getDatasetTag()
            print(f'Evaluating fixed library: {purple}{self.datasetTag}{resetColor}')
        else:
            print(f'Evaluating library: {purple}{self.enzymeName}{resetColor}')


        def printDNA(printedSeqs):
            # Select full DNA seq
            DNA = str(datapoint.seq)
            if not forwardRead:
                DNA = Seq(DNA).reverse_complement()

            # Get: Quality score
            QS = []
            if printQS and 'phred_quality' in datapoint.letter_annotations:
                QS = datapoint.letter_annotations['phred_quality']
            print(f'DNA sequence: {DNA}')

            # Inspect full DNA seq
            if startSeq in DNA and endSeq in DNA:
                # Find: Substrate indices
                start = DNA.find(startSeq) + len(startSeq)
                end = DNA.find(endSeq)
                if not startSeq:
                    start = 0
                if not endSeq:
                    end = len(DNA)

                # Extract substrate DNA seq
                substrateDNA = DNA[start:end].strip()
                print(f'     Inspected substrate: '
                      f'{greenLightB}{substrateDNA}{resetColor}')

                if len(substrateDNA) == self.substrateLength * 3:
                    # Translate substrate
                    substrate = str(Seq.translate(substrateDNA))
                    print(f'     Inspected Substrate:'
                          f'{greenLightB} {substrate}{resetColor}')

                    # Inspect substrate seq: PRINT ONLY
                    if 'X' not in substrate and '*' not in substrate:
                        print(f'     Extracted substrate: '
                              f'{pink}{substrate}{resetColor}')
                        if printQS:
                            QS = QS[start:end]
                            print(f'     QS Substrate: {QS}')
                        printedSeqs += 1
            print()

            return printedSeqs


        def inspectDNA():
            # Inspect full DNA seq
            if startSeq in DNA and endSeq in DNA:
                # Find: Substrate indices
                start = DNA.find(startSeq) + len(startSeq)
                end = DNA.find(endSeq)
                if not startSeq:
                    start = 0
                if not endSeq:
                    end = len(DNA)

                # Get: Quality score
                QS = []
                if printQS and 'phred_quality' in datapoint.letter_annotations:
                    QS = datapoint.letter_annotations['phred_quality']

                # Extract substrate DNA seq
                substrate = DNA[start:end].strip()
                if len(substrate) == self.substrateLength * 3:
                    # Express substrate
                    substrate = str(Seq.translate(substrate))

                    # Inspect substrate seq: Keep good fixed datapoints
                    if 'X' not in substrate and '*' not in substrate:
                        # Inspect quality score
                        if printQS and QS:
                            QS = datapoint.letter_annotations['phred_quality']
                            QS = QS[start:end]
                            if any(score < self.minQS for score in QS):
                                return

                        # Record datapoint
                        if substrate in subSequence:
                            subSequence[substrate] += 1
                        else:
                            subSequence[substrate] = 1


        def inspectDNAFixed():
            # Inspect full DNA seq
            if startSeq in DNA and endSeq in DNA:
                # Find: Substrate indices
                start = DNA.find(startSeq) + len(startSeq)
                end = DNA.find(endSeq)
                if not startSeq:
                    start = 0
                if not endSeq:
                    end = len(DNA)

                # Extract substrate DNA seq
                substrate = DNA[start:end].strip()
                if len(substrate) == self.substrateLength * 3:
                    # Express substrate
                    substrate = str(Seq.translate(substrate))

                    # Inspect substrate seq: Keep good fixed datapoints
                    if 'X' not in substrate and '*' not in substrate:
                        # Inspect quality score
                        if printQS and 'phred_quality' in datapoint.letter_annotations:
                            QS = datapoint.letter_annotations['phred_quality']
                            QS = QS[start:end]
                            if any(score < self.minQS for score in QS):
                                return

                        # Fix AAs
                        keep = True
                        for idx, AA in enumerate(self.fixedAA):
                            if isinstance(AA, list):
                                if substrate[self.fixedPos[idx] - 1] not in AA:
                                    keep = False
                                    break
                            else:
                                if substrate[self.fixedPos[idx] - 1] != AA:
                                    keep = False
                                    break
                        if keep:
                            if substrate in subSequence:
                                subSequence[substrate] += 1
                            else:
                                subSequence[substrate] = 1



        def evaluateDNAQuality(totalSeqsDNA, read):
            throwaway = totalSeqsDNA - self.printNumber
            throwawayPercent = (throwaway / totalSeqsDNA) * 100
            print(f'\nExtraction Efficiency: '
                  f'{greenLight}{fileName}{resetColor}\n'
                  f'     Number of discarded sequences until {red}'
                  f'{self.printNumber} substrates{resetColor} were found in '
                  f'{purple}{read}{resetColor}: {red}{throwaway:,}{resetColor}\n'
                  f'     {yellow}Percent throwaway{resetColor}:'
                  f'{red} {round(throwawayPercent, self.roundVal)} %'
                  f'{resetColor}')


        # Load the data
        if gZipped:
            # Open the file
            with (gzip.open(path, 'rt', encoding='utf-8') as file):
                data = SeqIO.parse(file, fileType)
                warnings.simplefilter('ignore', BiopythonWarning)

                for datapoint in data:
                    totalSeqsDNA += 1
                    printedSeqs = printDNA(printedSeqs)
                    if printedSeqs == self.printNumber:
                        evaluateDNAQuality(totalSeqsDNA, read)
                        break
                if fixData:
                    print(f'\nNote: The displayed substrates were not filtered for '
                          f'{purple}{self.datasetTag}{resetColor}\n'
                          f'      The filter will be applied to the extracted substrates')
                print('')


                # Extract the substrates
                totalSeqsDNA = 0
                if forwardRead:
                    if fixData:
                        for datapoint in data:
                            # Select full DNA seq
                            DNA = str(datapoint.seq)
                            totalSeqsDNA += 1
                            inspectDNAFixed()
                    else:
                        for datapoint in data:
                            # Select full DNA seq
                            DNA = str(datapoint.seq)
                            totalSeqsDNA += 1
                            inspectDNA()
                else:
                    if fixData:
                        for datapoint in data:
                            # Select full DNA seq
                            DNA = str(datapoint.seq)
                            DNA = Seq(DNA).reverse_complement()
                            totalSeqsDNA += 1
                            inspectDNAFixed()
                    else:
                        for datapoint in data:
                            # Select full DNA seq
                            DNA = str(datapoint.seq)
                            DNA = Seq(DNA).reverse_complement()
                            totalSeqsDNA += 1
                            inspectDNA()
        else:
            # Open the file
            with open(path, 'r') as file:
                data = SeqIO.parse(file, fileType)
                warnings.simplefilter('ignore', BiopythonWarning)

                for datapoint in data:
                    totalSeqsDNA += 1
                    printedSeqs = printDNA(printedSeqs)
                    if printedSeqs >= self.printNumber:
                        evaluateDNAQuality(totalSeqsDNA, read)
                        break
                if fixData:
                    print(f'\nNote: The displayed substrates were not filtered for '
                          f'{purple}{self.datasetTag}{resetColor}\n'
                          f'      The filter will be applied to the extracted substrates')
                print('')

                # Extract the substrates
                totalSeqsDNA = 0
                if forwardRead:
                    if fixData:
                        for datapoint in data:
                            # Select full DNA seq
                            DNA = str(datapoint.seq)
                            totalSeqsDNA += 1
                            inspectDNAFixed()
                    else:
                        for datapoint in data:
                            # Select full DNA seq
                            DNA = str(datapoint.seq)
                            totalSeqsDNA += 1
                            inspectDNA()
                else:
                    if fixData:
                        for datapoint in data:
                            # Select full DNA seq
                            DNA = str(datapoint.seq)
                            DNA = Seq(DNA).reverse_complement()
                            totalSeqsDNA += 1
                            inspectDNAFixed()
                    else:
                        for datapoint in data:
                            # Select full DNA seq
                            DNA = str(datapoint.seq)
                            DNA = Seq(DNA).reverse_complement()
                            totalSeqsDNA += 1
                            inspectDNA()


        # Verify if substrates have been extracted
        if len(subSequence) == 0:
            print(f'\nNo substrates were extracted from file at:\n{path}\n\n'
                  f'Recommend: adjust variables\n'
                  f'     startSeq: {red}{startSeq}{resetColor}\n'
                  f'     endSeq: {red}{endSeq}{resetColor}')
            sys.exit(1)
        else:
            extractionCount = sum(subSequence.values())
            throwaway = (totalSeqsDNA - extractionCount)
            throwawayPercent = (throwaway / totalSeqsDNA) * 100
            self.fileSize.append(totalSeqsDNA)
            self.countExtractedSubs.append(extractionCount)
            self.percentUnusableDNASeqs.append(throwawayPercent)
            print(f'Evaluate All DNA Sequences in {purple}{read}{resetColor}: '
                  f'{greenLight}{fileName}{resetColor}\n'
                  f'     Total DNA sequences in the file: '
                  f'{red}{totalSeqsDNA:,}{resetColor}\n'
                  f'     Number of extracted Substrates: '
                  f'{red}{extractionCount:,}{resetColor}\n'
                  f'     {yellow}Percent throwaway{resetColor}:'
                  f'{red} {round(throwawayPercent, self.roundVal)} %{resetColor}\n\n')

        # Rank the substrates
        subSequence = dict(sorted(
            subSequence.items(), key=lambda x: x[1], reverse=True)
        )

        return subSequence



    def extractionEfficiency(self, files):
        print('======================== Substrate Extraction Efficiency '
              '========================')
        for index, file in enumerate(files):
            print(f'{pink}Evaluate file{resetColor}:{greenLight} {file}{resetColor}\n'
                  f'     Total DNA sequences in the file: '
                  f'{red}{self.fileSize[index]:,}{resetColor}\n'
                  f'     Number of extracted Substrates: '
                  f'{red}{self.countExtractedSubs[index]:,}{resetColor}\n'
                  f'     {yellow}Percent throwaway{resetColor}: '
                  f'{red}{round(self.percentUnusableDNASeqs[index], self.roundVal)} %'
                  f'{resetColor}\n')
        print('')



    def countResidues(self, substrates, datasetType):
        print('============================= Calculate: AA Counts '
              '==============================')
        print(f'Dataset: {purple}{datasetType}{resetColor}\n'
              f'Total substrate count: {red}{sum(substrates.values()):,}{resetColor}\n'
              f'Unique substrate count: {red}{len(substrates):,}{resetColor}')

        # Determine substrate length
        firstSub, lengthSubstrate = None, None
        countMotif = False
        for substrate in substrates.keys():
            firstSub = substrate
            lengthSubstrate = len(substrate)
            if lengthSubstrate != len(self.xAxisLabels):
                countMotif = True
            break

        # Initialize the count matrix
        if countMotif and self.xAxisLabelsMotif is not None:
            countedData = pd.DataFrame(0,
                                       index=self.letters,
                                       columns=self.xAxisLabelsMotif,
                                       dtype=int)
        else:
            countedData = pd.DataFrame(0,
                                       index=self.letters,
                                       columns=self.xAxisLabels,
                                       dtype=int)

        # Verify consistent substrate lengths
        countModulus = 100 # Check the Nth substrate length
        for index, substrate in enumerate(substrates.keys()):
            if index % countModulus == 0:
                if len(substrate) != lengthSubstrate:
                    print(f'{orange}ERROR: The substrate lengths do not match\n'
                          f'     {firstSub}: {lengthSubstrate}AA\n'
                          f'     {substrate}: {len(substrate)}AA{resetColor}\n')
                    sys.exit(1)

        # Count the occurrences of each residue
        totalSubs = 0
        if countMotif:
            # Count the AAs
            print(f'Axis: {countedData.columns}')

            for substrate, counts in substrates.items():
                totalSubs += counts
                indicesResidue = [self.letters.index(AA) for AA in substrate]
                for position, indexResidue in enumerate(indicesResidue):
                    countedData.iloc[indexResidue, position] += counts

            # Sum all columns
            columnSums = pd.DataFrame(np.sum(countedData, axis=0),
                                      columns=['Total Counts'],
                                      index=self.xAxisLabelsMotif)
        else:
            # Count the AAs
            for substrate, counts in substrates.items():
                totalSubs += counts
                indicesResidue = [self.letters.index(AA) for AA in substrate]
                for position, indexResidue in enumerate(indicesResidue):
                    countedData.iloc[indexResidue, position] += counts


            # Sum all columns
            columnSums = pd.DataFrame(np.sum(countedData, axis=0),
                                      columns=['Total Counts'],
                                      index=self.xAxisLabels)

        # Print: Counts
        columnSumsFormated = columnSums.apply(
            lambda col: col.map(includeCommas)
        ).copy()
        print(f'{columnSumsFormated}\n')
        countedDataPrint = countedData.apply(lambda col: col.map(includeCommas))
        print(f'Counted data: {purple}{self.enzymeName}{resetColor}\n'
              f'{countedDataPrint}\n\n')


        # Sanity Check: Do the sums of each column match the total number of substrates?
        #totalSubs = sum(countedData.iloc[:, 0])
        for indexColumn in countedData.columns:
            columnSum = sum(countedData.loc[:, indexColumn])
            if columnSum != totalSubs:
                print(
                    f'Counted data: {purple}{self.enzymeName}{resetColor}\n'
                    f'{countedData}\n\n'
                    f'{orange}ERROR: The total number of substrates '
                    f'({cyan}{totalSubs:,}{orange}) =/= the sum of column '
                    f'{indexColumn} ({cyan}{columnSum:,}{orange}){resetColor}\n')
                sys.exit(1)

        print(f'Total substrates: {red}{totalSubs:,}{resetColor}\n'
              f'Total Unique Substrates: {red}{len(substrates):,}{resetColor}\n\n')

        return countedData, totalSubs



    def getFilePath(self, datasetTag, sortType, motifPath=False,
                    customTag=None):
        print('============================== Define: File Paths '
              '===============================')
        # Define: File path
        enzName = self.enzymeName.replace(' - ', '-').replace(' ', '_')
        datasetTag = datasetTag.replace(' - ', '-').replace(' ', '_')
        if motifPath:
            if customTag is None:
                file = (f'{enzName}-Register_{self.datasetTagMotif}-'
                        f'{sortType}-MinCounts_{self.minSubCount}').replace(
                    '/', '_').replace(' ', '_')
                pathSubs = os.path.join(self.pathData,
                                        f'fixedMotifSubs-{file}.json')
                pathCounts = os.path.join(self.pathData,
                                          f'fixedMotifCounts-{file}.csv')
                pathCountsReleased = os.path.join(self.pathData,
                                                  f'fixedMotifCountsRel-{file}.csv')
                paths = [pathSubs, pathCounts, pathCountsReleased]
            else:
                print(f'Tag: {enzName}-{customTag}-{sortType}-'
                        f'MinCounts_{self.minSubCount}')
                file = (
                    f'{enzName}-{customTag}-{sortType}-MinCounts_{self.minSubCount}'
                ).replace('/', '_').replace(' ', '_')
                print(f'File: {file}')
                pathSubs = os.path.join(self.pathData, f'fixedMotifSubs-{file}.json')
                pathCounts = os.path.join(self.pathData, f'fixedMotifCounts-{file}.csv')
                pathCountsReleased = os.path.join(self.pathData,
                                                  f'fixedMotifCountsRel-{file}.csv')
                paths = [pathSubs, pathCounts, pathCountsReleased]
        else:
            file = (
                f'{enzName}-{datasetTag}-{sortType}-MinCounts_{self.minSubCount}'
            ).replace('/', '_')
            pathSubs = os.path.join(
                self.pathData, f'fixedSubs-{file}.json')
            pathCounts = os.path.join(
                self.pathData, f'counts-{file}.csv')
            self.pathFilteredSubs = pathSubs
            self.pathFilteredCounts = pathCounts
            paths = [pathSubs, pathCounts]

        print(f'File: {greenLight}{file}{resetColor}\n')
        print(f'File paths:{greenDark}')
        for path in paths:
            print(f'     {path}')
        print(f'{resetColor}\n')

        return paths



    def getFilePathCombined(self, sortType, loadSubs=False,
                            loadCountsRel=False):
        print('============================== Define: File Paths '
              '===============================')
        print(f'Dataset: {purple}{self.datasetTag}{resetColor}\n')
        self.motifTags = []
        paths = []

        # Define: Dataset
        dataset = None
        if loadSubs:
            dataset = 'fixedMotifSubs'
        elif loadCountsRel:
            dataset = 'fixedMotifCountsRel'
        else:
            dataset = 'fixedMotifCounts'
        if not self.releasedCounts:
            dataset = dataset.replace('Motif', '')
            if 'Counts' in dataset:
                dataset = dataset.replace('fixed', '')


        # Define: File tags
        if len(self.fixedAA) == len(self.fixedPos):
            tags = [[] for _ in range(len(self.fixedAA))]
            for index, AA in enumerate(self.fixedAA):
                if isinstance(AA, list):
                    if len(AA) == 1:
                        AA = AA[0]
                    else:
                        AA = f'[{','.join(AA)}]'
                position = self.fixedPos[index]
                if isinstance(position, int):
                    tag = f'{AA}@R{position}'
                    tags[index].append(tag)
                else:
                    if len(position) > 1:
                        for pos in position:
                            tag = f'{AA}@R{pos}'
                            tags[index].append(tag)
                    else:
                        tag = f'{AA}@R{position[0]}'
                        tags[index].append(tag)

            # Define: Motif tags
            for combo in zip(*tags):
                self.motifTags.append(" ".join(combo))
        else:
            self.motifTags = []
            if len(self.fixedAA) == 1:
                if len(self.fixedAA[0]) == 1:
                    AA = self.fixedAA[0]
                else:
                    AA = f'[{','.join(self.fixedAA[0])}]'
                for pos in self.fixedPos:
                    tag = f'{AA}@R{pos}'
                    self.motifTags.append(tag)
            else:
                print(f'{orange}ERROR: Unable to handle this condition')
                sys.exit(1)


        # Exclude residues
        excludeTag = None
        if self.excludeAAs:
            fixResidueList = []
            for index, removedAA in enumerate(self.excludeAA):
                fixResidueList.append(
                    f'{removedAA}@R'
                    f'{self.excludePos[index]}'.replace(' ', ''))
            excludeTag = f'Exclude {' '.join(fixResidueList)}'

        # Define: File path
        enzName = self.enzymeName.replace(' - ', '-').replace(' ', '_')
        for motifTag in self.motifTags:
            if self.excludeAAs:
                file = (
                    f'{enzName}-Register_{excludeTag}_Fixed_{motifTag}-{sortType}-'
                    f'MinCounts_{self.minSubCount}'
                ).replace('/', '_')
            else:
                file = (
                    f'{enzName}-Register_{motifTag}-{sortType}-'
                    f'MinCounts_{self.minSubCount}'
                ).replace('/', '_')
            paths.append(
                os.path.join(
                    self.pathData, f'{dataset}-{file.replace(' ','_')}'
                )
            )
        if loadSubs :
            paths = [f'{p}.json' for p in paths]
        else:
            paths = [f'{p}.csv' for p in paths]

        print(f'File paths 1:{greenDark}')
        for path in paths:
            print(f'     {path}')
        print(f'{resetColor}\n')

        return paths



    def plotFig(self, plt, fig):
        fig.canvas.mpl_connect('key_press_event', pressKey)
        plt.tight_layout()
        if self.setFigureTimer:
            plt.ion()
            plt.show()
            plt.pause(self.figureTimerDuration)
            plt.close(fig)
            plt.ioff()
        else:
            plt.show()



    def loadCounts(self, filter, fileType, datasetTag=None, dropColumn=False):
        print('================================== Load Counts '
              '==================================')
        if filter:
            labelFile = f'{self.enzymeName} {fileType} - Filter {datasetTag}'
        else:
            labelFile = f'{self.enzymeName} {fileType}'
        print(f'Loading Counts: {purple}{labelFile}{resetColor}\n')
        files = []
        totalCounts = 0

        # Define: File paths
        if filter:
            if self.pathFilteredCounts is None:
                print(f'{orange}ERROR: {cyan}self.pathFilteredCounts {orange}needs '
                      f'to be defined before you can load the counts.{resetColor}\n')
                sys.exit(1)
            files = [self.pathFilteredCounts]
        else:
            if 'initial' in fileType.lower():
                fileNames = self.filesInit
            else:
                fileNames = self.filesFinal

            for fileName in fileNames:
                files.append(os.path.join(self.pathData, f'counts_{fileName}.csv'))

        print(f'Loading data:')
        for filePath in files:
            # Verify if the file exists at its specified path
            if not os.path.exists(filePath):
                print(f'{orange}ERROR: File not found\n'
                      f'     {filePath}\n')
                sys.exit(1)
            print(f'     {greenDark}{filePath}{resetColor}')
        print('\n')


        #  Load: AA counts
        countedData = None
        for index, filePath in enumerate(files):
            fileName = filePath.replace(self.pathData, '')

            # Load: File
            if index == 0:
                data = pd.read_csv(filePath, index_col=0)
                data = data.astype(int)
                countedData = data
            else:
                data = pd.read_csv(filePath, index_col=0)
                data = data.astype(int)
                countedData += data

            # Format values to have commas
            formattedCounts = data.to_string(
                formatters={column: '{:,.0f}'.format for column in
                            data.select_dtypes(include='number').columns})

            # Print: Counts
            substrateCounts = sum(data.iloc[:, 0])
            totalCounts += substrateCounts
            print(f'Counts: {greenLightB}{fileName}{resetColor}\n'
                  f'{formattedCounts}\n'
                  f'Substrate Count: {red}'
                  f'{substrateCounts:,}{resetColor}\n\n')

        # Drop columns
        if dropColumn:
            countedData = self.dropColumnsFromMatrix(countMatrix=countedData,
                                                     datasetType=fileType,
                                                     dropColumn=dropColumn)

        # Sum each column
        columnSums = pd.DataFrame(np.sum(countedData, axis=0), columns=['Total Counts'])
        columnSumsFormat = columnSums.apply(lambda x: x.map('{:,}'.format))
        print(f'{columnSumsFormat}')

        # Sanity Check: Do the sums of each column match the total number of substrates?
        for indexColumn, columnSum in enumerate(columnSums.iloc[:, 0]):
            if columnSum != totalCounts:
                columnSums = columnSums.apply(lambda x: x.map('{:,}'.format))
                print(f'{orange}ERROR: The total number of substrates '
                      f'({cyan}{totalCounts:,}{orange}) =/= '
                      f'the sum of column {pink}{columnSums.index[indexColumn]}{orange} '
                      f'({cyan}{columnSum:,}{orange})\n')
                sys.exit(1)
        print('\n')

        return countedData, totalCounts



    def loadSubstrates(self, fileNames, fileType, dropColumn=False):
        print('============================= Load: Substrate Files '
              '=============================')
        substrates = {}
        substrateTotal = 0

        print(f'Loading dataset: {purple}{self.enzymeName} {fileType}{resetColor}')
        for fileName in fileNames:
            print(f'     {greenLightB}{fileName}{resetColor}')
        print()

        # Function to load each file
        def loadFile(fileName):
            fileLocation = os.path.join(self.pathData, f'substrates_{fileName}.json')
            print(f'File path:\n     {greenDark}{fileLocation}{resetColor}\n')
            with open(fileLocation, 'r') as openedFile:  # Open file
                data = json.load(openedFile) # Access the data
                dataTotalSubs = sum(data.values())
                print(f'     Total substrates in {greenLightB}{fileName}{resetColor}: '
                      f'{red}{dataTotalSubs:,}{resetColor}\n')

                # Combine the loaded dictionary into the main substrates
                nonlocal substrates
                for key, value in data.items():
                    if key in substrates:
                        substrates[key] += value
                    else:
                        substrates[key] = value

                nonlocal substrateTotal
                substrateTotal += dataTotalSubs

        threads = []
        for fileName in fileNames:
            thread = threading.Thread(target=loadFile, args=(fileName,))
            threads.append(thread)
            thread.start()

        # Wait for all threads to complete
        for thread in threads:
            thread.join()

        # Sort loaded data
        substrates = dict(sorted(
            substrates.items(), key=lambda x: x[1], reverse=True)
        )

        # Print: Loaded data substrates
        print(f'Loaded data: {purple}{fileType}{resetColor}')
        iteration = 0
        for substrate, count in substrates.items():
            iteration += 1
            print(f'     {pink}{substrate}{resetColor}, Counts: {red}{count:,}'
                  f'{resetColor}')
            if iteration >= self.printNumber:
                break
        print(f'\nTotal substrates: {purple}{fileType}\n'
              f'     {red} {substrateTotal:,}{resetColor}\n\n')

        if dropColumn:
            if isinstance(dropColumn, list) and dropColumn[0]:
                dropAA = True
            elif isinstance(dropColumn, str):
                dropAA = True
            else:
                dropAA = False

            if dropAA:
                substrates = self.truncateSubs(substrates=substrates,
                                               dropColumn=dropColumn)

        return substrates, substrateTotal



    def loadUnfilteredSubs(self, loadInitial=False, loadFinal=False,
                           dropColumn=False):
        def loadSubsThread(fileNames, fileType, result):
            subsLoaded, totalSubs = self.loadSubstrates(fileNames=fileNames,
                                                        fileType=fileType,
                                                        dropColumn=dropColumn)
            result[fileType] = (subsLoaded, totalSubs)


        # Initialize result dictionary
        loadedResults = {}
        if loadInitial:
            # Create threads for loading initial and final substrates
            threadInitial = threading.Thread(target=loadSubsThread,
                                             args=(self.filesInit, 'Initial Sort',
                                                   loadedResults))
            threadInitial.start() # Start the threads
            threadInitial.join() # Wait for the threads to complete

            # Retrieve the loaded substrates
            substratesInitial, totalSubsInitial = loadedResults['Initial Sort']

            return substratesInitial, totalSubsInitial
        elif loadFinal:
            # Create thread for the final substrates
            threadFinal = threading.Thread(target=loadSubsThread,
                                           args=(self.filesFinal, 'Final Sort',
                                                 loadedResults))
            threadFinal.start() # Start the thread
            threadFinal.join() # Wait for the thread to complete

            # Retrieve the loaded substrates
            substrates, totalSubs = loadedResults['Final Sort']

            return substrates, totalSubs

        else:
            # Create threads for loading initial and final substrates
            threadInitial = threading.Thread(target=loadSubsThread,
                                             args=(self.filesInit, 'Initial Sort',
                                                   loadedResults))
            threadFinal = threading.Thread(target=loadSubsThread,
                                           args=(self.filesFinal, 'Final Sort',
                                                 loadedResults))
            threadInitial.start()  # Start the threads
            threadFinal.start()
            threadInitial.join() # Wait for the threads to complete
            threadFinal.join()

            # Retrieve the loaded substrates
            substratesInitial, totalSubsInitial = loadedResults['Initial Sort']
            substratesFinal, totalSubsFinal = loadedResults['Final Sort']

            return substratesInitial, totalSubsInitial, substratesFinal, totalSubsFinal



    def loadSubstratesFiltered(self):
        print('=========================== Load: Filtered Substrate '
              '============================')
        print(f'Loading substrates: {purple}{self.enzymeName} Fixed {self.datasetTag}\n'
              f'     {greenDark}{self.pathFilteredSubs}{resetColor}\n\n')
        with open(self.pathFilteredSubs, 'r') as file:
            substrates = json.load(file)

        iteration = 0
        print(f'Substrates:')
        for substrate, count in substrates.items():
            print(f'     {pink}{substrate}{resetColor}, Count:{red} {count:,}'
                  f'{resetColor}')
            iteration += 1
            if iteration >= self.printNumber:
                print('')
                break

        totalSubsFinal = sum(substrates.values())
        print(f'Total substrates: {red}{totalSubsFinal:,}{resetColor}\n'
              f'Unique Substrates: {red}{len(substrates.keys()):,}{resetColor}\n\n')

        return substrates, totalSubsFinal



    def loadMotifCounts(self, motifLabel, motifIndex,
                        sortType='FinalSort', returnList=False,
                        loadCountsRel=True, dropColumn=False):
        print('================================ Combine Motifs '
              '=================================')
        initialMotifFrame = self.xAxisLabels[motifIndex[0]:motifIndex[1]]
        print('Combine Matrices: Substrate Profile')
        print(f'Datasets: {purple}{self.datasetTag}{resetColor}\n'
              f'Motif Label: {blue}{", ".join(motifLabel)}{resetColor}\n'
              f'Motif Frame: {blue}{", ".join(initialMotifFrame)}{resetColor}\n\n')

        frameLength = len(motifLabel)
        countsMotifsAll = []
        totalCountsFixedFrame = None


        # Define: File paths
        paths = self.getFilePathCombined(sortType=sortType, loadCountsRel=loadCountsRel)

        # Load the counts
        for index, pathFixedMotifRelCounts in enumerate(paths):
            # Look for the file
            if os.path.exists(pathFixedMotifRelCounts):
                print(f'Loading ({index}): {greenLightB}Substrate Profile\n{greenDark}'
                      f'     {pathFixedMotifRelCounts}{resetColor}\n')

                # Load file
                countsLoaded = pd.read_csv(pathFixedMotifRelCounts, index_col=0)

                if dropColumn and dropColumn in countsLoaded.columns:
                    countsLoaded = self.dropColumnsFromMatrix(
                        countMatrix=countsLoaded, datasetType=fileType,
                        dropColumn=dropColumn
                    )

                # Define motif positions and extract the sequence
                startPosition = motifIndex[0]
                startSubPrevious = startPosition
                if index != 0:
                    # Evaluate previous motif index
                    # print(f'Pos ({index}): {self.fixedPos}')
                    if isinstance(self.fixedPos[0], list):
                        frame = self.fixedPos[0]
                        # print(f'Frame: {frame}')
                        pos = frame[index]
                        pos2 = frame[index - 1]
                        fixedPosDifference = pos - pos2
                    else:
                        pos = self.fixedPos[index]
                        pos2 = self.fixedPos[index - 1]
                        fixedPosDifference = pos - pos2
                    # print(f'  Pos: {pos} - {pos2}')
                    # print(f' Diff: {fixedPosDifference}\n')

                    # Define: Frame indices
                    startSubPrevious += fixedPosDifference
                    startSub = index + startSubPrevious - 1
                    endSub = startSub + frameLength
                else:
                    startSub = startPosition
                    endSub = motifIndex[-1]

                fixedFramePos = countsLoaded.columns[startSub:endSub]
                countsFixedFrame = countsLoaded.loc[:, fixedFramePos]
                countsFixedFrame.columns = motifLabel
                countsMotifsAll.append(countsFixedFrame.values)

                formattedCounts = countsFixedFrame.to_string(
                    formatters={column: '{:,.0f}'.format for column in
                                countsFixedFrame.select_dtypes(include='number').columns})
                print(f'Selecting Positions: {purple}Fixed Motif {self.motifTags[index]}'
                      f'{resetColor}\n'
                      f'Counts: {blue}{fixedFramePos}{resetColor}\n'
                      f'{formattedCounts}\n\n')

                # Track totals
                if index == 0:
                    totalCountsFixedFrame = countsFixedFrame
                else:
                    totalCountsFixedFrame += countsFixedFrame
            else:
                print(f'{orange}ERROR: The file was not found\n'
                      f'     {pathFixedMotifRelCounts}\n')
                sys.exit(1)

        # Sum the columns
        posSumsCombinedMotif = []
        for column in motifLabel:
            posSumsCombinedMotif.append(sum(totalCountsFixedFrame.loc[:, column]))

        # Format the DataFrame with commas
        formattedCounts = totalCountsFixedFrame.to_string(
            formatters={column: '{:,.0f}'.format for column in
                        totalCountsFixedFrame.select_dtypes(include='number').columns})

        # Print: Counts
        print(f'{pink}Combined Counts{resetColor}: {purple}{self.enzymeName} - '
              f'{self.datasetTag}{resetColor}\n{formattedCounts}\n')
        print('Total Counts:')
        for index, position in enumerate(motifLabel):
            print(f'     {position}: {red}{posSumsCombinedMotif[index]:,}{resetColor}')
        print('\n')

        if returnList:
            # Convert a list into 3D array
            countsMotifsAll = np.stack(countsMotifsAll, axis=0)

            return countsMotifsAll, totalCountsFixedFrame, posSumsCombinedMotif
        else:
            return totalCountsFixedFrame, posSumsCombinedMotif



    def loadMotifSeqs(self, motifLabel, motifIndex,
                      loadCountsRel=True, sortType='FinalSort'):
        print('============================ Load: Substrate Motifs '
              '=============================')
        print(f'Dataset: {purple}{self.datasetTag}{resetColor}\n'
              f'Motif Label: {blue}{", ".join(motifLabel)}{resetColor}')
        frameLength = len(motifLabel)
        totalMotifs = 0
        motifs = {}
        substrates = {}
        motifTag = None
        foundMotifs = {}

        # Assign: Motif parameters
        self.motifIndex = motifIndex
        print(f'Motif Index: {blue}{self.motifIndex}{resetColor}\n\n')


        # Define: File paths
        paths = self.getFilePathCombined(sortType=sortType, loadSubs=True,
                                         loadCountsRel=loadCountsRel)

        # Load the substrates
        for index, pathFixedMotifSubs in enumerate(paths):
            # Look for the file
            if os.path.exists(pathFixedMotifSubs):
                print(f'Loading ({index}): {greenLightB}Filtered Substrates\n{greenDark}'
                      f'     {pathFixedMotifSubs}{resetColor}\n')

                # Load file
                with open(pathFixedMotifSubs, 'r') as file:
                    loadedSubs = json.load(file)
                    if substrates == {}:
                        substrates = loadedSubs
                    else:
                        for substrate, count in loadedSubs.items():
                            if substrate not in substrates.keys():
                            #     substrates[substrate] += count
                            # else:
                                substrates[substrate] = count

                # Define motif positions and extract the sequence
                startPosition = motifIndex[0]
                startSubPrevious = startPosition
                if index != 0:
                    # Evaluate previous motif index
                    # print(f'Pos ({index}): {self.fixedPos}')
                    if isinstance(self.fixedPos[0], list):
                        frame = self.fixedPos[0]
                        # print(f'Frame: {frame}')
                        pos = frame[index]
                        pos2 = frame[index - 1]
                        fixedPosDifference = pos - pos2
                    else:
                        pos = self.fixedPos[index]
                        pos2 = self.fixedPos[index - 1]
                        fixedPosDifference = pos - pos2
                    # print(f'  Pos: {pos} - {pos2}')
                    # print(f' Diff: {fixedPosDifference}\n')

                    # Define: Frame indices
                    startSubPrevious += fixedPosDifference
                    startSub = index + startSubPrevious - 1
                    endSub = startSub + frameLength
                else:
                    startSub = startPosition
                    endSub = motifIndex[-1]

                # Print: Loaded data
                iteration = 0
                print(f'* Motif Indices: {blue}{self.xAxisLabels[startSub]}-'
                      f'{self.xAxisLabels[endSub - 1]}{resetColor}')
                for substrate, count in loadedSubs.items():
                    print(f'    {pink}{substrate}{resetColor}, '
                          f'{blue}{substrate[startSub:endSub]}{resetColor}, '
                          f'Counts: {red}{count:,}'
                          f'{resetColor}')
                    iteration += 1
                    if iteration >= self.printNumber:
                        break

                # Record motifs
                totalCounts = 0
                for substrate, count in loadedSubs.items():
                    if substrate in foundMotifs.keys():
                        continue # Don't double count substrate
                    else:
                        foundMotifs[substrate] = True
                    totalMotifs += count
                    totalCounts += count
                    motif = substrate[startSub:endSub]
                    if motif in motifs.keys():
                        motifs[motif] += count
                    else:
                        motifs[motif] = count
                print(f'* Total Counts: {red}{totalCounts:,}{resetColor}\n\n')
            else:
                print(f'{orange}ERROR: The file was not found\n'
                      f'     {pathFixedMotifSubs}\n')
                sys.exit(1)

            self.motifIndexExtracted.append((startSub, endSub))

        # Evaluate: Loaded data
        self.motifLen = len(next(iter(motifs)))

        # Sort the substrate dictionary by counts
        substrates = dict(sorted(
            substrates.items(), key=lambda x: x[1], reverse=True)
        )
        motifs = dict(sorted(
            motifs.items(), key=lambda x: x[1], reverse=True)
        )

        iteration = 0
        print(f'Top Substrates:')
        for substrate, count, in substrates.items():
            print(f'     {pink}{substrate}{resetColor}, '
                  f'Counts: {red}{count:,}{resetColor}')
            iteration += 1
            if iteration >= self.printNumber:
                break
        iteration = 0
        print(f'\nTop Motifs:')
        for motif, count, in motifs.items():
            print(f'     {blue}{motif}{resetColor}, Counts: {red}{count:,}{resetColor}')
            iteration += 1
            if iteration >= self.printNumber:
                break
        print(f'\nTotal Motifs: {red}{totalMotifs:,}{resetColor}\n'
              f'Unique Motifs: {red}{len(motifs.keys()):,}'
              f'{resetColor}\n\n')

        return motifs, totalMotifs, substrates



    def dropColumnsFromMatrix(self, countMatrix, datasetType, dropColumn):
        if isinstance(dropColumn, str):
            dropColumn = [dropColumn]
        if dropColumn[0] not in countMatrix.columns:
            return countMatrix
        print('================================ Truncate Matrix '
              '================================')
        idx = countMatrix.columns.get_loc(dropColumn[0])
        dropPos =  list(countMatrix.columns[idx:])
        countMatrix = countMatrix.iloc[:,0:idx]
        if len(self.xAxisLabels) > idx:
            # print(self.xAxisLabels[:idx])
            self.xAxisLabels = self.xAxisLabels[:idx]
        if self.xAxisLabelsMotif is not None and len(self.xAxisLabelsMotif) > idx:
            # print(self.xAxisLabelsMotif[:idx])
            self.xAxisLabelsMotif = self.xAxisLabelsMotif[:idx]

        countsFormatted = countMatrix.to_string(
            formatters={column: '{:,.0f}'.format for column in
                        countMatrix.select_dtypes(include='number').columns})
        print(f'Dropping: {blue}{dropPos}{resetColor}')
        print(f'Counts: {purple}{datasetType}{resetColor}\n{countsFormatted}\n\n')

        return countMatrix



    def truncateSubs(self, substrates, dropColumn):
        if isinstance(dropColumn, str):
            dropColumn = [dropColumn]
        idx = int(dropColumn[0].replace('R', '')) - 1
        dropPos = list(substrates.keys())
        subLen = len(next(iter(substrates)))
        if subLen <= idx:
            return substrates
        print('============================== Truncate Substrates '
              '==============================')
        print(f'Dropping: {blue}{' '.join(dropColumn)}{resetColor}\n'
              f'Substrate len: {pink}{subLen}{resetColor}\n'
              f'Truncated len: {blue}{idx}{resetColor}\n')
        subs = {}
        for substrate, count in substrates.items():
            sub = substrate[0:idx]
            if sub in subs.keys():
                subs[sub] += count
            else:
                subs[sub] = count
        subs = dict(sorted(subs.items(), key=lambda x: x[1], reverse=True))
        print('Substrates:')
        for i, (sub, count) in enumerate(subs.items()):
            print(f'    {blue}{sub}{resetColor}: Counts: {red}{count:,}{resetColor}')
            if i >= self.printNumber:
                break
        print('\n')

        if len(self.xAxisLabels) > idx:
            print(self.xAxisLabels[:idx])
            self.xAxisLabels = self.xAxisLabels[:idx]
        if self.xAxisLabelsMotif is not None and len(self.xAxisLabelsMotif) > idx:
            print(self.xAxisLabelsMotif[:idx])
            self.xAxisLabelsMotif = self.xAxisLabelsMotif[:idx]

        return subs



    def saveSubstrateCSV(self, seqs, initialRF, finalRF, minCounts=100,
                         seqsBg=False, excludeAA='', maxCountsBg=5,
                         mod=25, modScale=False, combinedMotifs=False,
                         chopSeq=False):
        print('============================== Save Substrate CSV '
              '===============================')
        print(f'Minimum substrate count: {red}{minCounts:,}{resetColor}')
        excValid = False
        predScores = True


        # Chop off the C-terminal AAs
        print(f'Chop: {chopSeq}')
        subsCounts = {}  # Limit substrates by counts
        if chopSeq:
            sub = next(iter(seqs))
            subLen = len(sub)
            print(f'Cutting substrates from {pink}{subLen}AA{resetColor} to '
                  f'{blue}{chopSeq}AA{resetColor}')
            x = subLen - chopSeq
            if self.filterSubs:
                print(f'  {pink}{sub.replace(
                    sub[0:chopSeq], f'{blue}{sub[0:chopSeq]}{pink}')}'
                      f'{resetColor}')
                for seq, counts in seqs.items():
                    if counts >= minCounts:
                        subsCounts[seq[0:chopSeq]] = counts
            else:
                for i in range(x+1):
                    print(f'  {pink}{sub.replace(
                        sub[i:chopSeq+i], f'{blue}{sub[i:chopSeq+i]}{pink}')}'
                          f'{resetColor}')
                    for seq, counts in seqs.items():
                        if counts >= minCounts:
                            subsCounts[seq[i:chopSeq+i]] = counts
        else:
            for seq, counts in seqs.items():
                if counts >= minCounts:
                    subsCounts[seq] = counts
        N = len(subsCounts.keys())
        if N == 0:
            print(f'{orange}ERROR: There were no substrates with '
                  f'{cyan}counts{orange} >= {cyan}{minCounts}{resetColor}\n\n')
            return
        subLen = len(next(iter(subsCounts)))
        subsCounts = dict(
            sorted(subsCounts.items(), key=lambda x: x[1], reverse=True)
        )
        print()

        # Add bg substrates
        tagMod = ''
        if seqsBg:
            s = {}
            print(f'N: {len(subsCounts.keys())}\nN: {len(seqsBg.keys())}')
            if isinstance(chopSeq, bool) or chopSeq > 1:
                for seq, count in seqsBg.items():
                    if (count <= maxCountsBg and seq not in subsCounts.keys()
                            and excludeAA not in seq):
                        s[seq] = count
                print(f'Chop: {chopSeq}, {isinstance(chopSeq, int)}')
            else:
                for seq, count in seqsBg.items():
                    if (count <= maxCountsBg and seq not in subsCounts.keys()
                            and excludeAA not in seq):
                        seq = seq[0:chopSeq]
                        s[seq] = count
            bg, Nexp, i = {}, N, 0
            print(f'Add Background substrates:')
            if modScale:
                tagMod = f'{mod}'
                mag = np.ceil(np.log10(len(s.keys())) - 1)
                divisor = int(5*10**(mag - 1))
                if divisor == 0:
                    divisor = 1
                print('D:', divisor)
                for seq, count in s.items():
                    i += 1
                    if i % divisor == 0:
                        # print(f'i: {red}{i}{resetColor} / {red}{f}{resetColor} = '
                        #       f'{i // divisor}')
                        mod += mod / (i // divisor)
                        # print(f'Mod: {blue}{mod:,}{resetColor}')
                    if i % mod == 0:
                        bg[seq] = count
            else:
                print('Mod:', mod)
                print(f'Seqs: {len(s.keys())}')
                for seq, count in s.items():
                        i += 1
                        if i % mod == 0:
                            bg[seq] = count
            bg = dict(sorted(bg.items(), key=lambda x: x[1], reverse=True))
            for i, (seq, count) in enumerate(bg.items()):
                if i >= self.printNumber:
                    break
                print(f'    {pink}{seq}{resetColor}: Count {red}{count:,}{resetColor}')
            for seq, count in bg.items():
                subsCounts[seq] = count
            N = len(subsCounts)
            print(f'\nN Experimental Substrates: {red}{Nexp:,}{resetColor}\n'
                  f'N Background Substrates:   {red}{len(bg.keys()):,}{resetColor}\n'
                  f'Total Unique Substrates:   {red}{N:,}{resetColor}')
        else:
            print(f'Total Unique Substrates: {red}{N:,}{resetColor}')
        print(f'Total Substrates: {red}{sum(subsCounts.values()):,}{resetColor}\n')

        # Normalize counts
        subsCountsNorm = {}
        maxCount = max(list(subsCounts.values()))
        for seq, count in subsCounts.items():
            subsCountsNorm[seq] = count / maxCount

        # Get data paths
        t = 'Counts'
        tag = (f'{self.enzyme}_{t}_{self.datasetTag}_{subLen}AA_'
               f'MinCounts{minCounts}.csv').replace(' ', '-')
        if seqsBg:
            if tagMod:
                tag = tag.replace('.csv',
                                  f'_MinBgCounts{maxCountsBg}_scaledMod{tagMod}.csv')
            else:
                tag = tag.replace('.csv',
                                  f'_MinBgCounts{maxCountsBg}_mod{mod}.csv')
        if excValid:
            tag = tag.replace(f'_{t}', f'_{t}_ExtValid')
        pathCSV = os.path.join(self.pathFolder, 'CSV')
        if not os.path.exists(pathCSV):
            os.makedirs(pathCSV, exist_ok=True)
        pathCSVFig = os.path.join(pathCSV, 'Figures')
        if not os.path.exists(pathCSVFig):
            os.makedirs(pathCSVFig, exist_ok=True)
        paths = [
            os.path.join(pathCSV, tag),
            os.path.join(pathCSV, tag.replace(f'_{t}', '_CountsNorm')),
            os.path.join(pathCSV, tag.replace(f'_{t}', '_ZCounts')),
        ]
        if predScores:
            paths.append(
                os.path.join(pathCSV, tag.replace(f'_{t}', '_ZPred'))
            )

        # Define: Save location
        figLabel = (f'{self.enzyme}-Bars-Counts-{self.datasetTag}-{subLen}AA-'
                    f'N_{N}-MinCounts_{minCounts}.png')
        if excValid:
            figLabel = figLabel.replace('Counts','CountsExtValid')
        if combinedMotifs:
            figLabel = figLabel.replace(self.datasetTag,
                                        f'SubstrateProfile_{self.datasetTag}.png')
        if seqsBg:
            if tagMod:
                figLabel = figLabel.replace(
                    '.png', f'-MinBgCounts_{maxCountsBg}-scaledMod_{tagMod}.png'
                )
            else:
                figLabel = figLabel.replace(
                    '.png', f'-MinBgCounts_{maxCountsBg}-mod_{mod}.png'
                )
        pathFigs = [
            os.path.join(pathCSVFig, figLabel),
            os.path.join(pathCSVFig, figLabel.replace('-Counts', '-CountsNorm')),
            os.path.join(pathCSVFig, figLabel.replace('-Counts', '-ZCounts')),
            os.path.join(pathCSVFig, figLabel.replace('-Counts', '-Pred')),
            os.path.join(pathCSVFig, figLabel.replace('-Counts', '-ZPred'))
        ]


        # Get prediction matrix
        matrix = self.normalizeProbRatios(
            finalRF=finalRF, initialRF=initialRF, pData=False
        )
        if chopSeq:
            dropCol = []
            for i in range(len(matrix.columns)):
                if i >= chopSeq:
                    dropCol.append(matrix.columns[i])
            matrix = matrix.drop(dropCol, axis=1)

        # Calculate: Z-scores Counts
        subsZCounts = {}
        mu = np.average(list(subsCounts.values()))
        sigma = np.std(list(subsCounts.values()))
        print(f'Z-Scores: Counts\n'
              f'* µ: {red}{round(mu, 3):,}{resetColor}\n'
              f'* σ: {red}{round(sigma, 3):,}{resetColor}\n')
        for seq, count in subsCounts.items():
            subsZCounts[seq] = (count - mu) / sigma

        # Predict substrate activity
        subsZ = {}
        subsZPred = {}
        if predScores:
            if subLen <= len(matrix.columns):
                print(f'Scoring Matrix: {purple}{self.datasetTag}{resetColor}\n'
                      f'{matrix}\n')
                for seq, count in subsCounts.items():
                    subsZ[seq] = self.scoreSubstrate(seq, matrix)

                # Normalize pred scores
                maxVal = max(subsZ.values())
                for seq, score in subsZ.items():
                    subsZ[seq] = score / maxVal

                # Calculate: Z-scores Predicted
                mu = np.average(list(subsZ.values()))
                sigma = np.std(list(subsZ.values()))
                print(f'Z-Scores: Predicted Activity\n'
                      f'* µ: {red}{round(mu, 3):,}{resetColor}\n'
                      f'* σ: {red}{round(sigma, 3):,}{resetColor}\n')
                for seq, score in subsZ.items():
                    subsZPred[seq] = (score - mu) / sigma
            else:
                print(f'{orange}Activity scores were not predicted:\n'
                      f'* Substrate length: {cyan}{subLen}{orange}\n'
                      f'* N matrix columns: {cyan}{len(matrix.columns)}{resetColor}\n')

        # Pre-format values
        pSeqs = list(subsCounts)[:self.printNumber]
        frmtCounts = [f'{subsCounts[seq]:,}' for seq in pSeqs]
        frmtZCounts = [f'{subsZCounts[seq]:.3f}' for seq in pSeqs]
        pSeqsLow = list(subsCounts)[-self.printNumber:]
        frmtCountsLow = [f'{subsCounts[seq]:,}' for seq in pSeqsLow]
        frmtZCountsLow = [f'{subsZCounts[seq]:.3f}' for seq in pSeqsLow]

        # Compute column widths
        widthCount = max(len(val) for val in frmtCounts)
        widthZCount = max(
            max(len(val) for val in frmtZCounts), max(len(val) for val in frmtZCountsLow)
        )

        # Print data
        print('Substrate Scores:')
        if subsZPred:
            frmtPred = [f'{subsZ[seq]:.3f}' for seq in pSeqs]
            frmtZPred = [f'{subsZPred[seq]:,.3f}' for seq in pSeqs]
            frmtPredLow = [f'{subsZ[seq]:.3f}' for seq in pSeqsLow]
            frmtZPredLow = [f'{subsZPred[seq]:,.3f}' for seq in pSeqsLow]
            widthProd = max(len(val) for val in frmtPred)
            widthZ = max(max(len(val) for val in frmtZPred),
                         max(len(val) for val in frmtZPredLow))
            for seq, countStr, zCountStr, prodStr, zStr in zip(
                    pSeqs, frmtCounts, frmtZCounts, frmtPred, frmtZPred):
                print(
                    f'    {pink}{seq}{resetColor}, '
                    f'Count: {red}{countStr:>{widthCount}}{resetColor}, '
                    f'Z-Score Counts: {red}{zCountStr:>{widthZCount}}{resetColor}, '
                    f'Pred: {red}{prodStr:>{widthProd}}{resetColor}, '
                    f'Z-Score Pred: {red}{zStr:>{widthZ}}{resetColor}'
                )
            print('    ...')
            for seq, countStr, zCountStr, prodStr, zStr in zip(
                    pSeqsLow, frmtCountsLow, frmtZCountsLow, frmtPredLow, frmtZPredLow):
                print(
                    f'    {pink}{seq}{resetColor}, '
                    f'Count: {red}{countStr:>{widthCount}}{resetColor}, '
                    f'Z-Score Counts: {red}{zCountStr:>{widthZCount}}{resetColor}, '
                    f'Pred: {red}{prodStr:>{widthProd}}{resetColor}, '
                    f'Z-Score Pred: {red}{zStr:>{widthZ}}{resetColor}'
                )
        else:
            for seq, countStr, zCountStr in zip(pSeqs, frmtCounts, frmtZCounts):
                print(f'    {pink}{seq}{resetColor}, '
                      f'Count: {red}{countStr:>{widthCount}}{resetColor}, '
                      f'Z Score Counts: {red}{zCountStr:>{widthZCount}}{resetColor}'
                )
            print('    ...')
            for seq, countStr, zCountStr in zip(
                    pSeqsLow, frmtCountsLow, frmtZCountsLow):
                print(f'    {pink}{seq}{resetColor}, '
                      f'Count: {red}{countStr:>{widthCount}}{resetColor}, '
                      f'Z Score Counts: {red}{zCountStr:>{widthZCount}}{resetColor}\n'
                )
        plotData = [False for _ in range(len(paths))]
        savedData = False
        print()


        # CSV: Scores
        savePath = paths[0]
        if not os.path.exists(savePath):
            if not savedData:
                savedData = True
                print(f'Saving {red}{N:,}{resetColor} substrates in a CSV file\n'
                      f'Save path:')
            print( f'    {greenDark}{savePath}{resetColor}')
            plotData[0] = True
            with open(savePath, 'w', newline='') as c:
                writer = csv.writer(c)
                writer.writerow(['sequence', self.enzyme])
                for seq, count in subsCounts.items():
                    writer.writerow([seq, count])

        # CSV: Normalized scores
        savePath = paths[1]
        if not os.path.exists(savePath):
            if not savedData:
                savedData = True
                print(f'Saving {red}{N:,}{resetColor} substrates in a CSV file\n'
                      f'Save path:')
            print(f'    {greenDark}{savePath}{resetColor}')
            plotData[1] = True
            with open(savePath, 'w', newline='') as c:
                writer = csv.writer(c)
                writer.writerow(['sequence', self.enzyme])
                for seq, value in subsCountsNorm.items():
                    writer.writerow([seq, value])

        # CSV: Z Scores
        savePath = paths[2]
        if not os.path.exists(savePath):
            if not savedData:
                savedData = True
                print(f'Saving {red}{N:,}{resetColor} substrates in a CSV file\n'
                      f'Save path:')
            print(f'    {greenDark}{savePath}{resetColor}')
            plotData[2] = True
            with open(savePath, 'w', newline='') as c:
                writer = csv.writer(c)
                writer.writerow(['sequence', self.enzyme])
                for seq, score in subsZCounts.items():
                    writer.writerow([seq, score])

        # CSV: Z-Pred
        if len(paths) > 3 and subsZPred and not os.path.exists(paths[3]):
            savePath = paths[3]
            if not savedData:
                savedData = True
                print(f'Saving {red}{N:,}{resetColor} substrates in a CSV file\n'
                      f'Save path:')
            print(f'    {greenDark}{savePath}{resetColor}')
            plotData[3] = True
            with open(savePath, 'w', newline='') as c:
                writer = csv.writer(c)
                writer.writerow(['sequence', self.enzyme])
                for seq, score in subsZPred.items():
                    writer.writerow([seq, f'{score:.2f}'])
        if not savedData:
            print(f'\nPrint no data was saved. Files already saved at:')
            for path in paths:
                print(f'    {greenDark}{path}{resetColor}')
        print()


        # Plot bar graphs
        self.plotBarGraphCSV(
            substrates=subsCounts, dataType='Counts',
            combinedMotifs=combinedMotifs, minCounts=minCounts,
            saveLocation=pathFigs[0]
        )
        self.plotBarGraphCSV(
            substrates=subsCountsNorm, dataType='Normalized Counts',
            combinedMotifs=combinedMotifs, minCounts=minCounts,
            saveLocation=pathFigs[1]
        )
        self.plotBarGraphCSV(
            substrates=subsZCounts, dataType='Z Counts',
            combinedMotifs=combinedMotifs, minCounts=minCounts,
            saveLocation=pathFigs[2]
        )
        if subsZ:
            self.plotBarGraphCSV(
                substrates=subsZ, dataType='Pred',
                combinedMotifs=combinedMotifs, minCounts=minCounts,
                saveLocation=pathFigs[3]
            )
        if subsZPred:
            self.plotBarGraphCSV(
                substrates=subsZPred, dataType='Z Pred',
                combinedMotifs=combinedMotifs, minCounts=minCounts,
                saveLocation=pathFigs[4]
            )



    def plotBarGraphCSV(self, substrates, dataType,
                        barColor='#BF5700', barWidth=3,
                        combinedMotifs=False, minCounts=False,
                        saveLocation=False):
        print('================================ Plot: Bar Graph '
              '================================')
        if minCounts:
            print(f'Minimum substrate count: {red}{minCounts:,}{resetColor}')
        print(f'Substrates:')

        # Evaluate data
        NSubs = len(substrates.keys())
        iteration = 0
        if ('counts' in dataType.lower() and 'z ' not in dataType.lower()
                and 'norm' not in dataType.lower()):
            for substrate, count in substrates.items():
                print(f'     {blue}{substrate}{resetColor}, '
                      f'Counts: {red}{count:,}{resetColor}')
                iteration += 1
                if iteration >= self.printNumber:
                    break
            countTotal = 0
            for count in substrates.values():
                countTotal += count
            print(f'Total Counts: {red}{countTotal:,}{resetColor}')
        else:
            for substrate, score in substrates.items():
                print(f'     {blue}{substrate}{resetColor}, '
                      f'Score: {red}{np.round(score, self.roundVal):,}{resetColor}')
                iteration += 1
                if iteration >= self.printNumber:
                    break
            scoreTotal = 0.0
            for score in substrates.values():
                scoreTotal += score
            scoreTotal = np.round(scoreTotal, self.roundVal)
            print(f'Total Values: {red}{scoreTotal:,}{resetColor}')
        print(f'Number of plotted sequences: {red}{NSubs:,}{resetColor}\n')

        # Collect substrates
        x = []
        y = []
        for substrate, count in substrates.items():
            x.append(str(substrate))
            y.append(count)
        N, val = len(x), 5000
        if N < 200:
            val = 50
            barWidth = 0.5
            while val < N:
                barWidth += 0.5
                val += 50
        else:
            while val < N:
                barWidth += 1
                val += 5000

        # Evaluate: Y axis
        maxValue = math.ceil(max(y))
        minValue = math.floor(min(y))
        yTicks = []
        if maxValue == 1:
            step = 0.1
            yMax = math.ceil(max(y)) + step
            yMin = math.floor(min(y))
            yTicks = np.arange(yMin, yMax+(step/2), 0.2)
            # print(f'Max Value: {blue}{maxValue:,}{resetColor}\n'
            #       f'Y Max: {blue}{yMax:,}{resetColor}\n'
            #       f'Y Min: {blue}{yMin:,}{resetColor}\n'
            #       f'Y Ticks: {pink}{yTicks}{resetColor}\n')
        else:
            if minValue < 0:
                diff = (maxValue - minValue)
                step = diff / 10
                ymax, ymin = max(y), min(y)
                if ymax + step / 2 > maxValue:
                    ymax = maxValue
                if ymin - step / 2 < minValue:
                    ymin = minValue

                yMax, yMin = 0, 0
                while yMax <= ymax:
                    yMax += step
                    # print(ymax, yMax)
                while yMin >= ymin:
                    yMin -= step
                yTicks = np.arange(yMin, yMax + step, step)
                # print(f'\nDiff: {blue}{diff:,}{resetColor}\n'
                #       f'Max Value: {blue}{maxValue:,}{resetColor}\n'
                #       f'Min Value: {blue}{minValue:,}{resetColor}\n'
                #       f'Step: {blue}{step:,}{resetColor}\n'
                #       f'Y Max: {blue}{yMax:,}{resetColor}\n'
                #       f'Y Min: {blue}{yMin:,}{resetColor}\n'
                #       f'Y Ticks: {pink}{yTicks}{resetColor}\n')
            else:
                magnitude = math.floor(math.log10(maxValue))
                step = 10 ** (magnitude) / 5
                yMax = 0
                while yMax <= maxValue:
                    yMax += step
                    yTicks = np.arange(0, yMax + 1, step, dtype=int)
                yMin = 0
            # print(f'Max Value: {blue}{maxValue:,}{resetColor}\n'
            #       f'Magnitude: {blue}{magnitude:,}{resetColor}\n'
            #       f'Step: {blue}{step:,}{resetColor}\n'
            #       f'Y Max: {blue}{yMax:,}{resetColor}\n'
            #       f'Y Min: {blue}{yMin:,}{resetColor}\n')


        # Evaluate: X axis
        magnitude = math.floor(math.log10(NSubs))
        step = 10 ** (magnitude) / 2
        xMax = 0
        while xMax < NSubs:
            xMax += step
        xTicks = np.arange(0, xMax+1, step, dtype=int)
        # print(f'N Subs: {blue}{NSubs:,}{resetColor}\n'
        #       f'Magnitude: {blue}{magnitude:,}{resetColor}\n'
        #       f'Step: {blue}{step:,}{resetColor}\n'
        #       f'X Max: {blue}{xMax:,}{resetColor}\n'
        #       f'X Labels: {pink}{", ".join([str(x) for x in xTicks])}{resetColor}\n')

        # Define: Figure title
        enzName = self.enzymeName.replace(' - ', '\n')
        title = f'{enzName}\n{self.datasetTag}'
        if minCounts:
            title += f'\nMinimum Substrate Count: {minCounts:,}'
        else:
            title += f'\nTop {NSubs:,} Substrates'


        # Plot the data
        fig, ax = plt.subplots(figsize=self.figSize)
        bars = plt.bar(x, y, color=barColor, width=barWidth)
        plt.xlabel('Substrates', fontsize=self.labelSizeAxis)
        plt.ylabel(dataType, fontsize=self.labelSizeAxis)
        plt.title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        plt.axhline(y=0, color='black', linewidth=self.lineThickness)

        # Set: x ticks
        ax.set_xticks(xTicks)
        ax.set_xticklabels(xTicks, ha='center')
        ax.set_xlim(-step/10, xTicks[-1])

        # Set: y ticks
        plt.ylim(yMin, yMax)
        plt.yticks(yTicks)


        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks, width=self.lineThickness)

        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(self.lineThickness)

        self.plotFig(plt=plt, fig=fig)

        # Save the figure
        if self.saveFigures:
            self.saveFig(fig=fig, path=saveLocation)
        else:
            print()



    def saveData(self, substrates: object, counts: object, sortType: object = 'FinalSort',
                 saveTag: object = None, countsReleased: object = None) -> None:
        if not isinstance(substrates, dict):
            print(f'{orange}ERROR: The substrates need to be stored in a dictionary\n'
                  f'     Current data structure: {type(substrates)}')
            sys.exit(1)
        if not isinstance(counts, pd.DataFrame):
            print(f'{orange}ERROR: The counts need to be stored in a Pandas DataFrame\n'
                  f'     Current data structure: {type(counts)}')
            sys.exit(1)

        # Define: Save path
        filePathCountsReleased = None
        if self.translateDNA:
            folder = os.path.join(self.pathFolder, 'Data')

            if not os.path.exists(folder):
                os.makedirs(folder, exist_ok=True)

            self.pathSaveFigs = os.path.join(self.pathFolder, 'Figures')

            filePathSubs = os.path.join(folder, f'substrates_{saveTag}.json')
            filePathCounts = os.path.join(folder, f'counts_{saveTag}.csv')
        else:
            if self.datasetTag == 'Unfiltered':
                return

            if countsReleased is None:
                filePathSubs, filePathCounts = (
                    self.getFilePath(datasetTag=self.datasetTagSave, sortType=sortType)
                )
            else:
                # Define: File paths
                filePathSubs, filePathCounts, filePathCountsReleased = (
                    self.getFilePath(datasetTag=self.datasetTagSave, sortType=sortType,
                                     motifPath=True)
                )

        if not os.path.exists(filePathSubs) or not os.path.exists(filePathCounts):
            print('================================= Save The data '
                  '=================================')
            if countsReleased is None:
                print(f'Substrate data saved at:\n'
                      f'     {greenDark}{filePathSubs}\n'
                      f'     {filePathCounts}{resetColor}\n\n')

                # Save the substrates
                with open(filePathSubs, 'w') as file:
                    json.dump(substrates, file)

                # Save the counts
                counts.to_csv(filePathCounts)
            else:
                print(f'Substrate data saved at:\n'
                      f'     {greenDark}{filePathSubs}\n'
                      f'     {filePathCounts}\n'
                      f'     {filePathCountsReleased}{resetColor}\n\n')

                # Save the substrates
                with open(filePathSubs, 'w') as file:
                    json.dump(substrates, file)

                # Save the counts
                counts.to_csv(filePathCounts)
                countsReleased.to_csv(filePathCountsReleased)



    def saveFig(self, fig, path):
        if os.path.exists(path):
            print(f'{yellow}The figure was not saved\n\n'
                  f'File was already found at path:\n'
                  f'     {path}{resetColor}\n\n')
        else:
            print(f'Saving figure at path:\n'
                  f'     {greenDark}{path}{resetColor}\n\n')
            fig.savefig(path, dpi=self.figureResolution)



    def saveFigure(self, fig, figType, seqLen, N=False, combinedMotifs=False):
        # Define: Save location
        enzName = self.enzymeName.replace(' - ', '-').replace(' ', '_')
        if self.motifFilter and not self.releasedCounts:
            datasetTag = self.datasetTagMotif.replace(' ', '_')
            figLabel = (f'{enzName}-{figType}-{self.saveFigureIteration}-'
                        f'Register_{datasetTag}-{seqLen}AA-'
                        f'MinCounts_{self.minSubCount}.png')
        elif self.releasedCounts:
            datasetTag = self.datasetTagMotif.replace(' ', '_')
            figLabel = (f'{enzName}-{figType}-SubstrateProfile-{datasetTag}-'
                        f'{seqLen}AA-MinCounts_{self.minSubCount}.png')
        elif combinedMotifs:
            figLabel = (f'{enzName}-{figType}-Combined-{self.datasetTagSave}-'
                        f'{seqLen}AA-MinCounts_{self.minSubCount}.png')
        else:
            figLabel = (f'{enzName}-{figType}-{self.datasetTagSave}-N_{self.nSubsFinal}-'
                        f'{seqLen}AA-MinCounts_{self.minSubCount}.png')
        if '/' in figLabel:
            figLabel = figLabel.replace('/', '_')
        if N:
            figLabel = figLabel.replace(f'N_{self.nSubsFinal}', f'N_{N}')
        saveLocation = os.path.join(self.pathSaveFigs, figLabel)
        saveLocation = saveLocation.replace(' ', '')

        # Save figure
        self.saveFig(fig=fig, path=saveLocation)



    def recordSampleSize(self, NInitial, NFinal, NFinalUnique):
        print('============================== Current Sample Size '
              '==============================')
        self.nSubsFinalUniqueSeqs = NFinalUnique
        if isinstance(NInitial, int) and isinstance(NFinal, int):
            # Update: Sample size
            self.nSubsInitial = NInitial
            self.nSubsFinal = NFinal
            print(f'Dataset: {purple}{self.datasetTag}{resetColor}\n'
                  f'Initial Sort: {red}{self.nSubsInitial:,}{resetColor}\n'
                  f'Final Sort: {red}{self.nSubsFinal:,}{resetColor}')
        else:
            print(f'Sample size was not recorded.\n'
                  f'N values must be integers:\n'
                  f'     N Initial: {pink}{type(NInitial)}{resetColor}\n'
                  f'     N Final: {pink}{type(NFinal)}{resetColor}')
        print(f'Unique Final Seqs: '
              f'{red}{self.nSubsFinalUniqueSeqs:,}{resetColor}\n')
        if NInitial == 0 and NFinal == 0:
            print(f'\n{orange}ERROR: No substrates were loaded in the '
                  f'{cyan}Background{orange} and the {cyan}Experimental{orange} sets\n')
            sys.exit(1)
        elif NInitial == 0:
            print(f'\n{orange}ERROR: No substrates were loaded in the '
                  f'{cyan}Background{orange} set\n')
            sys.exit(1)
        elif NFinal == 0:
            print(f'\n{orange}ERROR: No substrates were loaded in the '
                  f'{cyan}Experimental{orange} set\n')
            sys.exit(1)
        print()

        # Set figure titles
        enzName = self.enzymeName.replace(' - ', '\n')
        self.titleCombined = (f'{enzName}\n'
                              f'Combined Filter {self.datasetTag}')
        self.titleMotif = f'{enzName}\n{self.datasetTagMotif}'
        self.titleReleased = (f'{enzName}\n'
                              f'Substrate Profile\n'
                              f'{self.datasetTagMotif}')
        if self.showSampleSize:
            self.title = (f'{enzName}\n'
                          f'N Unsorted = {self.nSubsInitial:,}\n'
                          f'N Sorted = {self.nSubsFinal:,}')
            self.titleWeblogo = f'{enzName}\nN = {self.nSubsFinal:,}'

        else:
            self.title = f'{enzName}'
            self.titleWeblogo = f'{enzName}'
        if self.filterSubs:
            if self.motifFilter:
                self.titleWords = f'{enzName}\n{self.motifTag}'
            else:
                self.titleWords = f'{enzName}\n{self.datasetTagMotif}'
        elif self.excludeAAs:
            self.titleWords = f'{enzName}\n{self.datasetTag}'
        else:
            self.titleWords = f'{enzName}\nUnfiltered'



    def calculateRFCombinedMotif(self, countsCombinedMotifs):
        print('====================== Calculate: Fixed Motif Probability '
              '=======================')
        # Sum each column
        columnSums = np.sum(countsCombinedMotifs, axis=0)
        columnSums = pd.DataFrame(columnSums, index=countsCombinedMotifs.columns,
                                  columns=['Sum'])

        rf = pd.DataFrame(0.0, index=countsCombinedMotifs.index,
                            columns=countsCombinedMotifs.columns)

        # Calculate: Residue Probability
        for position in countsCombinedMotifs.columns:
            rf.loc[:, position] = (countsCombinedMotifs.loc[:, position] /
                                     columnSums.loc[position, 'Sum'])
        print(f'Dataset: {purple}{self.datasetTag}\n'
              f'{green}{rf}{resetColor}\n\n')

        return rf



    def calculateProbCodon(self, codonSeq):
        print('======================= Calculate: Residue Probabilities '
              '========================')
        print(f'Possible codons for: {codonSeq}')
        nucleotides = ['A', 'C', 'G', 'T']
        S = ['C', 'G']
        K = ['G', 'T']

        # Define what codons are associated with each residue
        codonsAA = {
            'A': ['GCT', 'GCC', 'GCA', 'GCG'],
            'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
            'N': ['AAT', 'AAC'],
            'D': ['GAT', 'GAC'],
            'C': ['TGT', 'TGC'],
            'E': ['GAA', 'GAG'],
            'Q': ['CAA', 'CAG'],
            'G': ['GGT', 'GGC', 'GGA', 'GGG'],
            'H': ['CAT', 'CAC'],
            'I': ['ATT', 'ATC', 'ATA'],
            'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
            'K': ['AAA', 'AAG'],
            'M': ['ATG'],
            'F': ['TTT', 'TTC'],
            'P': ['CCT', 'CCC', 'CCA', 'CCG'],
            'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
            'T': ['ACT', 'ACC', 'ACA', 'ACG'], # Count occurrences
            'W': ['TGG'],
            'Y': ['TAT', 'TAC'],
            'V': ['GTT', 'GTC', 'GTA', 'GTG']
        }

        codonsRes = {aa: [] for aa in codonsAA.keys()}

        # Initialize a list to store all possible combinations
        codons = []
        codonToAA = {codon: aa for aa, codons in codonsAA.items() for codon in codons}

        # Generate all possible combinations
        for combination in product(nucleotides, repeat=len(codonSeq)):
            # Check if the combination satisfies the conditions
            if all((c == 'N') or (c == 'S' and s in S) or (c == 'K' and s in K)
                   for c, s in zip(codonSeq, combination)):
                codon = ''.join(combination)
                if codon in codonToAA.keys():
                    AA = codonToAA[codon]
                codons.append(codon)
                codonsRes[AA].append(codon)

        # Print: All possible codon combinations
        for index, codon in enumerate(codons, 1):
            print(f'Codon {index}: {codon}')
        print('')

        print('Available Codons:')
        for AA, cod in codonsRes.items():
            print(f'     {AA}: {", ".join(cod)}')
        print('')

        # Count the possible codon combinations for each AA
        codonCounts = pd.DataFrame(index=self.letters, columns=['Counts'], data=0)
        for sequence in codons:
            for residue, codonsResidue in codonsAA.items():
                if sequence in codonsResidue:
                    if residue in codonCounts.index:
                        codonCounts.loc[residue, 'Counts'] += 1
                    break
        codonProbability = pd.DataFrame(index=self.letters,
                                        columns=['Probability'],
                                        data=0)
        codonProbability['Probability'] = codonCounts['Counts'] / len(codons)

        print('Amino Acid Probabilities:')
        for index, row in codonProbability.iterrows():
            print(f'     {index}    {round(row["Probability"] * 100, 2)} %')
        codonProb = round(sum(codonProbability["Probability"]) * 100, 2)
        print(f'Total probability of AA with {codonSeq}: {codonProb} %')
        print(f'Stop codon probability: {round(100 - codonProb, 2)} %\n\n')

        return codonProbability



    def calculateRF(self, counts, N, fileType, calcAvg=False):
        print('======================== Calculate: Residue Probability '
              '=========================')
        if self.filterSubs and 'initial' not in fileType.lower():
            print(f'Dataset: {purple}{self.enzymeName} {fileType} - Filter '
                  f'{self.datasetTag}{resetColor}\n')
        else:
            print(f'Dataset: {purple}{self.enzymeName} - {fileType}{resetColor}\n')

        # Calculate: Probability
        rf = counts / N
        print(f'{np.round(rf, 4)}\n\n')

        if calcAvg:
            rfAvg = np.sum(rf, axis=1) / len(rf.columns)
            rf = pd.DataFrame(rfAvg, index=rfAvg.index, columns=['Average RF'])
            print(f'{np.round(rf, 4)}\n\n')

        return rf



    def getDatasetTag(self, useCodonProb=False,
                      codon=None, combinedMotifs=False):
        if combinedMotifs:

            continuous, multiCombinedFrames = True, False
            if len(self.fixedPos) == 1:
                # tags = []
                # for idxList in range(len(self.fixedAA)):
                #     tag = ''
                #     # print(self.fixedPos, idxList)
                #     for idx, AA in enumerate(self.fixedAA[idxList]):
                #         tag += f'{AA}@R{self.fixedPos[idxList][idx]}_'
                #     tag = tag[:-1] # Remove the last char: "_"
                #     tags.append(tag)
                # print('tags', ','.join(tags))
                if isinstance(self.fixedAA[0], list):
                    self.datasetTag = f'[{",".join(self.fixedAA[0])}]@R{self.fixedPos[0]}'
                else:
                    self.datasetTag = f'{self.fixedAA[0]}@R{self.fixedPos[0]}'
            else:
                fixedPos, fixedAA = zip(*sorted(zip(self.fixedPos, self.fixedAA)))
                fixedPos = list(fixedPos)
                fixedAA = list(fixedAA)
                # print(f'Fixed Pos: {fixedPos}\n'
                #       f'Fixed AA:  {fixedAA}\n')
                for index in range(len(fixedPos) - 1):
                    # print(f'Idx: {index}')
                    pos1, pos2 = fixedPos[index], fixedPos[index + 1]
                    # print(f'Pos:\n'
                    #       f'     {pos1}\n'
                    #       f'     {pos2}\n')
                    if isinstance(pos1, int) and isinstance(pos2, int):
                        if pos1 == pos2 - 1 or pos1 == pos2 + 1:
                            continue
                        else:
                            if isinstance(pos1, list):
                                for indexPos in range(len(pos1)-1):
                                    if (pos1[indexPos] == pos1[indexPos + 1] - 1 or
                                            pos1[indexPos] == pos1[indexPos + 1] + 1):
                                        continue
                                    else:
                                        continuous = False
                            else:
                                continuous = False
                    elif isinstance(pos1, list) or isinstance(pos2, list):
                        multiCombinedFrames = True
                        continue
                        # Evaluate combined frames with multiple fixed positons
                        print(f'Pos1: {pos1}')
                        if isinstance(pos1, list):
                            for indexPos, posA in enumerate(pos1[:-1]):
                                posB = pos1[indexPos + 1]
                                if posB - posA != 1:
                                    continuous = False
                                    break
                        else:
                            for indexPos, posA in enumerate(pos2[:-1]):
                                posB = pos2[indexPos + 1]
                                if posB - posA != 1:
                                    continuous = False
                                    break
                    else:
                        continuous = False
                        break

                # Define the tag
                if continuous:
                    if multiCombinedFrames:
                        # Define subtags
                        fixedAA1 = self.fixedAA[0]
                        if isinstance(fixedAA1, list):
                            fixedAA1 = f'[{','.join(fixedAA1)}]'
                        fixedPos1 = self.fixedPos[0]
                        if isinstance(fixedPos1, list):
                            fixedPos1 = f'[{','.join(map(str, fixedPos1))}]'
                        fixedAA2 = self.fixedAA[-1]
                        if isinstance(fixedAA2, list):
                            fixedAA2 = f'[{','.join(fixedAA2)}]'
                        fixedPos2 = self.fixedPos[-1]
                        if isinstance(fixedPos2, list):
                            fixedPos2 = f'[{','.join(map(str, fixedPos2))}]'
                        # print(f'Fix: {fixedAA1}@R{fixedPos1}')
                        # print(f'Fix: {fixedAA2}@R{fixedPos2}\n')
                        self.datasetTag = (f'{fixedAA1}@R{fixedPos1}-'
                                           f'{fixedAA2}@R{fixedPos2}')
                    else:
                        tag = ''
                        for i in range(len(self.fixedAA)):
                            fixAA = self.fixedAA[i]
                            if isinstance(fixAA, list):
                                fixAA = f'[{','.join(fixAA)}]'
                            fixPos = self.fixedPos[i]
                            if tag:
                                tag += f' {fixAA}@R{fixPos}'
                            else:
                                tag = f'{fixAA}@R{fixPos}'
                        self.datasetTag = tag
                else:
                    tags = []
                    for idx, pos in enumerate(fixedPos):
                        # print(idx, pos)
                        aa = fixedAA[idx]
                        tag = f'{','.join(fixedAA[idx])}@R{pos}'
                        if isinstance(aa, list):
                            tag = f'[{','.join(fixedAA[idx])}]@R{pos}'
                        # print(f'Tag: {tag}')
                        tags.append(tag)
                    self.datasetTag = ' '.join(tags)
                    # print(f'Dataset: {self.datasetTag}\n')

            # Exclude residues
            if self.excludeAAs:
                fixResidueList = []
                for index, removedAA in enumerate(self.excludeAA):
                    fixResidueList.append(
                        f'{removedAA}@R'
                        f'{self.excludePos[index]}'.replace(' ', '')
                    )
                excludeTag = f'Exclude {' '.join(fixResidueList)}'
                self.datasetTag = f'{excludeTag} Fixed {self.datasetTag}'
        else:
            if self.filterSubs:
                fixResidueList = []
                if self.excludeAAs:
                    # Exclude residues
                    for index, removedAA in enumerate(self.excludeAA):
                        if index == 0:
                            fixResidueList.append(
                                f'Exclude_{removedAA}@R'
                                f'{self.excludePos[index]}'.replace(
                                    ' ', '')
                            )
                        else:
                            fixResidueList.append(
                                f'{removedAA}@R{self.excludePos[index]}'.replace(
                                    ' ', '')
                            )
    
                    # Fix residues
                    for index in range(len(self.fixedAA)):
                        tagFixedAA = self.fixedAA[index]
                        if isinstance(tagFixedAA, list) and len(tagFixedAA) == 1:
                            tagFixedAA = tagFixedAA[0]

                        if index == 0:
                            fixResidueList.append(
                                f'Fixed_{tagFixedAA}@R{self.fixedPos[index]}'.replace(
                                    ' ', '')
                            )
                        else:
                            fixResidueList.append(
                                f' {tagFixedAA}@R{self.fixedPos[index]}'.replace(
                                    ' ', '')
                            )
                    self.datasetTag = '_'.join(fixResidueList)
                else:
                    # Fix residues
                    try:
                        for index in range(len(self.fixedAA)):
                            fixResidueList.append(
                                f' {self.fixedAA[index]}@R'
                                f'{self.fixedPos[index]}'.replace(' ', ''))
                    except Exception as exc:
                        print(f'{orange}ERROR: {exc}\n')
                        print(f'Does the length of Fixed AA '
                              f'({cyan}{len(self.fixedAA)}{orange}) '
                              f'match the length of Fixed Pos '
                              f'({cyan}{len(self.fixedPos)}{orange})\n'
                              f'Fixed AA: {cyan}{self.fixedAA}{orange}\n'
                              f'Fixed Pos : {cyan}{self.fixedPos}{orange}\n')
                        sys.exit()
                    self.datasetTag = ' '.join(fixResidueList)
            else:
                if self.excludeAAs:
                    # Exclude residues
                    fixResidueList = []
                    for index, removedAA in enumerate(self.excludeAA):
                        if index == 0:
                            fixResidueList.append(
                                f'Exclude_{removedAA}@R'
                                f'{self.excludePos[index]}'.replace(
                                    ' ', '')
                            )
                        else:
                            fixResidueList.append(
                                f'{removedAA}@R{self.excludePos[index]}'.replace(
                                    ' ', '')
                            )
                    self.datasetTag = ' '.join(fixResidueList)
                else:
                    self.datasetTag = 'Unfiltered'
            # Condense the string
            if "'" in self.datasetTag:
                self.datasetTag = self.datasetTag.replace("'", '')
            self.datasetTag = self.datasetTag.replace("_", ' ')

            if self.excludeAAs:
                # Clean up fixed sequence tag
                if self.substrateLength == 9:
                    removeTag = 'Excl-Y@R1_Y@R2_Y@R3_Y@R4_Y@R6_Y@R7_Y@R8_Y@R9'
                    if removeTag in self.datasetTag:
                        # This should be reserved for simplifying dataset tags
                        self.datasetTag = self.datasetTag.replace(removeTag, '')
                        self.datasetTag = f'Exclude Y{self.datasetTag}'
                self.datasetTag = self.datasetTag.replace('_', ' ')

        # Define: Dataset Label
        if useCodonProb:
            self.datasetTag = f'{codon} Enrichment - {self.datasetTag}'
        if self.initialize:
            self.datasetTagMotif = self.datasetTag
            self.initialize = False

        # Define: Figure tag
        tag = self.datasetTag
        tag = tag.replace(' ', '_')
        self.datasetTagSave = tag.replace(' ', '_')
        # print(f'Dataset Tag: {purple}{self.datasetTag}{resetColor}')
        # print(f'Motif Tag: {purple}{self.datasetTagMotif}{resetColor}')
        # print(f'Save Tag: {purple}{self.datasetTagSave}{resetColor}\n\n')
        # sys.exit()

        return self.datasetTag



    def fixResidue(self, substrates, fixedString, printRankedSubs, sortType):
        print('=============================== Filter Substrates '
              '===============================')
        fixedSubs = {}
        fixedSubsTotal = 0
        print(f'Selecting {purple}{sortType} {resetColor}substrates with: '
              f'{purple}{fixedString}{resetColor}\n')

        # Sort the substrate dictionary by counts
        substrates = dict(sorted(
            substrates.items(), key=lambda x: x[1], reverse=True)
        )


        # Select substrates that contain selected AA at a specified position
        if self.excludeAAs:
            # Verify if the substrates contain the residue(s) you wish to remove
            for substrate, count in substrates.items():
                # Inspect substrate count
                if count < self.minSubCount:
                    break

                keepSub = True
                for indexExclude, AAExclude in enumerate(self.excludeAA):
                    if len(AAExclude) == 1:
                        indexRemoveAA = self.excludePos[indexExclude] - 1

                        # Is the AA acceptable?
                        if substrate[indexRemoveAA] == AAExclude:
                            keepSub = False
                            continue
                    else:
                        # Remove Multiple AA at a specific position
                        for AAExcludeMulti in AAExclude:
                            indexRemoveAA = self.excludePos[indexExclude] - 1
                            for AAExclude in AAExcludeMulti:

                                # Is the AA acceptable?
                                if substrate[indexRemoveAA] == AAExclude:
                                    keepSub = False
                                    continue

                # If the substrate has not been blacklisted, look for the desired AA
                if keepSub:
                    if len(self.fixedAA) == 1 and len(self.fixedAA[0]) == 1:
                        # Fix only one AA
                        if substrate[self.fixedPos[0] - 1] != self.fixedAA[0]:
                            keepSub = False
                    else:
                        for indexFixed, fixedAA in enumerate(self.fixedAA):
                            indexFixAA = self.fixedPos[indexFixed] - 1

                            if len(fixedAA) == 1:
                                # Fix one AA at a given position
                                if substrate[indexFixAA] != fixedAA:
                                    keepSub = False
                                    break
                            else:
                                # Fix multiple AAs at a given position
                                if substrate[indexFixAA] not in fixedAA:
                                    keepSub = False
                                    break
                # Extract the substrate
                if keepSub:
                    fixedSubs[substrate] = count
                    fixedSubsTotal += count
        else:
            # Fix AAs, and don't exclude any AAs
            if len(self.fixedAA) == 1 and len(self.fixedAA[0]) == 1:
                for substrate, count in substrates.items():
                    # Inspect substrate count
                    if count < self.minSubCount:
                        break

                    subAA = substrate[self.fixedPos[0] - 1]
                    if subAA in self.fixedAA[0]:
                        fixedSubs[substrate] = count
                        fixedSubsTotal += count
                        continue
            else:
                for substrate, count in substrates.items():
                    # Inspect substrate count
                    if count < self.minSubCount:
                        break

                    keepSub = []
                    for index in range(len(self.fixedAA)):
                        fixIndex = self.fixedPos[index] - 1
                        subAA = substrate[fixIndex]
                        selectAA = self.fixedAA[index]

                        if subAA in selectAA:
                            keepSub.append(True)
                        else:
                            keepSub.append(False)

                    if False not in keepSub:
                        fixedSubs[substrate] = count
                        fixedSubsTotal += count


        # Rank fixed substrates
        rankedFixedSubstrates = dict(sorted(
            fixedSubs.items(), key=lambda x: x[1], reverse=True)
        )

        # Print: Fixed substrates
        if printRankedSubs:
            iteration = 0
            fixedUniqueSubsTotal = len(rankedFixedSubstrates)
            print('Ranked Fixed Substrates:')
            if fixedUniqueSubsTotal == 0:
                print('')
                print(f'{orange}ERROR:\n'
                      f'     No substrates in {purple}{sortType}{orange} contained: '
                      f'{red}{fixedString}{resetColor}\n')
                sys.exit(1)
            else:
                for substrate, count in rankedFixedSubstrates.items():
                    print(f'     {pink}{substrate}{resetColor}, Counts: {red}{count:,}'
                          f'{resetColor}')
                    iteration += 1
                    if iteration >= self.printNumber:
                        break

        print(f'\nNumber of substrates with fixed {purple}{fixedString}{resetColor}: '
              f'{red}{fixedSubsTotal:,}{resetColor}\n\n')

        return rankedFixedSubstrates, fixedSubsTotal




    def exclResidue(self, substrates, fixedString, printRankedSubs, sortType):
        print('=============================== Filter Substrates '
              '===============================')
        fixedSubs = {}
        fixedSubsTotal = 0
        print(f'Selecting {purple}{sortType} {resetColor}substrates with: '
              f'{purple}{fixedString}{resetColor}\n')

        # Sort the substrate dictionary by counts
        substrates = dict(sorted(substrates.items(), key=lambda x: x[1], reverse=True))


        # Select substrates that contain selected AA at a specified position in the substrate
        if self.excludeAAs:
            # Verify if the substrates contain the residue(s) you wish to remove
            for substrate, count in substrates.items():
                # Inspect substrate count
                if count < self.minSubCount:
                    break

                keepSub = True
                for indexExclude, AAExclude in enumerate(self.excludeAA):
                    if len(AAExclude) == 1:
                        indexRemoveAA = self.excludePos[indexExclude] - 1
                        # Is the AA acceptable?
                        if substrate[indexRemoveAA] == AAExclude:
                            keepSub = False
                            continue
                    else:
                        # Remove Multiple AA at a specific position
                        for AAExcludeMulti in AAExclude:
                            indexRemoveAA = self.excludePos[indexExclude] - 1
                            for AAExclude in AAExcludeMulti:

                                # Is the AA acceptable?
                                if substrate[indexRemoveAA] == AAExclude:
                                    keepSub = False
                                    continue

                # Extract the substrate
                if keepSub:
                    fixedSubs[substrate] = count
                    fixedSubsTotal += count

        # Rank fixed substrates
        rankedFixedSubstrates = dict(sorted(
            fixedSubs.items(), key=lambda x: x[1], reverse=True)
        )

        # Print: Fixed substrates
        if printRankedSubs:
            iteration = 0
            fixedUniqueSubsTotal = len(rankedFixedSubstrates)
            print('Ranked Fixed Substrates:')
            if fixedUniqueSubsTotal == 0:
                print('')
                print(f'{orange}ERROR:\n'
                      f'     No substrates in {purple}{sortType}{orange} contained: '
                      f'{red}{fixedString}{resetColor}\n')
                sys.exit(1)
            else:
                for substrate, count in rankedFixedSubstrates.items():
                    print(f'     {pink}{substrate}{resetColor}, Counts: {red}{count:,}'
                          f'{resetColor}')
                    iteration += 1
                    if iteration >= self.printNumber:
                        break


        print(f'\nNumber of substrates with fixed {purple}{fixedString}{resetColor}: '
              f'{red}{fixedSubsTotal:,}{resetColor}\n\n')

        return rankedFixedSubstrates, fixedSubsTotal



    def identifyMotif(self, fixFullFrame):
        print('================================ Identify Motif '
              '=================================')
        if fixFullFrame:
            print(f'Selecting continuous motif')
        else:
            print(f'Selecting non-continuous motif')
        print(f'Minimum ΔS Value:{red} {self.minEntropy} Bits{resetColor}\n\n'
              f'Entropy:\n'
              f'{self.entropy}{resetColor}\n')
        subFrame = self.entropy.copy()
        lastPosition = len(self.entropy) - 1

        # Determine Substrate Frame
        if fixFullFrame:
            for indexPos, position in enumerate(self.entropy.index):
                if indexPos == 0 or indexPos == lastPosition:
                    if self.entropy.loc[position, 'ΔS'] < self.minEntropy:
                        subFrame.drop(position, inplace=True)
                else:
                    if self.entropy.loc[position, 'ΔS'] < self.minEntropy:
                        if (self.entropy.iloc[indexPos - 1, 0] > self.minEntropy and
                                self.entropy.iloc[indexPos + 1, 0] > self.minEntropy):
                            pass
                        else:
                            subFrame.drop(position, inplace=True)
        else:
            for indexPos, position in enumerate(self.entropy.index):
                if self.entropy.loc[position, 'ΔS'] < self.minEntropy:
                    subFrame.drop(position, inplace=True)

        # Sort the frame
        self.subFrame = subFrame.sort_values(by='ΔS', ascending=False).copy()
        print(f'Motif:\n'
              f'{blue}{subFrame}{resetColor}\n\n'
              f'Ranked Motif:\n'
              f'{blue}{self.subFrame}{resetColor}\n\n')

        # Define motif label
        indexSubFrameList = list(self.entropy.index)
        self.motifIndex = [indexSubFrameList.index(idx) for idx in self.subFrame.index]
        self.xAxisLabelsMotif = self.xAxisLabels[
                                min(self.motifIndex):max(self.motifIndex)+1]

        return self.subFrame



    def getMotif(self, substrates):
        print('================================= Extract Motif '
              '=================================')
        if self.motifIndex is None:
            print(f'The motif index is {cyan}{self.motifIndex}{resetColor}\n'
                  f'NGS.getMotif(substrates) will return the full length substrate '
                  f'sequences\n\n')
            return substrates
        motifs = {}
        indexStart = min(self.motifIndex)
        indexEnd = max(self.motifIndex) + 1
        print(f'Frame: {purple}{self.datasetTag}{resetColor}\n')

        # Print: data
        iteration = 0
        print(f'Substrates:')
        for substrate, count in substrates.items():
            print(f'     {pink}{substrate}{resetColor}, '
                  f'{blue}{substrate[indexStart:indexEnd]}{resetColor}, '
                  f'Count:{red} {count:,}{resetColor}')
            iteration += 1
            if iteration >= self.printNumber:
                break
        print(f'\nUnique Substrates: {red}{len(substrates):,}{resetColor}\n\n')

        # Extract motif
        countTotalSubstrates = 0
        for substrate, count in substrates.items():
            countTotalSubstrates += count
            motif = substrate[indexStart:indexEnd]
            if motif in motifs.keys():
                motifs[motif] += count
            else:
                motifs[motif] = count
        self.motifLen = len(next(iter(motifs)))
        motifs = dict(sorted(motifs.items(), key=lambda x: x[1], reverse=True))

        # Print: data
        iteration = 0
        print(f'Motifs:')
        for motif, count in motifs.items():
            print(f'     {blue}{motif}{resetColor}, Count:{red} {count:,}{resetColor}')
            iteration += 1
            if iteration >= self.printNumber:
                print()
                break

        print(f'Unique Motifs: {red}{len(motifs):,}{resetColor}\n\n')

        return motifs



    def calculateEnrichment(self, rfInitial, rfFinal, combinedMotifs=False,
                            posFilter=False, relFilter=False, relIteration=False,
                            relCounts=False, plotFigures=True):
        print('========================== Calculate: Enrichment Score '
              '==========================')
        print(f'Enrichment Scores:\n'
              f'     {magenta}log₂(RF FinalAA / RF InitialAA){resetColor}\n')
        print(f'RF Final:\n{rfFinal}\n\n'
              f'RF Initial:\n{rfInitial}\n\n')
        # Calculate: Enrichment scores
        if len(rfInitial.columns) == 1:
            matrix = pd.DataFrame(0.0, index=rfFinal.index,
                                  columns=rfFinal.columns)
            ln = False
            if ln:
                print(f'Using {purple}Natural Logarithmic{resetColor} Matrix')
                for position in rfFinal.columns:
                    # for AA in rfFinal.index:
                    matrix.loc[:, position] = np.log(
                        rfFinal.loc[:, position] / rfInitial.iloc[:, 0]
                    )
            else:
                for position in rfFinal.columns:
                    matrix.loc[:, position] = np.log2(rfFinal.loc[:, position] /
                                                      rfInitial.iloc[:, 0])
        else:
            if len(rfInitial.columns) != len(rfFinal.columns):
                print(f'{orange}ERROR: The number of columns in the Initial Sort '
                      f'({cyan}{len(rfInitial.columns)}{orange}) needs to equal to the '
                      f'number of columns in the Final Sort '
                      f'({cyan}{len(rfFinal.columns)}{orange})\n'
                      f'     Initial: {cyan}{rfInitial.columns}{orange}\n'
                      f'       Final: {cyan}{rfFinal.columns}\n\n')
                sys.exit(1)

            rfInitial.columns = rfFinal.columns
            matrix = np.log2(rfFinal / rfInitial)
            # matrix = rfFinal
        if self.releasedCounts:
            print(f'Enrichment Score: {purple}Substrate Profile{resetColor}\n'
                  f'{matrix.round(self.roundVal)}\n\n')
            print(f'Prob Initial:\n{rfInitial}\n\n'
                  f'Prob Final:\n{rfFinal}\n\n')
        else:
            print(f'Enrichment Score: {purple}{self.datasetTag}{resetColor}\n'
                  f'{matrix.round(self.roundVal)}\n\n')

        print('====================== Calculate: Scaled Enrichment Score '
              '=======================')
        if self.releasedCounts:
            print(f'Scale Enrichment Scores: {purple}Substrate Profile{resetColor}\n'
                  f'     {magenta}Enrichment Scores * ΔS{resetColor}')
        else:
            print(f'Scale Enrichment Scores:\n'
                  f'     {magenta}Enrichment Scores * ΔS{resetColor}')

        # if self.releasedCounts:
        #     print(f'Set Values: Fraud!!!')
        #     matrix.loc['G', 'R4'] = -0.556
        #     matrix.loc['S', 'R4'] = -0.762
        #     matrix.loc['R', 'R5'] = -0.352

        # Calculate: Letter heights
        heights = pd.DataFrame(0, index=matrix.index,
                               columns=matrix.columns, dtype=float)
        for indexColumn in heights.columns:
            heights.loc[:, indexColumn] = (matrix.loc[:, indexColumn] *
                                           self.entropy.loc[indexColumn, 'ΔS'])

        # Record values
        if self.releasedCounts:
            self.eMapReleased = matrix
            self.eMapReleasedScaled = heights.copy()
        else:
            self.eMap = matrix
            self.eMapScaled = heights.copy()

        # Calculate: Max positive
        columnTotals = []
        for indexColumn in heights.columns:
            totalPos = 0
            for value in heights.loc[:, indexColumn]:
                if value > 0:
                    totalPos += value
            columnTotals.append(totalPos)
        yMax = max(columnTotals)

        # Adjust values
        for pos in heights.columns:
            if heights.loc[:, pos].isna().any():
                nValues = heights[pos].notna().sum()
                print(f'\nNumber non NaN values in {purple}{pos}{resetColor}:\n'
                      f'{red}{heights.loc[heights[pos].notna(), pos]}'
                      f'{resetColor}\n')
                for aa in heights.index:
                    val = heights.loc[aa, pos]
                    if val == np.inf or val == -np.inf:
                        print(f'{aa}@{pos}: {val} -> 0')
                        heights.loc[aa, pos] = 0
                # heights.loc[heights[pos].notna(), pos] = yMax / nValues
                heights.loc[:, pos] = heights.loc[:, pos].fillna(0)

        heights = heights.replace([np.inf, -np.inf], 0)
        self.heights = heights
        print(f'\nResidue Heights: {purple}{self.datasetTag}{resetColor}\n'
              f'{heights}\n\n')


        if plotFigures:
            # Plot: Enrichment Map
            if self.plotFigEM:
                self.plotEnrichmentScores(dataType='Enrichment',
                                          combinedMotifs=combinedMotifs,
                                          posFilter=posFilter,
                                          relFilter=relFilter,
                                          relIteration=relIteration,
                                          relCounts=relCounts)
            if self.plotFigEMScaled:
                self.plotEnrichmentScores(dataType='Scaled Enrichment',
                                          combinedMotifs=combinedMotifs,
                                          posFilter=posFilter,
                                          relFilter=relFilter,
                                          relIteration=relIteration,
                                          relCounts=relCounts)

            # Plot: Enrichment Logo
            if self.plotFigLogo:
                self.plotEnrichmentLogo(combinedMotifs=combinedMotifs,
                                        posFilter=posFilter,
                                        relFilter=relFilter,
                                        relIteration=relIteration,
                                        relCounts=relCounts)

            # Calculate & Plot: Weblogo
            if self.plotFigWebLogo:
                self.calculateWeblogo(probability=rfFinal,
                                      combinedMotifs=combinedMotifs,
                                      relIteration=relIteration,
                                      relCounts=relCounts)

        return self.eMap



    def plotEnrichmentScores(self, dataType, combinedMotifs=False, posFilter=False,
                             relFilter=False, relIteration=False, relCounts=False,
                             pca=False):
        print('============================ Plot: Enrichment Score '
              '=============================')
        # Select: Dataset
        if 'scaled' in dataType.lower():
            if self.releasedCounts:
                scores = self.eMapReleasedScaled
            else:
                scores = self.eMapScaled
        else:
            if self.releasedCounts:
                scores = self.eMapReleased
            else:
                scores = self.eMap

        # Define: Figure title
        if self.releasedCounts or relCounts: # and len(self.motifIndexExtracted) > 1:
            title = self.titleReleased
        elif combinedMotifs:
            title = self.titleCombined
        elif self.motifFilter:
            title = self.titleMotif
        else:
            title = self.title
        if len(self.datasetTag.replace('[', '').replace(']', '').replace('-', '')) > 40:
            title = title.replace('Register ', 'Register\n')
        if pca:
            title = f'{self.enzyme}\n{pca}'

        print(f'Dataset: {purple}{self.datasetTag}{resetColor}\n'
             f'Unique Substrates: {red}{self.nSubsFinalUniqueSeqs:,}{resetColor}')
        if self.motifFilter:
            print(f'Figure Number: '
                  f'{magenta}{self.saveFigureIteration}{resetColor}')
        if posFilter:
            if relFilter:
                print(f'Releasing Filter: {magenta}{posFilter}{resetColor}')
            else:
                print(f'Applying Filter: {magenta}{posFilter}{resetColor}')

        print(f'\nEnrichment Scores:\n'
              f'{scores}\n\n')

        # Create heatmap
        cMapCustom = self.createCustomColorMap(colorType='EM')

        # Define the yLabel
        if self.residueLabelType == 0:
            scores.index = [residue[0] for residue in self.residues]
        elif self.residueLabelType == 1:
            scores.index = [residue[1] for residue in self.residues]
        elif self.residueLabelType == 2:
            scores.index = [residue[2] for residue in self.residues]

        # Define color bar limits
        if np.max(scores) >= np.min(scores):
            cBarMax = np.max(scores)
            cBarMin = -1 * cBarMax
        else:
            cBarMin = np.min(scores)
            cBarMax = -1 * cBarMin
        fig, ax = plt.subplots(figsize=self.figSizeEM)

        # Plot the heatmap
        if self.figEMSquares:
            sns.heatmap(scores, annot=False, cmap=cMapCustom, cbar=False,
                        linewidths=self.lineThickness - 1, linecolor='black',
                        square=self.figEMSquares, center=None, vmax=cBarMax, vmin=cBarMin)
        else:
            # Plot enrichment scores
            sns.heatmap(scores, annot=True, fmt='.3f', cmap=cMapCustom, cbar=False,
                        linewidths=self.lineThickness - 1, linecolor='black',
                        square=self.figEMSquares, center=None, vmax=cBarMax, vmin=cBarMin,
                        annot_kws={'fontweight': 'bold', 'size': self.labelSizeEM})
        ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        ax.set_xlabel('Substrate Position', fontsize=self.labelSizeAxis)
        ax.set_ylabel('Residue', fontsize=self.labelSizeAxis)


        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(self.lineThickness)

        # Set tick parameters
        ax.tick_params(axis='both', which='major', rotation=0, length=self.tickLength,
                       labelsize=self.labelSizeTicks, width=self.lineThickness)

        # Set x-ticks
        xTicks = np.arange(len(scores.columns)) + 0.5
        ax.set_xticks(xTicks)
        ax.set_xticklabels(scores.columns)

        # Set y-ticks
        yTicks = np.arange(len(scores.index)) + 0.5
        ax.set_yticks(yTicks)
        ax.set_yticklabels(scores.index)

        # Colormap
        if isinstance(cMapCustom, Colormap):
            cmap = cMapCustom
        else:
            # Look up premade colormap
            cmap = matplotlib.colormaps[cMapCustom]
        cmap.set_bad(color='lightgrey') # Set invalid values to gray

        # Modify the colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.1)
        norm = plt.Normalize(vmin=cBarMin, vmax=cBarMax)
        cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cMapCustom), cax=cax)
        cbar.ax.tick_params(axis='y', which='major', labelsize=self.labelSizeTicks,
                            length=self.tickLength, width=self.lineThickness)
        cbar.outline.set_linewidth(self.lineThickness)
        cbar.outline.set_edgecolor('black')

        self.plotFig(plt=plt, fig=fig)


        # Save the figure
        if self.saveFigures:
            if 'Scaled' in dataType:
                figType = 'EM_Scaled'
            elif 'Enrichment' in dataType:
                figType = 'EM'
            else:
                print(f'{orange}ERROR: What do I do with this dataset type -'
                      f'{cyan} {dataType}{resetColor}\n')
                sys.exit(1)
            if pca:
                pca = (pca.replace(f'{self.datasetTag}', '').replace('\n', '')
                       .replace(' ', '_').replace('Population', 'Pop'))
                figType = f'{pca}-{figType}'
            if not isinstance(relIteration, bool):
                figType += f'_{relIteration}'
            self.saveFigure(fig=fig, figType=figType, seqLen=len(xTicks),
                            combinedMotifs=combinedMotifs)



    def plotEnrichmentLogo(self, combinedMotifs=False, posFilter=False, relFilter=False,
                           relIteration=False, relCounts=False, pca=False):
        print('============================= Plot: Enrichment Logo '
              '=============================')
        # Define: Figure title
        if self.releasedCounts or relCounts:
            title = self.titleReleased
        elif combinedMotifs:
            title = self.titleCombined
        elif self.motifFilter:
            title = self.titleMotif
        else:
            title = self.title
        if len(self.datasetTag.replace('[', '').replace(']', '').replace('-', '')) > 40:
            title = title.replace('Register ', 'Register\n')
        if pca:
            title = f'{self.enzyme}\n{pca}'

        # Print: data
        print(f'Dataset: {purple}{self.datasetTag}{resetColor}\n'
              f'Unique Substrates: {red}{self.nSubsFinalUniqueSeqs:,}{resetColor}')
        if self.motifFilter:
            print(f'Figure Number: '
                  f'{magenta}{self.saveFigureIteration}{resetColor}')
        if posFilter:
            if relFilter:
                print(f'Releasing Filter: {magenta}{posFilter}{resetColor}')
            else:
                print(f'Applying Filter: {magenta}{posFilter}{resetColor}')
        print(f'\nResidue heights:\n'
              f'{self.heights}\n')

        # Set local parameters
        if self.bigAAonTop:
            stackOrder = 'big_on_top'
        else:
            stackOrder = 'small_on_top'


        # Calculate: Max and min
        columnTotals = [[], []]
        for indexColumn in self.heights.columns:
            totalPos = 0
            totalNeg = 0
            for value in self.heights.loc[:, indexColumn]:
                if value > 0:
                    totalPos += value
                elif value < 0:
                    totalNeg += value
            columnTotals[0].append(totalPos)
            columnTotals[1].append(totalNeg)
        yMax = max(columnTotals[0])
        yMin = min(columnTotals[1])
        print(f'y Max: {red}{np.round(yMax, 4)}{resetColor}\n'
              f'y Min: {red}{np.round(yMin, 4)}{resetColor}\n')

        # Rename columns for logomaker script
        data = self.heights.copy()
        data.columns = range(len(data.columns))


        def plotLogo(limitYAxis=False, pca=False):
            # Plot the sequence motif
            fig, ax = plt.subplots(figsize=self.figSize)
            motif = logomaker.Logo(data.transpose(), ax=ax, color_scheme=self.colorsAA,
                                   width=0.95, stack_order=stackOrder)
            motif.ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')


            # Set tick parameters
            ax.tick_params(axis='both', which='major', length=self.tickLength,
                           labelsize=self.labelSizeTicks)

            # Set borders
            motif.style_spines(visible=False)
            motif.style_spines(spines=['left', 'bottom'], visible=True)
            for spine in motif.ax.spines.values():
                spine.set_linewidth(self.lineThickness)

            # Set x-ticks
            motif.ax.set_xticks([pos for pos in range(len(self.heights.columns))])
            motif.ax.set_xticklabels(self.heights.columns, fontsize=self.labelSizeTicks,
                                     rotation=0, ha='center')

            # Set y-ticks
            yTicks = [yMin, 0, yMax]
            yTickLabels = [f'{tick:.2f}' if tick != 0 else f'{int(tick)}'
                           for tick in yTicks]
            motif.ax.set_yticks(yTicks)
            motif.ax.set_yticklabels(yTickLabels, fontsize=self.labelSizeTicks)
            motif.ax.set_ylim(yMin, yMax)

            # Set tick width
            for tick in motif.ax.xaxis.get_major_ticks():
                tick.tick1line.set_markeredgewidth(self.lineThickness)
            for tick in motif.ax.yaxis.get_major_ticks():
                tick.tick1line.set_markeredgewidth(self.lineThickness)

            # Label the axes
            motif.ax.set_xlabel('Substrate Position', fontsize=self.labelSizeAxis)
            motif.ax.set_ylabel('Scaled Enrichment', fontsize=self.labelSizeAxis)

            # Set horizontal line
            motif.ax.axhline(y=0, color='black', linestyle='-',
                             linewidth=self.lineThickness)

            # Evaluate dataset for fixed residues
            spacer = np.diff(motif.ax.get_xticks())  # Find the space between each tick
            spacer = spacer[0] / 2

            # Use the spacer to set a gray background to fixed residues
            if self.filterSubs:
                for index, position in enumerate(self.xAxisLabels):
                    if position in self.fixedPos:
                        # Plot gray boxes on each side of the xtick
                        motif.ax.axvspan(index - spacer, index + spacer,
                                         facecolor='darkgrey', alpha=0.2)

            self.plotFig(plt=plt, fig=fig)

            # Save the figure
            if self.saveFigures:
                figType = 'Logo'
                if limitYAxis:
                    figType += '_yMin'
                if not isinstance(relIteration, bool):
                    figType += f'_{relIteration}'
                if pca:
                    pca = (pca.replace(f'{self.datasetTag}', '').replace('\n', '')
                           .replace(' ', '_').replace('Population', 'Pop'))
                    figType = f'{pca}-{figType}'
                self.saveFigure(fig=fig, figType=figType, seqLen=len(data.columns),
                                combinedMotifs=combinedMotifs)

        # Plot figure
        plotLogo(pca=pca) # Full y-axis

        # Adjust yMin to fit the largest negative AA
        yMin = 0
        for col in self.heights.columns:
            for row in self.heights.index:
                if self.heights.loc[row, col] < yMin:
                    yMin = self.heights.loc[row, col]
        print(f'Adjusting Y Min:\n'
              f'y Max: {red}{np.round(yMax, 4)}{resetColor}\n'
              f'y Min: {red}{np.round(yMin, 4)}{resetColor}\n')
        plotLogo(limitYAxis=True, pca=pca) # Limited y-axis



    def calculateWeblogo(self, probability, combinedMotifs=False,
                         relIteration=False, relCounts=False):
        print('============================= Calculate: Weblogo '
              '================================')
        print(f'Probability: {self.datasetTag}\n{probability}\n\n'
              f'Max Entropy: {red}{self.entropyMax.round(6)}{resetColor}\n\n')

        # Calculate: Weblogo
        self.weblogo = pd.DataFrame(0, index=probability.index,
                                    columns=probability.columns, dtype=float)
        for indexColumn in self.weblogo.columns:
            self.weblogo.loc[:, indexColumn] = (probability.loc[:, indexColumn] *
                                                self.entropy.loc[indexColumn, 'ΔS'])

        if self.plotFigWebLogo:
            self.plotWeblogo(
                combinedMotifs=combinedMotifs, relIteration=False, relCounts=relCounts
            )



    def plotWeblogo(self, combinedMotifs=False,
                    relIteration=False, relCounts=False):
        print('================================= Plot: WebLogo '
              '=================================')
        if self.motifFilter:
            print(f'Figure Number: '
                  f'{magenta}{self.saveFigureIteration}{resetColor}')
        print(f'Letter Heights: {purple}{self.datasetTag}{resetColor}\n'
              f'{self.weblogo}\n\n')

        # Define: Figure title
        if self.releasedCounts or relCounts:
            title = self.titleReleased
        elif combinedMotifs:
            title = self.titleCombined
        else:
            title = self.titleWeblogo

        # Set local parameters
        if self.bigAAonTop:
            stackOrder = 'big_on_top'
        else:
            stackOrder = 'small_on_top'

        # Rename column headers
        dataColumnsReceived = self.weblogo.columns
        self.weblogo.columns = range(len(self.weblogo.columns))

        # Set -inf to zero
        if self.weblogo.isin([np.inf, -np.inf]).any().any():
            self.weblogo.replace([np.inf, -np.inf], 0, inplace=True)

        # Plot the sequence motif
        fig, ax = plt.subplots(figsize=self.figSize)
        motif = logomaker.Logo(self.weblogo.transpose(), ax=ax,
                               color_scheme=self.colorsAA, width=0.95,
                               stack_order=stackOrder)
        motif.ax.set_title(title, fontsize=self.labelSizeTitle,
                           fontweight='bold')


        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks)

        # Set borders
        motif.style_spines(visible=False)
        motif.style_spines(spines=['left', 'bottom'], visible=True)
        for spine in motif.ax.spines.values():
            spine.set_linewidth(self.lineThickness)

        # Set x-ticks
        if len(dataColumnsReceived) == len(self.xAxisLabels):
            motif.ax.set_xticks([pos for pos in range(len(self.xAxisLabels))])
            motif.ax.set_xticklabels(self.xAxisLabels, fontsize=self.labelSizeTicks,
                                     rotation=0, ha='center')
        else:
            motif.ax.set_xticks([pos for pos in range(len(dataColumnsReceived))])
            motif.ax.set_xticklabels(dataColumnsReceived, fontsize=self.labelSizeTicks,
                                     rotation=0, ha='center')

        for tick in motif.ax.xaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness)  # Set tick width

        # Set y-ticks
        yTicks = range(0, 5)
        yTickLabels = [f'{tick:.0f}' if tick != 0 else f'{int(tick)}'
                       for tick in yTicks]
        motif.ax.set_yticks(yTicks)
        motif.ax.set_yticklabels(yTickLabels, fontsize=self.labelSizeTicks)
        motif.ax.set_ylim(0, self.entropyMax)
        for tick in motif.ax.yaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness) # Set tick width

        # Label the axes
        motif.ax.set_xlabel('Substrate Position', fontsize=self.labelSizeAxis)
        motif.ax.set_ylabel('Bits', fontsize=self.labelSizeAxis)

        # Set horizontal line
        motif.ax.axhline(y=0, color='black', linestyle='-', linewidth=self.lineThickness)

        # Evaluate dataset for fixed residues
        spacer = np.diff(motif.ax.get_xticks()) # Find the space between each tick
        spacer = spacer[0] / 2

        # Use the spacer to set a gray background to fixed residues
        if self.filterSubs:
            for index, position in enumerate(self.xAxisLabels):
                if position in self.fixedPos:
                    # Plot gray boxes on each side of the xtick
                    motif.ax.axvspan(index - spacer, index + spacer,
                                     facecolor='darkgrey', alpha=0.2)

        self.plotFig(plt=plt, fig=fig)


        # Save the figure
        if self.saveFigures:
            datasetType = 'Weblogo'
            if not isinstance(relIteration, bool):
                datasetType += f'_{relIteration}'
            self.saveFigure(fig=fig, figType=datasetType,
                            seqLen=len(self.weblogo.columns),
                            combinedMotifs=combinedMotifs)



    def fixedMotifStats(self, countsList, initialRF,
                        motifFrame, datasetTag):
        print('================== Statistical Evaluation: Fixed Motif Counts '
              '===================')
        print(f'Count Lists: {len(countsList)}\n{countsList}\n\n')
        print(f'Evaluate: {purple}{datasetTag}{resetColor}\n')
        countsFrameTotal = pd.DataFrame(0, index=range(0, len(self.fixedPos)),
                                        columns=motifFrame)
        frameRF = pd.DataFrame(0.0, index=initialRF.index,
                                 columns=motifFrame)
        frameES = frameRF.copy()
        frameESList = []

        for index, countsFrame in enumerate(countsList):
            countsFrame = pd.DataFrame(countsFrame, index=initialRF.index,
                                       columns=motifFrame)

            # Format values to have commas
            formattedCounts = countsFrame.to_string(
                formatters={column: '{:,.0f}'.format for column in
                            countsFrame.select_dtypes(include='number').columns})
            print(f'Counts: {purple}{self.motifTags[index]}{resetColor}\n'
                  f'{formattedCounts}\n')

            for column in countsFrame.columns:
                countsFrameTotal.loc[index, column] = countsFrame[column].sum()

                # Calculate: RF
                frameRF.loc[:, column] = (countsFrame.loc[:, column] /
                                            countsFrameTotal.loc[index, column])

                # Calculate: ES
                frameES.loc[:, column] = np.log2(
                    frameRF.loc[:, column] / initialRF['Average RF'])

            print(f'Enrichment Score: {purple}{self.motifTags[index]}\n'
                  f'{greenLight}{frameES}{resetColor}\n\n')
            frameESList.append(frameES.copy())

        # Combine the ES DFs
        frameESCombined = pd.concat(frameESList, axis=0, keys=range(len(frameESList)))

        def calcAverage(x):
            # Function to calculate standard deviation ignoring -inf values

            # Remove -inf values
            xFiltered = x.replace([-np.inf, np.inf], np.nan).dropna()

            return np.average(xFiltered)

        def calcStDev(x):
            # Function to calculate standard deviation ignoring -inf values

            # Remove -inf values
            xFiltered = x.replace([-np.inf, np.inf], np.nan).dropna()

            return np.std(xFiltered)

        # Calculate standard deviation for each value across the corresponding positions
        frameESAvg = frameESCombined.groupby(level=1, sort=False).agg(calcAverage)
        frameESStDev = frameESCombined.groupby(level=1, sort=False).agg(calcStDev)

        print(f'Average: {purple}Enrichment Score{resetColor}\n'
              f'{frameESAvg}\n\n'
              f'Standard Deviation: {purple}Enrichment Score{resetColor}\n'
              f'{frameESStDev}{resetColor}\n\n')


        # Plot: Standard deviation
        self.plotStats(
            data=frameESStDev, totalCounts=None,
            dataType='Standard Deviation', combinedMotifs=True
        )



    def plotStats(self, data, totalCounts, dataType,
                  combinedMotifs=False):
        print('========================= Plot: Statistical Evaluation '
              '==========================')
        print(f'{dataType}: {purple}{self.datasetTag}{resetColor}\n{data}\n\n')

        # Set figure title
        enzName = self.enzymeName.replace(' - ', '\n')
        if totalCounts is not None and self.showSampleSize:
            title = f'{enzName}\n{self.datasetTag}\nN={totalCounts:,}'
            title = title.replace('N= ', f'{dataType}\n')
        else:
            if self.releasedCounts:
                title = self.titleReleased
                if len(self.datasetTag.replace('[', '').replace(
                        ']', '').replace('-', '')) > 40:
                    title = title.replace('Register ', 'Register\n')
            elif combinedMotifs:
                title = f'{self.titleCombined}'
            else:
                title = f'{enzName}\n{self.datasetTag}'
            # title += f'\n{dataType}'


        # Create heatmap
        cMapCustom = self.createCustomColorMap(colorType=dataType)

        # Convert the counts to a data frame for Seaborn heatmap
        if self.residueLabelType == 0:
            data.index = [residue[0] for residue in self.residues]
        elif self.residueLabelType == 1:
            data.index = [residue[1] for residue in self.residues]
        elif self.residueLabelType == 2:
            data.index = [residue[2] for residue in self.residues]

        # Define color bar limits
        cBarMax = np.max(data)
        cBarMax = math.ceil(cBarMax * 10) / 10 # Round up the 1st decimal
        if int(cBarMax * 10) % 2 != 0:
            # Set max to an even value
            cBarMax = (cBarMax * 10 + 1) / 10


        # Plot the heatmap with numbers centered inside the squares
        fig, ax = plt.subplots(figsize=self.figSizeEM)
        sns.heatmap(data, annot=True, fmt='.3f', cmap=cMapCustom, cbar=False,
                    linewidths=self.lineThickness-1, linecolor='black', square=False,
                    center=None, vmax = cBarMax, vmin = 0,
                    annot_kws={'fontweight': 'bold', 'size': self.labelSizeEM})
        ax.set_xlabel('Substrate Position', fontsize=self.labelSizeAxis)
        ax.set_ylabel('Residue', fontsize=self.labelSizeAxis)
        ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')


        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(self.lineThickness)

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks, width=self.lineThickness)
        ax.tick_params(axis='y', labelrotation=0)

        # Set x-ticks
        xTicks = np.arange(len(data.columns)) + 0.5
        ax.set_xticks(xTicks)
        ax.set_xticklabels(data.columns)

        # Set y-ticks
        yTicks = np.arange(len(data.index)) + 0.5
        ax.set_yticks(yTicks)
        ax.set_yticklabels(data.index)

        # Set the edge thickness
        for _, spine in ax.spines.items():
            spine.set_visible(True)

        # Modify the colorbar
        # cbar = heatmap.collections[0].colorbar
        # cbar.ax.tick_params(axis='y', which='major', labelsize=self.labelSizeTicks,
        #                     length=self.tickLength, width=self.lineThickness)
        # cbar.outline.set_linewidth(self.lineThickness)
        # cbar.outline.set_edgecolor('black')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.1)
        norm = plt.Normalize(vmin=0, vmax=cBarMax)
        cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cMapCustom), cax=cax)
        cbar.ax.tick_params(axis='y', which='major', labelsize=self.labelSizeTicks,
                            length=self.tickLength, width=self.lineThickness)
        cbar.outline.set_linewidth(self.lineThickness)
        cbar.outline.set_edgecolor('black')

        self.plotFig(plt=plt, fig=fig)

        # Save the figure
        figLabel = dataType
        if dataType == 'Standard Deviation':
            figLabel = 'StDev'
        if self.saveFigures:
            self.saveFigure(fig=fig, figType=figLabel, seqLen=len(xTicks),
                            combinedMotifs=combinedMotifs)



    def processSubstrates(self, motifs, subLabel, combinedMotifs=False,
                          predActivity=False, predModel=False):
        if predActivity:
            # Calculate: Motif enrichment
            for predType, predictions in motifs.items():
                print(f'Evaluating Predictions: {purple}{predType}{resetColor}\n'
                      f'Predictions: {type(predictions)}')
            # Plot: Word cloud
            self.plotWordCloud(substrates=motifs, combinedMotifs=combinedMotifs,
                               predActivity=predActivity, predModel=predModel)
        else:
            # Plot: Word cloud
            if self.plotFigWords:
                self.plotWordCloud(substrates=motifs, combinedMotifs=combinedMotifs)

        # Plot: Bar graphs
        if self.plotFigBars:
            self.plotBarGraph(substrates=motifs, dataType='Counts',
                              combinedMotifs=combinedMotifs)
            self.plotBarGraph(substrates=motifs, dataType='Counts',
                              combinedMotifs=combinedMotifs, plotAllSubs=True)
            # self.plotBarGraph(substrates=motifs, dataType='RF',
            #                   combinedMotifs=combinedMotifs)

        # PCA
        if self.plotFigPCA:
            # Convert substrate data to numerical
            tokensESM, subsESM, subCountsESM = self.ESM(
                substrates=motifs, subLabel=subLabel)

            # Cluster substrates
            subPopulations = self.plotPCA(substrates=motifs, data=tokensESM,
                                          indices=subsESM, N=subCountsESM,
                                          combinedMotifs=combinedMotifs)
            for NCluster, motifCluster in enumerate(subPopulations, start=1):


                # Plot: Bar graphs
                self.plotBarGraph(
                    substrates=motifCluster, dataType=f'Counts PCA {NCluster}',
                    combinedMotifs=combinedMotifs
                )
                self.plotBarGraph(
                    substrates=motifCluster, dataType=f'Counts PCA {NCluster}',
                    combinedMotifs=combinedMotifs, plotAllSubs=True
                )

                # Plot: Word cloud
                self.plotWordCloud(substrates=motifCluster, clusterNumPCA=NCluster,
                                   combinedMotifs=combinedMotifs)



    def predMotifEnrichment(self, motifES, scaledMatrix=False):
        print('=========================== Predict: Motif Enrichment '
              '===========================')
        if scaledMatrix:
            matrixType = 'Released Count Scaled Enrichment Scores'
            matrix = self.eMapReleasedScaled
        else:
            matrixType = 'Released Count Enrichment Scores'
            matrix = self.eMapReleased
        print(f'Prediction Matrix: {magenta}{matrixType}{resetColor}\n'
              f'{self.eMap}\n\n')

        predScores = {}
        maxScore, minScore = 0, 0
        normalizeValues = True

        for motif in motifES:
            score = 0
            for index, AA in enumerate(motif):
                position = self.xAxisLabelsMotif[index]
                score += matrix.loc[AA, position]
            # score = np.log2(score)
            predScores[motif] = score
            if score > maxScore:
                maxScore = score
            if score < minScore:
                minScore = score

        if normalizeValues:
            # Set minimum predicted score = 0
            if minScore < 0:
                for substrate, score in predScores.items():
                    newScore = score - minScore
                    predScores[substrate] = newScore

                    # Update max score
                    if newScore > maxScore:
                        maxScore = newScore
                minScore -= minScore
            else:
                for substrate, score in predScores.items():
                    newScore = score + minScore
                    predScores[substrate] = newScore
                maxScore += minScore

            # Normalize Values
            for substrate, score in predScores.items():
                predScores[substrate] = score / maxScore

        iteration = 0
        print(f'Predicted Scores:')
        for substrate, score in predScores.items():
            print(f'     {pink}{substrate}{resetColor}, '
                  f'Score: {red}{np.round(score, 2)}{resetColor}')
            iteration += 1
            if iteration >= self.printNumber:
                print('\n')
                break

        self.plotScatter(valuesExp=motifES, valuesPred=predScores, matrixType=matrixType)



    def plotScatter(self, valuesExp, valuesPred, matrixType,
                    color='#BF5700'):
        x, y = valuesExp.values(), valuesPred.values()

        # Normalize predicted values
        xMax, xMin, xNorm = max(x), min(x), []
        xAdj = xMax - xMin
        for value in x:
            xNorm.append((value - xMin) / xAdj)
        print(f'Experimental Values:\n'
              f'     Max: {red}{max(x)}{resetColor}'
              f'     Min: {red}{min(x)}{resetColor}\n'
              f'Normalized Experimental Values:\n'
              f'     Max: {red}{max(xNorm)}{resetColor}'
              f'     Min: {red}{min(xNorm)}{resetColor}\n\n'
              f'Predicted Values:\n'
              f'     Max: {red}{max(y)}{resetColor}'
              f'     Min: {red}{min(y)}{resetColor}\n\n')


        accuracy = 1.0
        label = f'R² = {np.round(accuracy, 3)}'
        title = (f'{self.enzymeName.replace(' - ', '\n')}\n'
                 f'{self.datasetTag}\n{matrixType}')

        # Create scatter plot
        fig, ax = plt.subplots(figsize=self.figSize)
        ax.scatter(xNorm, y, color=color)
        plt.title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        plt.xlabel('Motif Enrichment', fontsize=self.labelSizeAxis)
        plt.ylabel('Predicted Scores', fontsize=self.labelSizeAxis)

        # Define: Axis ticks
        ticks = np.linspace(0, 1, 5)
        tickLabels = ['0', '0.25', '0.5', '0.75', '1']

        # Set x-ticks
        ax.set_xticks(ticks)
        ax.set_xticklabels(tickLabels)
        ax.set_xlim(-0.03, 1.03)

        # Set y-ticks
        ax.set_yticks(ticks)
        ax.set_yticklabels(tickLabels)
        ax.set_ylim(-0.03, 1.03)

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks, width=self.lineThickness)

        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(self.lineThickness)

        self.plotFig(plt=plt, fig=fig)



    def decayRate(self, y):
        print('========================= Calculate: Exponential Decay '
              '==========================')
        print(f'Dataset: {purple}{self.datasetTag}{resetColor}\n')

        # Set: xValues 
        x = np.arange(0, len(y))
        smoothx = np.linspace(x[0], x[-1], 20)

        guessA, guessB, guessC = 4000, -0.005, 100
        guess = [guessA, guessB, guessC]

        exp_decay = lambda x, A, t, y0: A * np.exp(x * t) + y0

        params, cov = curve_fit(exp_decay, x, y, p0=guess)

        A, t, y0 = params

        print("A = %s\nt = %s\ny0 = %s\n" % (A, t, y0))



    def plotBarGraph(self, substrates, dataType, barColor='#BF5700',
                     barWidth=0.75, combinedMotifs=False,
                     plotAllSubs=False):
        print('================================ Plot: Bar Graph '
              '================================')
        print(f'Dataset: {purple}{self.datasetTag}{resetColor}')
        if plotAllSubs:
            barWidth = 1
            limitNSubs = len(substrates.keys())
            print(f'Plotting {red}{limitNSubs:,}{resetColor} substrates\n')
        else:
            limitNSubs = self.NSubBars
            print(f'Collecting the top {red}{limitNSubs:,}{resetColor} substrates\n')

        # Evaluate data
        totalCounts, x, y = sum(substrates.values()), [], []
        print(f'Substrates:')
        if 'counts' in dataType.lower():
            for i, (substrate, count) in enumerate(substrates.items(), start=1):
                print(f'     {blue}{substrate}{resetColor}, '
                      f'Counts: {red}{count:,}{resetColor}')
                if i >= self.printNumber:
                    break
            labelY = 'Counts'

            # Evaluate: Substrates
            if plotAllSubs:
                for i, (substrate, count) in enumerate(substrates.items(), start=1):
                    x.append(str(substrate))
                    y.append(count)
            else:
                for i, (substrate, count) in enumerate(substrates.items(), start=1):
                    x.append(str(substrate))
                    y.append(count)
                    if i == limitNSubs:
                        break

            # Evaluate: Y axis
            maxValue, yMin, mag = max(y), 0, 10
            magnitude = math.floor(math.log10(maxValue))
            if magnitude > 1:
                mag = 10 ** (magnitude - 1)
            yMax = math.ceil(maxValue / mag) * mag
        elif 'rf' in dataType.lower():
            for i, (substrate, count) in enumerate(substrates.items()):
                print(f'     {blue}{substrate}{resetColor}, '
                      f'RF: {red}{count / totalCounts:.{self.roundVal}f}{resetColor}')
                if i >= self.printNumber:
                    break
            labelY = 'Relative Frequency'

            # Evaluate: Substrates
            if plotAllSubs:
                for i, (substrate, count) in enumerate(substrates.items(), start=1):
                    x.append(str(substrate))
                    y.append(count / totalCounts)
            else:
                for i, (substrate, count) in enumerate(substrates.items(), start=1):
                    x.append(str(substrate))
                    y.append(count / totalCounts)
                    if i == limitNSubs:
                        break

            # Evaluate: Y axis
            maxValue = max(y)
            magnitude = math.floor(math.log10(maxValue))
            adjustedMax = maxValue * 10 ** abs(magnitude)
            yMax = math.ceil(adjustedMax) * 10 ** magnitude
            adjVal = 5 * 10 ** (magnitude - 1)
            yMaxAdjusted = yMax - adjVal
            if yMaxAdjusted > maxValue:
                yMax = yMaxAdjusted
            yMin = 0
        else:
            print(f'{orange}ERROR: What data type is: {cyan}{dataType}{resetColor}\n\n')
            return

        NSubs = len(x)
        print(f'Total Counts: {red}{totalCounts:,}{resetColor}\n'
              f'Unique Sequences: {red}{len(substrates.values()):,}{resetColor}\n')
        print(f'Number of plotted sequences: {red}{NSubs:,}{resetColor}\n\n')
        # print(f'Y Axis:\n'
        #       f'* Max: {yMax:,}\n'
        #       f'* Min: {yMin:,}\n')
        # print()

        # Define: Figure title
        enzName = self.enzymeName.replace(' - ', '\n')
        if self.releasedCounts:
            title = self.titleReleased
        elif combinedMotifs:
            title = self.titleCombined
        else:
            title = self.titleWords
            title += f'\nTop {NSubs} Substrates'
        title = title.replace(self.enzymeName, enzName)


        # Plot the data
        fig, ax = plt.subplots(figsize=self.figSize)
        bars = plt.bar(x, y, color=barColor, width=barWidth)
        plt.ylabel(labelY, fontsize=self.labelSizeAxis)
        plt.title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        plt.axhline(y=0, color='black', linewidth=self.lineThickness)
        plt.ylim(yMin, yMax)


        # Set: x ticks
        if plotAllSubs:
            # Determine x values
            magnitude = math.floor(math.log10(NSubs))
            step = 10 ** (magnitude) / 2
            xMax = 0
            while xMax < NSubs:
                xMax += step
            xTicks = np.arange(0, xMax + 1, step, dtype=int)

            ax.set_xticks(xTicks)
            ax.set_xticklabels(xTicks, ha='center')
            ax.set_xlim(-step / 10, xTicks[-1])
        else:
            xTicks = np.arange(0, NSubs)
            ax.set_xticks(xTicks)
            ax.set_xticklabels(x, rotation=0, ha='center')
            ax.set_xlim(left=xTicks[0] - barWidth, right=xTicks[-1] + barWidth)
        plt.xticks(rotation=90, ha='center')

        # Set: y ticks
        plotYTicks = True
        yTicks = []
        if max(y) == 1.0:
            yMax = 1.0
            if yMin == 0:
                yTicks = np.linspace(yMin, yMax, 6)
            else:
                plotYTicks = False
                dist = yMax - yMin
                vals = [0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55]
                for val in vals:
                    ratio = dist / val
                    if ratio == int(ratio) and ratio < 10:
                        plotYTicks = True
                        yTicks = np.arange(yMin, yMax + val, val)
                        break
            yMax += 0.1
        else:
            step = yMax / 10
            yTicks = np.linspace(yMin, yMax, 11)
        plt.ylim(yMin, yMax)
        if plotYTicks:
            plt.yticks(yTicks)

        # Set the edge color
        if not plotAllSubs:
            for bar in bars:
                bar.set_edgecolor('black')
                bar.set_linewidth(self.lineThickness)

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks, width=self.lineThickness)
        if not plotAllSubs and self.labelSizeTicks > 13:
            ax.tick_params(axis='x', which='major', labelsize=13)

        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(self.lineThickness)

        self.plotFig(plt=plt, fig=fig)

        # Save the figure
        if self.saveFigures:
            dataType = dataType.replace(' ', '_')
            if self.motifLen == None:
                seqLength = len(self.xAxisLabels)
            else:
                seqLength = self.motifLen

            # Define: Save location
            enzName = self.enzymeName.replace(' - ', '-').replace(' ', '_')
            enzName = enzName.replace('/', '-')
            datasetTag = self.datasetTagMotif.replace(' ', '_')
            figLabel = (f'{enzName}-Bars-{dataType}-{datasetTag}-{seqLength}AA-'
                        f'NBars_{NSubs}-N_{totalCounts}-MinCounts_{self.minSubCount}.png')
            if self.releasedCounts or combinedMotifs:
                figLabel = figLabel.replace(datasetTag,
                                            f'SubstrateProfile-{datasetTag}')
            if 'rf' in dataType.lower():
                figLabel = figLabel.replace(dataType, 'RF')
            saveLocation = os.path.join(self.pathSaveFigs, figLabel)
            saveLocation = saveLocation.replace(' ', '')

            # Save figure
            self.saveFig(fig=fig, path=saveLocation)



    def ESM(self, substrates, subLabel, useSubCounts=True):
        print('================== Convert Substrates To Numerical Values: ESM '
              '==================')
        print(f'Dataset: {purple}{self.datasetTag}{resetColor}\n\n'
              f'Collecting up to {red}{self.NSubsPCA:,}{resetColor} substrates\n'
              f'Total unique substrates: {red}{len(substrates):,}{resetColor}')

        # Extract: Datapoints
        iteration = 0
        collectedTotalValues = 0
        evaluateSubs = {}
        for substrate, value in substrates.items():
            evaluateSubs[str(substrate)] = value
            iteration += 1
            collectedTotalValues += value
            if iteration >= self.NSubsPCA:
                break
        sampleSize = len(evaluateSubs)
        print(f'Collected substrates:{red} {sampleSize:,}{resetColor}')
        if isinstance(collectedTotalValues, float):
              print(f'Total Values:{red} {round(collectedTotalValues, 1):,}'
                    f'{resetColor}\n\n')
        else:
            print(f'Total Values:{red} {collectedTotalValues:,}{resetColor}\n\n')

        # Step 1: Convert substrates to ESM model format and generate embeddings
        subs = []
        counts = []
        if useSubCounts:
            for index, (seq, count) in enumerate(evaluateSubs.items()):
                subs.append((f'Sub{index}', seq))
                counts.append(count)
        else:
            for index, seq in enumerate(evaluateSubs.keys()):
                subs.append((f'Sub{index}', seq))


        # Step 2: Load the ESM model and batch converter
        model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
        # model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()

        batch_converter = alphabet.get_batch_converter()


        # Step 3: Convert substrates to ESM model format and generate embeddings
        try:
            batchLabels, batchSubs, batchTokens = batch_converter(subs)
        except Exception as exc:
            print(f'{orange}ERROR: The ESM has failed to evaluate your substrates\n\n'
                  f'Exception:\n{exc}\n\n'
                  f'Suggestion:'
                  f'     Try replacing: {cyan}esm.pretrained.esm2_t36_3B_UR50D()'
                  f'{orange}\n'
                  f'     With: {cyan}esm.pretrained.esm2_t33_650M_UR50D()'
                  f'{resetColor}\n')
            sys.exit(1)

        print(f'Batch Tokens:{greenLight} {batchTokens.shape}{resetColor}\n'
              f'{greenLight}{batchTokens}{resetColor}\n\n')
        slicedTokens = pd.DataFrame(batchTokens[:, 1:-1],
                                    index=batchSubs,
                                    columns=subLabel)
        if useSubCounts:
            slicedTokens['Counts'] = counts
        print(f'Sliced Tokens:\n'
              f'{greenLight}{slicedTokens}{resetColor}\n\n')

        return slicedTokens, batchSubs, sampleSize



    def plotPCA(self, substrates, data, indices, N, combinedMotifs=False):
        print('====================================== PCA '
              '======================================')
        print(f'Dataset: {purple}{self.datasetTag}{resetColor}\n')
        print(f'Data:\n{data}\n\n')

        # Initialize lists for the clustered substrates
        self.selectedSubstrates = []
        self.selectedDatapoints = []
        rectangles = []

        # Define component labels
        pcaHeaders = []
        for componentNumber in range(1, self.numPCs + 1):
            pcaHeaders.append(f'PC{componentNumber}')
        headerCombinations = list(combinations(pcaHeaders, 2))

        # # Cluster the datapoints
        # Step 1: Apply PCA on the standardized data
        print(f'NPC: {self.numPCs}\n\n')
        pca = PCA(n_components=self.numPCs) # Adjust the number of components as needed
        scaler = StandardScaler()
        data = scaler.fit_transform(data)
        dataPCA = pca.fit_transform(data)
        # loadings = pca.components_.T

        # Step 2: Create a DataFrame for PCA results
        dataPCA = pd.DataFrame(dataPCA, columns=pcaHeaders, index=indices)
        print(f'PCA data:{red} # of components = {self.numPCs}\n'
              f'{greenLight}{dataPCA}{resetColor}\n\n')

        # Step 3: Print explained variance ratio
        varRatio = pca.explained_variance_ratio_ * 100
        print(f'Explained Variance Ratio: '
              f'{red}{" ".join([f"{x:.3f}" for x in varRatio])}{resetColor} %\n\n')


        # Define: Figure title
        if combinedMotifs:
            title = (f'{self.enzymeName}\n'
                     f'{self.datasetTagMotif}\n'
                     f'{N:,} Unique Substrates')
            if 'Excl' in self.datasetTagMotif:
                title = title.replace('Excl', 'Exclude')
        else:
            title = (f'{self.enzymeName}\n'
                     f'{self.datasetTag}\n'
                     f'{N:,} Unique Substrates')
            if 'Excl' in self.datasetTag:
                title = title.replace('Excl', 'Exclude')

        # Plot the data
        for components in headerCombinations:
            fig, ax = plt.subplots(figsize=self.figSize)

            def selectDatapoints(eClick, eRelease):
                # Function to update selection with a rectangle
                nonlocal ax, rectangles

                # Define x, y coordinates
                x1, y1 = eClick.xdata, eClick.ydata  # Start of the rectangle
                x2, y2 = eRelease.xdata, eRelease.ydata  # End of the rectangle

                # Collect selected datapoints
                selection = []
                selectedSubs = []
                for index, (x, y) in enumerate(zip(dataPCA.loc[:, 'PC1'],
                                                   dataPCA.loc[:, 'PC2'])):
                    if (min(x1, x2) <= x <= max(x1, x2) and
                            min(y1, y2) <= y <= max(y1, y2)):
                        selection.append((x, y))
                        selectedSubs.append(dataPCA.index[index])
                if selection:
                    self.selectedDatapoints.append(selection)
                    self.selectedSubstrates.append(selectedSubs)

                # Draw the boxes
                if self.selectedDatapoints:
                    for index, box in enumerate(self.selectedDatapoints):
                        # Calculate the bounding box for the selected points
                        padding = 0.05
                        xMinBox = min(x for x, y in box) - padding
                        xMaxBox = max(x for x, y in box) + padding
                        yMinBox = min(y for x, y in box) - padding
                        yMaxBox = max(y for x, y in box) + padding

                        # Draw a single rectangle around the bounding box
                        boundingRect = plt.Rectangle((xMinBox, yMinBox),
                                                     width=xMaxBox - xMinBox,
                                                     height=yMaxBox - yMinBox,
                                                     linewidth=2,
                                                     edgecolor='black',
                                                     facecolor='none')
                        ax.add_patch(boundingRect)
                        self.rectangles.append(boundingRect)

                        # Add text only if there are multiple boxes
                        if len(self.selectedDatapoints) > 1:
                            # Calculate the center of the rectangle for text positioning
                            centerX = (xMinBox + xMaxBox) / 2
                            centerY = (yMinBox + yMaxBox) / 2

                            # Number the boxes
                            text = ax.text(centerX, centerY, f'{index + 1}',
                                           horizontalalignment='center',
                                           verticalalignment='center',
                                           fontsize=25,
                                           color='#F8971F',
                                           fontweight='bold')
                            text.set_path_effects(
                                [path_effects.Stroke(linewidth=2, foreground='black'),
                                 path_effects.Normal()])
                plt.draw()
            plt.scatter(dataPCA[components[0]], dataPCA[components[1]],
                        c='#BF5700', edgecolor='black')
            plt.title(title, fontsize=self.labelSizeTitle, fontweight='bold')
            plt.xlabel(f'Principal Component {components[0][-1]} '
                       f'({np.round(varRatio[0], self.roundVal)} %)',
                       fontsize=self.labelSizeAxis)
            plt.ylabel(f'Principal Component {components[1][-1]} '
                       f'({np.round(varRatio[1], self.roundVal)} %)',
                       fontsize=self.labelSizeAxis)


            # Set tick parameters
            ax.tick_params(axis='both', which='major', length=self.tickLength,
                           labelsize=self.labelSizeTicks, width=self.lineThickness)

            # Set the thickness of the figure border
            for _, spine in ax.spines.items():
                spine.set_visible(True)
                spine.set_linewidth(self.lineThickness)


            # Create a RectangleSelector
            selector = RectangleSelector(ax,
                                         selectDatapoints,
                                         useblit=True,
                                         minspanx=5,
                                         minspany=5,
                                         spancoords='pixels',
                                         interactive=True)

            # Change rubber band color
            selector.set_props(facecolor='none', edgecolor='green', linewidth=3)

            self.plotFig(plt=plt, fig=fig)


        # Save the Figure
        if self.saveFigures:
            # Define: Save location
            figLabel = (f'{self.enzymeName}-PCA-{self.datasetTag}-'
                        f'{N}-MinCounts_{self.minSubCount}.png')
            if combinedMotifs :
                figLabel = figLabel.replace(self.datasetTag,
                                            f'Combined {self.datasetTag}')
            elif not combinedMotifs and len(self.motifIndexExtracted) > 1:
                figLabel = figLabel.replace(self.datasetTag,
                                            f'SubstrateProfile_{self.datasetTag}')
            saveLocation = os.path.join(self.pathSaveFigs, figLabel)
            saveLocation = saveLocation.replace(' ', '')

            # Save figure
            self.saveFig(fig=fig, path=saveLocation)


        # Create a list of collected substrate dictionaries
        if self.selectedSubstrates:
            collectedSubs = []
            for index, substrateSet in enumerate(self.selectedSubstrates):
                print(f'PCA Cluster: {greenLightB}{index + 1}{resetColor}')
                iteration = 0
                collectionSet = {}
                for substrate in substrateSet:
                    collectionSet[substrate] = substrates[substrate]

                # Sort collected substrates and add to the list
                collectionSet = dict(sorted(
                    collectionSet.items(), key=lambda x: x[1], reverse=True)
                )
                collectedSubs.append(collectionSet)

                # Print: Collected substrates
                for substrate, score in collectionSet.items():
                    if isinstance(score, float):
                        print(f'     {blue}{substrate}{resetColor}, '
                              f'Counts: {red}{round(score, 1):,}{resetColor}')
                    else:
                        print(f'     {blue}{substrate}{resetColor}, '
                              f'Counts: {red}{score:,}{resetColor}')
                    iteration += 1
                    if iteration >= self.printNumber:
                        print('\n')
                        break

            return collectedSubs
        else:
            return None



    def KLDivergence(self, P, Q, printProb, scaler):
        print('================================= KL Divergence '
              '=================================')
        P.columns = Q.columns
        if printProb:
            print(f'Baseline Probability Distribution:\n{Q}\n\n\n'
                  f'Probability Distribution:\n{P}\n\n')

        divergence = pd.DataFrame(0, columns=Q.columns,
                                  index=[self.datasetTag], dtype=float)
        divergenceMatrix = pd.DataFrame(0, columns=Q.columns,
                                        index=Q.index, dtype=float)

        for position in Q.columns:
            p = P.loc[:, position]
            q = Q.loc[:, position]
            divergence.loc[self.datasetTag, position] = (
                np.sum(np.where(p != 0, p * np.log2(p / q), 0)))

            for residue in Q.index:
                initial = Q.loc[residue, position]
                final = P.loc[residue, position]
                if initial == 0 or final == 0:
                    divergenceMatrix.loc[residue, position] = 0
                else:
                    divergenceMatrix.loc[residue, position] = (final *
                                                               np.log2(final / initial))

        # Scale the values
        if scaler is not None:
            for position in Q.columns:
                divergenceMatrix.loc[:, position] = (divergenceMatrix.loc[:, position] *
                                                     scaler.loc[position, 'ΔS'])

        print(f'{greyDark}KL Divergence:'
              f'{pink} Fixed Final Sort - {self.datasetTag}{resetColor}\n'
              f'{divergence}\n\n\n{greyDark}Divergency Matrix:'
              f'{pink} Fixed Final Sort - {self.datasetTag}'
              f'{resetColor}\n{divergenceMatrix.round(4)}\n\n')

        return divergenceMatrix, divergence



    def optimalWord(self, matrix, matrixType, maxResidues, dropPos,
                    printOptimalAA, normalizeValues):
        print('========================= Synthesize Optimal Substrates '
              '=========================')

        # Determine the OS
        combinations = 1
        optimalAA = []
        substratesOS = {}

        # Drop irrelevant positions
        if dropPos:
            dropColumns = []
            for indexDropPosition in dropPos:
                dropColumns.append(matrix.columns[indexDropPosition])
            # print(f'Dropping positions: {dropColumns}\n')
            matrix = matrix.drop(columns=dropColumns)
        # print(f'Values used to synthesize optimal substrates:\n{matrix}\n\n')



        for indexColumn, column in enumerate(matrix.columns):
            # Find the best residues at this position
            optimalAAPos = matrix[column].nlargest(maxResidues)

            # Filter the data
            for rank, (AA, score) in (enumerate(
                    zip(optimalAAPos.index, optimalAAPos.values), start=1)):
                if score <= 0:
                    optimalAAPos = optimalAAPos.drop(index=AA)
            optimalAA.append(optimalAAPos)


        if printOptimalAA:
            print(f'Optimal Residues: {purple}{self.enzymeName} - '
                  f'Fixed {self.datasetTag}{resetColor}')
            for index, data in enumerate(optimalAA, start=1):
                # Determine the number of variable residues at this position
                numberAA = len(data)
                combinations *= numberAA

                # Define substrate position
                positionSub = self.xAxisLabels[index-1]
                print(f'Position: {purple}{positionSub}{resetColor}')

                for AA, datapoint in data.items():
                    print(f'     {AA}:{red} {datapoint:.6f}{resetColor}')
                print('\n')
        else:
            for index, data in enumerate(optimalAA, start=1):
                # Determine the number of variable residues at this position
                combinations *= len(data)

        print(f'Possible Substrate Combinations:{pink} {combinations:,}{resetColor}\n')

        # Use the optimal residues to determine OS
        substrate = ''
        score = 0
        for index, data in enumerate(optimalAA, start=1):
            # Select the top-ranked AA for each position
            topAA = data.idxmax()
            topES = data.max()

            # Construct the OS
            substrate += ''.join(topAA)
            score += topES

        # Update OS dictionary
        substratesOS[substrate] = score

        # Create additional substrates
        for indexColumn, column in enumerate(matrix.columns):

            # Collect new substrates to add after the iteration
            newSubstratesList = []

            for substrate, ESMax in list(substratesOS.items()):
                AAOS = substrate[indexColumn]
                scoreSubstrate = optimalAA[indexColumn][AAOS]

                # Access the correct filtered data for the column
                optimalAAPos = optimalAA[indexColumn]
                # print(f'optimalAA:\n{optimalAA}\n\n'
                #       f'     optimalAAPos:\n{optimalAAPos}\n\n')

                for AA, score in optimalAAPos.items():
                    # print(f'AA: {AA}\nScore: {score}\n\n')
                    if AA != AAOS:

                        # Replace AAOS with AA
                        newSubstrate = (substrate[:indexColumn] + AA +
                                        substrate[indexColumn + 1:])
                        newES = ESMax + (score - scoreSubstrate)
                        # print(f'{greyDark}New Substrate{resetColor}:'
                        #       f'{greenLightB} {newSubstrate}{resetColor}, '
                        #       f'ES:{red} {newES}{resetColor}\n'
                        #       f'     Residue Score New:{red} {score}{resetColor}\n'
                        #       f'     Residue Score Old:{red} {scoreSubstrate}'
                        #       f'{resetColor}\n\n')

                        # Collect new substrate and ES to add later
                        newSubstratesList.append((newSubstrate, newES))
            # Update substratesOS with new substrates after the iteration
            for newSubstrate, newES in newSubstratesList:
                substratesOS[newSubstrate] = newES

        substratesOS = dict(sorted(substratesOS.items(),
                                   key=lambda x: x[1], reverse=True))
        if normalizeValues:
            _, topScore = substratesOS.popitem(last=False) # Top score
            print(f'Top Score:{red} {topScore}{resetColor}\n')

            # Normalize the values
            substratesOS = {key: value / topScore for key, value in substratesOS}
            # substratesOS = sorted(substratesOS.items(), key=lambda x: x[1], reverse=True)

        print(f'Top {self.printNumber} Optimal Substrates:{greyDark} '
              f'{matrixType}{resetColor}')
        iteration = 0
        for substrate, ES in substratesOS.items():
            print(f'     Substrate:{red} {substrate}{resetColor}, '
                  f'ES:{red} {ES:.6f}{resetColor}')
            iteration += 1
            if iteration == self.printNumber:
                break

        print(f'\nNumber of synthesized substrates:'
              f'{pink} {len(substratesOS):,}{resetColor}\n\n')

        return substratesOS



    def substrateEnrichment(self, initialSubs, finalSubs,
                            NSubs, saveData):
        print('========================= Evaluate Substrate Enrichment '
              '=========================')
        if self.datasetTag == None:
            datasetType = 'NNS'
        else:
            datasetType = f'{self.datasetTag}'

        # Define headers
        headersInitial = ['Initial Subs', 'Counts']
        headersFinal = ['Final Subs', 'Counts']
        headersEnriched = ['Enriched Subs', 'log₂(RF Final / RF Initial)']
        # headerWidth = {1: 14.6,
        #                4: 14.6,
        #                7: 14.6,
        #                8: 24} # Adjust column width at these indices for an excel sheet


        iteration = 0
        totalSubstratesInitial = 0
        totalUniqueSubstratesInitial = len(initialSubs)
        totalsubstrates = 0
        totalUniquesubstrates = len(finalSubs)

        # Evaluate the substrates
        if self.datasetTag == None:
            # Process: initial sort
            print(f'Ranked Substrates: {purple}Initial Sort{resetColor} -'
                  f'{red} NNS{resetColor}')
            if totalUniqueSubstratesInitial >= self.printNumber:
                for substrate, count in initialSubs.items():
                    if iteration <= self.printNumber:
                        print(f'     {pink}{substrate}{resetColor}, '
                              f'Counts: {red}{count:,}{resetColor}')
                        iteration += 1
                        totalSubstratesInitial += count
                    else:
                        totalSubstratesInitial += count
            else:
                print(f'{orange}The number of unique substrates '
                      f'({red}{totalUniqueSubstratesInitial}'
                      f'{orange}) is less than the number you requested to be see '
                      f'({red}{self.printNumber}{orange}){resetColor}')
                for substrate, count in initialSubs.items():
                    print(f'     {pink}{substrate}{resetColor}, Counts: {red}{count:,}'
                          f'{resetColor}')
                    totalSubstratesInitial += count
            iteration = 0
            # Print: Dataset totals
            print(f'\n     Total substrates {purple}Initial Sort{resetColor}:'
                  f'{red} {totalSubstratesInitial:,}{resetColor}\n'
                  f'     Unique substrates {purple}Initial Sort{resetColor}:'
                  f'{red} {totalUniqueSubstratesInitial:,}{resetColor}\n\n')


            # Process: final sort
            print(f'Ranked Substrates: {purple}Final Sort{resetColor} -'
                  f'{red} NNS{resetColor}')
            if totalUniquesubstrates >= self.printNumber:
                for substrate, count in finalSubs.items():
                    if iteration <= self.printNumber:
                        print(f'     {pink}{substrate}{resetColor}, '
                              f'Counts: {red}{count:,}{resetColor}')
                        iteration += 1
                        totalsubstrates += count
                    else:
                        totalsubstrates += count
            else:
                print(f'{orange}The number of unique substrates '
                      f'({red}{totalUniquesubstrates}'
                      f'{orange}) is less than the number you requested to be see '
                      f'({red}{self.printNumber}{orange}){resetColor}\n')
                for substrate, count in finalSubs.items():
                    print(f'     {pink}{substrate}{resetColor}, '
                          f'Counts: {red}{count:,}{resetColor}')
                    totalsubstrates += count
            iteration = 0
            # Print: Dataset totals
            print(f'\n     Total substrates {purple}Final Sort{resetColor}:'
                  f'{red} {totalsubstrates:,}{resetColor}\n'
                  f'     Unique substrates {purple}Final Sort{resetColor}:'
                  f'{red} {totalUniquesubstrates:,}{resetColor}\n\n')
        else:
            fixedSort = True

            # Print: Initial sort
            print(f'Ranked Substrates: {purple}Initial Sort{resetColor}{resetColor}')
            if totalUniqueSubstratesInitial >= self.printNumber:
                for substrate, count in initialSubs.items():
                    if iteration < self.printNumber:
                        print(f'     {pink}{substrate}{resetColor}, '
                              f'Counts: {red}{count:,}{resetColor}')
                        iteration += 1
                        totalSubstratesInitial += count
                    else:
                        totalSubstratesInitial += count
            else:
                print(f'{orange}The number of unique substrates '
                      f'({red}{totalUniqueSubstratesInitial}'
                      f'{orange}) is less than the number you requested to be see '
                      f'({red}{self.printNumber}{orange}){resetColor}')
                for substrate, count in initialSubs.items():
                    print(f'     {pink}{substrate}{resetColor}, '
                          f'Counts: {red}{count:,}{resetColor}')
                    totalSubstratesInitial += count
            print(f'\n     Total substrates {purple}Initial Sort{resetColor}:'
                  f'{red} {totalSubstratesInitial:,}{resetColor}\n'
                  f'     Unique substrates {purple}Initial Sort{resetColor}:'
                  f'{red} {totalUniqueSubstratesInitial:,}{resetColor}\n\n')
            iteration = 0


            # Process: Final sort
            print(f'Ranked Substrates: {purple}Final Sort{resetColor} -'
                  f'{red} Fixed {self.datasetTag}{resetColor}')
            if totalUniquesubstrates >= self.printNumber:
                for substrate, count in finalSubs.items():
                    if iteration < self.printNumber:
                        print(f'     {pink}{substrate}{resetColor}, '
                              f'Counts: {red}{count:,}{resetColor}')
                        iteration += 1
                        totalsubstrates += count
                    else:
                        totalsubstrates += count
            else:
                print(f'{orange}The number of unique substrates '
                      f'({red}{totalUniquesubstrates}'
                      f'{orange}) is less than the number you requested to be see '
                      f'({red}{self.printNumber}{orange}){resetColor}')
                for substrate, count in finalSubs.items():
                    print(f'     {pink}{substrate}{resetColor}, '
                              f'Counts: {red}{count:,}{resetColor}')
                    totalsubstrates += count
            print(f'\n     Total substrates {purple}Final Sort{resetColor}:'
                  f'{red} {totalsubstrates:,}{resetColor}\n'
                  f'     Unique substrates {purple}Final Sort{resetColor}:'
                  f'{red} {totalUniquesubstrates:,}{resetColor}\n\n')
            iteration = 0


        # Calculate: Substrate enrichment
        enrichedSubs = {}
        setMinCountFinal = False
        if setMinCountFinal:
            print(f'Mininum Substrate Count:{red} {self.minSubCount}{resetColor}')
            for substrate, count in finalSubs.items():
                if count < self.minSubCount:
                    continue

                if substrate in initialSubs.keys():
                    countInitial = initialSubs[substrate]
                else:
                    countInitial = 1
                probFinal = count / totalsubstrates
                probInitial = countInitial / totalSubstratesInitial
                enrichment = np.log2(probFinal/probInitial)
                enrichedSubs[substrate] = enrichment
        else:
            for substrate, count in finalSubs.items():
                if substrate in initialSubs.keys():
                    countInitial = initialSubs[substrate]
                else:
                    countInitial = 1
                probFinal = count / totalsubstrates
                probInitial = countInitial / totalSubstratesInitial
                enrichment = np.log2(probFinal/probInitial)
                enrichedSubs[substrate] = enrichment
        # sys.exit(1)

        # Sort enrichment dictionary
        enrichedSubs = dict(sorted(enrichedSubs.items(),
                                   key=lambda x: x[1], reverse=True))

        # Print: top enriched substrates
        print(f'{purple}Enriched substrates{resetColor}:')
        for substrate, score in enrichedSubs.items():
            if iteration < self.printNumber:
                print(f'     {magenta}{substrate}{resetColor}, '
                      f'ES: {red}{score:.3f}{resetColor}')
                iteration += 1
            else:
                break
        print('\n')


        if saveData:
            # Define file path
            filePathCSV = os.path.join(
                self.pathData,
                f'{self.enzymeName}-EnrichedSubs-{self.datasetTag}-'
                f'MinCounts_{self.minSubCount}.csv')

            # Convert dictionaries to a data frame and save as an Excel file
            clipDataset = False
            if totalSubstratesInitial > NSubs:
                clipDataset = True

                # Did you ask for too many substrates?
                if NSubs == 10**6:
                    NSubs -= 1
                elif NSubs > 10**6:
                    print(f'     The list of substrates in the  {purple}Initial Sort'
                          f'{resetColor} that you attempted to save '
                          f'({red}N = {NSubs:,}{resetColor}) '
                          f'is to large to save in an Excel file.\n'
                          f'          Extracting the first{red} {10**6:,}{resetColor} '
                          f'substrates and discarding the rest.\n')
                    NSubs = 10**6

                initialSubsDF = pd.DataFrame.from_dict(dict(
                    list(initialSubs.items())[:NSubs]).items())
                initialSubsDF.columns = headersInitial
            else:
                initialSubsDF = pd.DataFrame.from_dict(initialSubs.items())
                initialSubsDF.columns = headersInitial
            if totalsubstrates > NSubs:
                clipDataset = True

                # Did you ask for too many substrates?
                if NSubs > 10**6:
                    print(f'     The list of substrates in the  {purple}Final Sort'
                          f'{resetColor} that you attempted to save '
                          f'({red}N = {NSubs:,}{resetColor}) '
                          f'is to large to save in an Excel file.\n'
                          f'          Extracting the first{red} {10**6:,}{resetColor} '
                          f'substrates and discarding the rest.\n')
                    NSubs = 10**6

                finalSubsDF = pd.DataFrame.from_dict(dict(
                    list(finalSubs.items())[:NSubs]).items())
                finalSubsDF.columns = headersFinal

                enrichedSubsDF = pd.DataFrame.from_dict(dict(
                    list(enrichedSubs.items())[:NSubs]).items())
                enrichedSubsDF.columns = headersEnriched
            else:
                finalSubsDF = pd.DataFrame.from_dict(finalSubs.items())
                finalSubsDF.columns = headersFinal

                enrichedSubsDF = pd.DataFrame.from_dict(enrichedSubs.items())
                enrichedSubsDF.columns = headersEnriched
            if clipDataset:
                print(f'Saving dataset with a maximum of{red} {NSubs:,}'
                      f'{resetColor} substrates.\n\n'
                      f'Any data returned from this function that will be used for further '
                      f'analysis\nwill include the complete dataset\n\n')

            # Print: Sample sizes
            print(f'{greyDark}Dataset Size{resetColor}:\n'
                  f'     Initial Substrates: {red}{len(initialSubsDF):,}{resetColor}\n'
                  f'     Final Substrates: {red}{len(finalSubsDF):,}{resetColor}\n'
                  f'     Enriched Substrates: {red}{len(enrichedSubsDF):,}{resetColor}\n')

            if os.path.exists(filePathCSV):
                print(f'{orange}The{red} {datasetType}{orange} dataset at was '
                      f'found at the path:'
                      f'\n     {resetColor}{filePathCSV}\n\n'
                      f'{orange}The file was not overwritten{resetColor}\n\n')
            else:
                print(f'Saving the{red} {datasetType}{resetColor} dataset at the path:'
                      f'\n     {filePathCSV}{resetColor}\n\n')

                # Combine the data frames
                if datasetType == 'NNS':
                    enrichedSubsCSV = pd.concat([initialSubsDF,
                                                 finalSubsDF,
                                                 enrichedSubsDF],
                                                axis=1)
                    enrichedSubsCSV.to_csv(filePathCSV, index=True)
                else:
                    enrichedSubsCSV = pd.concat([finalSubsDF, enrichedSubsDF],
                                                axis=1)
                    enrichedSubsCSV.to_csv(filePathCSV, index=True)

        return enrichedSubs



    def plotSubstratePopulations(self, substrates, probInit, clusterIndex, numClusters,
                                 datasetTag, saveTag):
        print('=============================== Plot PCA Clusters '
              '===============================')
        print(f'Dataset: {purple}{datasetTag}{resetColor}\n'
              f'Save Tag: {purple}{saveTag}{resetColor}\n\n'
              f'Total Clusters:{red} {numClusters}{resetColor}\n'
              f'     Cluster Number:{red} {clusterIndex + 1}{resetColor}\n\n')

        # Define figure titles
        if numClusters == 1:
            figureTitle = f'PCA Population\n{self.datasetTag}'
        else:
            figureTitle = (f'PCA Population {clusterIndex + 1}\n'
                           f'{self.datasetTag}')

        # Count fixed substrates
        countFullSubstrate = False
        if countFullSubstrate:
            countsFinal, countsTotalFinal = self.countResidues(substrates=substrates,
                                                              datasetType='')
        else:
            countsFinal, countsTotalFinal = self.countResidues(substrates=substrates,
                                                              datasetType='')
        # self.sampleSizeUpdate(NSubs=countsTotalFinal, sortType='Final Sort',
        #                      datasetTag=self.datasetTag)

        # Adjust the zero counts at nonfixed positions
        countsFinalAdjusted = countsFinal.copy()
        # print(f'Counts:\n{countsFinalAdjusted}\n')
        # if inAdjustZeroCounts:
        #     for indexColumn in countsFinalAdjusted.columns:
        #         for AA in countsFinalAdjusted.index:
        #             if countsFinalAdjusted.loc[AA, indexColumn] == 0:
        #                 countsFinalAdjusted.loc[AA, indexColumn] = 1
        #     print(f'Adjusted Final Counts:{pink} {self.enzymeName}\n'
        #           f'{red}{countsFinalAdjusted}{resetColor}\n\n')

        # Calculate: RF
        probFinal = self.calculateRF(
            counts=countsFinal, N=countsTotalFinal, fileType='Final Sort'
        )
        probFinalAdjusted = self.calculateRF(
            counts=countsFinalAdjusted, N=countsTotalFinal, fileType='Final Sort'
        )

        # if 0 == 0:
        #     # Plot: Positional entropy
        #     self.plotEntropy(entropy=self.entropy)

        # Calculate: Enrichment scores
        fixedFramePopES = self.calculateEnrichment(
            rfInitial=probInit, rfFinal=probFinal, plotFigures=False
        )
        fixedFramePopESAdjusted = self.calculateEnrichment(
            rfInitial=probInit, rfFinal=probFinalAdjusted, plotFigures=False
        )

        # Calculate: Enrichment scores scaled
        fixedFramePCAESScaled = pd.DataFrame(0.0, index=fixedFramePopES.index,
                                             columns=fixedFramePopES.columns)
        fixedFramePCAESScaledAdjusted = pd.DataFrame(0.0, index=fixedFramePopES.index,
                                                     columns=fixedFramePopES.columns)

        # Scale enrichment scores with Shannon Entropy
        for positon in fixedFramePopES.columns:
            fixedFramePCAESScaled.loc[:, positon] = (fixedFramePopES.loc[:, positon] *
                                                     self.entropy.loc[
                                                         positon, 'ΔS'])
            fixedFramePCAESScaledAdjusted.loc[:, positon] = (
                    fixedFramePopESAdjusted.loc[:, positon] *
                    self.entropy.loc[positon, 'ΔS'])
        print(f'Motif:{greenLight} Scaled{resetColor}\n{fixedFramePCAESScaled}\n\n')
        yMax = max(fixedFramePCAESScaled[fixedFramePopES > 0].sum())
        yMin = min(fixedFramePCAESScaled[fixedFramePopES < 0].sum())

        # Plot data
        self.plotEnrichmentScores(dataType='Enrichment', pca=figureTitle)

        self.plotEnrichmentLogo(pca=figureTitle)

        # Plot: Word cloud
        self.plotWordCloud(substrates=substrates, pca=figureTitle)



    def plotMatrix(self, data, figLabel, totalCounts=False, printData=True):
        # Create heatmap
        cMapCustom = self.createCustomColorMap(colorType='Counts')

        # if 'prediction matrix' in figLabel.lower():
        #     cBarMax = 1
        cBarMax = np.ceil(data.values.max() * 10) / 10
        cBarMin = 0
        if printData:
            if totalCounts:
                print(f'Data:\n{data}\n\nMin: {int(cBarMin):,}\nMax: {int(cBarMax):,}')
            else:
                print(f'Data:\n{data}\n\nMin: {cBarMin:,}\nMax: {cBarMax:,}')

        # Convert the counts to a data frame for Seaborn heatmap
        if self.residueLabelType == 0:
            data.index = [residue[0] for residue in self.residues]
        elif self.residueLabelType == 1:
            data.index = [residue[1] for residue in self.residues]
        elif self.residueLabelType == 2:
            data.index = [residue[2] for residue in self.residues]

        # Set figure title
        enzName = self.enzymeName.replace(' - ', '\n')
        if totalCounts:
            title = f'{enzName}\n{figLabel}\nN={totalCounts:,}'
        else:
            title = f'{enzName}\n{figLabel}'


        # Plot the heatmap with numbers centered inside the squares
        fig, ax = plt.subplots(figsize=self.figSize)
        if totalCounts:
            sns.heatmap(data, annot=True, fmt=',d', cmap=cMapCustom, cbar=False,
                        linewidths=self.lineThickness-1, linecolor='black', square=False,
                        center=None, annot_kws={'fontweight': 'bold'})
        else:
            sns.heatmap(data, annot=True, fmt='.4f', cmap=cMapCustom,
                        cbar=False, linewidths=self.lineThickness - 1,
                        linecolor='black', square=False, center=None,
                        annot_kws={'fontweight': 'bold', 'size': self.labelSizeEM})
        ax.set_xlabel('Substrate Position', fontsize=self.labelSizeAxis)
        ax.set_ylabel('Residue', fontsize=self.labelSizeAxis)
        ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')

        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(self.lineThickness)

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks, width=self.lineThickness)
        ax.tick_params(axis='y', labelrotation=0)

        # Set x-ticks
        xTicks = np.arange(len(data.columns)) + 0.5
        ax.set_xticks(xTicks)
        ax.set_xticklabels(data.columns)

        # Set y-ticks
        yTicks = np.arange(len(data.index)) + 0.5
        ax.set_yticks(yTicks)
        ax.set_yticklabels(data.index)


        for _, spine in ax.spines.items():
            spine.set_visible(True)

        # Colormap
        if isinstance(cMapCustom, Colormap):
            cmap = cMapCustom
        else:
            # Look up premade colormap
            cmap = matplotlib.colormaps[cMapCustom]
        cmap.set_bad(color='lightgrey') # Set invalid values to gray

        # Modify the colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.1)
        norm = plt.Normalize(vmin=cBarMin, vmax=cBarMax)
        cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cMapCustom),
                            cax=cax)
        cbar.ax.tick_params(axis='y', which='major', labelsize=self.labelSizeTicks,
                            length=self.tickLength, width=self.lineThickness)
        cbar.outline.set_linewidth(self.lineThickness)
        cbar.outline.set_edgecolor('black')

        self.plotFig(plt=plt, fig=fig)

        # Save the figure
        if self.saveFigures:
            enzName = enzName.replace(' ', '_')
            # Define: Save location
            if totalCounts:
                figLabel = (f'{enzName}-Counts-{self.datasetTag}-'
                            f'MinCounts_{self.minSubCount}.png')
            else:
                figLabel = (f'{enzName}-Prediction Matrix-{self.datasetTag}-'
                            f'{len(xTicks)}AA-MinCounts_{self.minSubCount}.png')
            if self.releasedCounts:
                figLabel = figLabel.replace(self.datasetTag,
                                            f'SubstrateProfile-{self.datasetTag}')
            if '/' in figLabel:
                figLabel = figLabel.replace('/', '_')
            saveLocation = os.path.join(self.pathSaveFigs, figLabel)
            saveLocation = saveLocation.replace(' ', '')

            # Save figure
            self.saveFig(fig=fig, path=saveLocation)



    def calculateEntropy(self, rf, fixFullFrame=None,
                         combinedMotifs=False, manualEntropy=False,
                         manualFrame=None, plotFig=False):
        print('============================== Calculate: Entropy '
              '===============================')
        print(f'Dataset: {purple}{self.datasetTag}{resetColor}\n'
              f'Total substrate count: {red}{self.nSubsFinal:,}{resetColor}\n'
              f'Unique Substrates: {red}{self.nSubsFinalUniqueSeqs:,}{resetColor}\n')

        self.entropy = pd.DataFrame(0.0, index=rf.columns, columns=['ΔS'])
        self.entropyMax = np.log2(len(rf.index))
        for indexColumn in rf.columns:
            S = 0
            for indexRow, probRatio in rf.iterrows():
                prob = probRatio[indexColumn]
                if prob == 0:
                    continue
                else:
                    S += -prob * np.log2(prob)
            self.entropy.loc[indexColumn, 'ΔS'] = self.entropyMax - S
        print(f'{self.entropy}\n\nMax Entropy: {self.entropyMax.round(6)}\n\n')

        # Identify motif frame
        if not manualEntropy:
            if fixFullFrame is not None:
                self.identifyMotif(fixFullFrame=fixFullFrame)
        else:
            print('================================ Identify Motif '
                  '=================================')
            print(f'Selecting a {yellow}user defined motif{resetColor}')
            self.subFrame = pd.DataFrame(1, index=manualFrame, columns=['ΔS'])
            print(f'Ranked Substrate Frame:\n'
                  f'{blue}{self.subFrame}{resetColor}\n\n')

        if self.plotFigEntropy or plotFig:
            self.plotEntropy(entropy=self.entropy, combinedMotifs=combinedMotifs)

        return self.entropy



    def plotEntropy(self, entropy, combinedMotifs=False):
        enzName = self.enzymeName.replace(' - ', '\n')
        if self.filterSubs:
            title = f'{enzName}\n{self.datasetTag}'
        else:
            title = f'{enzName}'
        if self.releasedCounts:
            title = self.titleReleased
        if combinedMotifs and not self.releasedCounts:
            title = title.replace(self.datasetTag,
                                  f'Combined Filter\n{self.datasetTag}')
        if (len(self.datasetTag) > 50 and self.datasetTag in title and
                'exclude' in title.lower() and 'fixed' in title.lower()):
            title = f'{enzName}\n{self.datasetTag.replace('Fixed', '\nFixed')}'

        # Figure parameters
        yMax = self.entropyMax + 0.2

        # Map entropy values to colors using the colormap
        colors = [(0, 'navy'),
                  (0.3/self.entropyMax, 'navy'),
                  (0.7/self.entropyMax, 'dodgerblue'),
                  (0.97/self.entropyMax, 'white'),
                  (0.98/self.entropyMax, 'white'),
                  (1.0/self.entropyMax, 'white'),
                  (1.65/self.entropyMax, 'red'),
                  (3/self.entropyMax, 'firebrick'),
                  (1, 'darkred')]
        colorBar = LinearSegmentedColormap.from_list('custom_colormap', colors)


        # Map entropy values to colors using the colormap
        normalize = Normalize(vmin=0, vmax=yMax) # Normalize the entropy values
        cMap = [colorBar(normalize(value)) for value in entropy['ΔS'].astype(float)]

        # Plotting the entropy values as a bar graph
        fig, ax = plt.subplots(figsize=self.figSize)
        plt.bar(entropy.index, entropy['ΔS'], color=cMap,
                edgecolor='black', linewidth=self.lineThickness, width=0.8)
        plt.xlabel('Substrate Position', fontsize=self.labelSizeAxis)
        plt.ylabel('ΔS', fontsize=self.labelSizeAxis, rotation=0, labelpad=15)
        plt.title(title, fontsize=self.labelSizeTitle, fontweight='bold')


        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks)

        # Set x-ticks
        xTicks = np.arange(0, len(entropy.iloc[:, 0]), 1)
        ax.set_xticks(xTicks)
        ax.set_xticklabels(entropy.index, rotation=0, ha='center')
        for tick in ax.xaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness) # Set tick width

        # Set y-ticks
        yTicks = range(0, 5)
        yTickLabels = [f'{tick:.0f}' if tick != 0 else f'{int(tick)}' for tick in yTicks]
        ax.set_yticks(yTicks)
        ax.set_yticklabels(yTickLabels)
        for tick in ax.yaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness) # Set tick width

        # Set the edge thickness
        for spine in ax.spines.values():
            spine.set_linewidth(self.lineThickness)

        # Set axis limits
        ax.set_ylim(0, yMax)

        # Add color bar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.1)
        cbar = plt.colorbar(plt.cm.ScalarMappable(norm=normalize, cmap=colorBar), cax=cax)
        cbar.ax.tick_params(axis='y', which='major', labelsize=self.labelSizeTicks,
                            length=self.tickLength, width=self.lineThickness)
        for tick in cbar.ax.yaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness) # Set tick width
        cbar.outline.set_linewidth(self.lineThickness)

        self.plotFig(plt=plt, fig=fig)


        # Save the figure
        if self.saveFigures:
            self.saveFigure(fig=fig, figType='∆S', seqLen=len(xTicks),
                            combinedMotifs=combinedMotifs)



    def plotLibraryAADist(self, rfInitial, rfFinal, codonType, datasetTag,
                          skipInitial=False):
        # Inspect data
        if rfInitial is None and rfFinal is None:
            print(f'{orange}ERROR: both of the inputs for rfInitial and '
                  f'rfFinal cannot be None.{resetColor}\n')
            sys.exit(1)

        # Initialize parameters
        plotInitial, plotFinal = False, False
        maxInitial = 0
        maxInitialAdj = 0
        maxFinal = 0
        maxFinalAdj = 0

        # Determine yMax
        if rfInitial is not None:
            plotInitial = True
            numPos = rfInitial.shape[1]
            numAA = rfInitial.shape[0]
            maxInitial = rfInitial.values.max()
            maxInitialAdj = np.ceil(maxInitial * 10) / 10
        if rfFinal is not None:
            plotFinal = True
            numPos = rfFinal.shape[1]
            numAA = rfFinal.shape[0]
            maxFinal = rfFinal.values.max()
            maxFinalAdj= np.ceil(maxFinal * 10) / 10
        if maxFinalAdj > maxInitialAdj:
            yMax = maxFinalAdj
            maxY = maxFinal
        else:
            yMax = maxInitialAdj
            maxY = maxInitial
        if yMax > 0.6:
            tickStepSize = 0.1
        else:
            tickStepSize = 0.05
        if yMax != 1.0:
            if yMax < maxY:
                while yMax < maxY:
                    yMax += tickStepSize
            yMax = round(yMax, 1)
        #yMax = 0.15
        tickStepSize = 0.05


        def plotFig(rf, sortType):
            print('======================= Plot: AA Probability Distribution '
                  '=======================')
            if codonType == datasetTag:
                print(f'Plotting Probability Distribution:'
                      f' {purple}{codonType} codon{resetColor}')
            else:
                print(f'Plotting Probability Distribution:'
                      f' {purple}{self.enzymeName} {sortType}{resetColor}')
            print(f'{rf}\n')

            enzName = self.enzymeName.replace(' - ', '-')
            if codonType == datasetTag:
                title = f'\n\n{codonType} Codon'
            elif 'initial' in sortType.lower():
                title = f'\n{enzName}\nUnsorted Library'
            else:
                title = f'\n{enzName}\nSorted Library'

            fig, ax = plt.subplots(figsize=self.figSize)
            plt.ylabel('Relative Frequency', fontsize=self.labelSizeAxis)
            plt.title(title, fontsize=self.labelSizeTitle, fontweight='bold')


            # Set tick parameters
            ax.tick_params(axis='both', which='major', length=self.tickLength,
                           width=self.lineThickness)

            # Set x-ticks
            if codonType == datasetTag:
                widthBar = 9
            else:
                widthBar = 2
            spacing = widthBar * numPos
            widthCluster = spacing + 5
            indices = np.arange(numAA) * widthCluster
            midPoint = (numPos - 1) / 2 * widthBar
            xTicks = indices + midPoint
            ax.set_xticks(xTicks)
            ax.set_xticklabels(rf.index, rotation=0, ha='center',
                               fontsize=self.labelSizeTicks)

            # Set y-ticks
            yTicks = np.arange(0, yMax + tickStepSize, tickStepSize)
            yTickLabels = [f'{tick:.0f}' if tick == 0 or tick == 1
                           else f'{tick:.2f}' for tick in yTicks]
            ax.set_yticks(yTicks)
            ax.set_yticklabels(yTickLabels, fontsize=self.labelSizeTicks)
            for tick in ax.yaxis.get_major_ticks():
                tick.tick1line.set_markeredgewidth(self.lineThickness)

            # Set the edge color
            for index, AA in enumerate(rf.index):
                xPos = indices[index] + np.arange(numPos) * widthBar
                if AA == 'F' or AA == 'W' or AA == 'Y': # AA == 'N' or AA == 'Q' or
                    ax.bar(xPos, rf.loc[AA], widthBar, label=AA,
                           color=self.colorsAA[AA], edgecolor='dimgray')
                else:
                    ax.bar(xPos, rf.loc[AA], widthBar, label=AA,
                           color=self.colorsAA[AA], edgecolor='black')

            # Set the edge thickness
            for spine in ax.spines.values():
                spine.set_linewidth(self.lineThickness)

            # Set axis limits
            plt.ylim(0, yMax)

            self.plotFig(plt=plt, fig=fig)

            if self.saveFigures:
                # Define: Save location
                enzName = enzName.replace('/','_')
                if datasetTag is None:
                    figLabel = (f'AA Distribution-{enzName}-{sortType}-'
                                f'YMax_{yMax}-{codonType}-{len(rf.columns)}AA-'
                                f'MinCounts_{self.minSubCount}.png')
                else:
                    if codonType == datasetTag:
                        figLabel = (f'AA Distribution-{codonType}_Codon-'
                                    f'YMax_{yMax}.png')
                    else:
                        figLabel = (f'AA Distribution-{enzName}-{datasetTag}-{sortType}-'
                                    f'YMax_{yMax}-{codonType}-{len(rf.columns)}AA-'
                                    f'MinCounts_{self.minSubCount}.png')
                figLabel = figLabel.replace(' ', '_')
                saveLocation = os.path.join(self.pathSaveFigs, figLabel)

                # Save figure
                self.saveFig(fig=fig, path=saveLocation)

        # Plot the data
        if plotInitial and not skipInitial:
            plotFig(rf=rfInitial, sortType='Initial Sort')
        if plotFinal:
            if codonType == datasetTag:
                plotFig(rf=rfFinal, sortType=datasetTag)
            else:
                plotFig(rf=rfFinal, sortType='Final Sort')



    def plotPositionalProbDist(self, probability, entropyScores, sortType,
                               datasetTag):
        print('======================== Plot: Probability Distribution '
              '=========================')
        for position in entropyScores.index:
            yMax = 1

            # Extract values for plotting
            probabilities = list(probability.loc[:, position])

            # Calculate positional entropy
            shannonS = 0
            notNumber = False
            for AA in self.letters:
                prob = probability.loc[AA, position]
                shannonS += -prob * np.log2(prob)
            if math.isnan(shannonS):
                shannonS = 0
                notNumber = True

            # Plot the data
            setEdgeColor = True
            fig, ax = plt.subplots(figsize=self.figSizeMini)
            if setEdgeColor:
                widthBar = 0.8
                xPos = np.arange(len(probability.index))

                for index, AA in enumerate(probability.index):
                    edgeColor = 'dimgray' if AA in ['F', 'W', 'Y'] else 'black'
                    ax.bar(xPos[index], probability.loc[AA, position], widthBar, label=AA,
                           color=self.colorsAA[AA], edgecolor=edgeColor)
            else:
                plt.bar(self.letters, probabilities,
                        color=[self.colorsAA[AA] for AA in self.letters])
            # plt.xlabel('Amino Acids', fontsize=self.labelSizeAxis)
            plt.ylabel('Probability', fontsize=self.labelSizeAxis)
            if notNumber:
                plt.title(f'{self.enzymeName.replace(' - ', '\n')}\n'
                          f'Amino Acid Distribution at {position}\n'
                          f'ΔS = {entropyScores.loc[position, "ΔS"]:.3f}, '
                          f'Shannon Entropy = {shannonS:.0f}',
                          fontsize=self.labelSizeTitle, fontweight='bold')
            else:
                plt.title(f'{self.enzymeName.replace(' - ', '\n')}\n'
                          f'Amino Acid Distribution at {position}\n'
                          f'ΔS = {entropyScores.loc[position, "ΔS"]:.3f}, '
                          f'Shannon Entropy = {shannonS:.3f}',
                          fontsize=self.labelSizeTitle, fontweight='bold')
            plt.ylim(0, yMax)


            # Set tick parameters
            ax.tick_params(axis='both', which='major', length=self.tickLength,
                           labelsize=self.labelSizeTicks)

            # Set x-ticks
            xTicks = np.arange(0, len(self.letters), 1)
            ax.set_xticks(xTicks) # Set the positions of the ticks
            ax.set_xticklabels(self.letters, rotation=0, ha='center',
                               fontsize=self.labelSizeTicks)
            for tick in ax.xaxis.get_major_ticks():
                tick.tick1line.set_markeredgewidth(self.lineThickness) # Set tick width

            # Set y-ticks
            tickStepSize = 0.2
            yTicks = np.arange(0, yMax + tickStepSize, tickStepSize)
            yTickLabels = [f'{tick:.0f}' if tick == 0 or tick == 1 else f'{tick:.1f}'
                           for tick in yTicks]
            ax.set_yticks(yTicks) # Set the positions of the ticks
            ax.set_yticklabels(yTickLabels)
            for tick in ax.yaxis.get_major_ticks():
                tick.tick1line.set_markeredgewidth(self.lineThickness)

            # Set the edge thickness
            for spine in ax.spines.values():
                spine.set_linewidth(self.lineThickness)

            self.plotFig(plt=plt, fig=fig)

            # Save the figure
            if self.saveFigures:
                # Define: Save location
                if datasetTag is None:
                    figLabel = (f'AA_Distribution-{position}-{self.enzymeName}-'
                                f'{sortType}-Unfiltered-'
                                f'MinCounts_{self.minSubCount}.png')
                else:
                    figLabel = (f'AA_Distribution-{position}-{self.enzymeName}-'
                                f'{sortType}-{datasetTag}-'
                                f'MinCounts_{self.minSubCount}.png')
                saveLocation = os.path.join(self.pathSaveFigs, figLabel)

                # Save figure
                self.saveFig(fig=fig, path=saveLocation)



    def plotWordCloud(self, substrates, clusterNumPCA=None,
                      combinedMotifs=False, predActivity=False,
                      predModel=False, pca=False):
        print('=============================== Plot: Word Cloud '
              '================================')
        if clusterNumPCA is not None:
            print(f'Selecting PCA Population:{red} {clusterNumPCA}{resetColor}')
        else:
            print(f'Substrates: {purple}{self.datasetTag}{resetColor}')
        iteration = 0
        for substrate, count in substrates.items():
            print(f'     {blue}{substrate}{resetColor}, '
                  f'Count:{red} {round(count, 1):,}{resetColor}')
            iteration += 1
            if iteration == self.printNumber:
                break
        print('')

        # Limit the number of words
        if self.wordsLimit:
            print(f'Selecting: {red}{self.wordsTotal}{resetColor} words')
            subs = {}
            iteration = 0
            for substrate, count in substrates.items():
                subs[substrate] = count
                iteration += 1
                if iteration >= self.wordsTotal:
                    break
            substrates = subs
        totalWords = len(substrates)
        print(f'Plotting: {red}{totalWords:,}{resetColor} words\n\n')

        # Define: Figure title
        if predActivity:
            title = f'{self.enzymeName.replace(' - ', '\n')}\n{predModel}'
        elif self.releasedCounts:
            title = self.titleReleased
        elif combinedMotifs:
            title = self.titleCombined
        elif pca:
            title = f'{self.enzyme}\n{pca}'
        else:
            title = self.titleWords
            # title += f'\nTop {totalWords} Substrates'


        # Create word cloud
        cmap = self.createCustomColorMap(colorType='Word Cloud')
        wordcloud = (WordCloud(
            width=950,
            height=800,
            background_color='white',
            min_font_size=10, # Minimum font size
            max_font_size=120, # Maximum font size
            scale=5,  # Increase scale for larger words
            colormap=cmap
        ).generate_from_frequencies(substrates))


        # Display the word cloud
        fig = plt.figure(figsize=self.figSize, facecolor='white')
        plt.imshow(wordcloud, interpolation='bilinear')
        plt.title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        plt.axis('off')
        self.plotFig(plt=plt, fig=fig)


        # Save the Figure
        if self.saveFigures:
            if self.motifLen is None:
                seqLength = len(self.xAxisLabels)
            else:
                seqLength = self.motifLen
            figType = 'Words'
            if pca:
                pca = (pca.replace(f'{self.datasetTag}', '').replace('\n', '')
                       .replace(' ', '_').replace('Population', 'Pop'))
                figType = f'{pca}-{figType}'
            self.saveFigure(fig=fig, figType=figType, seqLen=seqLength,
                            combinedMotifs=combinedMotifs)



    def findSequence(self, substrates, sequence, sortType,
                     combinedMotifs=False):
        print('================================= Find Sequence '
              '=================================')
        if 'initial' in sortType.lower():
            print(f'Dataset: {purple}Unfiltered{resetColor}\n'
                  f' Enzyme: {purple}{self.enzymeName}{resetColor}\n'
                  f'   Sort: {purple}{sortType}{resetColor}')
        else:
            print(f'Dataset: {purple}{self.datasetTag}{resetColor}\n'
                  f' Enzyme: {purple}{self.enzymeName}{resetColor}\n'
                  f'   Sort: {purple}{sortType}{resetColor}')
        hitsTotal = {}
        for seq in sequence:
            hitsTotal[seq] = 0

        totalSubstrates = 0
        totalHits = 0
        hits = {}
        useStr = True
        if isinstance(sequence, str):
            for substrate, count in substrates.items():
                totalSubstrates += count
                if sequence in substrate:
                    totalHits += count
                    hits[substrate] = count
                    hitsTotal[sequence] += count
        else:
            useStr = False
            for substrate, count in substrates.items():
                totalSubstrates += count
                for seq in sequence:
                    if seq in substrate:
                        totalHits += count
                        hits[substrate] = count
                        hitsTotal[seq] += count
        if hits:
            print()
        if combinedMotifs:
            print(f'Unique Motifs: {red}{len(substrates.keys()):,}{resetColor}\n'
                  f' Total Motifs: {red}{totalSubstrates:,}{resetColor}\n')
        else:
            print(f'Unique Sequences: {red}{len(substrates.keys()):,}{resetColor}\n'
                  f' Total Sequences: {red}{totalSubstrates:,}{resetColor}\n')

        # Find matches
        if useStr:
            print(f'Finding Sequence: {purple}{sequence}{resetColor}')
        else:
            print(f'Finding Sequence: {purple}{", ".join(sequence)}{resetColor}')
        if len(hits.keys()) > 0:
            color = pink
            print(f'Hits:')
            if useStr:
                for index, (substrate, count) in enumerate(hits.items()):
                    substrate = substrate.replace(sequence, f'{blue}{sequence}{color}')
                    print(f'  {color}{substrate}{resetColor}, {red}{count:,}{resetColor}')
                    if index >= self.printNumber:
                        break
            else:
                for index, (substrate, count) in enumerate(hits.items()):
                    for seq in sequence:
                        substrate = substrate.replace(seq, f'{blue}{seq}{color}')
                    print(f'  {color}{substrate}{resetColor}, {red}{count:,}{resetColor}')
                    if index >= self.printNumber:
                        break
        else:
            print(f'Hits: {red}{totalHits}{resetColor}')
        print(f'Unique Sequences: {red}{len(hits.keys()):,}{resetColor}\n'
              f'   Total Matches: {red}{totalHits:,}{resetColor}\n')
        hitsPercent = (totalHits / totalSubstrates) * 100

        # Total occurrences
        print('Total occurrences:')
        hitsTotal = dict(sorted(
            hitsTotal.items(), key=lambda item: item[1], reverse=True)
        )
        for substrate, count in hitsTotal.items():
            print(f'  {blue}{substrate}{resetColor}, {red}{count:,}{resetColor}')
        print()
        print(f'Hit Population: {red}{totalHits:,}{resetColor} / {red}{totalSubstrates:,}'
              f'{resetColor} = {red}{round(hitsPercent,self.roundVal)} %'
              f'{resetColor}\n\n')


    def findAAInSequence(self, substrates, AA, idxPos, sortType,
                         combinedMotifs=False):
        print('================================= Find Sequence '
              '=================================')
        idxPos -= 1 # Adjust index
        if 'initial' in sortType.lower():
            print(f'Dataset: {purple}Unfiltered{resetColor}\n'
                  f' Enzyme: {purple}{self.enzymeName}{resetColor}\n'
                  f'   Sort: {purple}{sortType}{resetColor}')
        else:
            print(f'Dataset: {purple}{self.datasetTag}{resetColor}\n'
                  f' Enzyme: {purple}{self.enzymeName}{resetColor}\n'
                  f'   Sort: {purple}{sortType}{resetColor}')
        totalSubstrates = 0
        totalHits = 0
        hits = {}
        useStr = True
        if isinstance(AA, str):
            for substrate, count in substrates.items():
                totalSubstrates += count
                if substrate[idxPos] == AA:
                    totalHits += count
                    hits[substrate] = count
        else:
            useStr = False
            for substrate, count in substrates.items():
                totalSubstrates += count
                #print(f'Sub: {substrate}\n  I: {idxPos}, AA: {substrate[idxPos]}')
                if substrate[idxPos] in AA:
                    totalHits += count
                    hits[substrate] = count
        if combinedMotifs:
            print(f'Unique Motifs: {red}{len(substrates.keys()):,}{resetColor}\n'
                  f' Total Motifs: {red}{totalSubstrates:,}{resetColor}\n')
        else:
            print(f'Unique Sequences: {red}{len(substrates.keys()):,}{resetColor}\n'
                  f' Total Sequences: {red}{totalSubstrates:,}{resetColor}\n')

        # Find matches
        print(f'Finding: {purple}{AA}@R{idxPos + 1}{resetColor}')
        if len(hits.keys()) > 0:
            color = pink
            print(f'Hits:')
            if useStr:
                for index, (substrate, count) in enumerate(hits.items()):
                    substrate = substrate.replace(AA, f'{blue}{AA}{color}')
                    print(f'  {color}{substrate}{resetColor}, {red}{count:,}{resetColor}')
                    if index >= self.printNumber:
                        break
            else:
                for index, (substrate, count) in enumerate(hits.items()):
                    substrate = substrate.replace(
                        substrate[idxPos], f'{blue}{substrate[idxPos]}{color}')
                    print(f'  {color}{substrate}{resetColor}, {red}{count:,}{resetColor}')
                    if index >= self.printNumber:
                        break
        else:
            print(f'Hits: {red}{totalHits}{resetColor}')
        print(f'Unique Sequences: {red}{len(hits.keys()):,}{resetColor}\n'
              f'   Total Matches: {red}{totalHits:,}{resetColor}\n')
        hitsPercent = (totalHits / totalSubstrates) * 100
        print(f'Hit Population: {red}{totalHits:,}{resetColor} / {red}{totalSubstrates:,}'
              f'{resetColor} = {red}{round(hitsPercent,self.roundVal)} %'
              f'{resetColor}\n\n')



    def generateSubstrates(self, df, eMap, minES, dataType,
                           subsReq={}, filter={}):
        print('============================== Generate Substrates '
              '==============================')
        print(f'Dataset: {purple}{self.datasetTag}{resetColor}\n'
              f'DataFrame Type: {purple}{dataType}{resetColor}\n'
              f'{df}\n')

        # Select favorable AAs
        print(f'Preferred Residues:')
        preferredAAs = []
        addedAAsTag = ''
        for index, pos in enumerate(df.columns):
            if pos in filter.keys():
                nonzeroAAs = filter[pos]
            else:
                nonzeroAAs = df[df[pos] != 0].index.tolist()

                # Filter out AAs with low ES
                ES = eMap.loc[:, pos]
                nonzeroAAs = [aa for aa in nonzeroAAs if ES[aa] >= minES]

            # Add AAs for inSubsPred
            if subsReq != {}:
                addedAAs = []
                for key, substrates in subsReq.items():
                    for substrate in substrates:
                        AA = substrate[index]
                        if AA not in addedAAs and AA not in nonzeroAAs:
                            addedAAs.append(AA)
                            nonzeroAAs.append(AA)

                    # Record added AAs
                    if addedAAs:
                        if len(addedAAs) == 1:
                            label = f'{addedAAs[0]}@{pos}'
                        else:
                            label = f'[{",".join(addedAAs)}]@{pos}'
                        if index == len(df.columns) - 1:
                            addedAAsTag += f'{label}'
                        else:
                            addedAAsTag += f'{label}_'
            print(f'     Position {greenLight}{pos}{resetColor}: '
                  f'{pink}{", ".join(nonzeroAAs)}{resetColor}')
            preferredAAs.append(nonzeroAAs)
        print()

        if addedAAsTag != '':
            print(f'Added AAs: {pink}{addedAAsTag}{resetColor}')
        print(f'Minimum ES: {red}{minES}{resetColor}\n')
        # Generate all possible substrate combinations
        allDualModelss = list(product(*preferredAAs))
        # print(f'\n\nDualModelss:\n{allDualModelss}\n\n')

        # Convert tuples to strings (AA sequences)
        genSubstrates = [''.join(DualModels) for DualModels in allDualModelss]
        NSubs = len(genSubstrates)
        print(f'Generated substrates: {red}N = {NSubs:,}{resetColor}')
        for iteration in range(0, self.printNumber):
            index = random.randint(0, NSubs)
            print(f'     {pink}{genSubstrates[index]}{resetColor}')
        print('\n')

        return genSubstrates, addedAAsTag
    
    
    def normalizeValues(self, substrates, datasetTag):
        print(f'=============================== Normalize Values '
              f'===============================')
        print(f'Dataset: {purple}{datasetTag}{resetColor}')
        self.maxValue = max(substrates.values())
        print(f'Max Value: {red}{self.maxValue:,}{resetColor}\n')

        # Inspect datatype
        useIntegers = False
        for value in substrates.values():
            if isinstance(value, int):
                useIntegers = True
                break

        # Normalize the values
        substratesNormValues = {}
        for substrate, value in substrates.items():
            substratesNormValues[substrate] = (value / self.maxValue)
        print(f'Normalized Values: {purple}Top {self.printNumber} Substrates{resetColor}')
        if useIntegers:
            for iteration, (substrate, value) in enumerate(substratesNormValues.items()):
                print(f'     {pink}{substrate}, {red}{round(value, self.roundVal):,}'
                      f'{resetColor}')
                if iteration >= self.printNumber:
                    break
        else:
            for iteration, (substrate, value) in enumerate(substratesNormValues.items()):
                print(f'     {pink}{substrate}{resetColor}, {red}{round(value, 3):,}'
                      f'{resetColor}')
                if iteration >= self.printNumber:
                    break
        print('\n')

        return substratesNormValues



    def splitDataset(self, substrates, quantileSplit):
        print('================================= Split Dataset '
              '=================================')
        print(f'Dataset: {purple}{self.datasetTag}{resetColor}\n')

        # Convert data to a df
        pd.set_option('display.max_rows', 10)
        dfSubstrates = pd.DataFrame(substrates.values(),
                                            index=substrates.keys(),
                                            columns=['Activity'])

        # Split dataset
        cutoff = (100 - quantileSplit) / 100
        threshold = dfSubstrates['Activity'].quantile(cutoff)
        dfHigh = dfSubstrates[dfSubstrates['Activity'] > threshold]
        dfLow = dfSubstrates[dfSubstrates['Activity'] <= threshold]
        print(f'Substrates: {pink}Top {quantileSplit}% Activity Scores\n'
              f'{resetColor}{dfHigh}\n\n\n'
              f'Substrates: {pink}Remaining {100 - quantileSplit}% Activity Scores\n'
              f'{resetColor}{dfLow}\n\n')


        # Divide low substrate scores
        substratesFiltered = dfHigh
        numSections = 5
        maxScoreLow = max(dfLow.loc[:, 'Activity'])
        minScoreLow = min(dfLow.loc[:, 'Activity'])
        section = maxScoreLow / numSections
        activityRange = np.arange(minScoreLow, maxScoreLow + section, section)[::-1]
        activityRange[0] = maxScoreLow

        # Collect substrates
        print(f'Collecting substrates from the low activity dataset:')
        for index, valueA in enumerate(activityRange):
            if index == numSections - 1:
                break
            index = index + 1
            valueB = activityRange[index]

            # Select substrates
            selected = dfLow[(dfLow['Activity'] <= valueA)
                             & (dfLow['Activity'] > valueB)]
            selectedSubsN = len(selected.index)
            print(f'Low Activity Subset: {red}{index}{resetColor}\n'
                  f'   Total substrates: {red}{selectedSubsN}{resetColor}\n'
                  f'     Activity Range: {red}{round(valueA, 5)}{resetColor} - '
                  f'{red}{round(valueB, 5)}{resetColor}')

            if index > 2:
                # Select subset of the subset
                selectNthSubstrate = 10
                print(f'     Selecting every {red}{selectNthSubstrate}th{resetColor} '
                      f'substrate')
                selected = selected.iloc[::selectNthSubstrate]
            substratesFiltered = pd.concat([substratesFiltered, selected])
            print(f'{selected}\n\n')
        print(f'Filtered Substrates:\n{greenLight}{substratesFiltered}{resetColor}\n\n')

        collectedSubs = {}
        for substrate in substratesFiltered.index:
            collectedSubs[substrate] = substratesFiltered.loc[substrate, 'Activity']

        return collectedSubs



    def normalizeProbRatios(self, finalRF, initialRF,
                            pData=True, pHeader=True):
        if pData:
            if pHeader:
                print('========================= Normalize Probability Ratios '
                      '==========================')
            else:
                print(f'Normalize Probability Ratios: '
                      f'{purple}{self.datasetTag}{resetColor}')
            print(f'Relative Frequency: {purple}Final{resetColor}\n{finalRF}\n')
            print(f'Relative Frequency: {purple}Initial{resetColor}\n{initialRF}\n\n')

        # Calculate: RF ratios
        ratio = pd.DataFrame(0.0, index=finalRF.index, columns=finalRF.columns)
        for pos in ratio.columns:
            if len(initialRF.columns) == 1:
                ratio.loc[:, pos] = finalRF.loc[:, pos] / initialRF.iloc[:, 0]
            else:
                ratio.loc[:, pos] = finalRF.loc[:, pos] / initialRF.loc[:, pos]
        ratio = ratio.replace([-np.inf, np.inf], 0.0) # Remove inf values
        if pData:
            print(f'RF Ratios:\n{ratio}\n\n')


        # Calculate: Ratio probability
        norm = pd.DataFrame(0.0, index=finalRF.index, columns=finalRF.columns)
        for pos in norm.columns:
            for AA in norm.index:
                norm.loc[AA, pos] = ratio.loc[AA, pos] / sum(ratio.loc[:, pos])
        if pData:
            print(f'Normalized RF Ratios:\n{norm}\n\n')

        return norm



    def predictActivity(self, activityExp, finalRF, initialRF, predModel,
                        predLabel, errorBars=False, combinedMotifs=False,
                        barWidth=0.35, colorExp='#BF5700',
                        colorPred='#F8971F', rotateLabel=45):
        print('============================ Predict Substrate Activity '
              '=============================')
        N = len(activityExp.keys())
        print(f'Dataset: {purple}{predModel}{resetColor}\n'
              f'Evaluating: {purple}{predLabel}{resetColor}\n'
              f'N Substrate Sequences: {red}{N:,}{resetColor}\n')

        # Get prediction matrix
        matrix = self.normalizeProbRatios(finalRF=finalRF,
                                          initialRF=initialRF,
                                          pHeader=False)

        entropy = pd.DataFrame(0.0, index=matrix.columns, columns=['ΔS'])
        entropyMax = np.log2(len(matrix.index))
        for indexColumn in matrix.columns:
            S = 0
            for indexRow, probRatio in matrix.iterrows():
                prob = probRatio[indexColumn]
                if prob == 0:
                    continue
                else:
                    S += -prob * np.log2(prob)
            entropy.loc[indexColumn, 'ΔS'] = entropyMax - S
        maxS = max(entropy.iloc[:, 0])
        for pos in entropy.index:
            entropy.loc[pos, 'ΔS'] = entropy.loc[pos, 'ΔS'] / maxS
        print(f'{entropy}\n\n')

        # for pos in matrix.columns:
        #     matrix.loc[:, pos] = matrix.loc[:, pos] * entropy.loc[pos, 'ΔS']
        # print(f'Normalized RF Ratios:\n{matrix}\n\n')
        # self.calculateEntropy(rf=matrix, plotFig=True)


        def plotPredActivity(values, errorBars, tag):
            colorP, colorE = red, blue

            # Calculate: Activity
            activityPred = {}
            subLen = len(next(iter(activityExp)))
            neutral = 1 / len(values.index)
            spacerMax = max([len(pos) for pos in values.columns])
            for substrate in activityExp.keys():
                print(f'{pink}{substrate}{resetColor}')
                score = 0
                for index in range(subLen): ##
                    # Evaluate substrate
                    AA = substrate[index]
                    pos = values.columns[index]
                    value = values.loc[AA, pos]
                    S = entropy.loc[pos, 'ΔS']
                    x = S ** value
                    # y = value # R2: 0.622, 0.980
                    # y = value / (S ** value) # R2: 0.671, 0.978
                    y = value / x
                    spacer = (spacerMax - len(pos)) * ' '
                    print(f'* {blue}{AA}@{pos}{resetColor}:{spacer} '
                          f'ER: {purple}{value:.4f}{resetColor}, '
                          f'∆S: {S:.4f}, '
                          f'Adj: {greenLight}{x:.4f}{resetColor}, '
                          f'ER/Adj: {red}{y:.4f}{resetColor}')
                    # print(f'* {blue}{AA}@{pos}{resetColor} ({round(S, 4)}): '
                    #       f'{round(value, 4)} -> {round(x, 4)}')
                    # print(f'* {blue}{AA}@{pos}{resetColor}: {round(value,4)}')
                    if score == 0:
                        score = y
                    else:
                        score *= y
                print(f'Score: {colorP}{score:.3e}{resetColor}\n')
                activityPred[substrate] = float(score)
            ranked = pd.Series(activityPred.values()).rank(
                ascending=False, method='min').astype(int)


            dec = self.roundVal - 1
            maxActivity = max(activityPred.values())
            print(f'Matrix Type: {purple}{tag}{resetColor}\n'
                  f'Predicted Activity: '
                  f'(Max Score: {colorP}{maxActivity:.{dec}e}{resetColor})')
            for index, (substrate, activity) in enumerate(activityPred.items()):
                print(f'    {pink}{substrate}{resetColor} ({ranked[index]}), '
                      f'Score: {colorP}{activity:.{dec}e}{resetColor}')
                if index == self.printNumber - 1:
                    break
            print()

            def ZScores(data, tag=''):
                # Calculate: Z-scores
                z = {}
                mu = np.average(list(data.values()))
                sigma = np.std(list(data.values()))
                # print(f'Z-Scores: {tag}\n'
                #       f'* µ: {red}{mu:.3e}{resetColor}\n'
                #       f'* σ: {red}{sigma:.3e}{resetColor}\n')
                for seq, count in data.items():
                    z[seq] = float((count - mu) / sigma)
                return z

            # Normalize values
            activityExpNorm, maxExpActivity = {}, max(activityExp.values())
            activityPredNorm = {}
            for index, (substrate, activity) in enumerate(activityPred.items()):
                activityPredNorm[substrate] = activity / maxActivity
                activityExpNorm[substrate] = activityExp[substrate] / maxExpActivity
            activityPredZ = ZScores(activityPred, 'Predicted Activity')
            activityExpZ = ZScores(activityExp, 'Experimental Activity')
            # print(f'Predicted Normalized Activity:')
            scores = []
            for index, (substrate, activity) in enumerate(activityPred.items()):
                s = f'{activity:,.{self.roundVal}f}'
                l = np.log(activity)
                l = f'{l:,.{self.roundVal}f}'
                scores.append(float(s))

            # Summarize datasets
            expActivity = [round(a, self.roundVal) for a in activityExp.values()]
            expScores = [round(a, self.roundVal) for a in activityExpNorm.values()]
            expScoresZ = [round(a, self.roundVal) for a in activityExpZ.values()]
            predScores = [round(a, self.roundVal) for a in activityPredNorm.values()]
            predScoresZ = [round(a, self.roundVal) for a in activityPredZ.values()]

            expRank = self.rankScores(expScores)
            predRank = self.rankScores(predScores)

            # Compare predictions
            print(f'Normalized Activity:')
            data = pd.DataFrame(0.0, index=list(activityExp.keys()),
                                columns=['Experimental', 'Predicted'])
            data.loc[:, 'Experimental'] = expScores
            data.loc[:, 'Predicted'] = predScores
            print(data)

            # Spearman rank correlation
            from scipy.stats import pearsonr, spearmanr
            rho, pValue = spearmanr(expRank, predRank)
            rho, pValue = round(rho, 3), round(pValue, 3)
            print(f'\nSpearman ρ: {rho}, p={pValue}\n')

            print(f'Predicted Activity:\n'
                  f'* Normalized: {colorP}{predScores}{resetColor}\n'
                  f'* Z-Scores:   {colorP}{predScoresZ}{resetColor}\n')
            print(f'Experimental Activity:\n'
                  f'* Experimental: {colorE}{expActivity}{resetColor}\n'
                  f'* Normalized:   {colorE}{expScores}{resetColor}\n'
                  f'* Z-Scores:     {colorE}{expScoresZ}{resetColor}\n')

            # Plot prediction matrix
            self.plotMatrix(data=matrix, figLabel='Prediction Matrix', printData=False)

            # Set title
            enzName = self.enzymeName.replace(' - ', '\n')
            title = f'{enzName}\nSubstrate Activity'
            if predLabel:
                title += f'\n{predLabel}'

            # Plot bar graph
            labels = list(activityPred.keys())
            predVals = [activityPred[k] for k in labels]
            expVals = [activityExpNorm[k] for k in labels]
            xTicks = np.arange(len(labels))

            # Make figure
            fig, ax = plt.subplots(figsize=self.figSize)
            ax.bar(xTicks - barWidth / 2, predScoresZ, barWidth, label='Predicted',
                   color=colorPred, edgecolor='black', linewidth=self.lineThickness)
            ax.bar(xTicks + barWidth / 2, expScoresZ, barWidth, label='Experimental',
                   color=colorExp, edgecolor='black', linewidth=self.lineThickness)
            plt.title(f'{title}\nSpearman ρ: {rho}',
                      fontsize=self.labelSizeTitle, fontweight='bold')
            plt.axhline(y=0, color='black', linewidth=self.lineThickness)
            ax.set_ylabel('Z-Score', fontsize=self.labelSizeAxis)
            ax.legend(edgecolor='black', linewidth=self.lineThickness, loc='best',
                      prop=FontProperties(weight='bold', size=self.labelSizeTicks - 2))

            # Y-axis
            spacer = 0.2
            yMax = round(np.ceil(
                max(max(predScoresZ), max(expScoresZ)) * 10), 2) / 10
            yMin = round(np.floor(
                min(min(predScoresZ), min(expScoresZ)) * 10), 2) / 10
            yMax = np.ceil((round(yMax*10,2)/10)+spacer)
            yMin = np.floor((round(yMin*10,2)/10)-spacer)
            plt.ylim(yMin, yMax)

            # Set xticks
            ax.set_xticks(xTicks)
            ax.set_xticklabels(labels, rotation=rotateLabel)

            # Set tick parameters
            ax.tick_params(axis='both', which='major', length=self.tickLength,
                           labelsize=self.labelSizeTicks, width=self.lineThickness)

            # Set the edge thickness
            for tick in ax.xaxis.get_major_ticks():
                tick.tick1line.set_markeredgewidth(
                    self.lineThickness)  # Set tick width
            for tick in ax.yaxis.get_major_ticks():
                tick.tick1line.set_markeredgewidth(
                    self.lineThickness)  # Set tick width
            for spine in ax.spines.values():
                spine.set_linewidth(self.lineThickness)

            self.plotFig(plt=plt, fig=fig)

            # Save the Figure
            if self.saveFigures:
                # Define: Save location
                figTag = f'PredActivity-BarGraph-{predModel}'
                if predLabel:
                    figTag = figTag.replace(
                        'BarGraph',
                        f'BarGraph-{predLabel.replace(' ', '_')}'
                    )

                # Save figure
                self.saveFigure(fig=fig, figType=figTag, seqLen=subLen, N=N,
                                combinedMotifs=combinedMotifs)

            # ============================================================================
            # ============================================================================
            from scipy.optimize import curve_fit

            def fnExp(x, a, b, c):
                return a * np.exp(b * x) + c

            def fitData(x, y):
                x, y = np.array(x, dtype=float), np.array(y, dtype=float)

                # a and c are linear given b, so profile b on a grid to get a real p0
                def sse(b):
                    A = np.column_stack([np.exp(b * x), np.ones_like(x)])
                    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
                    return np.sum((y - A @ coef) ** 2)

                grid = np.linspace(0.05, 5.0, 200)
                b0 = grid[int(np.argmin([sse(b) for b in grid]))]
                A = np.column_stack([np.exp(b0 * x), np.ones_like(x)])
                (a0, c0), *_ = np.linalg.lstsq(A, y, rcond=None)

                # bound b away from 0 so the linear ramp is unreachable
                popt, pcov = curve_fit(fnExp, x, y, p0=[a0, b0, c0],
                                       bounds=([-np.inf, 1e-3, -np.inf],
                                               [np.inf, 10.0, np.inf]),
                                       maxfev=10000)

                if np.any(np.sqrt(np.diag(pcov)) > 100 * np.abs(popt)):
                    print(f'WARNING: unidentified fit, popt={popt}')

                yPred = fnExp(x, *popt)
                ss_res = np.sum((y - yPred) ** 2)
                ss_tot = np.sum((y - np.mean(y)) ** 2)
                r2 = 1 - (ss_res / ss_tot)

                xFit = np.linspace(x.min(), x.max(), 300)
                return xFit, fnExp(xFit, *popt), r2


            def fnLog(x, a, b):
                return a * np.log(x) - b

            def fitDataLog(x, y):
                popt, pcov = curve_fit(fnLog, x, y, maxfev=10000)
                a, b = popt

                # Generate smooth curve for plotting
                xFit = np.linspace(min(x), max(x), 300)
                yFit = fnLog(xFit, *popt)

                # R²
                yPred = fnLog(np.array(x), *popt)
                ss_res = np.sum((np.array(y) - yPred) ** 2)
                ss_tot = np.sum((np.array(y) - np.mean(y)) ** 2)
                r2 = 1 - (ss_res / ss_tot)
                return xFit, yFit, r2


            def fnPoly(x, a, b, c):
                return a * x ** 2 + b * x + c

            def fitDataPoly(x, y):
                popt, pcov = curve_fit(fnPoly, x, y, maxfev=10000)
                xFit = np.linspace(min(x), max(x), 300)
                yFit = fnPoly(xFit, *popt)
                yPred = fnPoly(np.array(x), *popt)
                ss_res = np.sum((np.array(y) - yPred) ** 2)
                ss_tot = np.sum((np.array(y) - np.mean(y)) ** 2)
                r2 = 1 - (ss_res / ss_tot)
                return xFit, yFit, r2
            # ============================================================================
            # ============================================================================
            # Fit curve: exponential
            x_fit, y_fit, r2 = fitData(x=np.array(expScoresZ), y=np.array(predScoresZ))
            print(f'R2 exp: {r2:.3f}')
            # x_fit, y_fit, r2 = fitDataPoly(x=np.array(expScoresZ), y=np.array(y))
            # print(f'R2 poly: {r2:.3f}')

            # Fit curve: linear
            coeffs = np.polyfit(expScoresZ, predScoresZ, 1)
            yPred = np.polyval(coeffs, expScoresZ)
            ss_res = np.sum((np.array(predScoresZ) - yPred) ** 2)
            ss_tot = np.sum((np.array(predScoresZ) - np.mean(predScoresZ)) ** 2)
            r2_linear = 1 - (ss_res / ss_tot)
            print(f'Linear R²: {r2_linear:.3f}')


            # Make figure
            fig, ax = plt.subplots(figsize=self.figSize)
            if errorBars:
                plt.errorbar(expScoresZ, predScoresZ, yerr=errorBars, fmt="o",
                             color=colorExp, ecolor='black')
            else:
                plt.scatter(expScoresZ, predScoresZ, color=colorExp, edgecolor='black')
            ax.plot(x_fit, y_fit, color='black', linestyle='-',
                    linewidth=self.lineThickness)
            plt.xlabel('Experimental Activity', fontsize=self.labelSizeAxis)
            plt.ylabel('Predicted Activity', fontsize=self.labelSizeAxis)
            plt.title(title, fontsize=self.labelSizeTitle, fontweight='bold')
            plt.axhline(y=0, color='black', linewidth=self.lineThickness)
            # plt.grid(True, linestyle='-', color='black')

            # Axis
            spacer = 0.2
            xMax = self.roundup(max(expScoresZ) + spacer)
            xMin = self.roundup(min(expScoresZ) - spacer, upperLim=False)
            yMax = self.roundup(max(predScoresZ) + spacer)
            yMin = self.roundup(min(predScoresZ) - spacer, upperLim=False)
            plt.xlim(xMin, xMax)
            plt.ylim(yMin, yMax)

            # Set tick parameters
            ax.tick_params(axis='both', which='major', length=self.tickLength,
                           labelsize=self.labelSizeTicks, width=self.lineThickness)

            # Set the edge thickness
            for tick in ax.xaxis.get_major_ticks():
                tick.tick1line.set_markeredgewidth(self.lineThickness) # Set tick width
            for tick in ax.yaxis.get_major_ticks():
                tick.tick1line.set_markeredgewidth(self.lineThickness) # Set tick width
            for spine in ax.spines.values():
                spine.set_linewidth(self.lineThickness)

            # Legend
            invisibleHandle = Line2D([], [], linestyle='None', marker='None',
                                      color='none')
            ax.legend(
                handles=[invisibleHandle],
                labels=[f'R² = {r2:.3f}\nρ: {rho}'],
                prop=FontProperties(size=self.labelSizeTicks - 2, weight='bold'),
                handlelength=0, handletextpad=0, edgecolor='black',
                linewidth=self.lineThickness, loc='upper left', framealpha=0.9
            )

            self.plotFig(plt=plt, fig=fig)

            # Save the Figure
            if self.saveFigures:
                # Define: Save location
                figTag = f'PredActivity-ScatterPlot-{predModel}'
                if predLabel:
                    figTag = figTag.replace(
                        'ScatterPlot',
                        f'ScatterPlot-{predLabel.replace(' ', '_')}'
                    )
                self.saveFigure(fig=fig, figType=figTag, seqLen=subLen, N=N,
                                combinedMotifs=combinedMotifs)
        plotPredActivity(values=matrix, errorBars=errorBars, tag='Probability Ratios')
