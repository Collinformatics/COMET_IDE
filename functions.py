import base64
from Bio import SeqIO
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor
from counter import counter
import gzip
import io
from itertools import batched
import json
import logomaker
import math
import matplotlib
from matplotlib.colors import Colormap, LinearSegmentedColormap, Normalize
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os.path
import pandas as pd
import queue
import seaborn as sns
import secrets
import sys
import threading
import time
from wordcloud import WordCloud


defaultResidues = (
    ('Alanine','Ala','A'), ('Arginine','Arg','R'), ('Asparagine','Asn','N'),
    ('Aspartic Acid','Asp','D'), ('Cysteine','Cys','C'), ('Glutamic Acid','Glu','E'),
    ('Glutamine','Gln','Q'), ('Glycine','Gly','G'), ('Histidine','His ','H'),
    ('Isoleucine','Ile','I'), ('Leucine','Leu','L'), ('Lysine','Lys','K'),
    ('Methionine','Met','M'), ('Phenylalanine','Phe','F'), ('Proline','Pro','P'),
    ('Serine','Ser','S'), ('Threonine','Thr','T'), ('Tryptophan','Typ','W'),
    ('Tyrosine','Tyr','Y'), ('Valine','Val','V')
)


# Generate figures entirely in memory without opening a window
matplotlib.use('Agg')  # Use a non-interactive backend for servers


class WebApp:
    def __init__(self):
        # Params: Job
        self.jobDone = False
        self.setS = False
        self.jobParams = {}
        self.datasetTag = ''
        self.datasetTagMotif = ''
        self.combineProfiles = False
        self.saveData = True

        # Params: Dataset
        self.useCustomAxis = False
        self.customAxis = ['P5','P4','P3','P2','P1','P1\'','P2\'','P3\'']
        self.enzymeName = ''
        self.seqLength = 0
        self.maxCounts = False
        self.minCounts = 1
        self.printN = 10
        self.roundVal = 3
        self.xAxisLabel = []
        self.xAxisSet = []
        self.dropPos = []
        self.subsExp = {}
        self.subsExpAll = {}
        self.subsPred = {}
        self.countsExp = pd.DataFrame()
        self.countsExpSets = []
        self.countExpTotal = 0
        self.countExpUnique = 0
        self.rfExp = None
        self.rfExpScaled = None
        self.subsBg = {}
        self.countsBg = pd.DataFrame()
        self.countBgTotal = 0
        self.countBgUnique = 0
        self.rfBg = None
        self.entropy = pd.DataFrame(0.0, index=[], columns=['∆S'])
        self.entropyMax = None
        self.eMap = None
        self.eMapScaled = None
        self.eMapReleased = None
        self.eMapReleasedScaled = None
        self.datasetTypes = {'Exp': 'Experimental', 'Bg': 'Background'}
        

        # Params: COMET
        self.subProfile = False
        self.motifFilter = False
        self.blockScript = False
        self.initCOMET = False
        self.initFixedPos = []
        self.iteration = 0
        self.minS = 0.6
        self.minES = 0
        self.minESRel = -1 # -0.5
        self.motifPos = {}
        self.substrateProfile = pd.DataFrame()
        self.substrateProfileScl = pd.DataFrame()

        # Params: Combined Profiles
        self.avgBgRF = False
        self.fileOrder = []
        self.idxMotif = {}
        self.profiles = []

        # Params: Files
        self.queueLog = queue.Queue()
        self.fileExp = []
        self.fileExpRev = []
        self.fileExpCounts = []
        self.fileBg = []
        self.fileBgRev = []
        self.pathDir = ''
        self.pathData = ''
        self.pathFigs = ''
        self.pathLog = ''
        self.errorLog = 'error'

        # Params: Process dna
        self.seq5Prime = ''
        self.seq3Prime = ''
        self.minPhred = False

        # Params: Filter Dataset
        self.filterPos = False
        self.fixAA = {}
        self.exclAA = {}

        # Params: Figures
        self.labelSizeTitle = 20 # Set fontsize
        self.labelSizeAxis = 18 # Set fontsize
        self.labelSizeTicks = 16 # Set fontsize
        self.labelSizeEM = 11  # Set fontsize
        self.lineThickness = 1.5
        self.tickLength = 4
        self.figTag = ''
        self.numSamples = 50
        self.figureResolution = 600
        self.titleCombined = ''
        self.titleReleased = ''
        self.titleWeblogo = ''
        self.titleWeblogoCombined = ''
        self.titleWeblogoReleased = ''
        self.titleWords = ''
        self.titleWordsCombined = ''
        self.figEMSquares = False
        self.figSize = (9.5, 8) # (width, height)
        self.figSizeSq = (5, 8)
        self.figSizeMini = (self.figSize[0], 6)
        self.residueLabelType = 2  # 0 = full AA name, 1 = 3-letter code, 2 = 1 letter
        self.colorsAA = self.residueColors()
        self.residues = defaultResidues
        self.AA = [residue[2] for residue in self.residues]
        self.bigAAonTop = False
        self.figures = {}

        # Colors
        self.orange = '#FA8128'
        self.orangeBurnt = '#BF5700'

        # Print params
        pd.options.display.float_format = '{:,.3f}'.format
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        # pd.option_context('display.width', 1000, 'display.max_columns', None)


    @staticmethod
    def createCustomColorMap(colorType):
        colorType = colorType.lower()
        if colorType == 'green':
            useGreen = True
            if useGreen:
                # Green
                colors = ['#FFFFFF','#ABFF9B','#39FF14','#2E9418','#2E9418','#005000']
            else:
                # Orange
                colors = ['white','white','#FF76FA','#FF50F9',
                          '#FF00F2','#CA00DF','#BD16FF']
        elif colorType == 'stdev':
            colors = ['white','white','#FF76FA','#FF50F9','#FF00F2','#CA00DF','#BD16FF']
        elif colorType == 'wordcloud':
            # ,'#F2A900','#2E8B57','black'
            colors = ['#CC5500','#CC5500','#F79620','#FAA338',
                      '#00C01E','#1D680D','#003000','black']
        elif colorType == 'em':
            colors = ['navy','royalblue','dodgerblue','lightskyblue','white',
                      'white','lightcoral','red','firebrick','darkred']
        else:
            self.logError(f'ERROR: Cannot create colormap. '
                          f'Unrecognized colorType parameter: {colorType}')
            sys.exit(1)

        # Create colormap
        if len(colors) == 1:
            colorList = [(0, colors[0]), (1, colors[0])]
        else:
            colorList = [(i / (len(colors) - 1), color) for i, color in enumerate(colors)]
        return LinearSegmentedColormap.from_list('custom_colormap', colorList)


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
    def getKey(app):  # required for CSRF
        app.config['SECRET_KEY'] = secrets.token_hex(nbytes=32)


    @staticmethod
    def pressButton(message):
        print(f'Received data: {message}')
        return {'key': 'Returned data'}


    @staticmethod
    def normValues(data):
        norm = {}
        if data.values():
            maxScore = max(data.values())
            for substrate in data.keys():
                norm[substrate] = data[substrate] / maxScore
        return norm


    @staticmethod
    def ZScores(data):
        # Calculate: Z-scores
        z = {}
        if data.values():
            mu = np.average(list(data.values()))
            sigma = np.std(list(data.values()))
            for seq, count in data.items():
                z[seq] = (count - mu) / sigma
        return z



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
        return val


    def encodeFig(self, fig):
        # Save figure to a memory buffer instead of disk
        buffer = io.BytesIO()
        plt.savefig(
            buffer, format='png', dpi=self.figureResolution, bbox_inches='tight'
        )
        buffer.seek(0)

        # Encode as base64 for embedding in HTML
        figBase64 = base64.b64encode(buffer.getvalue()).decode('utf-8')

        plt.close(fig) # or fig.clear() to release memory
        buffer.close()

        return figBase64


    def getDatasetTag(self, log=True):
        self.fixAA = dict(sorted(self.fixAA.items()))
        self.exclAA = dict(sorted(self.exclAA.items()))
        tagFix = 'Fix '
        tagExcl = 'Exclude '

        # Evaluate filters
        if self.exclAA:
            for index, (pos, AA) in enumerate(self.exclAA.items()):
                if len(AA) > 1:
                    tagExcl += f'[{','.join(AA)}]@{pos} '
                else:
                    tagExcl += f'{AA}@{pos} '
            tagExcl = tagExcl[:-1]
        if self.fixAA:
            for index, (pos, AA) in enumerate(self.fixAA.items()):
                if len(AA) > 1:
                    tagFix += f'[{','.join(AA)}]@{pos} '
                else:
                    tagFix += f'{AA}@{pos} '
            tagFix = tagFix[:-1]
        if tagExcl != 'Exclude ' and tagFix != 'Fix ':
            self.datasetTag = f'{tagExcl} {tagFix}'
        elif tagFix != 'Fix ':
            self.datasetTag = tagFix
        elif tagExcl != 'Exclude ':
            self.datasetTag = tagExcl
        elif self.combineProfiles:
            self.datasetTag = 'Combined'
        else:
            self.datasetTag = 'Unfiltered'
        if self.motifFilter and not self.datasetTagMotif:
            self.datasetTagMotif = f'Register {self.datasetTag}'
            if tagExcl == 'Exclude ':
                self.datasetTagMotif = self.datasetTagMotif.replace(
                    f'Fix ', ''
                )
            else:
                self.datasetTagMotif = self.datasetTagMotif.replace(
                    f'Register ', ''
                )
            self.jobParams['Filter'] = self.datasetTagMotif
            if log:
                self.log(f'Filter: {self.datasetTagMotif}')
        elif 'Dataset Tag' not in self.jobParams.keys():
            self.jobParams['Filter'] = self.datasetTag
            if log:
                self.log(f'Filter: {self.datasetTag}')


    def getFilter(self, data):
        self.fixAA = {}
        self.exclAA = {}
        if 'filterPos' in data.keys():
            self.filterPos = data['filterPos']
            for key, value in data.items(): # Get filter params
                if 'fix' in key:
                    self.fixAA[key.replace('fix', '')] = value
                if 'excl' in key:
                    self.exclAA[key.replace('excl', '')] = value
        if 'predSubs' not in data.keys():
            self.getDatasetTag()


    def getSaveTag(self):
        # Evaluate filters
        if self.initCOMET:
            self.initCOMET = False
            tag = f'{self.datasetTagMotif.replace(' ', '_')}'
        elif self.motifFilter:
            tag = f'SubProfile-{self.datasetTagMotif.replace(' ', '_')}'
        else:
            tag = self.datasetTag
            if 'Fix' in tag and 'Exclude' in tag:
                tag = tag.replace('Fix ', '-Fix')
            tag = tag.replace('Fix ', 'Fix_')
            tag = tag.replace('Exclude ', 'Exclude_')
            tag = tag.replace(' ', '_')
        return tag


    def getFileNameFig(self, tag, tag2=False):
        figName = (f'{tag}-{self.enzymeName}-{self.getSaveTag()}-'
                   f'{self.seqLength}AA-MinCounts_{self.minCounts}.png')
        if self.motifFilter and 'entropy' not in tag.lower():
            figName = figName.replace(tag, f'{tag}-{self.iteration}')
        if tag2:
            figName = figName.replace('.png', f'{tag2}.png')
        if self.dropPos:
            figName = figName.replace('.png', f'-Drop_{"_".join(self.dropPos)}.png')
        figName = figName.replace('--', '-')
        return figName
    
    
    def getFileName(self, fType='Subs', datasetType='Exp', subProfile=False):
        if subProfile:
            tag = f'SubProfile_{self.datasetTagMotif.replace(' ', '_')}'
        else:
            if not self.datasetTag:
                print(f'Dont save, dataset tag: {self.datasetTag}\n')
                sys.exit()
            tag = self.datasetTag.replace(' ', '_')
        if fType == 'Subs':
            fileName = (f'{self.enzymeName}-{fType}_{datasetType}-{tag}-'
                        f'{self.seqLength}AA-MinCounts_{self.minCounts}.json')
        else:
            fileName = (f'{self.enzymeName}-AA_{fType}_{datasetType}-{tag}-'
                        f'{self.seqLength}AA-MinCounts_{self.minCounts}.csv')
        return fileName


    def logError(self, msg):
        self.jobDone = True
        print(f'\n{msg}')
        d = os.path.join(self.errorLog, f'jobID-{self.jobParams['Job ID']}')
        if not os.path.exists(d):
            os.makedirs(d, exist_ok=True)

        # Record error message
        path = os.path.join(d, 'error.log')
        if os.path.exists(path):
            msg = f'\n\n{"-"*80}\n\n{msg}'
        with open(path, 'a') as log:
            log.write(f'{msg}')

        # Record job params
        path = os.path.join(d, 'jobParams.log')
        with open(path, 'w') as log:
            log.write(f'Job Params:')
            for key, value in self.jobParams.items():
                log.write(f'\n* {key}: {value}')


    def saveSubstrates(self, substrates, datasetType='Exp', subProfile=False):
        saveTag = self.getFileName(datasetType=datasetType, subProfile=subProfile)

        # Save the substrates
        path = os.path.join(self.pathData, saveTag)
        print(f'Saving Substrates: {path}')
        with open(path, "w") as f:
            json.dump(substrates, f)


    def saveCounts(self, counts, datasetType, subProfile=False):
        saveTag = self.getFileName(fType='Counts', datasetType=datasetType,
                                   subProfile=subProfile)
        # Save the counts
        path = os.path.join(self.pathData, saveTag)
        print(f'Saving Counts: {path}')
        counts.to_csv(path)


    def initParams(self):
        ## Manually set the axis
        # self.xAxisSet = ['P4','P3','P2','P1','P1\'','P2\'','P3\'','P4\''][:self.seqLength]
        self.xAxisLabel = [f'R{index}' for index in range(1, self.seqLength + 1)]
        if self.xAxisSet:
            self.countsExp = pd.DataFrame(0, index=self.AA, columns=self.xAxisSet)
        else:
            self.countsExp = pd.DataFrame(0, index=self.AA, columns=self.xAxisLabel)
        self.countsBg = pd.DataFrame(0, index=self.AA, columns=self.xAxisLabel)
        self.countsExpSets = []
        self.countExpTotal = 0
        self.countExpUnique = 0
        self.rfExp = None
        self.rfExpScaled = None
        self.subsBg = {}
        self.countBgTotal = 0
        self.countBgUnique = 0
        self.rfBg = None
        self.entropy = pd.DataFrame(0.0, index=[], columns=['∆S'])
        self.entropyMax = None
        self.eMap = None
        self.eMapScaled = None
        self.eMapReleased = None
        self.eMapReleasedScaled = None
        self.datasetTypes = {'Exp': 'Experimental', 'Bg': 'Background'}

        self.datasetTag = ''
        self.datasetTagMotif = ''
        self.maxCounts = False
        self.minCounts = 1
        self.fileExp = []
        self.fileExpRev = []
        self.fileExpCounts = []
        self.fileBg = []
        self.fileBgRev = []
        self.avgBgRF = False
        self.subsExp = {}
        self.subsBg = {}
        self.subsPred = {}
        self.fileOrder = []
        self.dropPos = []
        self.figures = {}
        self.figTag = ''
        self.subProfile = False
        self.motifFilter = False
        self.initCOMET = False
        self.initFixedPos = []
        self.iteration = 0
        self.minS = 0.6
        self.minES = 0
        self.minESRel = -1 # -0.5
        self.motifPos = {}
        self.substrateProfile = pd.DataFrame()
        self.substrateProfileScl = pd.DataFrame()


    def jobInit(self, form, job):
        self.figures, self.jobParams = {}, {}
        self.combineProfiles, self.jobDone = False, False
        self.jobParams['Job ID'] = form['jobID']
        self.jobParams['Job'] = job
        self.enzymeName = form['enzymeName']
        self.jobParams['Enzyme Name'] = self.enzymeName
        if self.enzymeName.lower() == 'mpro2':
            self.enzymeName = f'SARS-CoV-2 M{'ᵖʳᵒ'}'
        elif self.enzymeName.lower() == 'mpro':
            self.enzymeName = f'SARS-CoV M{'ᵖʳᵒ'}'

        # Initialize directories
        date = time.strftime("%Y_%m_%d-%H:%M:%S").replace(':','_')
        self.pathDir = os.path.join('datasets', self.enzymeName).replace(' ', '_')
        self.pathData = os.path.join(self.pathDir, 'data')
        self.pathFigs = os.path.join(self.pathDir,
                                     f'figures-{job}-{date}').replace(' ', '_')
        self.pathLog = os.path.join(self.pathFigs, 'log.txt')
        if self.pathDir is not None:
            if not os.path.exists(self.pathDir):
                os.makedirs(self.pathDir, exist_ok=True)
        if self.pathData is not None:
            if not os.path.exists(self.pathData):
                os.makedirs(self.pathData, exist_ok=True)
        if self.pathFigs is not None:
            if not os.path.exists(self.pathFigs):
                os.makedirs(self.pathFigs, exist_ok=True)

        self.log() # Clear the log
        self.log('================================ Job Summary '
                 '=================================')
        self.log(f'Job ID: {self.jobParams['Job ID']}')
        self.log(f'Job: {job}')
        self.log(f'Enzyme: {self.enzymeName}')
        self.seqLength = int(form['seqLength'])

        def addFile(data, val):
            if isinstance(val, list):
                for f in val:
                    data.append(f)
            else:
                data.append(val)

        # Job dependant parameters
        self.idxMotif = {}
        if job == 'Process DNA':
            self.seq5Prime = form['seq5Prime']
            self.seq3Prime = form['seq3Prime']
            self.minPhred = (
                int(round(float(form['minPhred'])))
            ) if form['minPhred'] != '' else 0
            self.log(f'5\' Sequence: {self.seq5Prime}\n'
                     f'3\' Sequence: {self.seq3Prime}\n'
                     f'Min Phred Score: {self.minPhred}')
        elif job == 'Filter Motif':
            self.saveData = False
            self.setS = True
            try:
                self.maxCounts = int(form['maxCounts'])
            except ValueError:
                self.maxCounts = False
        elif job == 'Combine Motifs' or job == 'Predict Substrate Activity':
            if job == 'Combine Motifs':
                self.combineProfiles = True
            self.seqLength = int(form['motifLength'])
            for key in form.keys():
                if 'idxStart' in key:
                    self.idxMotif[key] = int(form[key]) - 1
        elif job not in ['Filter AA']:
            print('ERROR: What Script Is Running')
            sys.exit()
        if job == 'Filter AA' or job == 'Filter Motif' or job == 'Combine Motifs':
            try:
                self.maxCounts = int(form['maxCounts'])
            except ValueError:
                self.maxCounts = False

        # Initialize params
        self.initParams()
        if job == 'Filter Motif':
            self.motifFilter = True
            self.initCOMET = True
        elif job == 'Predict Substrate Activity':
            # self.avgBgRF = True
            self.subsPred = [sub.strip() for sub in form['predSubs'].split(',')]

        if 'dropPos' in form.keys():
            self.dropPos = form['dropPos']
            if isinstance(self.dropPos, str):
                self.dropPos = [self.dropPos]
            self.jobParams['Drop position'] = ", ".join(self.dropPos)
            self.log(f'Drop position: {", ".join(self.dropPos)}')
            self.seqLength -= len(self.dropPos)
        self.jobParams['Substrate Length'] = self.seqLength
        self.log(f'Substrate Length: {self.seqLength}')

        # Conditional params
        if not job == 'Predict Substrate Activity':
            self.minCounts = int((round(float(form['minCounts']))))
            self.jobParams['Minimum Count'] = self.minCounts
            self.log(f'Minimum Count: {self.minCounts}')

        # Get the files
        for key, value in form.items():
            if 'fileBgRev' in key:
                addFile(self.fileBgRev, value)
            elif 'fileBg' in key:
                addFile(self.fileBg, value)
            elif 'fileExpCounts' in key:
                self.fileOrder.append(key)
                addFile(self.fileExpCounts, value)
            elif 'fileExpRev' in key:
                addFile(self.fileExpRev, value)
            elif 'fileExp' in key:
                addFile(self.fileExp, value)
        if len(self.fileExpCounts) > 1:
            self.avgBgRF = True

        # Get the filter and initialize the data structures
        self.getFilter(form)
        if job == 'Filter Motif':
            self.figTag = (f'Assign Register\n'
                           f'{self.datasetTagMotif.replace('Register', '')}')
        elif job == 'Combine Motifs':
            self.figTag = f'Specificity Landscape'


    def log(self, txt=None):
        if txt is None:
            with open(self.pathLog, 'w'):
                pass
        else:
            with open(self.pathLog, 'a') as log:
                log.write(f'{txt}\n')


    def logInQueue(self, logQueue):
        with open(self.pathLog, 'a') as log:
            while not logQueue.empty():
                log.write(logQueue.get() + '\n')


    def processSubs(self, substrates, datasetType, filteredAA):
        self.log('\n\n================================= Substrates '
                 '=================================')
        self.log(f'Dataset: {self.datasetTypes[datasetType]}')

        # Inspect sequences
        if not filteredAA:
            filteredSubs = {}
            for substrate, count in substrates.items():
                for AA in substrate:
                    if AA not in self.AA:
                        filteredSubs[substrate] = count
            if filteredSubs:
                self.log(f'\nFiltering Substrates:\n'
                         f'     If a substrate contains an '
                         f'unaccented AA it will be removed.\n'
                         f'     Accepted: {self.AA}\n\n'
                         f'     Removed Substrates:')
                for substrate, count in filteredSubs.items():
                    substrates.pop(substrate, count)
                    self.log(f'          {substrate}: {count}')

        # Sort data
        substrates = dict(
            sorted(substrates.items(),
                   key=lambda item: item[1], reverse=True)
        )


        # Count AAs
        countMatrix = pd.DataFrame(0, index=self.AA, columns=self.xAxisLabel)
        self.log('\nSubstrate Totals:')
        if datasetType == 'Exp':
            self.subsExp = substrates
            self.countExpTotal = sum(substrates.values())
            self.countExpUnique = len(substrates.keys())
            self.log(f'     Total Substrates: {self.countExpTotal:,}\n'
                     f'    Unique Substrates: {self.countExpUnique:,}\n')

            # Record job params
            self.jobParams['Total Experimental Substrates'] = f'{self.countExpTotal:,}'
            self.jobParams['Unique Experimental Substrates'] = f'{self.countExpUnique:,}'
        elif datasetType == 'Bg':
            self.subsBg = substrates
            self.countBgTotal = sum(substrates.values())
            self.countBgUnique = len(substrates.keys())
            self.log(f'     Total Substrates: {self.countBgTotal:,}\n'
                     f'    Unique Substrates: {self.countBgUnique:,}\n')

            # Record job params
            self.jobParams['Total Background Substrates'] = f'{self.countBgTotal:,}'
            self.jobParams['Unique Background Substrates'] = f'{self.countBgUnique:,}'
        else:
            self.logError(f'ERROR: processSubs()\n'
                          f'* Unknown dataset type: {datasetType}')

        self.log(f'Top {self.printN:,} {self.datasetTypes[datasetType]} Sequences')
        for index, (sub, count) in enumerate(substrates.items()):
            if index >= self.printN:
                break
            self.log(f'     {sub}: {count:,}')

        # Save data
        self.saveSubstrates(substrates=substrates, datasetType=datasetType)

        # Count AAs
        countMatrix = self.countAA(substrates=substrates, countMatrix=countMatrix,
                                   datasetType=datasetType)
        if datasetType == 'Bg':
            self.countsBg = countMatrix
        else:
            self.countsExp = countMatrix


    def loadSubstrates(self, data, queueData, queueLog):
        data.seek(0)  # Ensure at start
        substrates = json.load(data)
        queueData.put(substrates)
        queueLog.put(f'     {data.filename}')


    def loadCounts(self, data, queueData, queueLog):
        # Load file
        df = pd.read_csv(data, index_col=0)
        df = df.astype(int)
        queueData.put(df)
        queueLog.put(f'     {data.filename}\n')


    def countAA(self, substrates, countMatrix, datasetType, subProfile=False):
        self.log('\n\n================================== Count AA '
                 '==================================')
        self.log(f'Dataset: {self.datasetTypes[datasetType]}')
        countMatrix.loc[:, :] = 0
        totalCounts = pd.DataFrame(0, index=self.xAxisLabel, columns=['Sum'])

        def splitData(data, matrix):
            cores = os.cpu_count()
            size = max(1, int(np.ceil(len(data) / cores)))
            batches = list(batched(data.items(), size))
            #print(f'Cores: {cores}, Batches: {len(batches)}, Seq/Batch: {size:,}')
            args = [(dict(batch), matrix.columns.tolist(), matrix.index.tolist())
                    for batch in batches]

            with ProcessPoolExecutor(max_workers=cores) as executor:
                results = list(executor.map(counter, args))
            return sum(results[1:], results[0])

        start = time.time()
        countMatrix = splitData(substrates, countMatrix)
        end = time.time()
        runtime = (end - start) / 60
        runtime = round(runtime, 3)
        print(f'Counted {len(substrates.items()):,} unique substrates in: {runtime} min')
        self.log(f'\nCounts:\n{countMatrix}')

        for pos in countMatrix.columns:
            counts = sum(countMatrix.loc[:, pos])
            totalCounts.loc[pos, 'Sum'] = counts
        self.log(f'\nCount Totals:\n{totalCounts}')

        if self.saveData and not subProfile:
            self.saveCounts(
                counts=countMatrix, datasetType=datasetType, subProfile=subProfile
            )
        return countMatrix


    def loadDNA(self, path, seqs, datasetType, queueLog, reverseRead):
        fileName = path.filename if hasattr(path, 'filename') else path.name
        try:
            # Translate the dna
            self.translate(path, seqs, fileName, datasetType,
                           queueLog, reverseRead)
        except Exception as e:
            self.logError(f'ERROR: loadDNA()\n'
                          f'* File name: {path.filename}\n\n{e}')


    def translate(self, path, substrates, fileName, datasetType, queueLog, revRead):
        queueLog.put('\n\n================================ Translate DNA '
                     '===============================')
        queueLog.put(f'File Name: {fileName}')
        if revRead:
            queueLog.put('Read Type: Reverse Read')
        else:
            queueLog.put('Read Type: Forward Read')
        queueLog.put(f'  Dataset: {datasetType}')

        # Log tag
        if self.fixAA or self.exclAA:
            queueLog.put(f'   Filter: {self.datasetTag}')

        # Inspect the datasetType parameter
        if (datasetType != self.datasetTypes['Exp']
                and datasetType != self.datasetTypes['Bg']):
            self.logError(f'ERROR: translate()\n'
                          f'* Unknown dataset type: {datasetType}')


        def reverseComplement(seq):
            """Returns the reverse complement of a DNA sequence."""
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
            return ''.join(complement[base] for base in reversed(seq))


        def extractionEfficiency(nDataset, nExtracted, fullSet=False):
            if nDataset == 0:
                nDataset = 1
            perExtracted = (nExtracted / nDataset) * 100
            if fullSet:
                queueLog.put(f'\nExtraction Efficiency: - All Sequences')
            else:
                queueLog.put(f'\nExtraction Efficiency: {self.printN} Substrates')
            queueLog.put(f'     Evaluated DNA Sequences: {nDataset:,}\n'
                         f'        Extracted Substrates: {nExtracted:,}\n'
                         f'       Extraction Efficiency: {round(perExtracted, 3)} %')


        # Open the file
        if path.filename.endswith('.gz'):
            fileHandle = gzip.open(path, 'rt')
        else:
            path.seek(0) # Ensure at start
            fileHandle = io.StringIO(path.read().decode('utf-8'))
        data, useQS = None, False
        if path.filename.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
            data = SeqIO.parse(fileHandle, 'fastq')
            useQS = True
            queueLog.put(f'  Eval QS: {useQS}\n'
                         f'Min Phred: {self.minPhred}')
        elif path.filename.endswith(('.fasta', '.fa', '.fasta.gz', '.fa.gz')):
            data = SeqIO.parse(fileHandle, 'fasta')
        queueLog.put('\n')

        # Process file
        totalSeqs, totalSubsExtracted, logData = 0, 0, True
        for datapoint in data:
            if logData and totalSubsExtracted >= self.printN:
                extractionEfficiency(totalSeqs, totalSubsExtracted)
                logData = False
            totalSeqs += 1
            dna = datapoint.seq
            if revRead:
                dna = reverseComplement(dna)
            if logData:
                queueLog.put(f'DNA Seq: {dna}')
                queueLog.put(f'* 5\' Found: {'Yes' if self.seq5Prime in dna else 'No'}\n'
                             f'* 3\' Found: {'Yes' if self.seq3Prime in dna else 'No'}')
            if self.seq5Prime in dna or self.seq3Prime in dna:
                # Find: Substrate indices
                start = dna.find(self.seq5Prime) + len(self.seq5Prime)
                end = dna.find(self.seq3Prime)
                if end == -1:
                    end = False

                # Extract substrate dna seq
                if end:
                    substrateDNA = dna[start:end].strip()
                else:
                    substrateDNA = dna[start:].strip()
                if logData:
                    queueLog.put(f'* Sub DNA: {substrateDNA}')
                if len(substrateDNA) == self.seqLength * 3:
                    substrate = str(Seq.translate(substrateDNA))
                    if logData:
                        queueLog.put(f'* Sub Seq: {substrate}')

                    # Inspect substrate seq: log only
                    if 'X' not in substrate and '*' not in substrate:
                        keepSub = True
                        if useQS:
                            qs = datapoint.letter_annotations['phred_quality']
                            if revRead:
                                qs = qs[::-1]  # Reverse QS
                            if end:
                                qs = qs[start:end]
                            else:
                                qs = qs[start:]
                            if logData:
                                queueLog.put(f'* QS: {qs}')
                            if not all(score >= self.minPhred for score in qs):
                                keepSub = False

                        # Filter substrate
                        if self.fixAA or self.exclAA:
                            for posExcl, exclAA in self.exclAA.items():
                                if exclAA:
                                    if not isinstance(exclAA, list):
                                        exclAA = list(exclAA)
                                    idx = int(posExcl.replace('R', '')) - 1
                                    if substrate[idx] in exclAA:
                                        keepSub = False
                                        break
                            if not keepSub:
                                continue

                            for posFix, fixAA in self.fixAA.items():
                                if fixAA:
                                    if not isinstance(fixAA, list):
                                        fixAA = list(fixAA)
                                    idx = int(posFix.replace('R', '')) - 1
                                    if substrate[idx] not in fixAA:
                                        keepSub = False
                                        break
                        if keepSub:
                            if logData:
                                queueLog.put(f'Keep Substrate')
                            substrates[substrate] = (
                                    substrates.get(substrate, 0) + 1
                            )
                            totalSubsExtracted += 1
            if logData:
                queueLog.put('')
        extractionEfficiency(totalSeqs, totalSubsExtracted, fullSet=True)
        print(f'Finished Translating: {fileName}')


    def evalDNA(self, form):
        try:
            if form is None:
                self.jobParams['Invalid File'] = True
                return
            self.jobInit(form, job='Process DNA')

            # Load the data
            threads = []
            subsExp, queuesExpLog, queuesExpRevLog = [], [], []
            subsBg, queuesBgLog, queuesBgRevLog= [], [], []
            if self.fileExp:
                for file in self.fileExp:
                    subsLocal = {}
                    queueLog = queue.Queue() # maxsize=100)
                    queuesExpLog.append(queueLog)
                    thread = threading.Thread(
                        target=self.loadDNA,
                        args=(file, subsLocal, self.datasetTypes['Exp'],
                              queueLog, False,)
                    )
                    subsExp.append(subsLocal)
                    threads.append(thread)
                    thread.start()
            if self.fileExpRev:
                for file in self.fileExpRev:
                    subsLocal = {}
                    queueExpRevLog = queue.Queue() # maxsize=100)
                    queuesExpRevLog.append(queueExpRevLog)
                    thread = threading.Thread(
                        target=self.loadDNA,
                        args=(file, subsLocal, self.datasetTypes['Exp'],
                              queueExpRevLog, True,)
                    )
                    subsExp.append(subsLocal)
                    threads.append(thread)
                    thread.start()
            if self.fileBg:
                for file in self.fileBg:
                    subsLocal = {}
                    queueLog = queue.Queue() # maxsize=100)
                    queuesBgLog.append(queueLog)
                    thread = threading.Thread(
                        target=self.loadDNA,
                        args=(file, subsLocal, self.datasetTypes['Bg'],
                              queueLog, False,)
                    )
                    subsBg.append(subsLocal)
                    threads.append(thread)
                    thread.start()
            if self.fileBgRev:
                for file in self.fileBgRev:
                    subsLocal = {}
                    queueBgRevLog = queue.Queue() # maxsize=100)
                    queuesBgRevLog.append(queueBgRevLog)
                    thread = threading.Thread(
                        target=self.loadDNA,
                        args=(file, subsLocal, self.datasetTypes['Bg'],
                              queueBgRevLog, True,)
                    )
                    subsBg.append(subsLocal)
                    threads.append(thread)
                    thread.start()
            # Wait for all threads to finish
            for thread in threads:
                thread.join()

            # Combine dictionaries
            for subs in subsExp:
                for sub, count in subs.items():
                    self.subsExp[sub] = self.subsExp.get(sub, 0) + count
            for subs in subsBg:
                for sub, count in subs.items():
                    self.subsBg[sub] = self.subsBg.get(sub, 0) + count

            # Log the output
            if queuesExpLog:
                for log in queuesExpLog:
                    self.logInQueue(log)
            if queuesExpRevLog:
                for log in queuesExpRevLog:
                    self.logInQueue(log)
            if queuesBgLog:
                for log in queuesBgLog:
                    self.logInQueue(log)
            if queuesBgRevLog:
                for log in queuesBgRevLog:
                    self.logInQueue(log)
            if self.dropPos:
                self.removePos()

            # Make figures
            if self.subsExp:
                self.subsExp = dict(
                    sorted(self.subsExp.items(),
                           key=lambda item: item[1],
                           reverse=True)
                )

                # Sort substrates and count AA
                self.processSubs(
                    substrates=self.subsExp, datasetType='Exp', filteredAA=False
                )

                # Plot counts
                self.figures['exp_counts'] = (
                    self.plotMatrix(
                        data=self.countsExp, totalCounts=self.countExpTotal,
                        figLabel='Experimental Counts',
                        datasetType=self.datasetTypes['Exp']
                    )
                )

                self.figures['barCounts'] = self.plotBars(
                    self.subsExp, dataType='Counts'
                )
                self.figures['barRF'] = self.plotBars(
                    self.subsExp, dataType='Relative Frequency'
                )
            if self.subsBg:
                self.subsBg = dict(
                    sorted(self.subsBg.items(),
                           key=lambda item: item[1],
                           reverse=True)
                )

                # Sort substrates and count AA
                self.processSubs(
                    substrates=self.subsBg, datasetType='Bg', filteredAA=False
                )

                # Plot counts
                self.figures['bg_counts'] = (
                    self.plotMatrix(
                        data=self.countsBg, totalCounts=self.countBgTotal,
                        figLabel='Background Counts', datasetType=self.datasetTypes['Bg']
                    )
                )

            if self.subsExp and self.subsBg:
                self.calculateRF()
                self.calculateEntropy()
                self.evalEnrichment()
            self.jobDone = True
        except Exception as e:
            self.logError(f'ERROR: evalDNA()\n{e}')


    def removePos(self):
        self.log('\n\n================================ Drop Position '
                 '===============================')
        self.log(f'Removing: {", ".join(self.dropPos)}')

        # Organize drop indices
        indices = []
        for pos in reversed(self.dropPos):
            idx = int(pos.replace('R', '')) - 1
            self.xAxisLabel.remove(pos)
            indices.append(idx)

        def dropAA(seq):
            for i in indices:
                try:
                    seq = seq[:i] + seq[i + 1:]
                except:
                    seq = seq[:i]
            return seq

        # Remove AA
        if self.subsExp:
            self.log('\nExperimental Substrates:')
            for idx, (substrate, counts) in enumerate(self.subsExp.items(), start=1):
                sub = dropAA(substrate)
                self.log(f'* {substrate} -> {sub}')
                if idx >= self.printN:
                    break

            subs = {}
            for substrate, counts in self.subsExp.items():
                sub = dropAA(substrate)
                subs[sub] = subs.get(sub, 0) + counts
            self.subsExp = dict(sorted(
                subs.items(), key=lambda item: item[1], reverse=True)
            )

        if self.subsBg:
            self.log('\nBackground Substrates')
            for idx, (substrate, counts) in enumerate(self.subsBg.items(), start=1):
                sub = dropAA(substrate)
                self.log(f'* {substrate} -> {sub}')
                if idx >= self.printN:
                    break
            self.log('\n')

            subs = {}
            for substrate, counts in self.subsBg.items():
                sub = dropAA(substrate)
                subs[sub] = subs.get(sub, 0) + counts
            self.subsBg = dict(sorted(
                subs.items(), key=lambda item: item[1], reverse=True)
            )


    def evalData(self, form, filterMotifs=False, combineProfiles=False, predActivity=False):
        if form is None:
            self.jobParams['Invalid File'] = True
            return
        self.jobDone, self.blockScript = False, True
        if filterMotifs:
            self.jobInit(form, job='Filter Motif')
        elif combineProfiles:
            self.jobInit(form, job='Combine Motifs')
        elif predActivity:
            self.jobInit(form, job='Predict Substrate Activity')
        else:
            self.jobInit(form, job='Filter AA')
        self.log('\n\n================================== Load Data '
                 '=================================')

        # Load the data
        self.subsBg, self.subsExp, self.profiles = {}, {}, []
        threads = []
        queuesExp, queuesExpLog = [], []
        queuesExpCounts, queuesExpCountsLog = [], []
        queuesBg, queuesBgLog = [], []
        for file in self.fileExp:
            queueExp = queue.Queue()
            queueLog = queue.Queue()
            queuesExp.append(queueExp)
            queuesExpLog.append(queueLog)
            thread = threading.Thread(
                target=self.loadSubstrates, args=(file, queueExp, queueLog)
            )
            thread.start()
            threads.append(thread)
        for file in self.fileExpCounts:
            queueExp = queue.Queue()
            queueLog = queue.Queue()
            queuesExpCounts.append(queueExp)
            queuesExpCountsLog.append(queueLog)
            thread = threading.Thread(
                target=self.loadCounts, args=(file, queueExp, queueLog)
            )
            thread.start()
            threads.append(thread)
        for file in self.fileBg:
            queueBg = queue.Queue()
            queueLog = queue.Queue()
            queuesBg.append(queueBg)
            queuesBgLog.append(queueLog)
            thread = threading.Thread(
                target=self.loadCounts, args=(file, queueBg, queueLog)
            )
            thread.start()
            threads.append(thread)
        for thread in threads:
            thread.join() # Wait for all threads to finish

        # Process input
        if queuesExpLog:
            self.log('Loading Substrates: Experimental')
            # for log in queuesExpLog:
            #     self.logInQueue(log)
        if queuesExp:
            for q in queuesExp:
                substrates = q.get()
                if combineProfiles:
                    self.profiles.append(substrates)
                else:
                    for substrate, count in substrates.items():
                        if substrate in self.subsExp.keys():
                            self.subsExp[substrate] += count
                        else:
                            self.subsExp[substrate] = count
            if self.subsExp:
                self.log('\nSubstrates:')
                for i, (substrate, count) in enumerate(self.subsExp.items()):
                    self.log(f'    {substrate}, {count:,}')
                    if i >= self.printN:
                        break
                self.log(f'\nTotal Substrates: {sum(self.subsExp.values()):,}\n'
                         f'Unique Substrates: {len(self.subsExp.keys()):,}')
            elif self.profiles:
                for idx, profile in enumerate(self.profiles):
                    self.log(f'\nProfile: {idx}')
                    for i, (substrate, count) in enumerate(profile.items()):
                        self.log(f'    {substrate}, {count:,}')
                        if i >= self.printN:
                            break
                self.getMotifSeq()
            else:
                self.logError(f'ERROR: loadSubstrates()\n'
                              f'No experimental substrates were loaded.')

        if queuesExpCountsLog:
            self.log('Loading Dataset: Experimental')
        if queuesExpCounts:
            # Define columns
            columnLabels = self.xAxisLabel
            if self.xAxisSet:
                columnLabels = self.xAxisSet
            
            # Extract motif counts
            for idx, q in enumerate(queuesExpCounts):
                counts = q.get()
                idxN = self.idxMotif[self.fileOrder[idx].replace('fileExpCounts',
                                                                 'idxStart')]
                idxC = idxN + self.seqLength
                c = counts.iloc[:, idxN:idxC]
                if len(queuesExpCounts) > 1:
                    self.log(f'Counts: {idx}\n{counts}\n')
                    countMotif = c.copy()
                    countMotif.columns = columnLabels
                    self.countsExpSets.append(countMotif)
                else:
                    self.log(f'Counts:\n{counts}\n')
                if (len(c.columns) == len(counts.columns) and
                        (c.columns != counts.columns).any()):
                    self.log(f'Extracted Motif Counts:\n{c}')
                c.columns = columnLabels
                self.countsExp += c

            if len(queuesExpCounts) > 1:
                self.log(f'\nCombined Registers:\n{self.countsExp}\n')
            self.countExpTotal = self.countsExp.sum()
            for pos in columnLabels:
                if self.countExpTotal[pos] == 0:
                    self.logError(f'ERROR: loadCounts()\n'
                                  f'* Experimental count matrix contains an empty '
                                  f'column.\n\n* Total Counts:\n{self.countExpTotal}')
                    break
            self.log(f'Total Counts:\n{self.countExpTotal.to_string()}')
        if queuesBgLog:
            self.log('\n\nLoading Counts: Background')
        if queuesBg:
            for idx, q in enumerate(queuesBg):
                counts = q.get()
                self.countsBg += counts
                if len(queuesBg) > 1:
                    self.log(f'Counts: {idx}\n{counts}\n')
                else:
                    self.log(f'Counts:\n{counts}\n')
            self.countBgTotal = sum(self.countsBg.iloc[:, 0])
            if self.countBgTotal == 0:
                self.logError(f'ERROR: loadCounts()\n'
                              f'* No background counts were loaded.')
            if len(queuesBg) > 1:
                self.log(f'\nCombined Background Counts:\n{self.countsBg}\n')

        if self.dropPos:
            self.removePos()
            self.countsExp.drop(self.dropPos, axis=1, inplace=True)
            self.countsBg.drop(self.dropPos, axis=1, inplace=True)
            self.log(f'\nBackground Counts:\n{self.countsBg}\n')

        # Use data
        try:
            if filterMotifs:
                self.subsExpAll = self.subsExp.copy()
                self.filterSubs()
                self.initFixedPos = list(self.fixAA.keys())
                self.selectMotifPos()
                self.evalEnrichment(skipFigs=True)
            elif predActivity:
                self.predictActivity()
                try:
                    os.rmdir(self.pathData) # delete empty dir
                except OSError:
                    pass
            elif self.combineProfiles:
                self.evalProfiles()
            else:
                ###
                self.filterSubs()
                self.figTag = (f'N Experimental {sum(self.subsExp.values()):,}\n'
                               f'N Background {sum(self.countsBg.iloc[:, 1]):,}')
                self.evalEnrichment()
            if not filterMotifs and not predActivity:
                self.figures['barCounts'] = self.plotBars(
                    self.subsExp, dataType='Counts'
                )
                self.figures['barRF'] = self.plotBars(
                    self.subsExp, dataType='Relative Frequency'
                )
        except Exception as e:
            self.logError(f'ERROR: evalData()\n* Job: {self.jobParams['Job']}\n\n{e}')
        self.jobDone = True
        if filterMotifs:
            time.sleep(2)
            self.jobDone = False
        self.blockScript = False


    def filterSubs(self, allSubs=False, plotEntropy=True, subProfile=False,
                   filterPos=False, relPos=False):
        self.log('\n\n============================== Filter Substrates '
                 '=============================')
        if filterPos:
            self.log(f'Filtering: {filterPos}\n')
        elif relPos:
            self.log(f'Releasing: {relPos}\n')
        else:
            self.log(f'Filter: {self.datasetTag}\n')
        if self.fixAA:
            self.log('Filter AA:')
            for pos, aa in self.fixAA.items():
                self.log(f'* {pos}: {aa}')
            self.log('')
        if self.exclAA:
            self.log('Exclude AA:')
            for pos, aa in self.exclAA.items():
                self.log(f'* {pos}: {aa}')
            self.log('')

        # Select data
        if allSubs:
            substrates = self.subsExpAll
        else:
            substrates = self.subsExp
        totalSubs = 0
        for count in substrates.values():
            totalSubs += count
        totalSubsUnique = len(substrates.keys())
        self.log(f'Unfiltered Substrates:\n'
                 f'     Total: {totalSubs:,}\n'
                 f'    Unique: {totalSubsUnique:,}')

        if self.fixAA or self.exclAA:
            subs = {}
            for substrate, count in substrates.items():
                if count < self.minCounts:
                    continue

                keepSub = True
                for posExcl, exclAA in self.exclAA.items():
                    if exclAA:
                        if not isinstance(exclAA, list):
                            exclAA = list(exclAA)
                        idx = int(posExcl.replace('R', '')) - 1
                        if substrate[idx] in exclAA:
                            keepSub = False
                            break
                if not keepSub:
                    continue

                for posFix, fixAA in self.fixAA.items():
                    if fixAA:
                        if not isinstance(fixAA, list):
                            fixAA = list(fixAA)
                        idx = int(posFix.replace('R', '')) - 1
                        if substrate[idx] not in fixAA:
                            keepSub = False
                            break
                if keepSub:
                    subs[substrate] = count
            self.subsExp = subs
        self.countExpTotal = sum(self.subsExp.values())
        self.countExpUnique = len(self.subsExp.keys())
        self.subsExp = dict(sorted(
            self.subsExp.items(), key=lambda item: item[1], reverse=True)
        )

        # Log substrates
        self.log(f'\nFiltered Substrates:\n'
                 f'     Total: {self.countExpTotal:,}\n'
                 f'    Unique: {self.countExpUnique:,}\n')
        self.log('Substrates:')
        for i, (substrate, count) in enumerate(self.subsExp.items()):
            self.log(f'    {substrate}, {count:,}')
            i += 1
            if i >= self.printN:
                break

        # Save data
        if self.saveData:
            self.saveSubstrates(
                substrates=self.subsExp, datasetType='Exp', subProfile=subProfile
            )

        # Count AAs
        self.countsExp = self.countAA(
            substrates=self.subsExp, countMatrix=self.countsExp,
            datasetType='Exp', subProfile=subProfile
        )
        self.calculateRF()
        self.calculateEntropy(plotFig=plotEntropy)


    def comet(self, form):
        if form is None:
            self.jobParams['Invalid File'] = True
            return
        try:
            self.jobDone = False
            self.evalEnrichment(skipFigs=True)
            self.log('\n\n================================ Filter Motif '
                     '================================')
            self.minS = float(form['minS'])
            self.minES = float(form['minES'])
            self.minESRel = float(form['minESRel'])
            self.log(f'Minimum ∆S: {self.minS}\n'
                     f'Minimum ES Filter: {self.minES}\n'
                     f'Minimum ES Rel: {self.minESRel}\n')
            self.jobParams['Minimum ∆S'] = self.minS
            self.jobParams['Minimum ES Filter'] = self.minES
            self.jobParams['Minimum ES Release'] = self.minESRel
            self.selectMotifPos()
            self.log('\nRecognition Sites:')
            self.log(
                pd.DataFrame.from_dict(self.motifPos, orient='index', columns=['∆S'])
            )
            figTag = self.datasetTagMotif.replace('Register', '')

            def evalAAs(position, minES):
                self.fixAA[position] = []
                for aa in self.eMap.index:
                    if self.eMap.loc[aa, position] >= minES:
                        self.fixAA[position].append(aa)
                # print(f'Filter (minES={minES}): {position} - {self.fixAA[position]}')


            # Motif ID
            self.figTag = f'Motif ID\n{figTag}'
            for idx, pos in enumerate(self.motifPos.keys()):
                self.jobDone = False
                if pos in self.initFixedPos:
                    continue # Skip already fix pos
                if pos not in self.fixAA.keys():
                    evalAAs(pos, self.minES)
                    self.filterSubs(plotEntropy=False, filterPos=pos)
                    self.iteration += 1
                    self.evalEnrichment(skipFigs=True)


            # Refine Filter
            self.figTag = f'Refine Filter\n{figTag}'
            for posRel in self.motifPos.keys():
                self.iteration += 1
                self.jobDone = False
                self.fixAA, filterAA = {}, []

                # Release
                for pos in self.motifPos.keys():
                    if pos != posRel:
                        filterAA.append(pos)
                for posFix in filterAA:
                    evalAAs(posFix, self.minESRel)
                self.filterSubs(allSubs=True, plotEntropy=False, relPos=posRel)
                self.evalEnrichment(skipFigs=True)

                # Refilter
                for posFix in self.motifPos.keys():
                    evalAAs(posFix, self.minESRel)
                self.iteration += 1
                self.filterSubs(allSubs=True, plotEntropy=False, filterPos=posRel)
                self.evalEnrichment(skipFigs=True)


            # Release Filter
            self.figTag = f'Release Filter\n{figTag}'
            self.fixAA, fixAA = {}, self.fixAA
            exclAA = self.exclAA.copy()
            idxEnd = len(self.motifPos.keys()) - 1
            counts = self.countsExp.copy()
            self.substrateProfile = pd.DataFrame(0, index=self.eMap.index,
                                                 columns=self.eMap.columns)
            for idx, posRel in enumerate(self.motifPos.keys()):
                self.iteration += 1
                self.fixAA, self.exclAA = {}, exclAA
                if posRel in self.exclAA.keys():
                    self.exclAA.pop(posRel)
                for pos in self.motifPos.keys():
                    if pos != posRel:
                        # filter.append(pos)
                        self.fixAA[pos] = fixAA[pos]
                if idx == idxEnd:
                    self.exclAA = {}
                    self.saveData = True
                    self.filterSubs(allSubs=True, plotEntropy=False,
                                    subProfile=True, relPos=posRel)
                    self.substrateProfile.loc[:, posRel] = self.countsExp.loc[:, posRel]
                    self.countsExp = self.substrateProfile.copy()
                    self.calculateRF()
                    self.calculateEntropy(plotFig=False)
                    self.evalEnrichment(releasedCounts=True, skipFigs=True)

                    # Drop exclusions
                    if exclAA:
                        evalAAs(posRel, self.minESRel)
                        self.filterSubs(allSubs=True, plotEntropy=False,
                                        subProfile=True, filterPos=posRel)
                        newFig = False
                        for pos, aa in exclAA.items():
                            if pos not in self.motifPos.keys():
                                newFig = True
                                self.substrateProfile.loc[:, pos] = (
                                    self.countsExp.loc[:, pos])
                        if newFig:
                            self.countsExp = self.substrateProfile.copy()
                            self.calculateRF()
                            self.calculateEntropy(plotFig=False)
                            self.evalEnrichment(releasedCounts=True, skipFigs=True)
                else:
                    self.filterSubs(allSubs=True, plotEntropy=False, relPos=posRel)
                    self.substrateProfile.loc[:, posRel] = self.countsExp.loc[:, posRel]
                    self.countsExp = self.substrateProfile.copy()
                    self.calculateRF()
                    self.calculateEntropy(plotFig=False)
                    self.evalEnrichment(releasedCounts=True, skipFigs=True)

            # Populate non-motif pos
            for pos in counts.columns:
                if pos not in self.motifPos.keys() and pos not in exclAA.keys():
                    self.substrateProfile.loc[:, pos] = counts.loc[:, pos]
                    self.countsExp.loc[:, pos] = counts.loc[:, pos]
            self.subProfile = True
            self.iteration += 1
            self.figTag = f'Substrate Profile\n{figTag}'
            self.calculateRF()
            self.calculateEntropy()
            self.evalEnrichment(releasedCounts=True, noFigs=True)
            self.saveCounts(counts=self.substrateProfile,
                            datasetType='Exp', subProfile=True)

            # Plot profile
            self.figures['eMapProfile'] = (
                self.plotEnrichmentScores(dataType='Enrichment')
            )
            self.figures['eMapScProfile'] = (
                self.plotEnrichmentScores(dataType='Scaled Enrichment')
            )
            self.plotEnrichmentLogo()
            self.figures['wordsProfile'] = self.plotWordCloud(self.subsExp)
            self.figures['barCountsProfile'] = self.plotBars(
                self.subsExp, dataType='Counts'
            )
            self.figures['barRFProfile'] = self.plotBars(
                self.subsExp, dataType='Relative Frequency'
            )
        except Exception as e:
            self.logError(f'ERROR: comet()\n\n{e}')
        time.sleep(3)
        self.jobDone = True


    def getMotifSeq(self):
        motifs = {}
        for idx, profile in enumerate(self.profiles):
            # print(f'Profile {idx}, {list(self.idxMotif.keys())[idx]}')
            idxN = self.idxMotif[list(self.idxMotif.keys())[idx]]
            idxC = idxN + self.seqLength
            # print(f'Profile ({idx}): {idxN}-{idxC}')
            for i, (substrate, count) in enumerate(profile.items()):
                motif = substrate[idxN:idxC]
                if motif in motifs.keys():
                    motifs[motif] += count
                else:
                    motifs[motif] = count
        self.subsExp = dict(sorted(
            motifs.items(), key=lambda item: item[1], reverse=True)
        )

        self.log('\nMotif Sequences:')
        for i, (motif, count) in enumerate(self.subsExp.items()):
            self.log(f'* {motif}: {count:,}')
            if i >= self.printN:
                break
        self.log('\n')


    def evalProfiles(self):
        self.calculateRF()
        if len(self.countsExpSets) > 1:
            self.evalStats()
        self.calculateEntropy()
        self.evalEnrichment()


    def selectMotifPos(self):
        self.motifPos = {}
        entropy = self.entropy.sort_values(by=self.entropy.columns[0], ascending=False)
        for pos in entropy.index:
            S = round(float(self.entropy.loc[pos, self.entropy.columns[0]]), 2)
            if S >= self.minS:
                self.motifPos[pos] = S
        if not self.motifPos:
            self.motifPos = {False: 'No Positions Selected'}
        # print(f'Motif Pos:\n{self.motifPos}')


    def predictActivity(self):
        from scipy.stats import pearsonr, spearmanr

        self.calculateRF()

        # Initialize parameters
        self.jobDone = False
        self.log('\n\n========================= Predict Substrate Activity '
                 '=========================')

        # Process substrates
        activityExp = {}
        for idx, sub in enumerate(list(self.subsPred)):
            if ':' in sub:
                seq, score = sub.split(':')
                activityExp[seq] = float(score)
                self.subsPred[idx] = seq
            else:
                break

        # Define columns
        columnLabels = self.xAxisLabel
        if self.xAxisSet:
            columnLabels = self.xAxisSet

        # Evaluate prediction matrix
        predMatrix = pd.DataFrame(0.0, index=self.AA, columns=columnLabels)
        for idx in range(len(self.rfExp.columns)):
            if len(self.rfBg.columns) == 1:
                bgRF = self.rfBg['Average RF']
            else:
                bgRF = self.rfBg.iloc[:, idx]
            predMatrix.iloc[:, idx] = self.rfExp.iloc[:, idx] / bgRF
            predMatrix.iloc[:, idx] = (predMatrix.iloc[:, idx] /
                                       sum(predMatrix.iloc[:, idx]))
        self.log(f'Prediction Matrix:\n{predMatrix}\n')
        self.figures['predMatrix'] = (
            self.plotMatrix(
                data=predMatrix, figLabel='Prediction Matrix', datasetType=self.datasetTag
            )
        )

        # Predict substrate activity
        activityPred = {}
        for substrate in self.subsPred:
            score = 0
            for idx, aa in enumerate(substrate):
                val = predMatrix.loc[aa, predMatrix.columns[idx]]
                if score == 0:
                    score = val
                else:
                    score *= val
            activityPred[substrate] = score

        # Normalize activity scores
        activityExpNorm = self.normValues(activityExp)
        activityPredNorm = self.normValues(activityPred)

        # Evaluate Z-scores
        zExp = self.ZScores(activityExpNorm)
        zPred = self.ZScores(activityPredNorm)

        # Build table
        if zExp:
            columns = ['Exp Activity', 'Z-Score Exp', 'Ranked Exp',
                       'Pred Activity', 'Z-Score Pred', 'Ranked Pred']
        else:
            columns = ['Pred Activity', 'Z-Score Pred', 'Ranked Pred']
        data = pd.DataFrame(0.0, index=self.subsPred, columns=columns)
        for substrate in self.subsPred:
            if zExp:
                data.loc[substrate, 'Exp Activity'] = activityExpNorm[substrate]
                data.loc[substrate, 'Z-Score Exp'] = zExp[substrate]
            data.loc[substrate, 'Pred Activity'] = activityPredNorm[substrate]
            data.loc[substrate, 'Z-Score Pred'] = zPred[substrate]
        if zExp:
            data['Ranked Exp'] = (
                pd.Series(zExp).rank(ascending=False, method='min').astype(int)
            )
        data['Ranked Pred'] = (
            pd.Series(zPred).rank(ascending=False, method='min').astype(int)
        )
        self.log(data)
        
        # Log activities
        self.log(f'\nPredicted Activity:\n'
                 f'* Normalized: {list(data.loc[:, 'Pred Activity'])}\n'
                 f'* Z-Scores:   {list(data.loc[:, 'Z-Score Pred'])}\n')
        self.log(f'Experimental Activity:\n'
                 f'* Experimental: {list(activityExp.values())}\n'
                 f'* Normalized:   {list(data.loc[:, 'Exp Activity'])}\n'
                 f'* Z-Scores:     {list(data.loc[:, 'Z-Score Exp'])}\n')

        # Plot activity scores
        rho = False
        if zExp:
            rho, p = spearmanr(
                list(data.loc[:, 'Z-Score Exp']), list(data.loc[:, 'Z-Score Pred'])
            )
            self.log(f'\nStatistical Analysis:\n'
                     f'* Spearman correlation coefficient (ρ): '
                     f'{round(rho,3)}, p={round(p,3)}\n')
            self.figures['scatterActivity'] = self.plotActivityScatter(data, rho)
        self.figures['barPred'] = self.plotActivityBars(data, rho)


    def plotMatrix(self, data, figLabel, datasetType, totalCounts=False):
        # print(f'Data: {figLabel}\n{data}')

        # Create heatmap
        cMapCustom = self.createCustomColorMap(colorType='green')

        # if 'prediction matrix' in figLabel.lower():
        #     cBarMax = 1
        cBarMax = np.ceil(data.values.max() * 10) / 10
        cBarMin = 0

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
        fig.tight_layout()
        fig.set_size_inches(self.figSize)

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
        if self.useCustomAxis:
            ax.set_xticklabels(self.customAxis[0:len(xTicks)])

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


        # File path
        if totalCounts:
            figLabel = 'Counts'
        else:
            figLabel = figLabel.replace(' ', '')
        if '/' in figLabel:
            figLabel = figLabel.replace('/', '_')
        figName = f'{figLabel}-{self.enzymeName}.png'
        if 'PredictionMatrix' in figName:
            figName = figName.replace('PredictionMatrix', 'PredMatrix')
        if datasetType:
            figName = figName.replace('.png', f'-{datasetType}.png')
        figName = figName.replace(' ', '_').replace('/', '_')
        if totalCounts:
            figName = figName.replace('.png', f'-N_{totalCounts}.png')
        path = os.path.join(self.pathFigs, figName)
        print(f'Saving Figure: {path}')

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName


    def calculateRF(self):
        self.log('\n\n=============================== Calculate: RF '
                 '================================')
        if not self.subsPred:
            self.log(f'Filter: {self.datasetTag}\n')
        self.rfExp = pd.DataFrame(
            0.0, index=self.countsExp.index, columns=self.countsExp.columns
        )
        for pos in self.countsExp.columns:
            self.rfExp.loc[:, pos] = self.countsExp[pos] / sum(self.countsExp[pos])
        self.log(f'RF Experimental:\n{self.rfExp}')

        if self.rfBg is None:
            self.rfBg = pd.DataFrame(self.countsBg.copy(), dtype=float)
            if self.avgBgRF:
                self.rfBg = pd.DataFrame(self.rfBg.mean(axis=1), columns=['Average RF'])
                self.rfBg = self.rfBg / self.rfBg.sum()
            else:
                for pos in self.rfBg.columns:
                    totalCounts = sum(self.rfBg[pos])
                    for AA in self.rfBg.index:
                        count = self.rfBg.loc[AA, pos]
                        if count == 0:
                            count = 1
                        self.rfBg.loc[AA, pos] = count / totalCounts
            self.log(f'\nRF Background:\n{self.rfBg}')


    def evalStats(self):
        self.log('\n\n================================= Statistics '
                 '=================================')
        es = []
        columnLabels = self.xAxisLabel
        if self.xAxisSet:
            columnLabels = self.xAxisSet
        df = pd.DataFrame(0.0, index=self.countsExp.index, columns=columnLabels)

        # Evaluate: Counts
        for idx, counts in enumerate(self.countsExpSets):
            rf, eMap = df.copy(), df.copy()
            self.log(f'Counts: {idx}\n'
                     f'{counts}\n')

            # Calculate: RF
            for pos in counts.columns:
                rf.loc[:, pos] = counts.loc[:, pos] / sum(counts.loc[:, pos])
            self.log(f'RF: {idx}\n{rf}\n')

            # Calculate: Enrichment
            if 'Average RF' in self.rfBg.columns:
                for pos in rf.columns:
                    eMap.loc[:, pos] = (
                        np.log2(rf.loc[:, pos] / self.rfBg.loc[:, 'Average RF'])
                    )
            else:
                for pos in rf.columns:
                    eMap.loc[:, pos] = (
                        np.log2(rf.loc[:, pos] / self.rfBg.loc[:, pos])
                    )
            es.append(eMap)
            self.log(f'Enrichment Map: {idx}\n{eMap}\n\n')

        # Calculate: Standard Deviation
        stDevES = pd.DataFrame(np.std(es, axis=0),
                               index=self.countsExp.index, columns=columnLabels)

        # Plot: Standard deviation
        self.figures['stDev'] = self.plotStats(data=stDevES, dataType='stDev')


    def calculateEntropy(self, plotFig=True):
        self.log('\n\n============================= Calculate: Entropy '
                 '=============================')
        self.log(f'Filter: {self.datasetTag}\n')
        self.entropyMax = np.log2(len(self.rfExp.index))
        for indexColumn in self.rfExp.columns:
            S = 0
            for indexRow, probRatio in self.rfExp.iterrows():
                prob = probRatio[indexColumn]
                if prob == 0:
                    continue
                else:
                    S += -prob * np.log2(prob)
            self.entropy.loc[indexColumn, self.entropy.columns[0]] = self.entropyMax - S
        self.log(f'{self.entropy}\n\nMax Entropy: {self.entropyMax.round(6)}')

        if self.setS:
            self.selectMotifPos()

        if plotFig:
            if self.subProfile:
                self.figures['entropyProfile'] = self.plotEntropy()
            else:
                self.figures['entropy'] = self.plotEntropy()
            if self.setS:
                self.jobDone = True
                self.setS = False


    def plotEntropy(self):
        # Set figure title
        if self.figTag:
            title = f'{self.enzymeName}\n{self.figTag}'
        else:
            # tag = self.datasetTag
            # if 'Exclude' not in tag:
            #     tag = tag.replace('Fix', '')
            # title = f'{self.enzymeName}\n{tag}'
            title = self.enzymeName
        if 'Fix' in title and 'Exclude' in title:
            title = title.replace('Fix', '\nFix')

        # Figure parameters
        yMax = self.entropyMax + 0.2
        xMax = len(self.entropy.iloc[:, 0])

        # Map entropy values to colors using the colormap
        colors = [(0, 'navy'),
                  (0.3 / self.entropyMax, 'navy'),
                  (0.7 / self.entropyMax, 'dodgerblue'),
                  (0.97 / self.entropyMax, 'white'),
                  (0.98 / self.entropyMax, 'white'),
                  (1.0 / self.entropyMax, 'white'),
                  (1.65 / self.entropyMax, 'red'),
                  (3 / self.entropyMax, 'firebrick'),
                  (1, 'darkred')]
        colorBar = LinearSegmentedColormap.from_list('custom_colormap', colors)

        # Map entropy values to colors using the colormap
        normalize = Normalize(vmin=0, vmax=yMax) # Normalize the entropy values
        cMap = [
            colorBar(normalize(value))
            for value in self.entropy[self.entropy.columns[0]].astype(float)
        ]

        # Plotting the entropy values as a bar graph
        fig, ax = plt.subplots(figsize=self.figSize)
        plt.bar(self.entropy.index, self.entropy[self.entropy.columns[0]], color=cMap,
                edgecolor='black', linewidth=self.lineThickness, width=0.8)
        plt.xlabel('Substrate Position', fontsize=self.labelSizeAxis)
        plt.ylabel(self.entropy.columns[0], fontsize=self.labelSizeAxis,
                   rotation=0, labelpad=15)
        plt.title(title, fontsize=self.labelSizeTitle, fontweight='bold')

        # Set x-ticks
        xTicks = np.arange(0, xMax, 1)
        ax.set_xticks(xTicks)
        ax.set_xticklabels(self.entropy.index, rotation=0, ha='center',
                           fontsize=self.labelSizeTicks)
        if self.useCustomAxis:
            ax.set_xticklabels(self.customAxis[0:len(xTicks)], fontsize=self.labelSizeTicks)
        for tick in ax.xaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness) # Set tick width

        # Set y-ticks
        yTicks = range(0, 5)
        yTickLabels = [f'{tick:.0f}' if tick != 0 else f'{int(tick)}' for tick in yTicks]
        ax.set_yticks(yTicks)
        ax.set_yticklabels(yTickLabels, fontsize=self.labelSizeTicks)
        for tick in ax.yaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness) # Set tick width

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks)

        # Set the edge thickness
        for spine in ax.spines.values():
            spine.set_linewidth(self.lineThickness)

        # Set axis limits
        ax.set_xlim(-0.5, xMax-0.5)
        ax.set_ylim(0, yMax)

        # Set colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.1)
        norm = plt.Normalize(vmin=0, vmax=yMax)
        cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=colorBar), cax=cax)
        cbar.set_ticks(list(yTicks))
        cbar.set_ticklabels(yTickLabels)
        cbar.ax.tick_params(axis='y', which='major', labelsize=self.labelSizeTicks,
                            length=self.tickLength, width=self.lineThickness)
        for tick in cbar.ax.yaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness)
        for spine in cbar.ax.spines.values():
            spine.set_linewidth(self.lineThickness)

        plt.tight_layout()
        fig.set_size_inches(self.figSize)

        # File path
        figName = self.getFileNameFig('entropy')
        path = os.path.join(self.pathFigs, figName)
        print(f'Saving ∆S at:\n   {path}')

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName


    def calculateWebLogo(self):
        # Evaluate weblogo
        self.rfExpScaled = pd.DataFrame(
            0.0, index=self.rfExp.index, columns=self.rfExp.columns
        )
        for pos in self.entropy.index:
            self.rfExpScaled.loc[:, pos] = (
                    self.rfExp.loc[:, pos] *
                    self.entropy.loc[pos, self.entropy.columns[0]]
            )
        if self.subProfile:
            self.figures['wLogoProfile'] = self.plotWebLogo()
        else:
            self.figures['wLogo'] = self.plotWebLogo()


    def scaleMatrix(self, data):
        self.log(f'\n\nScale Enrichment Scores:\n'
                 f'     Enrichment Scores * ΔS\n')
        # Calculate: Letter heights
        heights = pd.DataFrame(0, index=data.index,
                               columns=data.columns, dtype=float)
        for idxColumn in heights.columns:
            heights.loc[:, idxColumn] = (
                    data.loc[:, idxColumn] * self.entropy.loc[idxColumn,
                    self.entropy.columns[0]]
            )

        # Calculate: Max positive
        columnTotals = []
        for pos in heights.columns:
            totalPos = 0
            for value in heights.loc[:, pos]:
                if value > 0:
                    totalPos += value
            columnTotals.append(totalPos)
        yMax = max(columnTotals)

        # Adjust values
        for column in heights.columns:
            if (heights[column] == 0).all():
                heights.loc[:, column] = -np.inf
                print(f'NaN: {column}:\n{heights[column]}')
            elif heights.loc[:, column].isna().any():
                nValues = heights[column].notna().sum()
                if nValues > 0:
                    self.log(
                        f'{len(heights[column]) - nValues} NaN values at: {column}')
                heights.loc[heights[column].notna(), column] = yMax / nValues
        if self.datasetTagMotif:
            self.log(f'Residue Heights: {self.datasetTagMotif}\n'
                     f'{heights}\n')
        else:
            self.log(f'Residue Heights: {self.datasetTag}\n'
                     f'{heights}\n')
        return heights


    def evalEnrichment(self, releasedCounts=False, skipFigs=False, noFigs=False):
        def evalMatrix(data):
            stacks = pd.DataFrame(0.0, index=data.columns,
                                  columns=['+Stack', '-Stack'])
            for pos in data.columns:
                totalPos = 0
                totalNeg = 0
                for value in data.loc[:, pos]:
                    if value > 0:
                        totalPos += value
                    if value < 0:
                        totalNeg += value
                stacks.loc[pos, '+Stack'] = totalPos
                stacks.loc[pos, '-Stack'] = totalNeg
            self.log(f'Stack Heights:\n{stacks}')

        matrix = pd.DataFrame(0.0, index=self.rfExp.index,
                              columns=self.rfExp.columns)
        if releasedCounts:
            self.log('\n\n============================= Substrate Profile '
                     '==============================')
        else:
            self.log('\n\n======================== Calculate: Enrichment Score '
                     '=========================')
        self.log(f'Enrichment Scores:\n'
                 f'     log₂(RF Experimental / RF Background)\n')

        # Calculate: Enrichment scores
        if len(self.rfBg.columns) == 1:
            # Eval: ES
            for pos in self.rfExp.columns:
                for AA in self.rfExp.index:
                    rf = self.rfExp.loc[AA, pos]
                    if rf == 0:
                        matrix.loc[AA, pos] = -np.inf
                    else:
                        rfBg = self.rfBg.loc[AA, self.rfBg.columns[0]]
                        matrix.loc[AA, pos] = np.log2(rf / rfBg)
        else:
            if len(self.rfBg.columns) != len(self.rfExp.columns):
                self.log(f'ERROR: The number of columns in the Initial Sort '
                      f'({len(self.rfBg.columns)}) needs to equal to the '
                      f'number of columns in the Final Sort '
                      f'({len(self.rfExp.columns)})\n'
                      f'     Initial: {self.rfBg.columns}\n'
                      f'       Final: {self.rfExp.columns}\n\n')
                sys.exit(1)

            # Eval: ES
            for pos in self.rfExp.columns:
                for AA in self.rfExp.index:
                    rf = self.rfExp.loc[AA, pos]
                    if rf == 0:
                        matrix.loc[AA, pos] = -np.inf
                    else:
                        matrix.loc[AA, pos] = np.log2(rf / self.rfBg.loc[AA, pos])
        self.eMap = matrix.copy()
        self.log(f'Enrichment Score: {self.datasetTag}\n'
              f'{matrix.round(self.roundVal)}\n')

        # Evaluate stack heights
        evalMatrix(matrix.replace([np.inf, -np.inf], 0))
        self.eMapScaled = self.scaleMatrix(self.eMap)
        evalMatrix(self.eMapScaled)
        if noFigs:
            return

        # Plot: Enrichment Map
        self.figures['eMap'] = (
            self.plotEnrichmentScores(dataType='Enrichment')
        )
        self.figures['eMapSc'] = (
            self.plotEnrichmentScores(dataType='Scaled Enrichment')
        )
        
        # Plot: Enrichment Logo
        self.plotEnrichmentLogo()

        if not skipFigs:
            # Plot: WebLogo
            self.calculateWebLogo()

            # Plot: WordCloud
            if self.subsExp:
                self.figures['words'] = self.plotWordCloud(self.subsExp)


    def plotEnrichmentScores(self, dataType):
        # Select: Dataset
        scaleData = False
        if 'scaled' in dataType.lower():
            scaleData = True
            scores = self.eMapScaled
        else:
            scores = self.eMap

        # Set figure title
        if self.figTag:
            title = f'{self.enzymeName}\n{self.figTag}'
        else:
            title = self.enzymeName
        if 'Fix' in title and 'Exclude' in title:
            title = title.replace('Fix', '\nFix')

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
            cBarMax = np.ceil(np.max(scores) * 10) / 10
            cBarMin = -1 * cBarMax
        else:
            cBarMin = np.floor(np.min(scores) * 10) / 10
            cBarMax = -1 * cBarMin

        # Plot the heatmap with numbers centered inside the squares
        if self.figEMSquares:
            fig, ax = plt.subplots(figsize=self.figSizeSq)
            sns.heatmap(
                scores, annot=False, cmap=cMapCustom, cbar=False,
                linewidths=self.lineThickness - 1, linecolor='black',
                square=True, center=None, vmax=cBarMax, vmin=cBarMin,
                cbar_kws={'pad': 0.02}
            )
        else:
            fig, ax = plt.subplots(figsize=self.figSize)
            sns.heatmap(
                scores, annot=True, fmt='.3f', cmap=cMapCustom, cbar=False,
                linewidths=self.lineThickness - 1, linecolor='black',
                square=False, center=None, vmax=cBarMax, vmin=cBarMin,
                annot_kws={'fontweight': 'bold', 'size': self.labelSizeEM},
                cbar_kws={'pad': 0.02}
            )
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
        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks)

        # Set x-ticks
        xTicks = np.arange(len(scores.columns)) + 0.5
        ax.set_xticks(xTicks)
        ax.set_xticklabels(scores.columns)
        if self.useCustomAxis:
            ax.set_xticklabels(self.customAxis[0:len(xTicks)])

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
        cmap.set_bad(color='lightgrey')  # Set invalid values to gray

        # Modify the colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.1)
        norm = plt.Normalize(vmin=cBarMin, vmax=cBarMax)
        cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cMapCustom), cax=cax)
        cbar.ax.tick_params(axis='y', which='major', labelsize=self.labelSizeTicks,
                            length=self.tickLength, width=self.lineThickness)
        cbar.outline.set_linewidth(self.lineThickness)
        cbar.outline.set_edgecolor('black')

        plt.tight_layout()
        fig.set_size_inches(self.figSize)

        # File path
        figName = self.getFileNameFig('eMap')
        if scaleData:
            figName = figName.replace('eMap', 'eMap_Scaled')
        path = os.path.join(self.pathFigs, figName)
        print(f'Saving EM at:\n   {path}')

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName


    def plotStats(self, data, dataType):
        self.log(f'Standard deviation:\n{data}\n\n')

        title = f'{self.enzymeName}\nStandard Deviation'

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
        cBarMax = np.ceil(np.max(data) * 10) / 10
        cBarMin = 0

        # Plot the heatmap with numbers centered inside the squares
        fig, ax = plt.subplots(figsize=self.figSize)
        sns.heatmap(
            data, annot=True, fmt='.3f', cmap=cMapCustom, cbar=False,
            linewidths=self.lineThickness - 1, linecolor='black',
            square=False, center=None, vmax=cBarMax, vmin=cBarMin,
            annot_kws={'fontweight': 'bold', 'size': self.labelSizeEM},
            cbar_kws={'pad': 0.02}
        )
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
        xTicks = np.arange(len(data.columns)) + 0.5
        ax.set_xticks(xTicks)
        ax.set_xticklabels(data.columns)
        if self.useCustomAxis:
            ax.set_xticklabels(self.customAxis[0:len(xTicks)])

        # Set y-ticks
        yTicks = np.arange(len(data.index)) + 0.5
        ax.set_yticks(yTicks)
        ax.set_yticklabels(data.index)

        # Colormap
        if isinstance(cMapCustom, Colormap):
            cmap = cMapCustom
        else:
            # Look up premade colormap
            cmap = matplotlib.colormaps[cMapCustom]
        cmap.set_bad(color='lightgrey')  # Set invalid values to gray

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

        plt.tight_layout()
        fig.set_size_inches(self.figSize)

        # File path
        figName = self.getFileNameFig(dataType)
        path = os.path.join(self.pathFigs, figName)
        print(f'Saving Standard Deviation at:\n   {path}')

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName


    def plotEnrichmentLogo(self):
        # Set figure title
        if self.figTag:
            title = f'{self.enzymeName}\n{self.figTag}'
        else:
            title = self.enzymeName
        if 'Fix' in title and 'Exclude' in title:
            title = title.replace('Fix', '\nFix')

        # Set parameters
        if self.bigAAonTop:
            stackOrder = 'big_on_top'
        else:
            stackOrder = 'small_on_top'

        # Rename columns for logomaker script
        # print(f'E Map Scaled:\n{self.eMapScaled}\n')
        scores = self.eMapScaled.copy().replace([np.inf, -np.inf], 0).replace(np.nan, 0)
        # print(f'Logo:\n{scores}')
        xTicks = scores.columns
        scores.columns = range(len(scores.columns))

        # Calculate: Max and min
        columnTotals = [[], []]
        for indexColumn in scores.columns:
            totalPos = 0
            totalNeg = 0
            for value in scores.loc[:, indexColumn]:
                if value > 0:
                    totalPos += value
                elif value < 0:
                    totalNeg += value
            columnTotals[0].append(totalPos)
            columnTotals[1].append(totalNeg)
        yMax = max(columnTotals[0])
        yMin = min(columnTotals[1])


        def plotLogo(matrix, limitYAxis=False):
            # Plot the sequence motif
            fig, ax = plt.subplots(figsize=self.figSize)
            motif = logomaker.Logo(matrix.transpose(), ax=ax, color_scheme=self.colorsAA,
                                   width=0.95, stack_order=stackOrder)
            motif.ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')

            # Set borders
            motif.style_spines(visible=False)
            motif.style_spines(spines=['left', 'bottom'], visible=True)
            for spine in motif.ax.spines.values():
                spine.set_linewidth(self.lineThickness)

            # Set x-ticks
            motif.ax.set_xticks([pos for pos in range(len(xTicks))])
            motif.ax.set_xticklabels(xTicks)
            if self.useCustomAxis:
                ax.set_xticklabels(self.customAxis[0:len(xTicks)])

            # Set y-ticks
            yTicks = [yMin, 0, yMax]
            yTickLabels = [f'{tick:.2f}' if tick != 0 else f'{int(tick)}'
                           for tick in yTicks]
            motif.ax.set_yticks(yTicks)
            motif.ax.set_yticklabels(yTickLabels, fontsize=self.labelSizeTicks)
            motif.ax.set_ylim(yMin, yMax)

            # Set tick parameters
            ax.tick_params(axis='both', which='major', length=self.tickLength,
                           labelsize=self.labelSizeTicks)

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

            # Scale fig
            plt.tight_layout()
            fig.set_size_inches(self.figSize)

            # File path
            figName = self.getFileNameFig('eLogo')
            if limitYAxis:
                figName = figName.replace('eLogo', 'eLogo_yMin')
            path = os.path.join(self.pathFigs, figName)
            print(f'Saving Enrichment Logo:\n     {path}')

            # Encode the figure
            figBase64 = self.encodeFig(fig)
            with open(path, "wb") as file:
                file.write(base64.b64decode(figBase64))

            # Close the figure to free memory
            plt.close(fig)

            return figName

        # Plot figure
        if self.subProfile: # Full y-axis
            self.figures['eLogoProfile'] = plotLogo(scores)
        else:
            self.figures['eLogo'] = plotLogo(scores)

        # Adjust yMin to fit the largest negative AA
        yMin = 0
        for col in scores.columns:
            for row in scores.index:
                if scores.loc[row, col] < yMin:
                    yMin = scores.loc[row, col]
        if self.subProfile: # Limited y-axis
            self.figures['eLogoMinProfile'] = plotLogo(scores, limitYAxis=True)
        else:
            self.figures['eLogoMin'] = plotLogo(scores, limitYAxis=True)


    def plotWordCloud(self, substrates):
        # Limit the number of words
        subs = {}
        iteration = 0
        for substrate, count in substrates.items():
            subs[substrate] = count
            iteration += 1
            if iteration >= self.numSamples:
                break
        substrates = subs
        totalWords = len(substrates)
        if not self.motifFilter:
            self.figTag = f'Top {totalWords} Substrates'

        # Set figure title
        if self.figTag:
            title = f'{self.enzymeName}\n{self.figTag}'
        else:
            title = self.enzymeName
        if 'Fix' in title and 'Exclude' in title:
            title = title.replace('Fix', '\nFix')


        # Create word cloud
        cmap = self.createCustomColorMap(colorType='WordCloud')
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

        # Scale fig
        plt.tight_layout()
        fig.set_size_inches(self.figSize)

        # File path
        figName = self.getFileNameFig('wordcloud', f'-{totalWords}_Words')
        path = os.path.join(self.pathFigs, figName)
        print(f'Saving WordCloud:\n     {path}')

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName


    def plotWebLogo(self):
        # Set figure title
        if self.figTag:
            title = f'{self.enzymeName}\n{self.figTag}'
        else:
            title = self.enzymeName
        if 'Fix' in title and 'Exclude' in title:
            title = title.replace('Fix', '\nFix')

        # Set parameters
        if self.bigAAonTop:
            stackOrder = 'big_on_top'
        else:
            stackOrder = 'small_on_top'

        # Rename columns for logomaker script
        data = self.rfExpScaled.copy().replace([np.inf, -np.inf], 0)
        xTicks = data.columns
        data.columns = range(len(data.columns))

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
        motif.ax.set_xticks([pos for pos in range(len(xTicks))])
        motif.ax.set_xticklabels(xTicks, fontsize=self.labelSizeTicks,
                                 rotation=0, ha='center')
        if self.useCustomAxis:
            ax.set_xticklabels(self.customAxis[0:len(xTicks)])

        # Set y-ticks
        yMax = self.entropyMax
        yTicks = range(0, 5)
        yTickLabels = [f'{tick:.0f}' if tick != yMax else f'{yMax:.2f}' for tick in
                       yTicks]
        # yTicks.append(4.32)
        # yTickLabels.append('')
        motif.ax.set_yticks(yTicks)
        motif.ax.set_yticklabels(yTickLabels, fontsize=self.labelSizeTicks)
        motif.ax.set_ylim(0, yMax)

        # Set tick width
        for tick in motif.ax.xaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness)
        for tick in motif.ax.yaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness)

        # Label the axes
        motif.ax.set_xlabel('Substrate Position', fontsize=self.labelSizeAxis)
        motif.ax.set_ylabel('Bits', fontsize=self.labelSizeAxis)

        # Set horizontal line
        motif.ax.axhline(y=0, color='black', linestyle='-',
                         linewidth=self.lineThickness)

        # Scale fig
        plt.tight_layout()
        fig.set_size_inches(self.figSize)

        # File path
        figName = self.getFileNameFig('webLogo')
        path = os.path.join(self.pathFigs, figName)
        print(f'Saving WebLogo:\n     {path}')

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName


    def plotBars(self, data, dataType, plotAll=False, barColor='#BF5700', barWidth=0.75):
        if not data:
            return
        x, y, totalCounts, limitNSubs = [], [], sum(data.values()), self.numSamples
        if plotAll:
            barWidth, limitNSubs = 1.5, len(data.keys())
            self.figTag = f'{limitNSubs:,} Substrates'
            for i, count in enumerate(data.values()):
                x.append(i), y.append(count)
        else:
            if not self.motifFilter:
                self.figTag = f'Top {limitNSubs:,} Substrates'
            for i, (substrate, count) in enumerate(data.items()):
                if i >= limitNSubs:
                    break
                x.append(substrate), y.append(count)

        # Evaluate data
        if 'counts' in dataType.lower():
            # Evaluate: Y axis
            maxValue, yMin, mag = max(y), 0, 10
            magnitude = math.floor(math.log10(maxValue))
            if magnitude > 1:
                mag = 10 ** (magnitude - 1)
            yMax = math.ceil(maxValue / mag) * mag
        elif 'relative' in dataType.lower():
            y = [v / totalCounts for v in y]

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
            print(f'ERROR: What data type is: {dataType}\n\n')
            return None
        NSubs = len(x)

        # Set figure title
        if self.figTag:
            title = f'{self.enzymeName}\n{self.figTag}'
        else:
            title = self.enzymeName
        if 'Fix' in title and 'Exclude' in title:
            title = title.replace('Fix', '\nFix')

        # Plot the data
        fig, ax = plt.subplots(figsize=self.figSize)
        bars = plt.bar(x, y, color=barColor, width=barWidth)
        plt.ylabel(dataType, fontsize=self.labelSizeAxis)
        plt.title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        plt.axhline(y=0, color='black', linewidth=self.lineThickness)
        plt.ylim(yMin, yMax)

        # Set: x ticks
        if plotAll:
            magnitude = math.floor(math.log10(NSubs))
            step = (10 ** magnitude) / 2
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
                vals = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55]
                for val in vals:
                    ratio = dist / val
                    if ratio == int(ratio) and ratio < 10:
                        plotYTicks = True
                        yTicks = np.arange(yMin, yMax + val, val)
                        break
            yMax += 0.1
        else:
            yTicks = np.linspace(yMin, yMax, 11)
        plt.ylim(yMin, yMax)
        if plotYTicks:
            plt.yticks(yTicks)

        # Set the edge color
        if not plotAll:
            for bar in bars:
                bar.set_edgecolor('black')
                bar.set_linewidth(self.lineThickness)

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks, width=self.lineThickness)
        if not plotAll:
            ax.tick_params(axis='x', which='major', labelsize=12)

        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(self.lineThickness)

        # Scale fig
        plt.tight_layout()
        fig.set_size_inches(self.figSize)

        # File path
        figName = self.getFileNameFig(
            f'bars-{dataType.replace(' ', '_')}-NBars_{len(xTicks)}'
        )
        path = os.path.join(self.pathFigs, figName)

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName


    def plotActivityScatter(self, data, rho=False, colorExp='#BF5700'):
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
            return xFit, fnExp(xFit, *popt), r2, popt

        def fnLinear(x, m, b):
            """Linear function: y = mx + b"""
            return m * x + b

        def fitDataLinear(x, y):
            # Fit the linear curve
            # p0=[1, 0] suggests initial slope=1, intercept=0
            popt, pcov = curve_fit(fnLinear, x, y, p0=[1, 0], maxfev=10000)
            # m, b = popt # y = mx + b

            # Generate smooth line for plotting
            xFit = np.linspace(min(x), max(x), 300)
            yFit = fnLinear(xFit, *popt)

            # R² calculation
            yPred = fnLinear(x, *popt)
            ss_res = np.sum((y - yPred) ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r2 = 1 - (ss_res / ss_tot)

            return xFit, yFit, r2, popt

        # Evaluate data
        x = list(data.loc[:, 'Z-Score Exp'])
        y = list(data.loc[:, 'Z-Score Pred'])
        xfit_exp, yFit_exp, r2_exp, popt_exp = fitData(x=np.array(x), y=np.array(y))
        xfit_lin, yFit_lin, r2_lin, popt_lin = fitDataLinear(x=np.array(x), y=np.array(y))
        if r2_exp >= r2_lin:
            xFit, yFit, r2 = xfit_exp, yFit_exp, r2_exp
            a, b, c = popt_exp # y = a · e^(b·x) + c
            self.log(f'Fitting data to an exponential curve:\n'
                     f'  y = {round(a,2)} * e^({round(b,2)}·x) + {round(c,2)}\n')
        else:
            xFit, yFit, r2 = xfit_lin, yFit_lin, r2_lin
            m, b = popt_lin # y = mx + b
            self.log(f'Fitting data to a linear equation:\n'
                     f'  y = {round(m,2)} * x + {round(b,2)}\n')

        
        # Set figure title
        title = f'{self.enzymeName}\nSubstrate Activity'

        # Make figure
        fig, ax = plt.subplots(figsize=self.figSize)
        plt.scatter(x, y, color=colorExp, edgecolor='black')
        ax.plot(xFit, yFit, color='black', linestyle='-',
                linewidth=self.lineThickness)
        plt.xlabel('Experimental Activity', fontsize=self.labelSizeAxis)
        plt.ylabel('Predicted Activity', fontsize=self.labelSizeAxis)
        plt.title(title, fontsize=self.labelSizeTitle, fontweight='bold')

        # Axis
        spacer = 0.2
        xMax = self.roundup(max(x) + spacer)
        xMin = self.roundup(min(x) - spacer, upperLim=False)
        yMax = self.roundup(max(y) + spacer)
        yMin = self.roundup(min(y) - spacer, upperLim=False)
        plt.xlim(xMin, xMax)
        plt.ylim(yMin, yMax)

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks, width=self.lineThickness)

        # Set the edge thickness
        for tick in ax.xaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness)  # Set tick width
        for tick in ax.yaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness)  # Set tick width
        for spine in ax.spines.values():
            spine.set_linewidth(self.lineThickness)

        # Legend
        invisibleHandle = Line2D([], [], linestyle='None',
                                 marker='None',color='none')
        ax.legend(
            handles=[invisibleHandle],
            labels=[f'R² = {round(r2, 3)}\nSpearman ρ: {round(rho, 3)}'],
            prop=FontProperties(size=self.labelSizeTicks - 2, weight='bold'),
            handlelength=0, handletextpad=0, edgecolor='black',
            linewidth=self.lineThickness, loc='upper left', framealpha=0.9
        )

        # Scale fig
        plt.tight_layout()
        fig.set_size_inches(self.figSize)

        # File path
        figName = self.getFileNameFig(f'PredActivity-Scatter-NSubs_{len(data.index)}')
        path = os.path.join(self.pathFigs, figName)
        print(f'Saving activity prediction scatterplot:\n     {path}')

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName


    def plotActivityBars(self, data, rho=False, colorExp='#F8971F', colorPred='#BF5700'):
        # Evaluate data
        actPred, NSubs = list(data['Z-Score Pred']), len(data.index)
        actExp = [0 for _ in range(NSubs)]
        title = f'{self.enzymeName}\nPredicted Activity' # Define: Figure title
        dualSets = False
        if 'Z-Score Exp' in data.columns:
            dualSets = True
            title = (f'{self.enzymeName}\nSubstrate Activity\n'
                     f'Spearman ρ: {round(rho, 3)}')
            actExp = list(data['Z-Score Exp'])
        print(f'Plotting {NSubs:,} substrates')
        xTicks = np.arange(0, NSubs)

        # Plot the data
        fig, ax = plt.subplots(figsize=self.figSize)
        spacer = 1 - 0.2
        if dualSets:
            barWidth = spacer / 2
            ax.bar(xTicks - barWidth / 2, actExp, barWidth, color=colorExp,
                   edgecolor='black', linewidth=self.lineThickness)
            ax.bar(xTicks + barWidth / 2, actPred, barWidth, color=colorPred,
                   edgecolor='black', linewidth=self.lineThickness)
        else:
            barWidth = spacer / 2
            plt.bar(xTicks, actPred, color=colorPred, width=barWidth,
                    edgecolor='black', linewidth=self.lineThickness)
        plt.ylabel('Z-Score', fontsize=self.labelSizeAxis)
        plt.title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        plt.axhline(y=0, color='black', linewidth=self.lineThickness)


        # Set: x ticks
        ax.set_xticks(xTicks)
        ax.set_xticklabels(data.index, rotation=0, ha='center')
        ax.set_xlim(left=xTicks[0] - (barWidth*1.5), right=xTicks[-1] + (barWidth*1.5))
        plt.xticks(rotation=90, ha='center')

        # Set: y ticks
        plotYTicks = True
        yTicks = []
        maxScore = max(max(actPred), max(actExp))
        minScore = min(min(actPred), min(actExp))
        if maxScore == 1.0:
            yMax = 1.0
            yMin = 0
            if minScore >= 0:
                yTicks = np.linspace(yMin, yMax, 6)
            else:
                plotYTicks = False
                dist = yMax - yMin
                vals = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55]
                for val in vals:
                    ratio = dist / val
                    if ratio == int(ratio) and ratio < 10:
                        plotYTicks = True
                        yTicks = np.arange(yMin, yMax + val, val)
                        break
            yMax += 0.1
        else:
            step, spacer = maxScore / 10, 0.2
            yMax = self.roundup(maxScore + spacer)
            yMin = self.roundup(minScore - spacer, upperLim=False)
            yTicks = np.linspace(yMin, yMax, 10)
        plt.ylim(yMin, yMax)
        if plotYTicks:
            plt.yticks(yTicks)
        # print(f'Y Axis:\n'
        #       f'* Max: {yMax:,}\n'
        #       f'* Min: {yMin:,}')

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks, width=self.lineThickness)

        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(self.lineThickness)

        # Legend
        if dualSets:
            colors = [colorExp, colorPred]
            handles = [Line2D([], [],
                              linestyle='None', marker='o', color='none',
                              markerfacecolor=color, markeredgewidth=self.lineThickness,
                              markersize=8) for color in colors]
            ax.legend(
                handles=handles, labels=['Experimental Activity', 'Predicted Activity'],
                prop=FontProperties(size=self.labelSizeTicks - 2, weight='bold'),
                handlelength=1, handletextpad=0.2, edgecolor='black',
                linewidth=self.lineThickness, loc='best', framealpha=0.9
            )

        # Scale fig
        plt.tight_layout()
        fig.set_size_inches(self.figSize)

        # File path
        figName = self.getFileNameFig(f'PredActivity-Bars-NSubs_{len(xTicks)}')
        path = os.path.join(self.pathFigs, figName)
        print(f'Saving activity prediction bar graph:\n     {path}')

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName
