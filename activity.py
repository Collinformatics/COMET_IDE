from functions import pressKey
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import pearsonr, spearmanr
import sys


# Set options
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)
pd.set_option('display.float_format', '{:,.3f}'.format)



# ========================================================================================
# Input: Data
inEnzyme = f'M{"ᵖʳᵒ"}2'
inEnzyme2 = f'M{"ᵖʳᵒ"}' # Secondary enzyme
inSigFigs = 0
inRoundVal = 3
inNatLog = False
inTableCols = ['% Product','Activity Z','Activity Rank','Predicted Z','Predicted Rank']
inSaveTables = True

# Input: Figure Params
inFigSize = (9.5, 8)
inTickLength = 4
inLinewidth = 1.5
inTitleSize = 20
inLabelSize = 18
inLabelTickSize = 16

# Input: Datasets
inSubstrates = ['AVLQSGFR', 'VILQSGFR', 'VILQTGFR', 'VILQSPFR',
                'VILHSGFR', 'VIMQSGFR', 'VPLQSGFR', 'NILQSGFR']
inExpActivity = [46.1, 49.5, 14.5, 0, 13.1, 37.0, 0, 16.1]
inExpActivity2 = [27.8, 34.8, 12.7, 0, 13.8, 31.3, 0, 15.2]
inExpActivity2 = [32.1, 39.1, 14.9, 0, 16.0, 36.5, 0, 15.6]
# inPredActivity = [0.595, 1.0, 0.008, 0.004, 0.055, 0.417, 0.003, 0.049] # 6 AA
# inPredActivity2 = [0.748, 1.0, 0.007, 0.009, 0.03, 0.453, 0.005, 0.027] # 6 AA
inPredActivity = [
    1.0, 0.99975, 0.83602, 0.75165, 0.83622, 0.93213, 0.75742, 0.80994
]# 8 AA
inPredActivity2 = [
    0.97208, 1.0, 0.75065, 0.65427, 0.80786, 0.9474, 0.68811, 0.82464
] # 8 AA
inStDev = [4.8, 2.6, 2.5, 0, 4.4, 3.2, 0, 2.9]
inStDev2 = [1.2, 5.5, 2.7, 0, 1.8, 4.7, 0, 1.9]
inSubstratesNat = ['AVLQSGFR', 'VTFQSAVK', 'ATVQSKMS', 'ATLQAIAS',
                   'VKLQNNEL', 'VRLQAGNA', 'PMLQSADA', 'TVLQAVGA',
                   'ATLQAENV', 'TRLQSLEN', 'PKLQSSQA']
inExpActivityNat = [53, 24, 2, 19,
                    3, 41, 17, 16,
                    10, 11, 21]
inPredActivityNat = [
    1.0, 0.86857, 0.76313, 0.84762, 0.84903,
    0.91229, 0.92121, 0.86465, 0.92956, 0.86279, 0.8557
]
inStDevNat = [2.6, 0.5, 0.5, 0.4, 0.2, 2.0, 1.2, 2.0, 1.0, 1.6, 0.6]
inDatasets = [ # Enzyme name, Substrates, Exp Activity, Exp StDev, Predicted Activity
    (f'M{"ᵖʳᵒ"}2', inSubstrates, inExpActivity, inStDev, inPredActivity),
    (f'M{"ᵖʳᵒ"}', inSubstrates, inExpActivity2, inStDev2, inPredActivity2),
    # (f'M{"ᵖʳᵒ"}2 pp1a/b', inSubstratesNat, inExpActivityNat, inStDevNat, inPredActivityNat)
] # Plot this data

# Input: Figures
inSaveTag = f'synthSubs_Mpros-MinES_neg0_5' # Add label to saved data
inPlotBarGraph = False
inPlotTable = True
inSavePath = 'Data/Figures/'
inFigTitle = f'Enzyme Activity'
inColor1 = '#BF5700'
inColor2 = '#F8971F'
inFigResolution = 600
inPlotNDatasets = 3 # Max num of plotted datasets
inPlotColors = ['#BF5700', '#101010', '#F8971F', '#2E9418',]
inPlotMarkers = ['D', 'o', '^', ]

# Input: Custom data
inPlotCleaveNet = False
inSubstratesCN, inPredActivityCN, inExpPredActivityCN, inStDevCN = [], [], [], []
if inPlotCleaveNet:
    inSubstrates = ['AVLQSGFR', 'VILQSGFR',
                    'VIMQSGFR',]
    inExpActivity = [
        1.588714482717045, 1.5637731129492973,
        0.5176798411319979
    ]
    inPredActivity = [
        2.2341681287344204, 2.133958673606005,
        -0.10741899173024665,
    ]
    inSubstratesCN = [
        'ARLQALFR', 'SRLQAGFR', 'VLLQALFR', 'VMLQALFR',
        'VNLQGLFR', 'VLLQGNFR', 'VALQCAFR',
        'IKLQCNFR', 'ILLQGFFR', 'ALLQGNFR',
        'VSLQAIFR', 'VLLQGYFR', 'VLLQSQFR', 'IWLQGSFR',
        'CLLQCVFR', 'VNLQSLFR', 'VLLQGRFR', 'CLLQALFR'
    ] # CleaveNet substrates
    inExpPredActivityCN = [
        0.9597151600034503, 0.3650697573358019, 0.9644517350670326, 1.0797913778321793,
        0.09705266047539371, 1.5495228973789, 0.7617482731121744,
        0.614607080696203, -0.39818537162702144, 1.4950726387522746,
        -0.013071362869361664, 0.08246156445857096, 0.9069049348566184, -0.5437405800340406,
        -0.1697590067954982, -0.2929163960274127, 0.8290939641364768, 0.8497066297693364
    ]
    inPredActivityCN = [
        0.48883328431314166, -0.2447286854435148, 0.49729298509131736, 0.7204989565593584,
        -0.4280625593065964, 2.078167824900689, 0.178651974732913,
        -0.0046997308775174745, -0.6353448494675541, 1.8743995174115495,
        -0.4865003551385704, -0.43630217773900515, 0.39800922050803417, -0.6746151958840482,
        -0.5562851208566394, -0.6018549567838597, 0.2752126683745324, 0.30652439474986753
    ]
    inStDev = [0 for _ in range(len(inSubstrates))]
    inStDevCN = [0 for _ in range(len(inSubstratesCN))]
    inDatasets = [ # Enzyme name, Substrates, Exp Activity, Exp StDev, Predicted Activity
        (f'M{"ᵖʳᵒ"}2', inSubstrates, inExpActivity, inStDev, inPredActivity),
        # (f'M{"ᵖʳᵒ"}', inSubstrates, inExpActivity2, inStDev2, inPredActivity2),
        # (f'M{"ᵖʳᵒ"}2 pp1a/b', inSubstratesNat, inExpActivityNat, inStDevNat, inPredActivityNat)
        (f'M{"ᵖʳᵒ"}2 CleaveNet', inSubstratesCN, inExpPredActivityCN, inStDevCN, inPredActivityCN)
    ] # Plot this data



# ========================================================================================
"""
    This takes the data given to "inDatasets" an builds the "data" dictionary.
    If you change a label, make sure it matches the labels in "dataTags" and "dataZTags".
"""
data = {}
dataTags = ['% Product', 'Predicted']
dataZTags = ['Activity Z', 'Predicted Z']
inNormalizeSub = '' # Optional: Normalize activity to a substrate
for idx, enzyme in enumerate(inDatasets, start=1):
    if idx <= inPlotNDatasets:
        emptyList = [0 for _ in range(len(enzyme[1]))]
        # print(enzyme)
        data[f'{enzyme[0]}'] = {
            'Substrates': enzyme[1],
            f'% Product': enzyme[2],
            f'% St Dev': enzyme[3],
            f'Activity Z': emptyList,
            f'Activity Rank': emptyList,
            f'Predicted': enzyme[4],
            f'Predicted Z': emptyList,
            f'Predicted Rank': emptyList,
        }


# ========================================================================================
if not os.path.exists(inSavePath):
    os.makedirs(inSavePath, exist_ok=True)


def normalizeData(data, tags):
    for enzyme in data.keys():
        for tag in tags:
            values = data[enzyme][tag]
            maxValue = max(values)
            if inNormalizeSub:
                idxMaxValue = data[enzyme]['Substrates'].index(inNormalizeSub)
                maxValue = values[idxMaxValue]
            data[enzyme][tag] = [v / maxValue for v in values]
    return data


def convertNum(data, key):
    keyStDev = '% St Dev'
    if key != '% St Dev':
        keyStDev = f'% St Dev{key.replace("% Product", "")}'

    activity = [] # Convert to %
    for idx, substrate in enumerate(data['Substrates']):
        act = data[key][idx] * 100
        stdev = int(data[keyStDev][idx] * 100)
        if inSigFigs == 0 or not inSigFigs:
            act = int(act)
        else:
            act = round(act, inSigFigs)
        if stdev == 0.0:
            activity.append(act)
        else:
            activity.append(f'{act} ± {stdev}')
    return activity


def zScore(data, tags, zTags):
    for enzyme in data.keys():
        # print(f'{enzyme}:')
        for i, tag in enumerate(tags):
            # print(f'* {tag}')
            values = data[enzyme][tag]
            avg = np.mean(values)
            stdev = np.std(values)
            z = []
            for x in values:
                z.append(float((x - avg) / stdev))
            zRank = pd.Series(z).rank(ascending=False, method='min')
            zRank = [int(x) for x in zRank]
            data[enzyme][f'{zTags[i]}'] = z
            data[enzyme][f'{zTags[i].replace(' Z', ' Rank')}'] = zRank
            # print(f'  * {values}\n'
            #       f'  * {z}\n')

    return data


def fnLinear(x, m, b):
    return m * x + b


def fitData(x, y, printOutputs=True):
    if printOutputs:
        print(f'Fitting data to an linear equation:')
    popt, pcov = curve_fit(fnLinear, x, y, maxfev=10000)
    xFit = np.linspace(min(x), max(x), 300)
    yFit = fnLinear(xFit, *popt)
    yPred = fnLinear(np.array(x), *popt)
    ss_res = np.sum((np.array(y) - yPred) ** 2)
    ss_tot = np.sum((np.array(y) - np.mean(y)) ** 2)
    r2 = 1 - (ss_res / ss_tot)
    if printOutputs:
        dec = 3
        m, b = popt # y = m · x + b
        if m > 10**-3:
            m = f'{m:,.{dec}f}'
        else:
            m = f'{m:,.{dec}e}'
        if b > 10**-3:
            b = f'{b:,.{dec}f}'
        else:
            b = f'{b:,.{dec}e}'
        print(f'\ty = {m} * x + {b}\n'
              f'\tR² = {round(r2, dec)}\n\n')
    return xFit, yFit, r2


def plotBars(dataset, barWidth=0.35):
    """
        Plot experimental and predicted activity of multiple enzymes.
    """
    columns = ['% Product', 'Predicted']
    for enzyme in dataset.keys():
        # Get data
        substrates, yValues = [], []
        if not yValues:
            substrates = dataset[list(dataset.keys())[0]]['Substrates']
        yValues.append(list(dataset[enzyme][columns[0]]))
        yValues.append(list(dataset[enzyme][columns[1]]))
        xTicks = np.arange(len(substrates))
        d = pd.DataFrame(0.0, index=substrates, columns=columns)
        for i in range(len(yValues)):
            d.loc[substrates, d.columns[i]] = yValues[i]

        # Spearman rank correlation
        rho, p = spearmanr(yValues[0], yValues[1])

        # Labels
        l1 = 'Experimental Activity'
        l2 = 'Predicted Activity'
        if f'M{"ᵖʳᵒ"}' in enzyme:
            enzyme = f'SARS-CoV'
        elif f'M{"ᵖʳᵒ"}' in enzyme:
            enzyme = f'SARS-CoV-2'
        if inFigTitle:
            title = f'{enzyme}\n{inFigTitle}\nSpearman ρ: {rho:.3f}'
        else:
            title = f'{enzyme}\nSpearman ρ: {rho:.3f}'

        # Plot bar graph
        fig, ax = plt.subplots(figsize=inFigSize)
        ax.bar(xTicks - barWidth / 2, yValues[0], barWidth, label=l1,
               color=inColor2, edgecolor='black', linewidth=inLinewidth)
        ax.bar(xTicks + barWidth / 2, yValues[1], barWidth, label=l2,
               color=inColor1, edgecolor='black', linewidth=inLinewidth)
        plt.title(title, fontsize=inTitleSize, fontweight='bold')
        ax.set_ylabel('Normalized Activity', fontsize=inLabelSize)

        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(inLinewidth)

        # Legend
        legend_props = {
            'size': inLabelTickSize-2,
            'weight': 'bold'
        }
        ax.legend(edgecolor='black', prop=legend_props, loc='best')

        # Set xticks
        ax.set_xticks(xTicks)
        ax.set_xticklabels(substrates, rotation=45)

        # Set yticks
        ax.set_ylim([0, 1.1])

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=inTickLength,
                       labelsize=inLabelTickSize, width=inLinewidth)

        fig.canvas.mpl_connect('key_press_event', pressKey)
        plt.tight_layout()
        plt.show()

        if inSavePath:
            figName = f'{enzyme}-enzActivity-bars.png'
            if inSaveTag:
                figName = figName.replace('.png', f'-{inSaveTag}.png')
            if inNormalizeSub:
                figName = figName.replace('.png',
                                          f'-Norm_{inNormalizeSub}.png')
            path = os.path.join(inSavePath, figName)
            fig.savefig(path, dpi=inFigResolution)
            print(f'Saving figure at path:\n'
                  f'     {path}\n\n')
        else:
            print(f'The figure was not saved\n\n')


def pdata(dataset):
    print(f'Data:')
    for k, v in dataset.items():
        print(f'{k}:')
        for t, d in v.items():
            print(f'  {t}: {d}')
        print()


def processData(dataset, tags, zTags):
    dataset = normalizeData(dataset, tags)
    dataset = zScore(dataset, tags, zTags)
    # pdata(dataset)
    if 'Mᵖʳᵒ2 CleaveNet' in dataset.keys():
        dataset['Mᵖʳᵒ2 CleaveNet']['Activity Z'] = inExpPredActivityCN
        dataset[f'Mᵖʳᵒ2 CleaveNet']['Predicted Z']= inPredActivityCN
        if f'M{"ᵖʳᵒ"}2' in dataset.keys():
            dataset[f'M{"ᵖʳᵒ"}2']['Activity Z'] = inExpActivity
            dataset[f'M{"ᵖʳᵒ"}2']['Predicted Z']= inPredActivity

    # Build tables
    tables = {}
    for enzyme in dataset.keys():
        columns = list(dataset[enzyme])
        df = pd.DataFrame(0.0, index=[], columns=[])
        for col in columns:
            df.loc[:, col] = dataset[enzyme][col]
        tables[enzyme] = df
    for enzyme, table in tables.items():
        print(f'{enzyme}:')
        print(f'{table.to_string(index=False)}\n')

        if inSavePath and inSaveTables:
            fileName = f'enzActivity-table-{enzyme}.csv'
            if inSaveTag:
                fileName = fileName.replace('.csv', f'-{inSaveTag}.csv')
            if inNormalizeSub:
                fileName = fileName.replace('.csv',
                                            f'-Norm_{inNormalizeSub}.csv')
            fileName = fileName.replace(' ', '_').replace('/', '')
            path = os.path.join(inSavePath, fileName)
            table.to_csv(path, index=False)
    print()

    return dataset, tables



# ========================================================================================
data, tables = processData(
    dataset=data, tags=dataTags, zTags=dataZTags
)

# Plot data
if inPlotBarGraph:
    plotBars(dataset=tables)



# ========================================================================================
def plotTrendline(axis, xVals, yVals, color):
    x_fit, y_fit, fitCurve = fitData(x=xVals, y=yVals, printOutputs=False)
    # color = '#101010'
    axis.plot(x_fit, y_fit, color=color,
            linestyle='-', linewidth=inLinewidth)


# Plot data
fig, ax = plt.subplots(figsize=inFigSize)
plt.title(inFigTitle, fontsize=inTitleSize, fontweight='bold')
x, y = f'Activity Z {inEnzyme}', f'Predicted Z {inEnzyme}'
for idx, enzyme in enumerate(data.keys()):
    xData = data[enzyme][dataZTags[0]] # Plot the line behind all datapoints
    yData = data[enzyme][dataZTags[1]]
    plotTrendline(axis=ax, xVals=xData, yVals=yData, color=inPlotColors[idx])

for idx, enzyme in enumerate(data.keys()):
    x = data[enzyme][dataZTags[0]]
    y = data[enzyme][dataZTags[1]]

    # Spearman rank correlation
    rho, p = spearmanr(x, y)
    print(f'Enzyme: {enzyme}\n* Spearman ρ: {rho:.3f}, p={p:.3f}\n')

    # Figure labels
    label = enzyme
    l1, l2, = f'M{"ᵖʳᵒ"}2', f'M{"ᵖʳᵒ"}'
    if l1 in label:
        label = label.replace(
            l1, f'SARS-CoV-2 {l1.replace("2", "")}'
        )
    elif l2 in label:
        label = label.replace(l2, f'SARS-CoV M{"ᵖʳᵒ"}')
    x_fit, y_fit, fitCurve = fitData(x=x, y=y)
    label += f' R² = {fitCurve:.3f}, ρ = {rho:.3f}'

    # Add data
    edgeWidth = 1
    ax.plot(x, y, color=inPlotColors[idx], marker=inPlotMarkers[idx],
            linestyle='none', markeredgecolor='black',
            markeredgewidth=edgeWidth, label=label)
ax.legend(prop=FontProperties(size=inLabelTickSize - 2, weight='bold'),
          edgecolor='black', linewidth=inLinewidth, loc='upper left', framealpha=0.9,
          handlelength=0.8,  handletextpad=0.2, borderpad=0.4, columnspacing=0.8)

# Set the thickness of the figure border
for _, spine in ax.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(inLinewidth)

# X Axis
ax.set_xlabel('Experimental Activity Z Scores', fontsize=inLabelSize)
yTicks = np.arange(2, -2, 0.5)
# ax.set_xticks(yTicks)
# ax.set_xlim([-2.0, 2.0]) # Set x-axis
ax.set_xticks(ax.get_xticks())
ax.set_xticklabels(ax.get_xticklabels(), ha='center', fontsize=inLabelSize-2)

# Y Axis
ax.set_ylabel('Predicted Activity Z Scores', fontsize=inLabelSize)
ax.set_yticks(ax.get_yticks())


# Set tick parameters
ax.tick_params(axis='both', which='major', length=inTickLength,
                   labelsize=inLabelTickSize, width=inLinewidth)

plt.tight_layout()
fig.canvas.mpl_connect('key_press_event', pressKey)
plt.show()

if inSavePath:
    figName = 'enzActivity-scatter.png'
    if inSaveTag:
        figName = figName.replace('.png', f'-{inSaveTag}.png')
    if inNormalizeSub:
        figName = figName.replace('.png',
                                  f'-Norm_{inNormalizeSub}.png')
    figName = figName.replace('.png', f'-N_{len(data.keys())}.png')
    path = os.path.join(inSavePath, figName)
    fig.savefig(path, dpi=inFigResolution)
    print(f'Saving figure at path:\n'
          f'     {path}\n\n')
else:
    print(f'The figure was not saved\n\n')
