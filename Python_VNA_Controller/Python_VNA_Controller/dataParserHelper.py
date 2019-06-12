from copy import deepcopy
import glob
import math
import csv
import statistics
from numpy import percentile
import matplotlib.pyplot as plt

#****************** Helper functions **********************
def parseRawCSV(directoryName):
    dataTable = []
    filenames = glob.glob(directoryName + "*.csv")
    if len(filenames) > 0:
        filename = filenames[0]
        try:
            with open(filename) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
        
                for row in csv_reader:
                    row_float = []
                    for i in range(0, len(row)):
                        try:
                            row_float.append(float(row[i]))
                        except ValueError:
                            continue
                    if len(row_float) > 0:
                        dataTable.append(row_float)
        except:
            print('error loading ' + filename + ', press Enter to continue...')
            input()
            return dataTable
    return dataTable

def areApproximatelyEqual(x, y, margin = 0.001, Abs_error = False):
    upperLim = y
    lowerLim = y
    if Abs_error:
        upperLim = y + abs(margin)
        lowerLim = y - abs(margin)
        return ((x >= lowerLim) and (x <= upperLim))
    else:
        upperLim = y * (1 + margin)
        lowerLim = y * (1 - margin)
        if y > 0:
            return ((x >= lowerLim) and (x <= upperLim))
        else:
            return ((x >= upperLim) and (x <= lowerLim))

def getColumn(dataTable, col):
    ret = []
    for i in range(0, len(dataTable)):
        ret.append(dataTable[i][col])
    return ret

def CompareDatasets(data, baseline, mag_err_margin, phase_err_margin, isStandardDeviationData = False):
    ret_indices = []
    for i in range(0, len(data)):
        mag_baseline = interpolate(data[i][0], baseline, ycol = 1)
        mag_phase = interpolate(data[i][0], baseline, ycol = 2)
        if isStandardDeviationData:
            if data[i][1] > mag_baseline * (1 + mag_err_margin):
                ret_indices.append(i)
            if data[i][2] > mag_phase * (1 + phase_err_margin):
                ret_indices.append(i)
        else:
            if not areApproximatelyEqual(data[i][1], mag_baseline, margin = mag_err_margin):
                ret_indices.append(i)
            if not areApproximatelyEqual(data[i][2], mag_phase, margin = phase_err_margin, Abs_error = True):
                ret_indices.append(i)
    if len(ret_indices) > 0:
        freqList = getColumn(data, 0)
        magList = getColumn(data, 1)
        phaseList = getColumn(data, 2)
        freqList_base = getColumn(baseline, 0)
        magList_base = getColumn(baseline, 1)
        phaseList_base = getColumn(baseline, 2)
        plt.plot(freqList, magList, 'r')
        plt.plot(freqList_base, magList_base, 'b')
        plt.xscale('log')
        plt.show()

    return ret_indices

def findOutliers(maglist, phaselist, margin = 1.5):
    #outliers are based on maglist data, since clipping doesn't affect the phase as much
    q25 = percentile(maglist, 25)
    q75 = percentile(maglist, 75)
    iqr = q75 - q25
    cut_off = iqr * margin
    lower, upper = q25 - cut_off, q75 + cut_off
    for i in range(len(maglist) - 1, -1, -1):    # loop thru indices in descending order
        if maglist[i] < lower or maglist[i] > upper:
            maglist.pop(i)
            phaselist.pop(i)

def interpolate(x, lookupTable, ycol = 1, semilog = True):
    #make sure table is sorted in descending order (by first column)
    if(lookupTable[0][0] < lookupTable[-1][0]):
        lookupTable.sort(key = sortByFirstElement)

    #if x is outside the bounds of the table, return one of the bounding values
    if x >= lookupTable[0][0]:
        return lookupTable[0][ycol]
    elif x <= lookupTable[-1][0]:
        return lookupTable[-1][ycol]
    #else interpolate
    else:
        indexR = 1
        for i in range(1, len(lookupTable)):
            if x > lookupTable[indexR][0]:
                break
            else:
                indexR += 1
        indexL = indexR - 1
        xL = 0.0; xR = 0.0; xi = x
        if (semilog):
            xL = math.log10(lookupTable[indexL][0])
            xR = math.log10(lookupTable[indexR][0])
            xi = math.log10(x)
        else:
            xL = lookupTable[indexL][0]
            xR = lookupTable[indexR][0]
        yL = lookupTable[indexL][ycol]
        yR = lookupTable[indexR][ycol]
        yi = yL + (xi - xL) / (xR - xL) * (yR - yL)
        return yi

def popSweep(masterTable, returnStats = False):
    ret = []
    sweepStatTable = []
    while True:
        mag = []
        phase = []
        freq = 0
        while True:
            row = masterTable.pop(0)
            freq = row[0]
            mag.append(row[1])
            phase.append(row[2])
            if (len(masterTable) == 0) or (masterTable[0][0] != freq):
                break
        findOutliers(mag, phase, margin = 5)
        avg_mag = statistics.mean(mag)
        ret.append([freq, avg_mag, statistics.mean(phase)])
        if returnStats:
            sweepStatTable.append([freq, statistics.stdev(mag) / avg_mag, (max(mag) - min(mag)) /avg_mag, statistics.stdev(phase), max(phase) - min(phase)])
        if (len(masterTable) == 0) or (masterTable[0][0] > freq):
            break
    return ret, sweepStatTable

#****************** "experimentDataSet" class **********************
class experimentDataSet:
    def __init__(self, rawdata, DatasetName = 'Dataset'):
        self.sweepParameters = []
        self.rawData = []
        self.sweepStats = []
        self.dataName = DatasetName
        self.loadComparisonDataset(DatasetName)
        temp = deepcopy(rawdata)
        while len(temp) > 0:
            sweep, sweep_stats = popSweep(temp, returnStats = True)
            self.rawData.append(sweep)
            self.sweepStats.append(sweep_stats)

    def loadComparisonDataset(self, dataset_name):
        self.comparisonData = []
        self.comparisonDataSweepStats = []
        raw_data = parseRawCSV('c:/potentiostat/squidstatcalibrator/bin/debug/ACdataStandards/' + dataset_name + '/')
        while len(raw_data) > 0:
            sweep, sweep_stats = popSweep(raw_data, returnStats = True)
            self.comparisonData.append(sweep)
            self.comparisonDataSweepStats.append(sweep_stats)

    def addSweep(self, FrequencyList = [], StartingFreq = 1e6, EndingFrequency = 1e4, PointsPerDecade = 10):
        if len(FrequencyList) == 0:
            FrequencyList = getFrequencyList(StartingFreq, EndingFrequency, PointsPerDecade)
        sweep = sweepParams()
        for i in range(0, len(FrequencyList)):
            sweep.appendFrequency(FrequencyList[i])
        self.sweepParameters.append(sweep)

    def checkForCompleteness(self, outStringList = []):
        sweepsMissing = len(self.sweepParameters) - len(self.rawData)
        if sweepsMissing > 0:
            outStringList.append(str(sweepsMissing) + ' sweeps missing from ' + self.dataName + '\n')
            return False
        elif sweepsMissing < 0:
            outStringList.append('Wrong dataset selected for ' + self.dataName + '\n')
            return False

        datasetComplete = True
        tempStringList = []
        for i in range(0, len(self.sweepParameters)):
            frequenciesPresent = []
            for j in range(0, len(self.rawData[i])):
                frequenciesPresent.append(self.rawData[i][j][0])
            missingFreqList = []
            if not self.sweepParameters[i].checkForCompleteness(frequenciesPresent, missingFreqList):
                datasetComplete = False
                tempStringList.append('\tSweep ' + str(i + 1) + ' incomplete:\n')
                tempStringList += missingFreqList

            # check that data falls within 5% of the "normal" data, and that noise is within expected bounds
            badIndices = CompareDatasets(self.rawData[i], self.comparisonData[i], 0.1, 5)
            badIndices += CompareDatasets(self.sweepStats[i], self.comparisonDataSweepStats[i], 1, 1, isStandardDeviationData = True)
            if len(badIndices) > 0:
                datasetComplete = False
                tempStringList.append('\tSweep ' + str(i + 1) + ' has bad data at the following frequencies: ')
                for badIndex in badIndices:
                    tempStringList.append('\t\t' + str(badIndex))

        if datasetComplete:
            outStringList.append(self.dataName + ' complete\n')
        else:
            outStringList.append(self.dataName + ' errors:\n')
            outStringList += tempStringList
        return datasetComplete
    def getStats(self):
        return self.sweepStats

class sweepParams:
    def __init__(self):
        self.frequencyList = []
    def appendFrequency(self, freq):
        self.frequencyList.append(freq)
    def checkForCompleteness(self, frequenciesPresent, outStringList = []):
        sweepComplete = True
        for frequency in self.frequencyList:
            frequencyPresentInList = False
            for i in range(0, len(frequenciesPresent)):
                if areApproximatelyEqual(frequency, frequenciesPresent[i]):
                    frequencyPresentInList = True
                    break
            if not frequencyPresentInList:
                sweepComplete = False
                outStringList.append('\t\tMissing frequency: ' + str(frequency) + 'Hz\n')
        return sweepComplete