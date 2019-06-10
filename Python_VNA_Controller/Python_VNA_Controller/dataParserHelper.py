from copy import deepcopy
import statistics
from numpy import percentile

def areApproximatelyEqual(x, y, margin = 0.001):
    upperLim = y * (1 + margin)
    lowerLim = y * (1 - margin)
    if y > 0:
        return ((x >= lowerLim) and (x <= upperLim))
    else:
        return ((x >= upperLim) and (x <= lowerLim))

def doesNumberApproximatelyMatchListEntry(number, list):
    for i in range(0, len(list)):
        if areApproximatelyEqual(number, list[i]):
            return True
    return False

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

class experimentDataSet:
    def __init__(self, rawdata, DatasetName = 'Dataset'):
        self.sweepParameters = []
        self.rawData = []
        self.dataName = DatasetName
        self.sweepStats = [[DatasetName]]
        temp = deepcopy(rawdata)
        while len(temp) > 0:
            self.rawData.append(self.popSweep(temp))
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
        if datasetComplete:
            outStringList.append(self.dataName + ' complete\n')
        else:
            outStringList.append(self.dataName + ' errors:\n')
            outStringList += tempStringList
        return datasetComplete

    def popSweep(self, masterTable):
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
            sweepStatTable.append([freq, statistics.stdev(mag) / avg_mag, (max(mag) - min(mag)) /avg_mag, statistics.stdev(phase), max(phase) - min(phase)])
            if (len(masterTable) == 0) or (masterTable[0][0] > freq):
                break
        self.sweepStats.append(sweepStatTable)
        return ret
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
            if not doesNumberApproximatelyMatchListEntry(frequency, frequenciesPresent):
                sweepComplete = False
                outStringList.append('\t\tMissing frequency: ' + str(frequency) + 'Hz\n')
        return sweepComplete
