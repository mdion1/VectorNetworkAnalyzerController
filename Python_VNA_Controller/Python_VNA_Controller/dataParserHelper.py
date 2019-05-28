from copy import deepcopy

class experimentDataSet:
    def __init__(self, rawdata, DatasetName = 'Dataset'):
        self.sweepParameters = []
        self.rawData = []
        self.dataName = DatasetName
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
            while True:
                mag = 0
                phase = 0
                freq = 0
                IsigMag = 0
                N = 0
                while True:
                    row = masterTable.pop(0)
                    freq = row[0]
                    mag += row[1]
                    phase += row[2]
                    if len(row) > 3:
                        IsigMag += row[3]
                    N += 1
                    if (len(masterTable) == 0) or (masterTable[0][0] != freq):
                        break
                ret.append([freq, mag / N, phase / N])
                if (len(masterTable) == 0) or (masterTable[0][0] > freq):
                    break
            return ret

class sweepParams:
    def __init__(self):
        self.frequencyList = []
    def appendFrequency(self, freq):
        self.frequencyList.append(freq)
    def checkForCompleteness(self, frequenciesPresent, outStringList = []):
        sweepComplete = True
        for frequency in self.frequencyList:
            if not (frequency in frequenciesPresent):
                sweepComplete = False
                outStringList.append('\t\tMissing frequency: ' + str(frequency) + 'Hz\n')
        return sweepComplete
