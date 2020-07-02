import glob
print('glob imported')
import math
import csv
import statistics
print('statistics imported')
#from numpy import percentile
import numpy as np
print('numpy imported')
import matplotlib.pyplot as plt
import PolMath

#****************** Helper functions **********************

FREQUENCY_COLUMN = 6
MAGNITUDE_COLUMN = 7
PHASE_COLUMN = 8
EWE_AMPLITUDE_COLUMN = 10
NUMBER_OF_SAVED_COLUMNS = 4

def condenseRawCSV(filename):
    dataTable = []
    with open(filename, mode = 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            if len(row) <= NUMBER_OF_SAVED_COLUMNS:
                return;
            row_float = []
            for i in [FREQUENCY_COLUMN, MAGNITUDE_COLUMN, PHASE_COLUMN, EWE_AMPLITUDE_COLUMN]:
                try:
                    row_float.append(float(row[i]))
                except ValueError:
                    continue
            if len(row_float) > 0:
                dataTable.append(row_float)
    with open(filename, 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', lineterminator = '\n')
        csv_writer.writerow(['Frequency', 'Magnitude', 'Phase(degrees)', 'AC amplitude (mV)'])
        for i in range(0, len(dataTable)):
            csv_writer.writerow(dataTable[i])
        csv_file.close()

def parseRawCSV(directoryName):
    dataTable = []
    filenames = glob.glob(directoryName + "*.csv")
    if len(filenames) > 0:
        filename = filenames[0]
        try:
            condenseRawCSV(filename)
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

def InvertPhase(baseTable):
    for i in range(0, len(baseTable)):
        baseTable[i][2] *= -1

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
                print('stdev mag')
            if data[i][2] > mag_phase * (1 + phase_err_margin):
                ret_indices.append(i)
                print('stdev phase')
        else:
            if not areApproximatelyEqual(data[i][1], mag_baseline, margin = mag_err_margin):
                ret_indices.append(i)
                print('mag')
            if not areApproximatelyEqual(data[i][2], mag_phase, margin = phase_err_margin, Abs_error = True):
                ret_indices.append(i)
                print('phase')
    #if len(ret_indices) > 0:
    #    freqList = getColumn(data, 0)
    #    magList = getColumn(data, 1)
    #    phaseList = getColumn(data, 2)
    #    freqList_base = getColumn(baseline, 0)
    #    magList_base = getColumn(baseline, 1)
    #    phaseList_base = getColumn(baseline, 2)
    #    plt.plot(freqList, magList, 'r')
    #    plt.plot(freqList_base, magList_base, 'b')
    #    plt.xscale('log')
    #    plt.show()

    return ret_indices

def findFurthestPoint(list):
    avg = statistics.mean(list)
    index = 0
    max_err = 0
    for i in range(0, len(list)):
        if abs(list[i] - avg) > max_err:
            index = i
            max_err = abs(list[i] - avg)
    return index

def findOutliers(maglist, phaselist, margin = 1.5):
    #outliers are based on maglist data, since clipping doesn't affect the phase as much
    q25 = np.percentile(maglist, 25)
    q75 = np.percentile(maglist, 75)
    iqr = q75 - q25
    cut_off = iqr * margin
    lower, upper = q25 - cut_off, q75 + cut_off
    for i in range(len(maglist) - 1, -1, -1):    # loop thru indices in descending order
        if maglist[i] < lower or maglist[i] > upper:
            maglist.pop(i)
            phaselist.pop(i)

def popSweepAvg(masterTable):
    ret = []
    sweepRaw = masterTable.pop(0)
    for i in range(0, len(sweepRaw)):
        ret.append([sweepRaw[i][0], statistics.mean(sweepRaw[i][1]), statistics.mean(sweepRaw[i][2])])
    return ret

def popSweepAvgWithStats(masterTable):
    ret = []
    retStats = []
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
        if len(mag) > 1:
            ret.append([freq, statistics.mean(mag), statistics.mean(phase)])
            retStats.append([freq, statistics.stdev(mag), statistics.stdev(phase)])
        else:
            ret.append([freq, mag[0], phase[0]])
            retStats.append([freq, mag[0], phase[0]])
        if (len(masterTable) == 0) or (masterTable[0][0] > freq):
            break
    return ret, retStats

def popSweepFull(masterTable):
    ret = []
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
        if len(mag) > 1:
            ret.append([freq, mag, phase])
        else:
            ret.append([freq, mag[0], phase[0]])
        if (len(masterTable) == 0) or (masterTable[0][0] > freq):
            break
    return ret

#****************** "experimentDataSet" class **********************
class experimentDataSet:
    def __init__(self, rawdata, DatasetName = 'Dataset'):
        self.sweepObjs = []
        self.sweepParameters = []
        self.rawData = []
        self.dataName = DatasetName
        self.loadComparisonDataset(DatasetName)
        i = 0
        while len(rawdata) > 0:
            sweepData = popSweepFull(rawdata)
            self.rawData.append(sweepData)
            self.sweepObjs.append(sweep(sweepData, self.comparisonData[i], self.comparisonDataSweepStats[i]))   #not the cleanest way to add comparisonData and comparisonDataSweepStats...
            i += 1

    def getRawSweepData(self):
        return self.rawData

    def loadComparisonDataset(self, dataset_name):
        self.comparisonData = []
        self.comparisonDataSweepStats = []
        dummy = []
        raw_data = parseRawCSV('c:/potentiostat/squidstatcalibrator/bin/debug/ACdataStandards/' + dataset_name + ' impedance/')
        while len(raw_data) > 0:
            sweep = popSweepFull(raw_data)
            self.comparisonData.append(sweep)
        raw_data = parseRawCSV('c:/potentiostat/squidstatcalibrator/bin/debug/ACdataStandards/' + dataset_name + ' stdev/')
        while len(raw_data) > 0:
            sweep = popSweepFull(raw_data)
            self.comparisonDataSweepStats.append(sweep)

    def verifyData(self, outStringList = []):
        dataVerified = True
        tempStringList = []

        if len(self.sweepObjs) < len(self.comparisonData):
            tempStringList.append('\tSweeps missing: expected ' + str(len(self.comparisonData)) + ' sweeps but read ' + str(len(self.sweepObjs)) + '.')
            dataVerified = False

        for i in range(0, len(self.sweepObjs)):
            sweep = self.sweepObjs[i]
            sweepVerified = sweep.checkForCompleteness()
            if not sweepVerified:
                tempStringList.append('\tSweep ' + str(i + 1) + ' is missing data points.')
            
            sweepVerified = sweep.cleanData()
            if not sweepVerified:
                tempStringList.append('\tSweep ' + str(i + 1) + ' has noisy or inaccurate data.')
            
            dataVerified &= sweepVerified

        if dataVerified:
            outStringList.append(self.dataName + ' complete\n')
        else:
            outStringList.append(self.dataName + ' errors:\n')
            outStringList += tempStringList
        return dataVerified
    def getStats(self):
        return self.sweepStats

class sweep:
    def __init__(self, experimentalData, comparisonData, stdevData):
        self.loadRawData(experimentalData)
        self.loadComparisonData(comparisonData)
        self.loadStDevData(stdevData)
    def loadRawData(self, rawData):
        self.allData = {}
        ZModData = []
        phaseData = []
        freqList = []
        for i in range(0, len(rawData)):
            freqList.append(rawData[i][0])
            ZModData.append(rawData[i][1])
            phaseData.append(rawData[i][2])
        self.allData["Zmod"] = ZModData
        self.allData["Phase"] = phaseData
        self.allData["Frequency list"] = freqList

    def loadComparisonData(self, rawData):
        ZModData = []
        phaseData = []
        freqData = []
        while len(rawData) > 0:
            row = rawData.pop(0)
            freqData.append(row[0])
            ZModData.append(row[1])
            phaseData.append(row[2])
        self.allData["Freq comparison data"] = freqData
        self.allData["Zmod comparison data"] = ZModData
        self.allData["Phase comparison data"] = phaseData

    def loadStDevData(self, rawData):
        ZModStDevData = []
        phaseStDevData = []
        while len(rawData) > 0:
            row = rawData.pop(0)
            ZModStDevData.append(row[1])
            phaseStDevData.append(row[2])
        self.allData["Zmod stdev comparison data"] = ZModStDevData
        self.allData["Phase stdev comparison data"] = phaseStDevData

    def checkForCompleteness(self):
        frequencyList = self.allData["Freq comparison data"]
        frequenciesPresent = self.allData["Frequency list"]
        sweepComplete = True
        for frequency in frequencyList:
            frequencyPresentInList = False
            for i in range(0, len(frequenciesPresent)):
                if areApproximatelyEqual(frequency, frequenciesPresent[i]):
                    frequencyPresentInList = True
                    break
            if not frequencyPresentInList:
                sweepComplete = False
        return sweepComplete
    
    def cleanData(self):
        success = True
        zmod_margin = 0.5
        phase_margin = 20
        zmod_stdev_margin = 2
        phase_stdev_margin = 3
        for i in range(0, len(self.allData["Frequency list"])):
            sampleSize = len(self.allData["Zmod"][i])
            Zmod_baseline = self.allData["Zmod comparison data"][i]
            Phase_baseline = self.allData["Phase comparison data"][i]
            stdevZmod_max = self.allData["Zmod stdev comparison data"][i] * zmod_stdev_margin
            stdevPhase_max = self.allData["Phase stdev comparison data"][i] * phase_stdev_margin
            while True:
                dataOk = True
                _Zmod = statistics.mean(self.allData["Zmod"][i])
                _phase = statistics.mean(self.allData["Phase"][i])
                stdevZmod = statistics.stdev(self.allData["Zmod"][i])
                stdevPhase = statistics.stdev(self.allData["Phase"][i])
                
                if (stdevZmod > stdevZmod_max * zmod_stdev_margin) or not areApproximatelyEqual(_Zmod, Zmod_baseline, zmod_margin):
                    dataOk = False
                    badIndex = findFurthestPoint(self.allData["Zmod"][i])
                    self.allData["Zmod"][i].pop(badIndex)
                    self.allData["Phase"][i].pop(badIndex)
                if (stdevPhase > stdevPhase_max * phase_stdev_margin) or not areApproximatelyEqual(_phase, Phase_baseline, phase_margin, Abs_error = True):
                    dataOk = False
                    badIndex = findFurthestPoint(self.allData["Phase"][i])
                    self.allData["Zmod"][i].pop(badIndex)
                    self.allData["Phase"][i].pop(badIndex)
                if dataOk:
                    break
                elif len(self.allData["Zmod"][i]) <= sampleSize * 0.6:
                    success = False
                    break
            #if not success:
                #break
        if not success:
            success = self.decideFromPlot()
        return success

    def decideFromPlot(self):
        freqList = []
        magList = []
        magStDevList = []
        phaseList = []
        phaseStDevList = []
        for i in range(0, len(self.allData["Zmod"])):
            for j in range(0, len(self.allData["Zmod"][i])):
                freqList.append(self.allData["Frequency list"][i])
                magList.append(self.allData["Zmod"][i][j])
                phaseList.append(self.allData["Phase"][i][j])
            magStDevList.append(statistics.stdev(self.allData["Zmod"][i]))
            phaseStDevList.append(statistics.stdev(self.allData["Phase"][i]))
        freqList_base = self.allData["Freq comparison data"]
        magList_base = self.allData["Zmod comparison data"]
        phaseList_base = self.allData["Phase comparison data"]
        magStDevList_base = self.allData["Zmod stdev comparison data"]
        phaseStDevList_base = self.allData["Phase stdev comparison data"]

        # Zmod subplot
        plt.subplot(2,2,1)
        plt.title('Zmod')
        plt.plot(freqList, magList, 'r')
        plt.plot(freqList_base, magList_base, 'b')
        plt.xscale('log')
        # Phase subplot
        plt.subplot(2,2,2)
        plt.title('Phase')
        plt.plot(freqList, phaseList, 'r')
        plt.plot(freqList_base, phaseList_base, 'b')
        plt.xscale('log')
        # Zmod stdev subplot
        plt.subplot(2,2,3)
        plt.title('Zmod stdev')
        plt.plot(freqList_base, magStDevList, 'r')
        plt.plot(freqList_base, magStDevList_base, 'b')
        plt.xscale('log')
        # Phase stdev subplot
        plt.subplot(2,2,4)
        plt.title('Phase stdev')
        plt.plot(freqList_base, phaseStDevList, 'r')
        plt.plot(freqList_base, phaseStDevList_base, 'b')
        plt.xscale('log')
        plt.show()

        response = input()
        return response[0] == 'Y' or response[0] == 'y'