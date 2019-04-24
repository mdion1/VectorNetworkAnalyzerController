import csv
import cmath
import math

class experimentParamsReader:
    def __init__(self, filename):
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            self.masterDictionary = {};
            for row in csv_reader:
                self.masterDictionary[row[0]] = row[1:]

    def getParam(self, param):
        if param in self.masterDictionary:
            if (param == 'sweepCenterFreq'):
                return self.masterDictionary[param]
            else:
                return self.masterDictionary[param][0]
        else:
            return ''

class dataWriter:
    def __init__(self, filename):
        self._filename = filename
    def writeData(self, dataList, headerList):
        with open(self._filename, mode = 'w') as csv_file:
            self.csv_writer = csv.writer(csv_file, delimiter=',')
            self.csv_writer.writerow(headerList)
            for i in range(0, len(dataList)):
                self.csv_writer.writerow(dataList[i])

class experiment:
    centerFrequenciesIndex = 0
    centerFrequencies = []
    def __init__(self, paramsFile, VNAobj, SquidstatObj):
        self.paramsReader = experimentParamsReader(paramsFile)
        self.VNA = VNAobj
        self.Squid = SquidstatObj
        self.sweepType = self.paramsReader.getParam('sweepType')
        if self.sweepType == 'power':
            self.centerFrequencies = self.paramsReader.getParam('sweepCenterFreq')
        self.isExperimentComplete = False

    def setup(self):
        self.Squid.ac_cal_mode(self.paramsReader.getParam('AC_CAL_MODE'))
        self.VNA.setup_basline_settings()
        if self.sweepType == 'power':
            self.VNA.setSweepType(sweeptype = 'power',
                                  start = self.paramsReader.getParam('sweepStart'),
                                  stop = self.paramsReader.getParam('sweepEnd'),
                                  numPoints = self.paramsReader.getParam('NumPoints'),
                                  centerFreq = self.centerFrequencies[self.centerFrequenciesIndex])
        else:
            self.VNA.setSweepType(sweeptype = 'frequency',
                                  start = self.paramsReader.getParam('sweepStart'),
                                  stop = self.paramsReader.getParam('sweepEnd'),
                                  numPoints = self.paramsReader.getParam('NumPoints'),
                                  signalStrength = self.paramsReader.getParam('VoltageAmplitude'))
        self.VNA.setAverNum(self.paramsReader.getParam('AveragingNum'))
        
    def runExperiment(self):
        #trigger A/B meas sweep
        self.VNA.trigSweeps_AB(self.paramsReader.getParam('AveragingNum'))
        self.VNA.waitForDataReady()
        #get data
        rawdata1 = self.VNA.downloadPolarData(self.paramsReader.getParam('NumPoints'))

        #trigger B meas sweep
        self.VNA.trigSweeps_B(self.paramsReader.getParam('AveragingNum'))
        self.VNA.waitForDataReady()
        #get data
        rawdata2 = self.VNA.downloadPolarData(self.paramsReader.getParam('NumPoints'))
        
        xdata = []
        ret = []
        start = self.paramsReader.getParam('sweepStart')
        stop = self.paramsReader.getParam('sweepEnd')
        numPoints = self.paramsReader.getParam('NumPoints')
        if (self.paramsReader.getParam('sweepType') == 'power'):
            xdata = getLinearList(float(start), float(stop), int(numPoints))
            ret = combineData(xdata, rawdata1, rawdata2, fourthColumn = '')
            #todo: finish handling power sweep data
        else:
            xdata = getLogList(float(start), float(stop), int(numPoints))
            ret = combineData(xdata, rawdata1, rawdata2, fourthColumn = 'power')

        #increment centerFrequenciesIndex, update isExperimentDone
        self.centerFrequenciesIndex += 1
        if self.centerFrequenciesIndex >= len(self.centerFrequencies):
            self.isExperimentComplete = True

        return ret

    def IsExperimentComplete(self):
        return self.isExperimentComplete

    def normalizeData(self, dataList):
        for i in range(0, len(dataList)):
            dataList[i][1] /= dataList[-1][1]
            dataList[i][2] -= dataList[-1][2]
            if dataList[i][2] > 180:
                dataList[i][2] -= 360
            if dataList[i][2] <= -180:
                dataList[i][2] += 360
        return dataList

def getLinearList(start, end, points):
    span = end - start
    delta = span / (points - 1)
    linList = []
    for i in range(0, points):
        linList.append(start + i * delta)
    return linList


def getLogList(start: float, end: float, points: int):
    span = end / start
    delta = span ** (1 / (points - 1))
    logList = []
    for i in range(0, points):
        logList.append(start * (delta ** i))
    return logList

def combineData(xdata, polarData, powerData, signalB_atten = 10, fourthColumn = 'power'):
    masterList = []
    phaseList = []
    magList = []
    powerList = []
    #parse polarData list for complex pairs
    for i in range(0, len(polarData), 2):
        a = complex(polarData[i], polarData[i+1])
        magList.append(abs(a))
        phaseList.append(cmath.phase(a) * 180 / cmath.pi)

    #parse powerData list for complex pairs
    for i in range(0, len(powerData), 2):
        a = complex(powerData[i], powerData[i+1])
        powerList.append(20 * math.log10(abs(a) * signalB_atten))

    #order sweep from lowest to highest frequence
    magList.reverse()
    phaseList.reverse()
    powerList.reverse()
    xdata.reverse()

    #combine xdata and ydata
    for i in range(0, len(xdata)):
        if (fourthColumn == 'power'):
            masterList.append([xdata[i], magList[i], phaseList[i], powerList[i]])
        else:
            masterList.append([powerList[i], magList[i], phaseList[i]])
    return masterList

