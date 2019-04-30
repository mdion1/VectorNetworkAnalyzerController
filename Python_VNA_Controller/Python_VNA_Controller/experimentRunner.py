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
        with open(self._filename, mode = 'w') as csv_file:
            self.csv_writer = csv.writer(csv_file, delimiter=',')
        csv_file.close()
    def writeHeader(self, headerList):
        with open(self._filename, mode = 'a') as csv_file:
            self.csv_writer = csv.writer(csv_file, delimiter=',')
            self.csv_writer.writerow(headerList)
        csv_file.close()
    def writeData(self, dataList):
        with open(self._filename, mode = 'a') as csv_file:
            self.csv_writer = csv.writer(csv_file, delimiter=',')
            for i in range(0, len(dataList)):
                self.csv_writer.writerow(dataList[i])
        csv_file.close()

class experiment:
    centerFrequenciesIndex = 0
    centerFrequencies = []
    pwrsweepTempStart = 0
    pwrsweepTempEnd = 0
    pwrsweepTempNumPoints = 0
    pwrsweepSegmentCount = 0
    sweepComplete = False
    def __init__(self, paramsFile, VNAobj, SquidstatObj):
        self.paramsReader = experimentParamsReader(paramsFile)
        self.VNA = VNAobj
        self.Squid = SquidstatObj
        self.sweepType = self.paramsReader.getParam('sweepType')
        if self.sweepType == 'power':
            self.centerFrequencies = self.paramsReader.getParam('sweepCenterFreq')
        self.isExperimentComplete = False

    def getCenterFrequencies(self):
        return self.centerFrequencies

    def initPwrSweep(self):
        start = float(self.paramsReader.getParam('sweepStart'))
        end = float(self.paramsReader.getParam('sweepEnd'))
        numPoints = int(self.paramsReader.getParam('NumPoints'))
        pwrsweepStep = (end - start) / (numPoints - 1)

        if self.pwrsweepSegmentCount == 0:
            self.pwrsweepTempStart = start
        
        self.pwrsweepTempNumPoints = int(min(20, end - self.pwrsweepTempStart) / pwrsweepStep) + 1
        self.pwrsweepTempEnd = self.pwrsweepTempStart + (self.pwrsweepTempNumPoints - 1) * pwrsweepStep
        self.VNA.setSweepType(sweeptype = 'power',
                                  start = str(self.pwrsweepTempStart),
                                  stop = str(self.pwrsweepTempEnd),
                                  numPoints = str(self.pwrsweepTempNumPoints),
                                  centerFreq = self.centerFrequencies[self.centerFrequenciesIndex])
        self.pwrsweepTempStart = self.pwrsweepTempEnd + pwrsweepStep
        self.pwrsweepSegmentCount += 1
        if self.pwrsweepTempStart > end:
            return True
        else:
            return False

    def initFreqSweep(self):
        start = float(self.paramsReader.getParam('sweepStart'))
        end = float(self.paramsReader.getParam('sweepEnd'))
        numPoints = int(self.paramsReader.getParam('NumPoints'))
        pwrsweepStep = math.pow(end / start, 1 / (numPoints - 1))

        if self.pwrsweepSegmentCount == 0:
            self.pwrsweepTempStart = start

        if self.pwrsweepTempStart < 1500:
            self.pwrsweepTempNumPoints = int(math.log(1500 / self.pwrsweepTempStart) / math.log(pwrsweepStep)) + 1
            self.pwrsweepTempEnd = math.pow(pwrsweepStep, self.pwrsweepTempNumPoints - 1) * self.pwrsweepTempStart
        else:
            self.pwrsweepTempEnd = end
            self.pwrsweepTempNumPoints = int(math.log(end / self.pwrsweepTempStart) / math.log(pwrsweepStep)) + 1

        self.VNA.setSweepType(sweeptype = 'frequency',
                                  start = str(self.pwrsweepTempStart),
                                  stop = str(self.pwrsweepTempEnd),
                                  numPoints = str(self.pwrsweepTempNumPoints),
                                  signalStrength = self.paramsReader.getParam('VoltageAmplitude'))
        self.pwrsweepTempStart = self.pwrsweepTempEnd * pwrsweepStep
        self.pwrsweepSegmentCount += 1
        if self.pwrsweepTempStart > end:
            return True
        else:
            return False

    def setup(self):
        self.Squid.ac_cal_mode(self.paramsReader.getParam('AC_CAL_MODE'))
        self.VNA.setup_basline_settings()
        if self.sweepType == 'power':
            self.sweepComplete = self.initPwrSweep()
        else:
            self.sweepComplete = True
            self.sweepComplete = self.initFreqSweep()
        self.VNA.setAverNum(self.paramsReader.getParam('AveragingNum'))
        
    def runExperiment(self):
        self.pwrsweepSegmentCount = 0
        xdata = []
        ret = []
        while True:
            #setup parameters, send to VNA
            self.setup()

            #trigger A/B meas sweep
            self.VNA.trigSweeps_AB(self.paramsReader.getParam('AveragingNum'))
            self.VNA.waitForDataReady()
            #get data
            rawdata1 = self.VNA.downloadPolarData(self.paramsReader.getParam('NumPoints'))

            #trigger B meas sweep
            self.VNA.trigSweeps_B('3')
            self.VNA.waitForDataReady()
            #get data
            rawdata2 = self.VNA.downloadPolarData(self.paramsReader.getParam('NumPoints'))
        
            if (self.paramsReader.getParam('sweepType') == 'power'):
                start = self.pwrsweepTempStart
                stop = self.pwrsweepTempEnd
                numPoints = self.pwrsweepTempNumPoints
                xdata = getLinearList(start, stop, numPoints)
                ret += combineData(xdata, rawdata1, rawdata2, fourthColumn = '')
            else:
                xdata = getLogList(float(self.paramsReader.getParam('sweepStart')),
                                   float(self.paramsReader.getParam('sweepEnd')),
                                   int(self.paramsReader.getParam('NumPoints')))
                ret += combineData(xdata, rawdata1, rawdata2, fourthColumn = 'power')

            if self.sweepComplete:
                break;

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
            dataList[i][2] = fixPhase(dataList[i][2])
        return dataList

def fixPhase(phase):
    if phase > 135:
        return phase - 360
    if phase < -225:
        return phase + 360

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
        phaseList.append(fixPhase(-cmath.phase(a) * 180 / cmath.pi))

    #parse powerData list for complex pairs
    for i in range(0, len(powerData), 2):
        a = complex(powerData[i], powerData[i+1])
        powerList.append(20 * math.log10(abs(a) * signalB_atten))

    #combine xdata and ydata
    for i in range(0, len(xdata)):
        if (fourthColumn == 'power'):
            masterList.append([xdata[i], magList[i], phaseList[i], powerList[i]])
        else:
            masterList.append([powerList[i], magList[i], phaseList[i]])
    return masterList

