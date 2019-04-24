import csv
import cmath

class experimentParamsReader:
    def __init__(self, filename):
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            self.masterDictionary = {};
            for row in csv_reader:
                self.masterDictionary[row[0]] = row[1]

    def getParam(self, param):
        if param in self.masterDictionary:
            return self.masterDictionary[param]
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
    def __init__(self, paramsFile, VNAobj, SquidstatObj):
        self.paramsReader = experimentParamsReader(paramsFile)
        self.VNA = VNAobj
        self.Squid = SquidstatObj
        self.setup()

    def setup(self):
        self.Squid.ac_cal_mode(self.paramsReader.getParam('AC_CAL_MODE'))
        self.VNA.setup_basline_settings() 
        self.VNA.setSweepType(sweeptype = self.paramsReader.getParam('sweepType'),
                              start = self.paramsReader.getParam('sweepStart'),
                              stop = self.paramsReader.getParam('sweepEnd'),
                              numPoints = self.paramsReader.getParam('NumPoints'),
                              signalStrength = self.paramsReader.getParam('VoltageAmplitude'),
                              centerFreq = self.paramsReader.getParam('sweepCenterFreq'))
        self.VNA.setAverNum(self.paramsReader.getParam('AveragingNum'))
        
    def runExperiment(self):
        #trigger sweeps
        self.VNA.trigSweeps(self.paramsReader.getParam('AveragingNum'))
        self.VNA.waitForDataReady()
        #get data
        rawdata = self.VNA.downloadData(self.paramsReader.getParam('NumPoints'))
        xdata = []
        ret = []
        start = self.paramsReader.getParam('sweepStart')
        stop = self.paramsReader.getParam('sweepEnd')
        numPoints = self.paramsReader.getParam('NumPoints')
        if (self.paramsReader.getParam('sweepType') == 'power'):
            xdata = getLinearList(float(start), float(stop), int(numPoints))
            ret = getNormalizedPolarData(xdata, rawdata)
            #todo: finish handling power sweep data
        else:
            xdata = getLogList(float(start), float(stop), int(numPoints))
            ret = getNormalizedPolarData(xdata, rawdata)
        return ret

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

def getNormalizedPolarData(xdata, ydata):
    masterList = []
    phaseList = []
    magList = []
    #parse ydata list for complex pairs
    for i in range(0, len(ydata), 2):
        a = complex(ydata[i], ydata[i+1])
        magList.append(abs(a))
        phaseList.append(cmath.phase(a) * 180 / cmath.pi)

    #order sweep from lowest to highest frequence
    magList.reverse()
    phaseList.reverse()
    xdata.reverse()

    #normalize complex mag and phase to the last
    for i in range(0, len(magList)):
        magList[i] /= magList[-1]
        phaseList[i] = phaseList[-1] - phaseList[i]
        if phaseList[i] > 180:
            phaseList[i] -= 360
        if phaseList[i] <= -180:
            phaseList[i] += 360

    #combine xdata and ydata
    for i in range(0, len(xdata)):
        masterList.append([xdata[i], magList[i], phaseList[i]])
    return masterList