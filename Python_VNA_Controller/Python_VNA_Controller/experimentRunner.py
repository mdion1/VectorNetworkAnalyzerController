import csv

class experimentParamsReader:
    def __init__(self, filename):
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            self.masterDictionary = {};
            for row in csv_reader:
                self.masterDictionary[row[0]] = row[1]

    def getParam(self, param):
        return self.masterDictionary[param]

class dataWriter:
    def __init__(self, filename):
        with open(filename) as csv_file:
            self.csv_writer = csv.writer(csv_file, delimiter=',')
    def writeData(self, dataList, headerList):
        self.csv_writer.writerow(headerList)
        for subList in dataList:
            self.csv_writer.writerow(subList)

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
                              stop = self.paramsReader.getParam('sweepStop'),
                              numPoints = self.paramsReader.getParam('sweepNumPoints'),
                              signalStrength = self.paramsReader.getParam('signalStrength'),
                              centerFreq = self.paramsReader.getParam('sweepCenterFreq'))
        self.VNA.setAverNum(self.paramsReader.getParam('AveragingNum'))
        
    def runExperiment(self):
        #trigger sweeps
        self.VNA.trigSweeps(self.paramsReader.getParam('AveragingNum'))
        self.VNA.waitForDataReady()
        ydata = self.VNA.downloadData()
        xdata = []
        if (self.paramsReader.getParam('sweepType') == 'frequency'):
            xdata = getLinearList()
        else:
            xdata = getLogList()
        return getDataMasterList(xdata, ydata)

def getLinearList(start, end, points):
    span = end - start
    delta = span / (points - 1)
    linList = []
    for i in range(0, points - 1):
        linList[i] = start + i * delta
    return linList


def getLogList(start, end, points):
    span = end / start
    delta = span ** (1 / (points - 1))
    logList = []
    for i in range(0, points - 1):
        logList[i] = start * (delta ** i)

def getDataMasterList(xlist, ylist : list[complex]):
    masterList = [[]]
    for i in range(0, len(xlist) - 1):
        masterList[i][0] = xlist[i]
        masterList[i][1] = abs(ylist[i])
        masterList[i][2] = cmath.phase(yList[i]) * 180 / cmath.pi
    return masterList

def getDataMasterList(xlist, ylist : list[list[float]]):
    masterList = [[]]
    for i in range(0, len(xlist) - 1):
        masterList[i][0] = xlist[i]
        for j in range(0, len(ylist[i]) - 1):
            masterList[i][j + 1] = ylist[i][j]
    return masterList