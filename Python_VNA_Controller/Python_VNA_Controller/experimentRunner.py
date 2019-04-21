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

class experiment:
    def __init__(self, paramsFile, VNAobj, SquidstatObj):
        self.paramsReader = experimentParamsReader(paramsFile)
        self.VNA = VNAobj
        self.Squid = SquidstatObj
        self.setup()

    def setup(self):
        self.Squid.ac_cal_mode(self.paramsReader.getParam('AC_CAL_MODE'))
        self.VNA.setup_basline_settings() 
        #self.VNA.setSweepType(sweeptype = self.paramsReader.getParam('sweepType'),
        #                      start = self.paramsReader.getParam('sweepStart'),
        #                      stop = self.paramsReader.getParam('sweepStop'),
        #                      numPoints = self.paramsReader.getParam('sweepNumPoints'),
        #                      signalStrength = self.paramsReader.getParam('signalStrength'),
        #                      centerFreq = self.paramsReader.getParam('sweepCenterFreq'))
        self.VNA.setAverNum(self.paramsReader.getParam('AveragingNum'))
        
    def runExperiment(self):
        #trigger sweeps
        self.VNA.trigSweeps('10') #self.VNA.trigSweeps(self.paramsReader.getParam('AveragingNum'))
        self.VNA.waitForDataReady()
        print("done")