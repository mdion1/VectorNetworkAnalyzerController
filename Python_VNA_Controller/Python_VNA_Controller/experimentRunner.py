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
        self.VNA.setup_basline_settings()
        self.VNA.setSweepType(sweeptype = 'frequency', start = '100', stop = '1000', numPoints = '10')
        
        #self.VNA.setSweepType(self.paramsReader.getParam('sweepType'))
        