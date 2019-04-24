import sys
import serial
from VNA_ctrl import VNA_ctrl
from Squidstat_ctrl import Squidstat_ctrl
from experimentRunner import experiment
from experimentRunner import dataWriter

#connect to Squidstat and VNA
if len(sys.argv) < 5:
    sys.exit()
Squidstat = Squidstat_ctrl(sys.argv[1])
Squidstat.ac_cal_mode(6)
VNA = VNA_ctrl(sys.argv[2])

#run experiment, write data
exp = experiment(sys.argv[3], VNA, Squidstat)
data = [];
dataheader = []
if exp.sweepType == 'frequency':
    dataheader = ['Frequency', 'Magnitude', 'Phase', 'Input power (dBm)']
else:
    dataheader = ['Input power (dBm)', 'Magnitude', 'Phase']
while True:
    exp.setup()
    data += exp.runExperiment()
    if exp.IsExperimentComplete():
        break;

data = exp.normalizeData(data)

#write data

data_writer = dataWriter(sys.argv[4])
data_writer.writeData(dataList = data, headerList = dataheader)