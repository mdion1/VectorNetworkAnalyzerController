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
data_writer = dataWriter(sys.argv[4])
exp = experiment(sys.argv[3], VNA, Squidstat)
dataTable = []
dataheader = []
experimentIndex = 0
if exp.sweepType == 'frequency':
    dataheader = ['Frequency', 'Magnitude', 'Phase', 'Input power (dBm)']

while True:
    dataTable.append(exp.runExperiment())
    if exp.IsExperimentComplete():
        break;

for i in range(0, len(dataTable)):
    if exp.sweepType == 'frequency':
        dataTable[i] = exp.normalizeData(dataTable[i])
        data_writer.writeHeader(dataheader)
        data_writer.writeData(dataTable[i])
    elif exp.sweepType == 'power':
        data_writer.writeHeader('frequency = ' + str(exp.getCenterFrequencies()[i]) + 'Hz')
        data_writer.writeData(dataTable[i])