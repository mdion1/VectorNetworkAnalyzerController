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
powerBaselineSweep = False
if len(sys.argv) == 6:
    if sys.argv[5] == 'power_sweep_baseline':
        powerBaselineSweep = True

#run experiment, write data
exp = experiment(sys.argv[3], VNA, Squidstat)
dataTable = []
dataheader = []
experimentIndex = 0
if exp.sweepType == 'frequency':
    if powerBaselineSweep:
        dataheader = ['Frequency', 'Magnitude', 'Phase', 'Output power (dB)']
    else:
        dataheader = ['Frequency', 'Magnitude', 'Phase']

while True:
    if (exp.sweepType == 'frequency'):
        if powerBaselineSweep:
            dataTable.append(exp.runPowerBaselineSweep())
        else:
            dataTable.append(exp.runFreqSweepExperiment())
    elif (exp.sweepType == 'power'):
        dataTable.append(exp.runPowerSweepExperiment())
    if exp.IsExperimentComplete():
        break;

data_writer = dataWriter(sys.argv[4])
for i in range(0, len(dataTable)):
    if exp.sweepType == 'frequency':
        dataTable[i].reverse()
        #if not powerBaselineSweep:
        #    dataTable[i] = exp.normalizeData(dataTable[i])
        data_writer.writeHeader(dataheader)
        data_writer.writeData(dataTable[i])
    elif exp.sweepType == 'power':
        data_writer.writeHeader(['frequency', str(exp.getCenterFrequencies()[i])])
        data_writer.writeData(dataTable[i])