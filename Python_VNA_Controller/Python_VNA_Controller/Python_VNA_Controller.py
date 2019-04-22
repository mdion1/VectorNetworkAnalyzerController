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

#run experiment
exp = experiment(sys.argv[3], VNA, Squidstat)
data = exp.runExperiment()

#write data
dataheader = {'Frequency', 'Magnitude', 'Phase'}
data_writer = dataWriter(sys.argv[4])
data_writer.writeData(dataList = data, headerList = dataheader)
