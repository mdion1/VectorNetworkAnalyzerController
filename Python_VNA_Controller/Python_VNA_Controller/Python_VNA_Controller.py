import sys
import serial
from VNA_ctrl import VNA_ctrl
from Squidstat_ctrl import Squidstat_ctrl
from experimentRunner import experiment

#connect to Squidstat and VNA
if len(sys.argv) < 4:
    sys.exit()
Squidstat = Squidstat_ctrl(sys.argv[1])
Squidstat.ac_cal_mode(6)
VNA = VNA_ctrl(sys.argv[2])
exp = experiment(sys.argv[3], VNA, Squidstat)
exp.runExperiment()
