import sys
import csv
import math
import cmath
import os
from copy import deepcopy
from scipy.signal import savgol_filter
from dataParserHelper import experimentDataSet
from dataParserHelper import popSweepAvg
from dataParserHelper import popSweepAvgWithStats
from dataParserHelper import parseRawCSV
from dataParserHelper import InvertPhase
from dataParserHelper import getColumn
from dataParser3 import writeCSV
from PolMath import *
import matplotlib.pyplot as plt

def monotonicallyDescendingOrder(dataList):
    max = 0
    for i in range(len(dataList) - 1, -1, -1):
        if dataList[i] > max:
            max = dataList[i]
        else:
            dataList[i] = max

def windowedRoundup(dataList, windowWidth):
    for i in range(0, len(dataList)):
        a = max(dataList[i:i+windowWidth])
        a = math.log10(a)
        a = math.ceil(a * 10) / 10
        a = math.pow(10, a)
        dataList[i] = a

def stdevSmoothing(sweepData):
    magList = getColumn(sweepData, 1)
    phaseList = getColumn(sweepData, 2)
    windowedRoundup(magList, 3)
    monotonicallyDescendingOrder(magList)
    windowedRoundup(phaseList, 3)
    monotonicallyDescendingOrder(phaseList)
    for i in range(0, len(sweepData)):
        sweepData[i][1] = magList[i]
        sweepData[i][2] = phaseList[i]

def parse_process_write(basefolder, experimentName):
    if not (basefolder[-1] == '\\' or basefolder[-1] == '/'):
        basefolder += '\\'
    rawData = parseRawCSV(basefolder + experimentName + '/')
    InvertPhase(rawData)
    ImpedanceDataToWrite = []
    StdevDataToWrite = []

    while len(rawData) > 0:
        sweep, sweep_stats = popSweepAvgWithStats(rawData)
        stdevSmoothing(sweep_stats)
        ImpedanceDataToWrite += sweep
        StdevDataToWrite += sweep_stats

    directoryName1 = 'c:/potentiostat/squidstatcalibrator/bin/debug/ACdataStandards/' + experimentName + ' impedance/'
    directoryName2 = 'c:/potentiostat/squidstatcalibrator/bin/debug/ACdataStandards/' + experimentName + ' stdev/'
    if not os.path.isdir('c:/potentiostat/squidstatcalibrator/bin/debug/ACdataStandards/'):
        os.mkdir('c:/potentiostat/squidstatcalibrator/bin/debug/ACdataStandards/')
    if not os.path.isdir(directoryName1):
        os.mkdir(directoryName1)
    if not os.path.isdir(directoryName2):
        os.mkdir(directoryName2)
    writeCSV(directoryName1 + 'data.csv', ImpedanceDataToWrite)
    writeCSV(directoryName2 + 'data.csv', StdevDataToWrite)

def main():
    basefolder = sys.argv[1]
    parse_process_write(basefolder, '1 Ohm run')
    parse_process_write(basefolder, '10 Ohm run')
    parse_process_write(basefolder, '100 Ohm run')
    parse_process_write(basefolder, '1kOhm run')
    parse_process_write(basefolder, '10kOhm run')
    parse_process_write(basefolder, '100kOhm run')
    parse_process_write(basefolder, '10MOhm run')
    print("Script complete, press enter to continue...")
    input()


if __name__ == '__main__':
    main()