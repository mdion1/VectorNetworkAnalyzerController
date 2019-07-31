import sys
import csv
import math
import cmath
from copy import deepcopy
import traceback
#from scipy.signal import savgol_filter
#print('savgol_filter imported')
from dataParserHelper import experimentDataSet
print('experimentDataSet imported')
from dataParserHelper import popSweepAvg
from dataParserHelper import interpolate
from dataParserHelper import parseRawCSV
from dataParserHelper import InvertPhase
from dataParserHelper import getColumn
print('helper functions imported')

#****************** Helper functions **********************
smoothing_on = False

def PolSmooth(datatable, window = 5, order = 1, re_normalize = True):
    if smoothing_on:
        real_list = []
        imag_list = []
        for i in range(0,len(datatable)):
            mag = datatable[i][1]
            phi = datatable[i][2] / 180 * math.pi
            real_list.append(mag * math.cos(phi))
            imag_list.append(mag * math.sin(phi))
        real_list = savgol_filter(real_list, window, order)
        imag_list = savgol_filter(imag_list, window, order)
        for i in range(0,len(datatable)):
            x = complex(real_list[i], imag_list[i])
            datatable[i][1] = abs(x)
            datatable[i][2] = cmath.phase(x) * 180 / math.pi
    if (re_normalize):
        normalize(datatable)

def normalize(table):
    BasePhase = table[-1][2]
    BaseMag = table[-1][1]
    for i in range(0,len(table)):
        table[i][1] /= BaseMag
        table[i][2] -= BasePhase

def writeRow(csvwriter, datarow):
    if isinstance(datarow, list):
        if isinstance(datarow[0], list):
            for i in range(0, len(datarow)):
                writeRow(csvwriter, datarow[i])
        else:
            csvwriter.writerow(datarow)
    else:
        csvwriter.writerow([datarow])

def writeCSV(filename, data, header = ['Frequency', 'Magnitude', 'Phase(degrees)']):
    with open(filename, mode = 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', lineterminator = '\n')
        writeRow(csv_writer, header)
        writeRow(csv_writer, data)
    csv_file.close()

def PolDivide(numerator, denominator):
    for i in range(0, len(numerator)):
        freq = numerator[i][0]
        numerator[i][1] /= interpolate(freq, denominator, 1)
        numerator[i][2] -= interpolate(freq, denominator, 2)

def PolInv(baseTable):
    ret = deepcopy(baseTable)
    for i in range(0, len(baseTable)):
        ret[i][1] = 1 / baseTable[i][1]
        ret[i][2] = -baseTable[i][2]
    return ret

def PolMult(baseTable, auxTable):
    for i in range(0, len(baseTable)):
        freq = baseTable[i][0]
        baseTable[i][1] *= interpolate(freq, auxTable, 1)
        baseTable[i][2] += interpolate(freq, auxTable, 2)

def PolAdd(baseTable, auxTable):
    for i in range(0, len(baseTable)):
        freq = baseTable[i][0]
        base_phi = baseTable[i][2] / 180 * math.pi
        base = complex(baseTable[i][1] * math.cos(base_phi), baseTable[i][1] * math.sin(base_phi))
        aux_phi = interpolate(freq, auxTable, 2) / 180 * math.pi
        aux_mag = interpolate(freq, auxTable, 1)
        aux = complex(aux_mag * math.cos(base_phi), aux_mag * math.sin(base_phi))
        sum = base + aux
        baseTable[i][1] = abs(sum)
        baseTable[i][2] = cmath.phase(sum) * 180 / math.pi

def PolSubtract(baseTable, auxTable):
    for i in range(0, len(baseTable)):
        freq = baseTable[i][0]
        base_phi = baseTable[i][2] / 180 * math.pi
        base = complex(baseTable[i][1] * math.cos(base_phi), baseTable[i][1] * math.sin(base_phi))
        aux_phi = interpolate(freq, auxTable, 2) / 180 * math.pi
        aux_mag = interpolate(freq, auxTable, 1)
        aux = complex(aux_mag * math.cos(base_phi), aux_mag * math.sin(base_phi))
        diff = base - aux
        baseTable[i][1] = abs(diff)
        baseTable[i][2] = cmath.phase(diff) * 180 / math.pi

def sortByFirstElement(val):
    return val[0]

def ScalarMult(baseTable, scalar):
    for i in range(0, len(baseTable)):
        baseTable[i][1] *= scalar

#*************** Main script ****************************

def main():
    try:
        if smoothing_on:
            from scipy.signal import savgol_filter
        basefolder = sys.argv[1]
        if not (basefolder[-1] == '\\' or basefolder[-1] == '/'):
            basefolder += '\\'

        DataComplete = True
        ErrorStrList = []
    
        # import raw data from "1 Ohm Run" experiment
        rawData1R = parseRawCSV(basefolder + '1 Ohm run/')
        InvertPhase(rawData1R)
        dataParseCheck = experimentDataSet(rawData1R, DatasetName = '1 Ohm Run')
        DataComplete &= dataParseCheck.verifyData(ErrorStrList)
        rawData1R = dataParseCheck.getRawSweepData()
        print('1 Ohm run parsed')

        # import raw data from "10 Ohm Run" experiment
        rawData10R = parseRawCSV(basefolder + '10 Ohm run/')
        InvertPhase(rawData10R)
        dataParseCheck = experimentDataSet(rawData10R, DatasetName = '10 Ohm Run')
        DataComplete &= dataParseCheck.verifyData(ErrorStrList)
        rawData10R = dataParseCheck.getRawSweepData()
        print('10 Ohm run parsed')

        # import raw data from "100 Ohm Run" experiment
        rawData100R = parseRawCSV(basefolder + '100 Ohm run/')
        InvertPhase(rawData100R)
        dataParseCheck = experimentDataSet(rawData100R, DatasetName = '100 Ohm Run')
        DataComplete &= dataParseCheck.verifyData(ErrorStrList)
        rawData100R = dataParseCheck.getRawSweepData()
        print('100 Ohm run parsed')

        # import raw data from "1kOhm Run" experiment
        rawData1k = parseRawCSV(basefolder + '1kOhm run/')
        InvertPhase(rawData1k)
        dataParseCheck = experimentDataSet(rawData1k, DatasetName = '1kOhm Run')
        DataComplete &= dataParseCheck.verifyData(ErrorStrList)
        rawData1k = dataParseCheck.getRawSweepData()
        print('1 kOhm run parsed')

        # import raw data from "10kOhm Run" experiment
        rawData10k = parseRawCSV(basefolder + '10kOhm run/')
        InvertPhase(rawData10k)
        dataParseCheck = experimentDataSet(rawData10k, DatasetName = '10kOhm Run')
        DataComplete &= dataParseCheck.verifyData(ErrorStrList)
        rawData10k = dataParseCheck.getRawSweepData()
        print('10 kOhm run parsed')

        # import raw data from "100kOhm Run" experiment
        rawData100k = parseRawCSV(basefolder + '100kOhm run/')
        InvertPhase(rawData100k)
        dataParseCheck = experimentDataSet(rawData100k, DatasetName = '100kOhm Run')
        DataComplete &= dataParseCheck.verifyData(ErrorStrList)
        rawData100k = dataParseCheck.getRawSweepData()
        print('100 kOhm run parsed')

        # import raw data from "10MOhm Run" experiment
        rawData10M = parseRawCSV(basefolder + '10MOhm run/')
        InvertPhase(rawData10M)
        dataParseCheck = experimentDataSet(rawData10M, DatasetName = '10MOhm Run')
        DataComplete &= dataParseCheck.verifyData(ErrorStrList)
        rawData10M = dataParseCheck.getRawSweepData()
        print('10 MOhm run parsed')

        #*****************************************************************************************
        if not DataComplete:
            logFile = open('C:\\potentiostat\\squidstatcalibrator\\bin\\debug\\ACcalibration_logfile.txt', 'w')
            for substr in ErrorStrList:
                print(substr)
                logFile.write(substr)
                logFile.write('\r\n')
            logFile.close()
        else:
            #****************************************** parse sweeps from each test***********************************************
            # 1 Ohm run
            R2_V100_I1_1R = popSweepAvg(rawData1R)
            R2_V500_I1_1R = popSweepAvg(rawData1R)
            R2_V1000_I1_1R = popSweepAvg(rawData1R)

            # 10 Ohm run -- Range0/Range1 sweeps
            R2_V1_I1_10R = popSweepAvg(rawData10R)
            R2_V10_I1_10R = popSweepAvg(rawData10R)
            R1_V10_I1_10R = popSweepAvg(rawData10R)
            R1_V10_I10_10R = popSweepAvg(rawData10R)
            R0_V10_I10_10R = popSweepAvg(rawData10R)

            # 10 Ohm run -- Vgain50/100/200 sweeps
            R2_V20_I1_10R = popSweepAvg(rawData10R)
            R2_V50_I1_10R = popSweepAvg(rawData10R)
            R2_V50_I2_10R = popSweepAvg(rawData10R)
            R2_V50_I5_10R = popSweepAvg(rawData10R)
            R2_V50_I10_10R = popSweepAvg(rawData10R)
            R2_V50_I20_10R = popSweepAvg(rawData10R)
            R2_V50_I50_10R = popSweepAvg(rawData10R)
            R2_V100_I1_10R = popSweepAvg(rawData10R)
            R2_V100_I2_10R = popSweepAvg(rawData10R)
            R2_V100_I5_10R = popSweepAvg(rawData10R)
            R2_V100_I10_10R = popSweepAvg(rawData10R)
            R2_V100_I20_10R = popSweepAvg(rawData10R)
            R2_V100_I50_10R = popSweepAvg(rawData10R)
            R2_V200_I1_10R = popSweepAvg(rawData10R)
            R2_V200_I2_10R = popSweepAvg(rawData10R)
            R2_V200_I5_10R = popSweepAvg(rawData10R)
            R2_V200_I10_10R = popSweepAvg(rawData10R)
            R2_V200_I20_10R = popSweepAvg(rawData10R)
            R2_V200_I50_10R = popSweepAvg(rawData10R)

            # 100 Ohm run
            R2_V1_I1_100R = popSweepAvg(rawData100R)
            R2_V2_I1_100R = popSweepAvg(rawData100R)
            R2_V5_I1_100R = popSweepAvg(rawData100R)
            R2_V10_I1_100R = popSweepAvg(rawData100R)
            R2_V1_I2_100R = popSweepAvg(rawData100R)
            R2_V1_I5_100R = popSweepAvg(rawData100R)
            R2_V1_I10_100R = popSweepAvg(rawData100R)
            #R2_V10_I10_100R = popSweep(rawData100R)
            R3_V10_I1_100R = popSweepAvg(rawData100R)

            # 1kOhm run
            R2_V1_I10_1k = popSweepAvg(rawData1k)
            R2_V1_I20_1k = popSweepAvg(rawData1k)
            R2_V1_I50_1k = popSweepAvg(rawData1k)
            R2_V1_I100_1k = popSweepAvg(rawData1k)
            R3_V10_I1_1k = popSweepAvg(rawData1k)  #R3_V10_I10_1k = popSweep(rawData1k)
            R4_V10_I1_1k = popSweepAvg(rawData1k)
            R5_V10_I1_1k = popSweepAvg(rawData1k)

            # 10kOhm run
            R2_V1_I100_10k = popSweepAvg(rawData10k)
            R2_V1_I200_10k = popSweepAvg(rawData10k)
            R2_V1_I500_10k = popSweepAvg(rawData10k)
            R2_V1_I1000_10k = popSweepAvg(rawData10k)

            # 100kOhm run
            R5_V10_I1_100k = popSweepAvg(rawData100k)
            R6_V10_I1_100k = popSweepAvg(rawData100k)
            R7_V100_I1_100k = popSweepAvg(rawData100k)

            # 10MOhm run
            R6_V10_I100_10M = popSweepAvg(rawData10M)
            R8_V10_I1_10M = popSweepAvg(rawData10M)
            R9_V100_I1_10M = popSweepAvg(rawData10M)

            #****************************************** calculations ***********************************************
            range2 = deepcopy(R2_V1_I1_100R)
            H_BaselineV = deepcopy(R2_V1_I1_100R)

            # I/V Gains 1, 2, 5, 10 (from 100 Ohm run)
            V2_I1 = deepcopy(R2_V2_I1_100R)
            PolDivide(V2_I1, H_BaselineV); PolSmooth(V2_I1)
            V5_I1 = deepcopy(R2_V5_I1_100R)
            PolDivide(V5_I1, H_BaselineV); PolSmooth(V5_I1)
            V10_I1 = deepcopy(R2_V10_I1_100R)
            PolDivide(V10_I1, H_BaselineV); PolSmooth(V10_I1)
            V1_I2 = deepcopy(R2_V1_I2_100R)
            PolDivide(V1_I2, H_BaselineV); PolSmooth(V1_I2)
            V1_I5 = deepcopy(R2_V1_I5_100R)
            PolDivide(V1_I5, H_BaselineV); PolSmooth(V1_I5)
            V1_I10 = deepcopy(R2_V1_I10_100R)
            PolDivide(V1_I10, H_BaselineV); # smooth later
            #V10_I10 = deepcopy(R2_V10_I10_100R)
            #PolDivide(V10_I10, H_BaselineV); PolSmooth(V10_I10)

            # VGains 20, 50, 100 (from 10 Ohm run)
            V20_I1 = deepcopy(R2_V20_I1_10R); PolDivide(V20_I1, R2_V1_I1_10R); PolSmooth(V20_I1)
            V50_I1 = deepcopy(R2_V50_I1_10R); PolDivide(V50_I1, R2_V1_I1_10R); PolSmooth(V50_I1)
            V50_I2 = deepcopy(R2_V50_I2_10R); PolDivide(V50_I2, R2_V1_I1_10R); PolSmooth(V50_I2)
            V50_I5 = deepcopy(R2_V50_I5_10R); PolDivide(V50_I5, R2_V1_I1_10R); PolSmooth(V50_I5)
            V50_I10 = deepcopy(R2_V50_I10_10R); PolDivide(V50_I10, R2_V1_I1_10R); PolSmooth(V50_I10)
            V50_I20 = deepcopy(R2_V50_I20_10R); PolDivide(V50_I20, R2_V1_I1_10R); PolSmooth(V50_I20)
            V50_I50 = deepcopy(R2_V50_I50_10R); PolDivide(V50_I50, R2_V1_I1_10R); PolSmooth(V50_I50)
            V100_I1 = deepcopy(R2_V100_I1_10R); PolDivide(V100_I1, R2_V1_I1_10R); # smooth later
            V100_I2 = deepcopy(R2_V100_I2_10R); PolDivide(V100_I2, R2_V1_I1_10R); PolSmooth(V100_I2)
            V100_I5 = deepcopy(R2_V100_I5_10R); PolDivide(V100_I5, R2_V1_I1_10R); PolSmooth(V100_I5)
            V100_I10 = deepcopy(R2_V100_I10_10R); PolDivide(V100_I10, R2_V1_I1_10R); PolSmooth(V100_I10)
            V100_I20 = deepcopy(R2_V100_I20_10R); PolDivide(V100_I20, R2_V1_I1_10R); PolSmooth(V100_I20)
            V100_I50 = deepcopy(R2_V100_I50_10R); PolDivide(V100_I50, R2_V1_I1_10R); PolSmooth(V100_I50)
            V200_I1 = deepcopy(R2_V200_I1_10R); PolDivide(V200_I1, R2_V1_I1_10R); PolSmooth(V200_I1)
            V200_I2 = deepcopy(R2_V200_I2_10R); PolDivide(V200_I2, R2_V1_I1_10R); PolSmooth(V200_I2)
            V200_I5 = deepcopy(R2_V200_I5_10R); PolDivide(V200_I5, R2_V1_I1_10R); PolSmooth(V200_I5)
            V200_I10 = deepcopy(R2_V200_I10_10R); PolDivide(V200_I10, R2_V1_I1_10R); PolSmooth(V200_I10)
            V200_I20 = deepcopy(R2_V200_I20_10R); PolDivide(V200_I20, R2_V1_I1_10R); PolSmooth(V200_I20)
            V200_I50 = deepcopy(R2_V200_I50_10R); PolDivide(V200_I50, R2_V1_I1_10R); PolSmooth(V200_I50)

            #VGains 500, 1000 (from 1 Ohm run)
            V500_I1 = deepcopy(R2_V500_I1_1R)
            PolDivide(V500_I1, R2_V100_I1_1R); PolMult(V500_I1, V100_I1); PolSmooth(V500_I1)
            V1000_I1 = deepcopy(R2_V1000_I1_1R)
            PolDivide(V1000_I1, R2_V100_I1_1R); PolMult(V1000_I1, V100_I1); PolSmooth(V1000_I1)
            PolSmooth(V100_I1)

            # IGains 20, 50, 100 (from 1kOhm run)
            V1_I20 = deepcopy(R2_V1_I20_1k)
            PolDivide(V1_I20, R2_V1_I10_1k); PolMult(V1_I20, V1_I10); PolSmooth(V1_I20)
            V1_I50 = deepcopy(R2_V1_I50_1k)
            PolDivide(V1_I50, R2_V1_I10_1k); PolMult(V1_I50, V1_I10); PolSmooth(V1_I50)
            V1_I100 = deepcopy(R2_V1_I100_1k)
            PolDivide(V1_I100, R2_V1_I10_1k); PolMult(V1_I100, V1_I10) #smooth later
            PolSmooth(V1_I10)

            # IGains 200, 500, 1000 (from 10kOhm run)
            V1_I200 = deepcopy(R2_V1_I200_10k)
            PolDivide(V1_I200, R2_V1_I100_10k); PolMult(V1_I200, V1_I100); PolSmooth(V1_I200)
            V1_I500 = deepcopy(R2_V1_I500_10k)
            PolDivide(V1_I500, R2_V1_I100_10k); PolMult(V1_I500, V1_I100); PolSmooth(V1_I500)
            V1_I1000 = deepcopy(R2_V1_I1000_10k)
            PolDivide(V1_I1000, R2_V1_I100_10k); PolMult(V1_I1000, V1_I100); PolSmooth(V1_I1000)
            PolSmooth(V1_I100)

            #range2
            normalize(range2)
            range2 = PolInv(range2)

            #range1
            range1 = deepcopy(R1_V10_I1_10R)
            PolDivide(range1, R2_V10_I1_10R)
            normalize(range1)
            range1 = PolInv(range1)
            PolMult(range1, range2)

            #range0
            range0 = deepcopy(R0_V10_I10_10R)
            PolDivide(range0, R1_V10_I10_10R)
            normalize(range0)
            range0 = PolInv(range0)
            PolMult(range0, range1)

            #range3
            range3 = deepcopy(R3_V10_I1_100R)
            PolDivide(range3, R2_V10_I1_100R)
            normalize(range3)
            range3 = PolInv(range3)
            PolMult(range3, range2)

            #range4
            range4 = deepcopy(R4_V10_I1_1k)
            PolDivide(range4, R3_V10_I1_1k)
            normalize(range4)
            range4 = PolInv(range4)
            PolMult(range4, range3)

            #range5
            range5 = deepcopy(R5_V10_I1_1k)
            PolDivide(range5, R3_V10_I1_1k)
            normalize(range5)
            range5 = PolInv(range5)
            PolMult(range5, range3)

            #range6
            range6 = deepcopy(R6_V10_I1_100k)
            normalize(R5_V10_I1_100k)
            PolSmooth(R5_V10_I1_100k)
            PolDivide(range6, R5_V10_I1_100k)
            normalize(range6)
            range6 = PolInv(range6)
            PolMult(range6, range5)

            #range7
            range7 = deepcopy(R7_V10_I10_10M)
            PolDivide(range7, R5_V10_I1_100k)
            normalize(range7)
            range7 = PolInv(range7)
            PolMult(range7, range5)

            #range8
            range8 = deepcopy(R8_V10_I1_10M)
            normalize(R6_V10_I100_10M)
            PolSmooth(R6_V10_I100_10M)
            PolDivide(range8, R6_V10_I100_10M)
            normalize(range8)
            range8 = PolInv(range8)
            PolMult(range8, range6)

            #range9
            range9 = deepcopy(R9_V100_I1_10M)
            PolDivide(range9, R6_V10_I100_10M)
            normalize(range9)
            range9 = PolInv(range9)
            PolMult(range9, range6)

            #****************************************** debugging section ***********************************************

            #r3_1k_05ac = popSweep(parseRawCSV(basefolder + 'r3_v10_i10_1k_05ac/'))
            #r3_1k_04ac = popSweep(parseRawCSV(basefolder + 'r3_v10_i10_1k_04ac/'))
            #r3_1k_03ac = popSweep(parseRawCSV(basefolder + 'r3_v10_i10_1k_03ac/'))
            #r3_1k_02ac = popSweep(parseRawCSV(basefolder + 'r3_v10_i10_1k_02ac/'))
            #r3_1k_01ac = popSweep(parseRawCSV(basefolder + 'r3_v10_i10_1k_01ac/'))
            #PolDivide(r3_1k_04ac, r3_1k_05ac)
            #PolDivide(r3_1k_03ac, r3_1k_05ac)
            #PolDivide(r3_1k_02ac, r3_1k_05ac)
            #PolDivide(r3_1k_01ac, r3_1k_05ac)
            #x = [[['40mV']], r3_1k_04ac, [['30mV']], r3_1k_03ac, [['20mV']], r3_1k_02ac, [['10mV']], r3_1k_01ac]
            #writeCSV(basefolder + 'amp_comparison.csv', x)

            #r3_1k = popSweep(parseRawCSV(basefolder + 'r3_v10_i10_1k_05ac/'))
            #r3_100R = popSweep(parseRawCSV(basefolder + 'r3_v10_i10_100R_01ac/'))
            #PolDivide(r3_1k, r3_100R)
            #writeCSV(basefolder + '1k_curve.csv', r3_1k)

            #****************************************** write CSV files ***********************************************
            writeCSV(basefolder + 'WEgain1_Igain2.csv', V1_I2)
            writeCSV(basefolder + 'WEgain1_Igain5.csv', V1_I5)
            writeCSV(basefolder + 'WEgain1_Igain10.csv', V1_I10)
            writeCSV(basefolder + 'WEgain1_Igain20.csv', V1_I20)
            writeCSV(basefolder + 'WEgain1_Igain50.csv', V1_I50)
            writeCSV(basefolder + 'WEgain1_Igain100.csv', V1_I100)
            writeCSV(basefolder + 'WEgain1_Igain200.csv', V1_I200)
            writeCSV(basefolder + 'WEgain1_Igain500.csv', V1_I500)
            writeCSV(basefolder + 'WEgain1_Igain1000.csv', V1_I1000)
            writeCSV(basefolder + 'WEgain2_Igain1.csv', V2_I1)
            writeCSV(basefolder + 'WEgain5_Igain1.csv', V5_I1)
            writeCSV(basefolder + 'WEgain10_Igain1.csv', V10_I1)
            writeCSV(basefolder + 'WEgain20_Igain1.csv', V20_I1)
            writeCSV(basefolder + 'WEgain50_Igain1.csv', V50_I1)
            writeCSV(basefolder + 'WEgain50_Igain2.csv', V50_I2)
            writeCSV(basefolder + 'WEgain50_Igain5.csv', V50_I5)
            writeCSV(basefolder + 'WEgain50_Igain10.csv', V50_I10)
            writeCSV(basefolder + 'WEgain50_Igain20.csv', V50_I20)
            writeCSV(basefolder + 'WEgain50_Igain50.csv', V50_I50)
            writeCSV(basefolder + 'WEgain100_Igain1.csv', V100_I1)
            writeCSV(basefolder + 'WEgain100_Igain2.csv', V100_I2)
            writeCSV(basefolder + 'WEgain100_Igain5.csv', V100_I5)
            writeCSV(basefolder + 'WEgain100_Igain10.csv', V100_I10)
            writeCSV(basefolder + 'WEgain100_Igain20.csv', V100_I20)
            writeCSV(basefolder + 'WEgain100_Igain50.csv', V100_I50)
            writeCSV(basefolder + 'WEgain200_Igain1.csv', V200_I1)
            writeCSV(basefolder + 'WEgain200_Igain2.csv', V200_I2)
            writeCSV(basefolder + 'WEgain200_Igain5.csv', V200_I5)
            writeCSV(basefolder + 'WEgain200_Igain10.csv', V200_I10)
            writeCSV(basefolder + 'WEgain200_Igain20.csv', V200_I20)
            writeCSV(basefolder + 'WEgain200_Igain50.csv', V200_I50)
            writeCSV(basefolder + 'WEgain500_Igain1.csv', V500_I1)
            writeCSV(basefolder + 'WEgain1000_Igain1.csv', V1000_I1)
            writeCSV(basefolder + 'Range0.csv', range0)
            writeCSV(basefolder + 'Range1.csv', range1)
            writeCSV(basefolder + 'Range2.csv', range2)
            writeCSV(basefolder + 'Range3.csv', range3)
            writeCSV(basefolder + 'Range4.csv', range4)
            writeCSV(basefolder + 'Range5.csv', range5)
            writeCSV(basefolder + 'Range6.csv', range6)
            writeCSV(basefolder + 'Range7.csv', range7)
            writeCSV(basefolder + 'Range8.csv', range8)
            writeCSV(basefolder + 'Range9.csv', range9)

            # write "dummy" files
            dummyData = [[1000000, 1, 0], [1, 1, 0]]
            V1_I1 = deepcopy(dummyData)
            writeCSV(basefolder + 'V_inputBuffer.csv', dummyData)
            writeCSV(basefolder + 'WEgain1_Igain1.csv', V1_I1)

            #****************************** data testing section*******************************************
            results = []

            # 1 Ohm Sweeps
            results += [['1 Ohm sweeps']]
            DataList = [R2_V100_I1_1R, R2_V500_I1_1R, R2_V1000_I1_1R]
            RangeList = [range2, range2, range2]
            GainList1 = [V100_I1, V500_I1, V1000_I1]
            for i in range(0, len(DataList)):
                data = deepcopy(DataList[i])
                PolMult(data, RangeList[i])
                PolDivide(data, GainList1[i])
                results += data

            # 10 Ohm Sweeps, Range0/Range1
            results += [['10 Ohm sweeps, Range0/1']]
            DataList = [R2_V1_I1_10R, R2_V10_I1_10R, R1_V10_I1_10R, R1_V10_I10_10R, R0_V10_I10_10R]
            RangeList = [range2, range2, range1, range1, range0]
            GainList1 = [V1_I1, V10_I1, V10_I1, V10_I1, V10_I1]
            GainList2 = [V1_I1, V10_I1, V10_I1, V1_I10, V1_I10]
            for i in range(0, len(DataList)):
                data = deepcopy(DataList[i])
                PolMult(data, RangeList[i])
                PolDivide(data, GainList1[i])
                PolDivide(data, GainList2[i])
                results += data

            # 10 Ohm Sweeps, Vgain50/100/200
            results += [['10 Ohm sweeps, Vgain50/100/200']]
            DataList = [R2_V20_I1_10R, R2_V50_I1_10R, R2_V50_I2_10R, R2_V50_I5_10R, R2_V50_I10_10R, R2_V50_I20_10R, R2_V50_I50_10R,
                        R2_V100_I1_10R, R2_V100_I2_10R, R2_V100_I5_10R, R2_V100_I10_10R, R2_V100_I20_10R, R2_V100_I50_10R,
                        R2_V200_I1_10R, R2_V200_I2_10R, R2_V200_I5_10R, R2_V200_I10_10R, R2_V200_I20_10R, R2_V200_I50_10R]
            RangeList = [range2, range2, range2, range2, range2, range2, range2,
                         range2, range2, range2, range2, range2, range2,
                         range2, range2, range2, range2, range2, range2]
            GainList1 = [V20_I1, V50_I1, V50_I2, V50_I5, V50_I10, V50_I20, V50_I50,
                         V100_I1, V100_I2, V100_I5, V100_I10, V100_I20, V100_I50,
                         V200_I1, V200_I2, V200_I5, V200_I10, V200_I20, V200_I50]
            for i in range(0, len(DataList)):
                data = deepcopy(DataList[i])
                PolMult(data, RangeList[i])
                PolDivide(data, GainList1[i])
                results += data

            # 100 Ohm Sweeps
            results += [['100 Ohm sweeps']]
            DataList = [R2_V1_I1_100R, R2_V2_I1_100R, R2_V5_I1_100R, R2_V10_I1_100R, R2_V1_I2_100R, R2_V1_I5_100R, R2_V1_I10_100R, R3_V10_I1_100R]
            RangeList = [range2, range2, range2, range2, range2, range2, range2, range3]
            GainList1 = [V1_I1, V2_I1, V5_I1, V10_I1, V1_I2, V1_I5, V1_I10, V10_I1]
            for i in range(0, len(DataList)):
                data = deepcopy(DataList[i])
                PolMult(data, RangeList[i])
                PolDivide(data, GainList1[i])
                results += data

            # 1kOhm Sweeps
            results += [['1kOhm sweeps']]
            DataList = [R2_V1_I10_1k, R2_V1_I20_1k, R2_V1_I50_1k, R2_V1_I100_1k, R3_V10_I1_1k, R4_V10_I1_1k, R5_V10_I1_1k]
            RangeList = [range2, range2, range2, range2, range3, range4, range5]
            GainList1 = [V1_I10, V1_I20, V1_I50, V1_I100, V10_I1, V10_I1, V10_I1]
            for i in range(0, len(DataList)):
                data = deepcopy(DataList[i])
                PolMult(data, RangeList[i])
                PolDivide(data, GainList1[i])
                results += data

            # 10kOhm Sweeps
            results += [['10kOhm sweeps']]
            DataList = [R2_V1_I100_10k, R2_V1_I200_10k, R2_V1_I500_10k, R2_V1_I1000_10k]
            RangeList = [range2, range2, range2, range2]
            GainList1 = [V1_I100, V1_I200, V1_I500, V1_I1000]
            for i in range(0, len(DataList)):
                data = deepcopy(DataList[i])
                PolMult(data, RangeList[i])
                PolDivide(data, GainList1[i])
                results += data

            # 100kOhm Sweeps
            results += [['100kOhm sweeps']]
            DataList = [R5_V10_I1_100k, R6_V10_I1_100k]
            RangeList = [range5, range6]
            GainList1 = [V10_I1, V10_I1]
            for i in range(0, len(DataList)):
                data = deepcopy(DataList[i])
                PolMult(data, RangeList[i])
                PolDivide(data, GainList1[i])
                results += data

            # 10MOhm Sweeps
            results += [['10MOhm sweeps']]
            DataList = [R6_V10_I100_10M, R8_V10_I1_10M]
            RangeList = [range6, range8]
            for i in range(0, len(DataList)):
                data = deepcopy(DataList[i])
                PolMult(data, RangeList[i])
                results += data

            writeCSV(basefolder + '/testresults.csv', results)

        print("Script complete, press enter to close this window...")
        input()
    except Exception as e:
        traceback.print_exc(file=sys.stdout)
        input()

if __name__ == '__main__':
    main()