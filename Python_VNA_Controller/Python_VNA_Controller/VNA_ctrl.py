import serial
import math
import array
import time

class VNA_ctrl:
    def __init__(self, comport):
        self.ser = serial.Serial(comport)
        self.ser.timeout = 0.2
        self.write('*IDN?')
        print(self.read())

    def setup_basline_settings(self):
        self.write('HOLD') 		#takes instrument out of free-run mode
        self.write('FORM3')        #sets data output mode to 64-bit double
        self.write('NA')           #set network analyzer mode
        self.write('ATTR 0DB')     #set channel R attenuation
        self.write('ATTA 0DB')     #set channel A attenuation
        self.write('ATTB 0DB')     #set channel B attenuation
        self.write('BWAUTO 0')     #turn off automatic IF BW selection
        self.write('AVER 1')       #turns on sample averaging

    def setSweepType(self, sweeptype, start, stop, numPoints, signalStrength = '0.01', centerFreq = '1000000'):   
        if sweeptype == 'power':
            self.write('CWFREQ ' + centerFreq)      #set center frequency
            self.write('SWPT POWE')                 #set sweep type to power
            self.write('STAR ' + start + 'DB')      #set sweep start
            self.write('STOP ' + stop + 'DB')       #set sweep end
        else:                                       #default sweeptype is 'frequency'
            self.write('POWE ' + convertVoltsToDBM(signalStrength))      #set the input signal strength
            self.write('SWPT LOGF')                 #set sweep type to log frequency
            self.write('STAR ' + start + 'HZ')      #set sweep start
            self.write('STOP ' + stop + 'HZ')       #set sweep end
        self.write('POIN ' + numPoints)
        self.setIFBW(sweeptype, start, stop, centerFreq)
    
    def setIFBW(self, sweeptype, start, stop, centerFreq):
        ifbw = ''
        if sweeptype == 'power':
            ifbw = self.get_IFBW_val_str((float(centerFreq)))
        else:
            ifbw = self.get_IFBW_val_str(min(float(start), float(stop)))
        self.write('BW ' + ifbw)
    
    def get_IFBW_val_str(self, freq):
        if freq >= 1500:
            return '300HZ'
        elif freq >=100:
            return '30HZ'
        elif freq >= 30:
            return '10HZ'
        else:
            return '2HZ'

    def setAverNum(self, averageNum):
        self.write('AVER 1')
        self.write('AVERFACT ' + averageNum)

    def trigSweeps_AB(self, numSweeps):
        self.write('MEAS AB')      #set A/B meas mode
        self.write('NUMG ' + numSweeps)

    def trigSweeps_B(self, numSweeps):
        self.write('MEAS B')       #set B meas mode
        self.write('NUMG ' + numSweeps)

    def waitForDataReady(self):
        while True:
            self.write('ESB?')
            resp = self.read()
            if len(resp) >= 3:
                if resp[1] == 0x31:     # 0x31 is ASCII '1'
                    break

    def downloadPolarData(self, pointsPerSweep):
        self.write('OUTPDATA?')
        minNumBytes = 8 * 2 * int(pointsPerSweep) + 8 + 1        # 8 bytes * 2 doubles per point + 8-byte header + 1-byte postscript
        rawBytes = self.read(minNumBytes)
        trimmedBytes = rawBytes[8:len(rawBytes) - 1]
        x = array.array('d', trimmedBytes)
        x.byteswap()
        return x

    def write(self, msg):
        if msg[-1] != '\n':
            msg += '\n'
        rawBytes = msg.encode('utf-8')
        self.ser.write(rawBytes)
        time.sleep(0.2)

    def read(self, maxLength = 256):
        x = self.ser.read(maxLength)
        return x

def convertVoltsToDBM(volts):
    volts_float = float(volts)
    dBm = 10 * math.log10(1000 * volts_float * volts_float / 50);
    dBm = round(dBm, 1)
    return str(dBm)