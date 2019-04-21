import serial
import math

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
        self.write('MEAS AB')      #set A/B meas mode
        self.write('ATTR 0DB')     #set channel R attenuation
        self.write('ATTA 0DB')     #set channel A attenuation
        self.write('ATTB 0DB')     #set channel B attenuation
        self.write('BWAUTO 0')     #turn off automatic IF BW selection
        self.write('AVER 1')       #turns on sample averaging

    def setSweepType(self, sweeptype, start, stop, numPoints, signalStrength, centerFreq):
        if sweeptype == 'frequency':
            self.write('POWE ' + convertVoltsToDBM(signalStrength))      #set the input signal strength
            self.write('SWPT LOGF')                 #set sweep type to log frequency
            self.write('STAR ' + start + 'HZ')      #set sweep start
            self.write('STOP ' + stop + 'HZ')       #set sweep end
        elif sweeptype == 'power':
            self.write('???')                       #todo: set center frequency
            self.write('???')                       #todo: set sweep type to power
            self.write('STAR ' + start + 'DB')      #set sweep start
            self.write('STOP ' + stop + 'DB')       #set sweep end
        else:
            return
        self.write('POIN ' + numPoints)

    def setAverNum(self, averageNum):
        self.write('AVERFACT ' + averageNum)

    def trigSweeps(self, numSweeps):
        self.write('NUMG ' + numSweeps)

    def waitForDataReady(self):
        while True:
            self.write('ESB?')
            resp = self.read()
            if len(resp) >= 3:
                if resp[1] == 0x31:     # 0x31 is ASCII '1'
                    break

    def write(self, msg):
        if msg[-1] != '\n':
            msg += '\n'
        rawBytes = msg.encode('utf-8')
        self.ser.write(rawBytes)

    def read(self, maxLength = 256):
        x = self.ser.read(maxLength)
        return x

def convertVoltsToDBM(volts):
    volts_float = float(volts)
    dBm = 10 * math.log10(1000 * volts_float * volts_float / 50);
    dBm = round(dBm, 1)
    return str(dBm)