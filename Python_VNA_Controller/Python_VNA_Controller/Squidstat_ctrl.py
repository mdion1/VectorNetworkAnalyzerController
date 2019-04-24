import serial

class Squidstat_ctrl:
    def __init__(self, comport):
        self.ser = serial.Serial(comport)
        self.ser.timeout = 0.2
        connectSuccess = self.handshake()
        if connectSuccess:
            print("Successfully connected to Squidstat on " + comport)
        else:
            print("Failed to connect to Squidstat on " + comport)
        
    def handshake(self):
        handshake = [0xEE, 0xFF, 0x41, 0x00, 0x00, 0x00]
        self.write(handshake)
        rsp = self.read()
        handshakeMatch = True
        if len(rsp) == 6:
            for i in range(0,5):
                if handshake[i] != rsp[i]:
                    handshakeMatch = False
        return handshakeMatch

    def write(self, msg):
        self.ser.write(msg)

    def read(self, maxLength = 256):
        x = self.ser.read(maxLength)
        return x

    def ac_cal_mode(self, calmode):
        calmode_int = int(calmode)
        command = [0xEE, 0xFF, 0x63, 0x00, 0x01, 0x00, calmode_int]
        self.write(command)