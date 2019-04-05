#ifndef SIMPLE_SQUIDSTAT_COMMS_H
#define SIMPLE_SQUIDSTAT_COMMS_H

#include <QSerialPort>
#include <QDebug>
#include "c:/Potentiostat/_arduino/retrofit_firmware/global_typedefs.h"

#pragma pack(push, 1)
struct CommandPacket {
	uint16_t frame;
	FramelessComPacketHeader_t hdr;
	char data[0];
};
#pragma pack(pop)


class simpleSquidstatComms : QObject
{
    Q_OBJECT
public:
	void setup(QString portname);
	void send_AC_cal_mode_cmd(uint8_t calmode);

private:
	void sendMessage(PCcommand_t command, int channel, const QByteArray &data);
    QSerialPort _serialPort;
	int port_timeout = 5000;
};

#endif // SIMPLE_SQUIDSTAT_COMMS_H
