#include "simpleSquidstatComms.h"

void simpleSquidstatComms::setup(QString portname)
{
	/* setup serial port connection */
	_serialPort.setPortName(portname);
	_serialPort.open(QIODevice::ReadWrite);

	/* Send handshake to Squidstat */
	sendMessage(HANDSHAKE, 0, QByteArray());
	_serialPort.waitForReadyRead(port_timeout);
	char handshakeResponse_arr[] = { 0xEE, 0xFF, HANDSHAKE_RESPONSE, 0, 0, 0 };
	if (_serialPort.readAll() == QByteArray::fromRawData(handshakeResponse_arr, 6))
		qDebug() << "Handshake received from Squidstat";
	else
		qDebug() << "Connection error with Squidstat";

	//connect(&_serialPort, &QSerialPort::bytesWritten, this, &simpleSquidstatComms::handleBytesWritten);
}

void simpleSquidstatComms::send_AC_cal_mode_cmd(uint8_t calmode)
{
	sendMessage(AC_CAL_MODE, 0, QByteArray(1, calmode));
}

void simpleSquidstatComms::sendMessage(PCcommand_t command, int channel, const QByteArray &data)
{
	QByteArray toSend(sizeof(CommandPacket) + data.size(), 0x00);
	CommandPacket *pack = (CommandPacket*)toSend.data();

	pack->frame = 0xFFEE;
	pack->hdr.command = command;
	pack->hdr.channelNum = channel;
	pack->hdr.dataLength = data.size();
	memcpy(pack->data, data.data(), data.size());

	_serialPort.write(toSend);
}

void simpleSquidstatComms::handleBytesWritten(qint64 bytes)
{
	qDebug() << "bytes written:" << bytes;
}