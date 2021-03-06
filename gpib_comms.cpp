#include "gpib_comms.h"

void GPIB_comms::setup(QString portname)
{
	/* setup serial port connection */
	_serialPort.setPortName(portname);
	_serialPort.open(QIODevice::ReadWrite);

	/* For asynchronous port operation, uncomment the following line */
	//connect(&_serialPort, &QSerialPort::readyRead, this, &GPIB_comms::Read_async, Qt::QueuedConnection);
	//connect(&_serialPort, &QSerialPort::bytesWritten, this, &GPIB_comms::handleBytesWritten, Qt::QueuedConnection);
	connect(&_serialPort, static_cast<void (QSerialPort::*)(QSerialPort::SerialPortError)>(&QSerialPort::error), this, &GPIB_comms::handleError);

	writeToInstr("++addr 16");
}

void GPIB_comms::handleBytesWritten(qint64 bytes)
{
	qDebug() << "bytes written: " << bytes;
}

void GPIB_comms::handleError(QSerialPort::SerialPortError serialPortError)
{
	qDebug() << "error: " << (int)serialPortError;
}

QByteArray GPIB_comms::Read_sync()
{
	_serialPort.waitForReadyRead(port_timeout);
	return _serialPort.readAll();
}

void GPIB_comms::writeToInstr(QString msg)
{
	msg.append("\r");
	_serialPort.write(msg.toUtf8());
	_serialPort.waitForBytesWritten(-1);
	_sleep(25);
}

void GPIB_comms::writeToInstr(const char * msg)
{
    int strLen = 0;
    while (msg[strLen])
    {
        strLen++;
    }
    _serialPort.write(msg, strLen);
	_serialPort.write("\r");
	_serialPort.waitForBytesWritten(-1);
}

void GPIB_comms::Read_async()
{
    qDebug(_serialPort.readAll());
}