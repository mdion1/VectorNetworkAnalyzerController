#include "gpib_comms.h"

void GPIB_comms::setup(QString portname)
{
	/* setup serial port connection */
	_serialPort.setPortName(portname);
	_serialPort.open(QIODevice::ReadWrite);

	/* For asynchronous port operation, uncomment the following line */
	//connect(_serialPort, &QSerialPort::readyRead, this, &GPIB_comms::DataArrived, Qt::QueuedConnection);

	/* setup VNA settings */

}

QByteArray GPIB_comms::Read_sync()
{
	_serialPort.waitForReadyRead(port_timeout);
	return _serialPort.readAll();
}

void GPIB_comms::writeToInstr(QString msg)
{
	//	msg.append("\r\n"); //todo: figure out if "\r\n" is needed
	_serialPort.write(msg.toUtf8());
}

void GPIB_comms::writeToInstr(const char * msg)
{
    int strLen = 0;
    while (msg[strLen])
    {
        strLen++;
    }
    _serialPort.write(msg, strLen);

	//_serialPort.write("\r\n");		//todo: determine if "\r\n" is needed
}

void GPIB_comms::Read_async()
{
    qDebug(_serialPort.readAll());
}