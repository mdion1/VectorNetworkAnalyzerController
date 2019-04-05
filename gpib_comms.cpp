#include "gpib_comms.h"

#include <QDebug>

GPIB_comms::GPIB_comms()
{
    //this->moveToThread(this);
    _serialPort = new QSerialPort(this);
}

void GPIB_comms::testFn()
{
	_serialPort->setPortName(_portname);
	_serialPort->open(QIODevice::ReadWrite);
	//connect(_serialPort, &QSerialPort::readyRead, this, &GPIB_comms::DataArrived, Qt::QueuedConnection);

	connect(_serialPort, &QSerialPort::bytesWritten, this, &GPIB_comms::testHandlerFn);
}

void GPIB_comms::syncReadTest()
{
	_serialPort->waitForReadyRead(5000);
	auto readresult = _serialPort->readAll();
	qDebug() << readresult;
}

void GPIB_comms::run()
{
    //_serialPort->setPortName(_portname);
    //connect(_serialPort, &QSerialPort::readyRead, this, &GPIB_comms::DataArrived, Qt::QueuedConnection);

    _serialPort->open(QIODevice::ReadWrite);
    exec();

    _serialPort->close();
    _serialPort->deleteLater();
}

void GPIB_comms::writeToInstr(const char * msg)
{
    int strLen = 0;
    while (msg[strLen])
    {
        strLen++;
    }
    _serialPort->write(msg, strLen);

	_serialPort->write("\r\n");
}

void GPIB_comms::DataArrived()
{
    qDebug(_serialPort->readAll());
}

void GPIB_comms::testHandlerFn()
{
	qDebug() << "written";
}
