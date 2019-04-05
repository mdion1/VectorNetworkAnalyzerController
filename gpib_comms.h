#ifndef GPIB_COMMS_H
#define GPIB_COMMS_H

#include <QSerialPort>
#include <QDebug>

class GPIB_comms : QObject
{
    Q_OBJECT
public:
	void setup(QString portname);
    void writeToInstr(const char * msg);
	QByteArray Read_sync();

private slots:
    void Read_async();

private:
    QSerialPort _serialPort;
	int port_timeout = 5000;
};

#endif // GPIB_COMMS_H
