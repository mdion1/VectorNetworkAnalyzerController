#ifndef GPIB_COMMS_H
#define GPIB_COMMS_H

#include <QSerialPort>
#include <QDebug>

class GPIB_comms : public QObject
{
    Q_OBJECT
public:
	void setup(QString portname);
    void writeToInstr(const char * msg);
	void writeToInstr(QString msg);
	QByteArray Read_sync();

private slots:
    void Read_async();
	void handleBytesWritten(qint64 bytes);
	void handleError(QSerialPort::SerialPortError serialPortError);

private:
    QSerialPort _serialPort;
	int port_timeout = 2000;
};

#endif // GPIB_COMMS_H
