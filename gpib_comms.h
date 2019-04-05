#ifndef GPIB_COMMS_H
#define GPIB_COMMS_H

#include <QSerialPort>
#include <QThread>

class GPIB_comms : public QThread
{
    Q_OBJECT
public:
    GPIB_comms();
    void run();
    void setPortName(QString portname) { _portname = portname; }
    void writeToInstr(const char * msg);
	void testFn();
	void syncReadTest();

private slots:
    void DataArrived();
	void testHandlerFn();

private:
    QSerialPort * _serialPort;
    QString _portname;
};

#endif // GPIB_COMMS_H
