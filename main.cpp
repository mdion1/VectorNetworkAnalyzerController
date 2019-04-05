#include <QCoreApplication>
#include "gpib_comms.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    GPIB_comms instr;
    instr.setPortName("COM15");
	instr.testFn();
    //instr.start();
	
    instr.writeToInstr("R0T1X");
	instr.syncReadTest();


    return a.exec();
}
