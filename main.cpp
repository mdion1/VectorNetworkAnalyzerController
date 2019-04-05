#include <QCoreApplication>
#include "gpib_comms.h"
#include "simpleSquidstatComms.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    GPIB_comms instr;
    instr.setup("COM15");
	instr.writeToInstr("R0T1X");
	qDebug() <<	instr.Read_sync();

	simpleSquidstatComms squidstat;
	squidstat.setup("COM8");

	system("pause");
	return 0;
    //return a.exec();
}
