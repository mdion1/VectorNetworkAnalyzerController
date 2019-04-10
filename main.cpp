#include <QCoreApplication>
#include "experimentRunner.h"

int main(int argc, char *argv[])
{
	//QCoreApplication a(argc, argv);

	experimentRunner experiment;
	if (argc == 4)
	{
		experiment.connectSquidstat(argv[1]);
		experiment.connectVNAnalyzer(argv[2]);
		experiment.setupExperiment(argv[3]);
	}
	else if (argc == 5)
	{
		experiment.connectSquidstat(argv[1]);
		experiment.connectVNAnalyzer(argv[2]);
		experiment.setupExperiment(argv[3], argv[4]);
	}
	else
		return -1;

	experiment.runExperiment();
	//experiment.testFunction();

	//system("pause");
	//a.exit(0);
    //return a.exec();
}
