#pragma once

#include "qtcsv/reader.h"
#include <qmap.h>
#include "gpib_comms.h"
#include "simpleSquidstatComms.h"
#include <qfile.h>
#include <cmath>

#define KEY_CAL_MODE_SETTING "AC_CAL_MODE"
#define KEY_VOLTAGE_AMPLITUDE "VoltageAmplitude"
#define KEY_STARTING_FREQ "StartingFreq"
#define KEY_ENDING_FREQ "EndingFreq"
#define KEY_NUM_POINTS "NumPoints"
#define KEY_AVERAGING_NUM	"AveragingNum"

class experimentRunner
{
	typedef struct
	{
		double real;
		double imag;
		double mag;
		double phase;
	}complexNum;

public:
	~experimentRunner();
	void setupExperiment(const char * inputParamsFile, const char * outputDest = nullptr);
	void runExperiment();
	void connectVNAnalyzer(const char * squidstatComport);
	void connectSquidstat(const char * squidstatComport);

	void testFunction();

private:
	QMap<QString, QString> experimentParams;
	QFile * outputFile = nullptr;
	simpleSquidstatComms squidstat;
	GPIB_comms VNAnalyzer;
	QList<double> getFreqList(int listMinLength);
	void normalizeListToLastPoint(QList<complexNum> &data);
	void reverseList(QList<complexNum> &list);
	QList<complexNum> averageList(QList<QList<complexNum>> &masterList);
	static QList<double> parseRawBytes(QByteArray raw, bool reverseEndianness = true);
	static QList<complexNum> NumListToComplexNumList(QList<double> raw);
	static double convertVoltsToDBM(double volts);
};

