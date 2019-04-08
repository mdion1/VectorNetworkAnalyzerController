#include "experimentRunner.h"

experimentRunner::~experimentRunner()
{
	delete outputFile;
}

void experimentRunner::connectSquidstat(const char * squidstatComport)
{
	squidstat.setup(QString(squidstatComport));
}

void experimentRunner::connectVNAnalyzer(const char * squidstatComport)
{
	VNAnalyzer.setup(QString(squidstatComport));
}

void experimentRunner::setupExperiment(const char * inputParamsFile, const char * outputDest)
{
	QList<QStringList> inputParams = QtCSV::Reader::readToList(QString(inputParamsFile), QChar(','));
	for (QStringList line : inputParams)
	{
		if (line.count() == 2)
			experimentParams.insert(line[0], line[1]);
	}

	if (outputDest)
		outputFile = new QFile(QString(outputDest));
}

void experimentRunner::runExperiment()
{
	squidstat.send_AC_cal_mode_cmd(experimentParams[KEY_CAL_MODE_SETTING].toInt());

	VNAnalyzer.writeToInstr("*IDN?");
	qDebug() << VNAnalyzer.Read_sync();

	VNAnalyzer.writeToInstr("NA");	//set network analyzer mode
	VNAnalyzer.writeToInstr("MEAS AB");		//set A/B meas mode

	/* set up sweep */
	VNAnalyzer.writeToInstr("SWPT LOGF ");																	//set sweep type (lin freq, log freq, list freq, power sweep)
	VNAnalyzer.writeToInstr(QString("STAR ") + experimentParams[KEY_STARTING_FREQ] + QString("HZ"));		//set sweep starting frequency in Hz
	VNAnalyzer.writeToInstr(QString("STOP ") + experimentParams[KEY_ENDING_FREQ] + QString("HZ"));			//set sweep ending frequency in Hz
	VNAnalyzer.writeToInstr(QString("POINT ") + experimentParams[KEY_NUM_POINTS]);							//set number of points swept
	
	/* set up sample averaging */
	//VNAnalyzer.writeToInstr("AVER ON");
	//VNAnalyzer.writeToInstr(QString("AVERFACT ") + experimentParams[KEY_AVERAGING_NUM]);

	/* set up signal strength */
	double amp_dBm = convertVoltsToDBM(experimentParams[KEY_VOLTAGE_AMPLITUDE].toDouble());
	QString amp_dBm_str;
	amp_dBm_str.setNum(amp_dBm, 'f', 1);
	VNAnalyzer.writeToInstr(QString("POWE ") + amp_dBm_str);

	/* trigger sweep */
	VNAnalyzer.writeToInstr("SING");										//trigger single sweep

	/* Wait for sweep to complete */
	do
	{
		/* Read Event Status Register B, bit 0, to detect sweep completion */
		VNAnalyzer.writeToInstr("ESB?");
		int statusReg = VNAnalyzer.Read_sync().toInt();
		if (statusReg & 1)
			break;
		_sleep(25);
	} while (true);

	/* Print out values */
	VNAnalyzer.writeToInstr("FORM4");
	VNAnalyzer.writeToInstr("OUTPRAW1?");
	qDebug() << VNAnalyzer.Read_sync();

	/* Read out values to array */
	VNAnalyzer.writeToInstr("FOMR3");
	VNAnalyzer.writeToInstr("OUTPRAW1?");
	auto numList = parseRawBytes(VNAnalyzer.Read_sync());
	auto pointList = NumListToComplexNumList(numList);
	auto freqList = getFreqList(pointList.count());

	/* Normalize values */
	pointList = normalizeListToLastPoint(pointList);

	/* Write values to output file */
	if (outputFile)
	{
		if (outputFile->open(QIODevice::WriteOnly))
		{
			QTextStream out(outputFile);
			out.setRealNumberPrecision(12);
			out << "Frequency,Mag,Phase\n";
			for (int i = 0; i < pointList.count(); i++)
			{
				out << freqList[i];			out << ",";
				out << pointList[i].mag;	out << ",";
				out << pointList[i].phase;	out << "\n";
			}
			outputFile->close();
		}
	}
}

QList<experimentRunner::complexNum> experimentRunner::normalizeListToLastPoint(QList<experimentRunner::complexNum> raw)
{
	QList<experimentRunner::complexNum> ret;
	experimentRunner::complexNum baseline;
	baseline.mag = raw.last().mag;
	baseline.phase = raw.last().phase;
	for (int i = 0; i < raw.count(); i++)
	{
		experimentRunner::complexNum x;
		x.mag = raw[i].mag / baseline.mag;
		x.phase = raw[i].phase - baseline.phase;
		x.real = x.mag * cos(x.phase * 3.14159265358979323846264338327950288 / 180);
		x.imag = x.mag * sin(x.phase * 3.14159265358979323846264338327950288 / 180);
		ret.append(x);
	}
	return ret;
}

QList<double> experimentRunner::parseRawBytes(QByteArray raw, bool reverseEndianness)
{
	QList<double> ret;

	/* If endianness is reversed */
	if (reverseEndianness)
	{
		for (int i = 0; i < raw.count() / sizeof(double); i += sizeof(double))
		{
			QByteArray reverse;
			for (int j = 0; j < sizeof(double); j++)
			{
				reverse.append(raw[i + sizeof(double) - j]);
			}
			double * dblPtr = (double *)reverse.data() + i;
			ret.append(*dblPtr);
		}
	}
	/* If endianness is ok */
	else
	{
		for (int i = 0; i < raw.count() / sizeof(double); i += sizeof(double))
		{
			double * dblPtr = (double *)raw.data() + i;
			ret.append(*dblPtr);
		}
	}

	return ret;
}

QList<experimentRunner::complexNum> experimentRunner::NumListToComplexNumList(QList<double> raw)
{
	QList<experimentRunner::complexNum> ret;
	for (int i = 0; i < raw.length() / 2; i += 2)
	{
		experimentRunner::complexNum x;
		x.real = raw[i];
		x.imag = raw[i + 1];
		x.mag = sqrt(x.real * x.real + x.imag * x.imag);
		x.phase = atan2(x.imag, x.real) * 180 / 3.14159265358979323846264338327950288;
		ret.append(x);
	}
	return ret;
}

QList<double> experimentRunner::getFreqList(int listMinLength)
{
	QList<double> ret;
	double fStart = experimentParams[KEY_STARTING_FREQ].toDouble();
	double fEnd = experimentParams[KEY_ENDING_FREQ].toDouble();
	int numPoints = experimentParams[KEY_NUM_POINTS].toInt();

	double fstep = pow(fEnd / fStart, 1. / numPoints);
	for (int i = 0; i < listMinLength; i++)
	{
		int j = i < numPoints ? i : numPoints;		//in case listMinLength > numPoints
		ret.append(fStart * pow(fstep, j));
	}
	return ret;
}

double experimentRunner::convertVoltsToDBM(double volts)
{
	double dBm = 10 * log10(volts * volts / 50);
	return dBm;
}