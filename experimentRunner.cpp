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

void experimentRunner::testFunction()
{
	int x = experimentParams[KEY_CAL_MODE_SETTING].toInt();
	squidstat.send_AC_cal_mode_cmd(x);
}

void experimentRunner::runExperiment()
{
	squidstat.send_AC_cal_mode_cmd(experimentParams[KEY_CAL_MODE_SETTING].toInt());

	//VNAnalyzer.writeToInstr(QString("PRES"));
	VNAnalyzer.writeToInstr(QString("*IDN?"));
	qDebug() << VNAnalyzer.Read_sync();

	VNAnalyzer.writeToInstr("FORM3");		//sets data output mode to 64-bit double
	VNAnalyzer.writeToInstr("NA");	//set network analyzer mode
	//VNAnalyzer.writeToInstr("NA?");
	//qDebug() << VNAnalyzer.Read_sync();
	VNAnalyzer.writeToInstr("MEAS AB");		//set A/B meas mode
	//VNAnalyzer.writeToInstr("MEAS?");
	//qDebug() << VNAnalyzer.Read_sync();

	/* set up sweep */
	VNAnalyzer.writeToInstr("SWPT LOGF");																	//set sweep type (lin freq, log freq, list freq, power sweep)
	//VNAnalyzer.writeToInstr("SWPT?");
	//qDebug() << VNAnalyzer.Read_sync();

	VNAnalyzer.writeToInstr("BWAUTO 1");

	VNAnalyzer.writeToInstr(QString("STAR ") + experimentParams[KEY_STARTING_FREQ] + QString("HZ"));		//set sweep starting frequency in Hz
	VNAnalyzer.writeToInstr(QString("STOP ") + experimentParams[KEY_ENDING_FREQ] + QString("HZ"));			//set sweep ending frequency in Hz
	VNAnalyzer.writeToInstr(QString("POINT ") + experimentParams[KEY_NUM_POINTS]);							//set number of points swept
	
	/* set up sample averaging */
	VNAnalyzer.writeToInstr("AVER 1");
	VNAnalyzer.writeToInstr(QString("AVERFACT ") + experimentParams[KEY_AVERAGING_NUM]);

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
		auto resp = VNAnalyzer.Read_sync();
		if (resp == "+1\n")
			break;
		//auto x = VNAnalyzer.Read_sync();
		//x.remove(0, 1);
		//int statusReg = x.toInt();//int statusReg = VNAnalyzer.Read_sync().toInt();
		//if (statusReg & 1)
		//	break;
		_sleep(25);
	} while (true);

	/* Print out values */
	//VNAnalyzer.writeToInstr("FORM4");
	//VNAnalyzer.writeToInstr("OUTPRAW1?");
	//qDebug() << VNAnalyzer.Read_sync();

	/* Read out values to array */
	VNAnalyzer.writeToInstr("OUTPRAW1?");

	auto resp = VNAnalyzer.Read_sync();
	//qDebug() << resp;
	resp.remove(0, 8); resp.remove(resp.count() - 1, 1);
	auto numList = parseRawBytes(resp);
	auto pointList = NumListToComplexNumList(numList);
	auto freqList = getFreqList(pointList.count());

	/* Reverse and normalize values */
	reverseList(pointList);
	pointList = normalizeListToLastPoint(pointList);

	/* Write values to output file */
	if (outputFile)
	{
		if (outputFile->open(QIODevice::WriteOnly))
		{
			QTextStream out(outputFile);
			out.setRealNumberPrecision(12);
			out << "Frequency,Mag,Phase\r\n";
			for (int i = 0; i < pointList.count(); i++)
			{
				out << freqList[i];			out << ",";
				out << pointList[i].mag;	out << ",";
				out << pointList[i].phase;	out << "\r\n";
			}
			outputFile->close();
		}
	}
}

void experimentRunner::reverseList(QList<experimentRunner::complexNum> &list)
{
	for (int k = 0; k < (list.size() / 2); k++) list.swap(k, list.size() - (1 + k));
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
		for (int i = 0; i < raw.count(); i += sizeof(double))
		{
			QByteArray reverse;
			for (int j = 7; j >=0; j--)
			{
				reverse.append(raw[i + j]);
			}
			double * dblPtr = (double *)reverse.data();
			ret.append(*dblPtr);
		}
	}
	/* If endianness is ok */
	else
	{
		for (int i = 0; i < raw.count(); i += sizeof(double))
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
	for (int i = 0; i < raw.length(); i += 2)
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
	double fStart = max(experimentParams[KEY_STARTING_FREQ].toDouble(), experimentParams[KEY_ENDING_FREQ].toDouble());
	double fEnd = min(experimentParams[KEY_STARTING_FREQ].toDouble(), experimentParams[KEY_ENDING_FREQ].toDouble());
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
	double dBm = 10 * log10(1000 * volts * volts / 50);
	return dBm;
}