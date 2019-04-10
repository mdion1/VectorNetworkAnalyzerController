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

	VNAnalyzer.writeToInstr("HOLD");		//takes instrument out of free-run mode


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
	VNAnalyzer.writeToInstr(QString("POIN ") + experimentParams[KEY_NUM_POINTS]);							//set number of points swept
	
	/* set up sample averaging */
	VNAnalyzer.writeToInstr("AVER 1");			//turns on sample averaging
	VNAnalyzer.writeToInstr("AVERREST");		//resets the averaging calculations
	VNAnalyzer.writeToInstr(QString("AVERFACT ") + experimentParams[KEY_AVERAGING_NUM]);		//specifies the number of samples to average

	/* set up signal strength */
	double amp_dBm = convertVoltsToDBM(experimentParams[KEY_VOLTAGE_AMPLITUDE].toDouble());
	QString amp_dBm_str;
	amp_dBm_str.setNum(amp_dBm, 'f', 1);
	VNAnalyzer.writeToInstr(QString("POWE ") + amp_dBm_str);

//#define MANUAL_AVERAGING
#ifdef MANUAL_AVERAGING
	QList<QList<experimentRunner::complexNum>> dataMasterList;
	for (int i = 0; i < experimentParams[KEY_AVERAGING_NUM].toInt(); i++)
	{
		/* trigger sweep */
		VNAnalyzer.writeToInstr("SING");										//trigger single sweep

		/* Wait for sweep to complete */
		do
		{
			/* Read Event Status Register B, bit 0, to detect sweep completion */
			VNAnalyzer.writeToInstr("ESB?");
			if (VNAnalyzer.Read_sync() == "+1\n")
				break;
			_sleep(25);
		} while (true);

		/* Read out values to array */
		VNAnalyzer.writeToInstr("OUTPRAW1?");
		auto resp = VNAnalyzer.Read_sync();
		resp.remove(0, 8); resp.remove(resp.count() - 1, 1);
		auto numList = parseRawBytes(resp);
		auto pointList = NumListToComplexNumList(numList);
		
		dataMasterList << pointList;
	}

	/* Average, reverse, and normalize values */
	auto avgData = averageList(dataMasterList);
	reverseList(avgData);
	normalizeListToLastPoint(avgData);
	auto freqList = getFreqList(avgData.count());

	/* Write values to output file */
	if (outputFile)
	{
		if (outputFile->open(QIODevice::WriteOnly))
		{
			QTextStream out(outputFile);
			out.setRealNumberPrecision(12);
			out << "Frequency,Mag,Phase\r\n";
			for (int i = 0; i < avgData.count(); i++)
			{
				out << freqList[i];			out << ",";
				out << avgData[i].mag;		out << ",";
				out << avgData[i].phase;	out << "\r\n";
			}
			outputFile->close();
		}
	}

#else
	/* trigger sweeps */
	VNAnalyzer.writeToInstr(QString("NUMG ") + experimentParams[KEY_AVERAGING_NUM]);										//triggers multiple sweeps before returning to hold mode
	
	/* Wait for sweep to complete */
	do
	{
		/* Read Event Status Register B, bit 0, to detect sweep completion */
		VNAnalyzer.writeToInstr("ESB?");
		if (VNAnalyzer.Read_sync() == "+1\n")
			break;
		_sleep(25);
	} while (true);

	/* Read out values to array */
	VNAnalyzer.writeToInstr("OUTPDATA?");
	auto resp = VNAnalyzer.Read_sync();
	resp.remove(0, 8); resp.remove(resp.count() - 1, 1);
	auto numList = parseRawBytes(resp);
	auto pointList = NumListToComplexNumList(numList);

	/* Average, reverse, and normalize values */
	reverseList(pointList);
	normalizeListToLastPoint(pointList);
	auto freqList = getFreqList(pointList.count());

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

#endif
}

QList<experimentRunner::complexNum> experimentRunner::averageList(QList<QList<experimentRunner::complexNum>> &masterList)
{
	QList<experimentRunner::complexNum> ret;
	int N = masterList.count();
	
	for (int i = 0; i < masterList.first().count(); i++)
	{
		experimentRunner::complexNum sum = { 0, 0, 0, 0 };
		for (int j = 0; j < N; j++)
		{
			sum.mag += masterList[j][i].mag;
			sum.phase += masterList[j][i].phase;
		}
		sum.mag /= N;
		sum.phase /= N;
		sum.real = sum.mag * cos(sum.phase * 3.14159265358979323846264338327950288 / 180);
		sum.imag = sum.mag * sin(sum.phase * 3.14159265358979323846264338327950288 / 180);
		ret << sum;
	}
	return ret;
}

void experimentRunner::reverseList(QList<experimentRunner::complexNum> &list)
{
	for (int k = 0; k < (list.size() / 2); k++) list.swap(k, list.size() - (1 + k));
}

void experimentRunner::normalizeListToLastPoint(QList<experimentRunner::complexNum> &data)
{
	experimentRunner::complexNum baseline;
	baseline.mag = data.last().mag;
	baseline.phase = data.last().phase;
	for (int i = 0; i < data.count(); i++)
	{
		experimentRunner::complexNum x;
		data[i].mag /= baseline.mag;
		data[i].phase -= baseline.phase;
		data[i].real = data[i].mag * cos(data[i].phase * 3.14159265358979323846264338327950288 / 180);
		data[i].imag = data[i].mag * sin(data[i].phase * 3.14159265358979323846264338327950288 / 180);
	}
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