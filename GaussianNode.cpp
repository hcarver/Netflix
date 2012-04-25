#include "GaussianNode.h"

namespace model {

// Constructor and Destructor
GaussianNode::GaussianNode(double mn, double prec) {
	mean = mn;
	precision = prec;
	initialiseMoments();
}

GaussianNode::~GaussianNode() {
}

// Accessors to prior and moments
double GaussianNode::getMean() {
	return mean;
}

double GaussianNode::getPrecision() {
	return precision;
}

// Functions inherited from EXandX2
double GaussianNode::getEX() {
	return moments[0];
}

double GaussianNode::getEX2() {
	return moments[1];
}

double* GaussianNode::getMoments() {
	return &moments[0];
}

// Functions inherited from Node
double GaussianNode::getBound() {
	double b = 0.0;

	// First, get standard bound calculation of (prior-parameter)*moment
	double prior[2] = { 0.0, 0.0 };
	updatePrior(prior);
	b += (prior[0] - parameters[0]) * moments[0]
			+ (prior[1] - parameters[1]) * moments[1];

	// Add PriorG component
	b += 0.5 * (log(precision) - precision * mean * mean - LOG2PI);

	// Add PostG component
	// Note that parameters[0] / (-2 * parameters[1]) is the mean
	// Note that -2 * parameters[1] is the precision
	b -= 0.5
			* (log(-2 * parameters[1])
					- parameters[0] * parameters[0] / (-2 * parameters[1])
					- LOG2PI);

	return b;
}

void GaussianNode::update(int index1, int index2, int index3,
		double message[]) {
	// First fill parameters with the correct values, using prior and message
	updatePrior(parameters);
	for (int j = 0; j < 2; j++) {
		parameters[j] += message[j];
	}

	// Then set moments using parameters. See notes above about values of parameters
	double variance = -0.5 / (parameters[1]);
	setMoments(parameters[0] * variance,
			parameters[0] * variance * parameters[0] * variance + variance);
}

void GaussianNode::IO(string description, ostream& file, bool loadPriors,
		bool debug) {
	file << description << "\n";
	file << "Prec" << "\n";
	file << precision << "\n";
	file << getEX() << "," << getEX2() << "\n";
}

void GaussianNode::IO(string description, istream& file, bool loadPriors,
		bool debug) {
	int stringSize = 5000;
	string input(stringSize, '\0');
	char * pos;
	double mom1;

	// Description
	file.getline(&input[0], stringSize - 1);
	if (debug) {
		cerr << "Read: " << input << " and ignored\n";
	}

	// "Prec" line
	input.assign(stringSize, '\0');
	file.getline(&input[0], stringSize - 1);
	if (debug) {
		cerr << "Read: " << input << " and ignored\n";
	}

	// Precision value line
	input.assign(stringSize, '\0');
	file.getline(&input[0], stringSize - 1);
	if (loadPriors) {
		precision = strtod(&input[0], &pos);
	}
	if (debug) {
		cerr << "Read: " << input << " and precision is " << precision << "\n";
	}

	// Moments value line
	input.assign(stringSize, '\0');
	file.getline(&input[0], stringSize - 1);
	mom1 = strtod(&input[0], &pos);
	setMoments(mom1, strtod(&pos[1], &pos));
	if (debug) {
		cerr << "Read: " << input << "\n and moments are " << getEX() << ","
				<< getEX2() << "\n";
	}
}

void GaussianNode::setMoments(double val1, double val2) {
	//	moments[inOne % maxInOneDim][(inOne / maxInOneDim) * dim2 + inTwo][0] = val1;
	//	moments[inOne % maxInOneDim][(inOne / maxInOneDim) * dim2 + inTwo][1] = val2;
	moments[0] = val1;
	moments[1] = val2;

	parameters[1] = 1.0 / ((moments[1] - moments[0] * moments[0]) * -2);

	parameters[0] = moments[0] * -2 * parameters[1];
}

bool GaussianNode::initialiseMoments() {
	srand(time(0));

	double variance = 1.0 / precision;

	double x1 = ((double) rand()) / (double) RAND_MAX;
	double x2 = ((double) rand()) / (double) RAND_MAX;

	double y1 = sqrt(-2.0 * log(x1)) * cos(2.0 * PI * x2);
	double y2 = sqrt(-2.0 * log(x1)) * sin(2.0 * PI * x2);

	double value1 = mean + y1 * sqrt(variance);
	setMoments(value1, value1 * value1 + variance);

	return true;
}

//double GaussianNode::getPriorG()
//{
//    return 0.5 * (log(precision) - precision * mean*mean - LOG2PI);
//}
//double GaussianNode::getPostG(double reuseable[])
//{
//    // first recover the parameters from the moments
//    // double[] theseMoms = getMoments(inOne, inTwo);
//    //    double theseParams[2] = {0,0};
//
//    reuseable[0] = parameters[0];
//    reuseable[1] = parameters[1];
//
//    double gamma = -2 * reuseable[1];
//    //    if (gamma <= 0)
//    //    {
//    //        cerr << "-ve gamma in getPostG";
//    //    }
//    return 0.5 * (log(gamma) - reuseable[0] * reuseable[0] / gamma - LOG2PI);
//}

// ************** UPDATES ***************
void GaussianNode::updatePrior(double priorParameters[]) {
	priorParameters[0] = precision * mean;
	priorParameters[1] = -0.5 * precision;
}
//void GaussianNode::setVal(double val)
//{
//    setMoments(val, val*val+exp(-log(precision)));
//}
//void GaussianNode::updateMoments()
//{
//
//    // double[] paras = getParams(pos1, pos2);
//    double variance = -0.5 / (parameters[1]);
//    // System.err.println("variance="+variance);
//    // if (variance<0) throw new IllegalStateException("Negative variance: "+variance);
//    setMoments(parameters[0] * variance, parameters[0] * variance * parameters[0]
//               * variance + variance);
//}

}
