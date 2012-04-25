#include "GammaNode.h"

namespace model {

// Constructor and destructor
GammaNode::GammaNode(double xiShape, double xiScale) {
	gamma_series[0] = 76.18009172947146;
	gamma_series[1] = -86.50532032941677;
	gamma_series[2] = 24.01409824083091;
	gamma_series[3] = -1.231739572450155;
	gamma_series[4] = 0.1208650973866179e-2;
	gamma_series[5] = -0.5395239384953e-5;

	shape = xiShape;
	scale = xiScale;
	initialiseMoments();
}

GammaNode::~GammaNode() {
}

// Moments-related functions. 
void GammaNode::setMoments(double EX, double ELnX) {
	moments[0] = EX;
	moments[1] = ELnX;
}

void GammaNode::initialiseMoments() {
	moments[0] = shape / scale;
	moments[1] = diGamma(shape) - log(scale);
	parameters[0] = -scale;
	parameters[1] = shape - 1;
}

void GammaNode::updatePrior(double priorParameters[]) {
	priorParameters[0] = -scale;
	priorParameters[1] = shape - 1;
}

// Functions inherited from EXandLnX
double GammaNode::getEX() {
	return moments[0];
}

double GammaNode::getELnX() {
	return moments[1];
}

// Functions inherited from Node
double GammaNode::getBound() {
	double b = 0.0;

	// First get the prior
	double prior[2] = { 0.0, 0.0 };
	updatePrior(prior);

	// Add the standard (prior-parameter)*moment expressions
	b += (prior[0] - parameters[0]) * moments[0]
			+ (prior[1] - parameters[1]) * moments[1];
	// Add PriorG
	b += shape * log(scale) - logGamma(shape);
	// Subtract PostG
	b -= (parameters[1] + 1) * log(-parameters[0])
			- logGamma(parameters[1] + 1);

	return b;
}

void GammaNode::update(int index1, int index2, int index3, double message[]) {
	// Fill in parameters from prior and message
	updatePrior(parameters);
	for (int j = 0; j < 2; j++) {
		parameters[j] += message[j];
	}

	// Set moments
	double shapeTemp = parameters[1] + 1;
	double scaleTemp = -parameters[0];
	moments[0] = shapeTemp / scaleTemp;
	moments[1] = diGamma(shapeTemp) - log(scaleTemp);
}

void GammaNode::IO(string description, ostream& file, bool loadPriors,
		bool debug) {
	file << description << "\n";
	file << "Prior" << "\n";
	file << shape << "," << scale << "\n";
	file << getEX() << "," << getELnX() << "\n";
}

void GammaNode::IO(string description, istream& file, bool loadPriors,
		bool debug) {
	int stringSize = 100;
	string input(stringSize, '\0');
	char * pos;
	double readValue;

	//description line
	file.getline(&input[0], stringSize - 1);
	if (debug) {
		cerr << "Read: " << input << " and ignored\n";
	}

	//"Prior" line
	input.assign(stringSize, '\0');
	file.getline(&input[0], stringSize - 1);
	if (debug) {
		cerr << "Read: " << input << " and ignored\n";
	}

	//Prior value line.
	// Note that priors are only used to describe constants - so the moments of
	// the gamma node aren't considered priors, but the shape and scale are.
	input.assign(stringSize, '\0');
	file.getline(&input[0], stringSize - 1);
	if (loadPriors) {
		shape = strtod(&input[0], &pos);
		scale = strtod(&pos[1], &pos);
	}
	if (debug) {
		cerr << "Read: " << input << " and shape is " << shape
				<< " and scale is " << scale << "\n";
	}

	// Moments line
	input.assign(stringSize, '\0');
	file.getline(&input[0], stringSize - 1);
	readValue = strtod(&input[0], &pos);
	setMoments(readValue, strtod(&pos[1], &pos));
	if (debug) {
		cerr << "Read: " << input << "\n and moments are " << getEX() << ","
				<< getELnX() << "\n";
	}
}

// Helper functions to compute gamma-related function values.
// These are adapted from functions in Tom Minka's Lightspeed library
double GammaNode::diGamma(double x) {
	double y = 0.0;

	if (x <= 0.0) {
		cerr << "x <= 0.0 in diGamma function\n";
	}
	// Use approximation if argument <= small.
	else if (x <= small) {
		y = d1 - 1.0 / x + d2 * x;
	}
	// Else use a recurrence formula to get x greater than the large bound
	// Reduce to digamma(X + N) where (X + N) >= large
	else {
		while (x < large) {
			y = y - 1.0 / x;
			x = x + 1.0;
		}
		// Now use deMoivre's expansion
		double r = 1.0 / x;
		y = y + log(x) - 0.5 * r;
		r = r * r;
		y = y - r * (s3 - r * (s4 - r * (s5 - r * (s6 - r * s7))));
	}

	return y;
}

double GammaNode::logGamma(double x) {
	if (x <= 0.0) {
		cerr << "x <= 0.0 in logGamma function\n";
	}

	double denom = x + 1;
	double x1 = x + 5.5;
	double series = 1.000000000190015;

	for (int i = 0; i < 6; i++) {
		series += gamma_series[i] / denom;
		denom += 1.0;
	}
	return (M_lnSqrt2PI + (x + 0.5) * log(x1) - x1 + log(series / x));
}

}
