#include "Constant.h"
#include <math.h>

namespace model {

// Constructor and destructor
Constant::Constant(double inD) {
	d[0] = inD;
	d[1] = inD * inD;
}

Constant::~Constant() {
}

// Get statistics functions

double Constant::getEX() {
	return d[0];
}
double Constant::getELnX() {
	return log(d[0]);
}
double Constant::getEX2() {
	return d[1];
}
double* Constant::getMoments() {
	return d;
}
// Set statistics function

void Constant::setVal(double newValue) {
	d[0] = newValue;
	d[1] = newValue * newValue;
}

// Functions inherited from Node

void Constant::IO(string description, istream& file, bool loadPriors,
		bool debug) {
	int stringSize = 100;
	string input(stringSize, '\0');
	char * pos;

	// Get the description line and ignore it
	file.getline(&input[0], stringSize - 1);
	if (debug) {
		cerr << "Read: " << input << " and ignored\n";
	}

	// Get the value line and use its value
	input.assign(stringSize, '\0');
	file.getline(&input[0], stringSize - 1);
	pos = &input[0];
	if (loadPriors) {
		setVal(strtod(pos, &pos));
	}

	// Output debug info if appropriate.
	if (debug) {
		cerr << "Read: " << input << "\n and value is " << getEX() << "\n";
	}
}

void Constant::IO(string description, ostream& file, bool loadPriors,
		bool debug) {
	file << description << "\n";
	file << getEX() << "\n";
}

double Constant::getBound() {
	return 0.0;
}

void Constant::update(int index1, int index2, int index3, double message[]) {
	return;
}

}
