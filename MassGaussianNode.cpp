#include "MassGaussianNode.h"
#include "GaussianNode.h"
#include "GammaNode.h"

namespace model {

// Constructor and Destructor
MassGaussianNode::MassGaussianNode(int size1, int size2, EXandX2* xiMean,
		EXandLnX* xiPrecision) {
	// Set the mean and precision, and create the moments arrays
	mean = xiMean;
	precision = xiPrecision;
	dim1 = size1;
	dim2 = size2;

	moments = new double**[size1];
	for (int i = 0; i < size1; i++) {
		moments[i] = new double *[size2];
		for (int j = 0; j < size2; j++) {
			moments[i][j] = new double[2];
		}
	}

	// Randomly initialise the moments
	initialiseMoments();
}

MassGaussianNode::~MassGaussianNode() {
	// Delete all the nested arrays and the mean and precision nodes
	for (int ii = 0; ii < dim1; ii++) {
		for (int jj = 0; jj < dim2; jj++) {
			delete[] (moments[ii][jj]);
		}
		delete[] (moments[ii]);
	}

	delete[] moments;
	delete mean;
	delete precision;
}

// Functions to update parents
void MassGaussianNode::updateMeanParent() {
	double message[2] = { 0.0, 0.0 };
	// Message[0] for each element is precision->getEX() * getEX()
	for (int ii = 0; ii < dim1; ii++) {
		for (int jj = 0; jj < dim2; jj++) {
			message[0] += getMoments(ii, jj)[0];
		}
	}
	message[0] *= precision->getEX();
	// Message[1] is constant for each element.
	message[1] = ((double) dim1) * ((double) dim2) * -0.5 * precision->getEX();

	// Now update the mean with the message
	mean->update(0, 0, 0, message);
}

void MassGaussianNode::updatePrecisionParent() {
	double message[2] = { 0.0, 0.0 };
	double meanMean = mean->getEX();
	double meanEX2 = mean->getEX2();
	double *lMoments;

	// Message[0] for one element is given by EX() * mean.getEX() - 0.5 * (EX2() + mean.getEX2())
	for (int ii = 0; ii < dim1; ii++) {
		for (int jj = 0; jj < dim2; jj++) {
			lMoments = getMoments(ii, jj);
			message[0] += lMoments[0] * meanMean
					- 0.5 * (lMoments[1] + meanEX2);
		}
	}

	// Message[1] for one element is 0.5
	message[1] = ((double) dim1) * ((double) dim2) * 0.5;

	// Now perform the update
	precision->update(0, 0, 0, message);
}

// Accessors for moments
double* MassGaussianNode::getMoments(int inOne, int inTwo) {
	return moments[inOne][inTwo];
}

double MassGaussianNode::getEX(int pos1, int pos2) {
	return moments[pos1][pos2][0];
}

double MassGaussianNode::getEX2(int pos1, int pos2) {
	return moments[pos1][pos2][1];
}

// Accessors for mean and precision
EXandX2* MassGaussianNode::getMean() {
	return mean;
}

EXandLnX* MassGaussianNode::getPrecision() {
	return precision;
}

// Function to get the vital statistics of the distribution
void MassGaussianNode::getVitalStatistics(double &xiPrior11, double &xiPrior12,
		double &xiPrior21, double &xiPrior22, double &xiMean,
		double &xiVarOfMean, double &xiAverageVariance) {
	// Priors first because they're easy
	xiPrior11 = mean->getEX();
	xiPrior12 = mean->getEX2();
	xiPrior21 = precision->getEX();
	xiPrior22 = precision->getELnX();

	// Then the computed statistics
	double obs = (double) (dim1 * dim2);
	double meanVal = 0.0;
	double meanValSquared = 0.0;
	double averageVariance = 0.0;

	for (int ii = 0; ii < dim1; ii++) {
		for (int jj = 0; jj < dim2; jj++) {
			double* thisMoms = getMoments(ii, jj);
			meanVal += thisMoms[0];
			meanValSquared += thisMoms[0] * thisMoms[0];
			averageVariance += thisMoms[1] - thisMoms[0] * thisMoms[0];
		}
	}

	xiMean = meanVal / obs;
	xiVarOfMean = (meanValSquared / obs) - xiMean * xiMean;
	xiAverageVariance = averageVariance / obs;
}

// Functions inherited from Node 
double MassGaussianNode::getBound() {
	double b = 0.0;

	// Get prior (constant over matrix)
	double prior[2] = { 0.0, 0.0 };
	updatePrior(prior);

	// For each matrix element, add the standard bound component and subtract the postG
	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim2; j++) {
			double params[2];
			double* moms;

			// Get parameters and moments
			getParams(i, j, params);
			moms = getMoments(i, j);

			// Add the standard (prior - parameter) * moment component
			b += (prior[0] - params[0]) * moms[0]
					+ (prior[1] - params[1]) * moms[1];

			// Subtract the PostG for this component
			double gamma = -2 * params[1];
			b -= 0.5 * (log(gamma) - params[0] * params[0] / gamma - LOG2PI);
		}
	}

	// Add priorG (constant over the matrix)
	b += ((double) dim1) * ((double) dim2) * 0.5
			* (precision->getELnX() - precision->getEX() * mean->getEX2()
					- LOG2PI);

	return b;
}

void MassGaussianNode::IO(string description, ostream& file, bool loadPriors,
		bool debug) {
	double * stats = 0;
	// Output a descriptive line
	file << description << " " << dim1 << ":" << dim2 << ":2\n";

	// Output precision info
	precision->IO("Prec", file, true, debug);

	// Output stats for each element of the matrix
	for (int ii = 0; ii < dim1; ii++) {
		for (int jj = 0; jj < dim2; jj++) {
			stats = getMoments(ii, jj);
			file << (stats[0]) << "," << (stats[1]) << ";";
		}
		file << "\n";
	}
}

void MassGaussianNode::IO(string description, istream& file, bool loadPriors,
		bool debug) {
	int stringSize = 5000;
	string input(stringSize, '\0');

	// Input descriptive line
	file.getline(&input[0], stringSize - 1);
	if (debug) {
		cerr << "Discarded: " << input << "\n";
	}

	// Input precision info
	input.assign(stringSize, '\0');
	precision->IO("Prec", file, loadPriors, debug);

	// Input stats for each matrix element
	char * pos = 0;
	double mom1, mom2;

	for (int ii = 0; ii < dim1; ii++) {
		// Reset input string and get next line
		input.assign(stringSize, '\0');
		file.getline(&input[0], stringSize - 1);
		pos = &input[0];

		// Get the stats for each inner element
		for (int jj = 0; jj < dim2; jj++) {
			mom1 = strtod(pos, &pos);
			pos = &pos[1];
			mom2 = strtod(pos, &pos);
			pos = &pos[1];

			// Correct for any errors in the save file. It *could* happen.
			if (mom2 < mom1 * mom1) {
				mom2 = (mom1 * mom1) + 0.1;
			}

			// Set moments accordingly
			setMoments(ii, jj, mom1, mom2);
			if (debug && (ii == 0 || ((dim1 - ii) == 1))
					&& (jj == 0 || ((dim2 - jj) == 1))) {
				cerr << "Read: " << input << "\n and set moments to: " << mom1
						<< "," << mom2 << "\n";
			}
		}
	}
}

void MassGaussianNode::update(int index1, int index2, int index3,
		double message[]) {
	double parameters[2];
	// Get the parameters from the prior and the message
	updatePrior(parameters);
	for (int j = 0; j < 2; j++) {
		parameters[j] += message[j];
	}

	// Update the moments accordingly.
	double variance = -0.5 / (parameters[1]);
	double mean = parameters[0] * variance;
	setMoments(index1, index2, mean, mean * mean + variance);
}

// Helper functions to get priors, moments and parameters
void MassGaussianNode::updatePrior(double priorParameters[]) {
	priorParameters[0] = precision->getEX() * mean->getEX();
	priorParameters[1] = -0.5 * precision->getEX();
}

void MassGaussianNode::setMoments(int inOne, int inTwo, double val1,
		double val2) {
	moments[inOne][inTwo][0] = val1;
	moments[inOne][inTwo][1] = val2;
}

void MassGaussianNode::getParams(int inOne, int inTwo, double* theseParams) {
	double* theseMoms = getMoments(inOne, inTwo);
	// Note:
	// 1. variance = -0.5 / (params[1]);
	// 2. moments = {params[0] * variance, params[0] * variance * params[0] * variance +  variance};

	// From 2: variance = theseMoms[1] - theseMoms[0]*theseMoms[0];
	// From 1: theseParams[1] = -0.5 / variance;
	// From 2: theseParams[0] = theseMoms[0] * -2 * theseParams[0]
	theseParams[1] = 1.0 / ((theseMoms[1] - theseMoms[0] * theseMoms[0]) * -2);

	theseParams[0] = theseMoms[0] * -2 * theseParams[1];
}

// Helper function to randomly initialise the moments array
void MassGaussianNode::initialiseMoments() {
	// Variables
	srand(time(0));

	double myMean = ((GaussianNode*) mean)->getMean();
	double variance = 1.0 / precision->getEX();

	// For each element, we obtain a Gaussian-distributed variable. We generat two at a time.
	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim2; j++) {
			double x1 = ((double) rand()) / (double) RAND_MAX;
			double x2 = ((double) rand()) / (double) RAND_MAX;

			double y1 = sqrt(-2.0 * log(x1)) * cos(2.0 * PI * x2);
			double y2 = sqrt(-2.0 * log(x1)) * sin(2.0 * PI * x2);

			// Use the first value, and use the second if it won't go over edge of array
			double value1 = myMean + y1 * sqrt(variance);
			setMoments(i, j, value1, value1 * value1 + variance);
			j++;
			if (j < dim2) {
				value1 = myMean + y2 * sqrt(variance);
				setMoments(i, j, value1, value1 * value1 + variance);
			}
		}
	}
}

}

