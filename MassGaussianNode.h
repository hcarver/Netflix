#ifndef MASSGAUSSIANNODE_H_
#define MASSGAUSSIANNODE_H_

#include "nfx.h"
#include "EXandX2.h"
#include "EXandLnX.h"
#include "Node.h"
#include <math.h>

namespace model {

class MassGaussianNode: public Node {

private:
	// Member variables for defining the prior and size of the array represented
	EXandX2* mean;
	EXandLnX* precision;
	int dim1, dim2;

	// Moments store
	double*** moments;

public:
	// Constructor and Destructor
	MassGaussianNode(int size1, int size2, EXandX2* mean, EXandLnX* precision);
	~MassGaussianNode();

	// Functions to update parents
	void updateMeanParent();
	void updatePrecisionParent();

	// Accessors for moments
	double* getMoments(int inOne, int inTwo);
	double getEX(int pos1, int pos2);
	double getEX2(int pos1, int pos2);

	// Accessors for mean and precision
	EXandX2* getMean();
	EXandLnX* getPrecision();

	// Function to get the vital statistics of the distribution
	void getVitalStatistics(double &xiPrior11, double &xiPrior12,
			double &xiPrior21, double &xiPrior22, double &xiMean,
			double &xiVarOfMean, double &xiAverageVariance);

	// Functions inherited from Node
	double getBound();
	virtual void IO(string description, istream& file, bool loadPriors,
			bool debug);
	virtual void IO(string description, ostream& file, bool loadPriors,
			bool debug);
	void update(int index1, int index2, int index3, double message[]);

private:
	// Helper functions to get priors, moments and parameters
	void updatePrior(double priorParameters[]);
	void setMoments(int inOne, int inTwo, double val1, double val2);
	void getParams(int inOne, int inTwo, double* params);

	// Helper function to randomly initialise the moments array
	void initialiseMoments();

};

}

#endif /*MASSGAUSSIANNODE_H_*/
