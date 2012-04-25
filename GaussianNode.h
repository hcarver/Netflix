#ifndef GAUSSIANNODE_H_
#define GAUSSIANNODE_H_

#include "EXandX2.h"
#include "nfx.h"
#include <math.h>

namespace model {

class GaussianNode: public EXandX2 {
private:
	// Member variables for the prior mean and precision, moments and parameters
	double mean, precision;
	double moments[2];
	double parameters[2];

public:
	// Constructor and Destructor
	GaussianNode(double mn, double prec);
	~GaussianNode();

	// Accessors to prior and moments
	double getMean();
	double getPrecision();

	// Functions inherited from EXandX2
	double getEX();
	double getEX2();
	double* getMoments();

	// Functions inherited from Node
	double getBound();
	void update(int index1, int index2, int index3, double message[]);
	virtual void IO(string description, istream& file, bool loadPriors,
			bool debug);
	virtual void IO(string description, ostream& file, bool loadPriors,
			bool debug);

private:
	// Private functions for controlling moments
	bool initialiseMoments();
	void setMoments(double val1, double val2);
	void updatePrior(double priorParameters[]);
};

}

#endif /*GAUSSIANNODE_H_*/
