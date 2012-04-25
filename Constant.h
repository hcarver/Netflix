#ifndef CONSTANT_H_
#define CONSTANT_H_

#include "EXandX2.h"
#include "EXandLnX.h"
#include "Node.h"

namespace model {

class Constant: public EXandX2, public EXandLnX {

private:
	double d[2];

public:
	// Constructor and destructor
	Constant(double d);
	~Constant();

	// Get statistics functions
	double getEX();
	double getELnX();
	double getEX2();
	double* getMoments();

	// Set statistics functions
	void setVal(double newValue);

	// Functions inherited from Node
	virtual void IO(string description, istream& file, bool loadPriors,
			bool debug);
	virtual void IO(string description, ostream& file, bool loadPriors,
			bool debug);
	double getBound();
	void update(int index1, int index2, int index3, double message[]);
};

}

#endif /*CONSTANT_H_*/
