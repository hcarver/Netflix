#ifndef GAMMANODE_H_
#define GAMMANODE_H_

#include"EXandX2.h"
#include"EXandLnX.h"
#include "Node.h"
#include "math.h"

namespace model {

// Define constants for calculating bound and moments
#define euler       0.57721566490153286060651209008240243104215933593992
#define M_lnSqrt2PI 0.91893853320467274178

#define large 9.5
#define d1    -0.5772156649015328606065121
#define d2    1.64493406684822643647
#define small 1e-5
//static const double s3    = 1.0 / 12.0;
#define s3    0.083333333333333
//static const double s4    = 1.0 / 120.0;
#define s4    0.0083333333333333
//static const double s5    = 1.0 / 252.0;
#define s5    0.003968253968254
//static const double s6    = 1.0 / 240.0;
#define s6    0.004166666666666
//static const double s7    = 1.0 / 132.0;
#define s7    0.007575757575758
//static const double s8    = 691.0 / 32760.0;
#define s8    0.021092796092796
//static const double s9    = 1.0 / 12.0;
#define s9    0.083333333333333
//static const double s10   = 3617.0 / 8160.0;
#define s10   0.443259803921569

class GammaNode: public EXandLnX {

private:

// Private member variables - those that define distribn + useful gamma_series
	double shape, scale;
	double moments[2], parameters[2];
	double gamma_series[6];

public:
// Constructor and Destructor
	GammaNode(double xiShape, double xiScale);
	virtual ~GammaNode();

// Functions inherited from EXandLnX
	double getEX();
	double getELnX();

// Functions inherited from Node    
	double getBound();
	void update(int index1, int index2, int index3, double message[]);
	virtual void IO(string description, istream& file, bool loadPriors,
			bool debug);
	virtual void IO(string description, ostream& file, bool loadPriors,
			bool debug);

private:
// Moments-related functions. 
	void setMoments(double EX, double ELnX);
	void initialiseMoments();
	void updatePrior(double inHere[]);

// Helper functions to compute gamma-related function values
	double diGamma(double x);
	double logGamma(double x);

};

}

#endif /*GAMMANODE_H_*/
