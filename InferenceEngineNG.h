#ifndef INFERENCEENGINENG_H_
#define INFERENCEENGINENG_H_

#include "NetflixDataNG.h"
#include "MassGaussianNode.h"
#include "Constant.h"
#include "DotProductGaussianAddNG.h"
#include "GaussianNode.h"
#include "PreProc.h"
#include "GammaNode.h"

#include <string>
#include <math.h>
#include <vector>

using namespace model;
using namespace utils;
using namespace std;

namespace infer {

class InferenceEngineNG {

private:
	// Data accessor
	NetflixDataNG data;

	// Model elements
	MassGaussianNode* DPusers;
	MassGaussianNode* DPmovies;
	MassGaussianNode* W;
	MassGaussianNode* userAddition;
	MassGaussianNode* movieAddition;
	MassGaussianNode* userHistoryMultiplier1;
	MassGaussianNode* userHistoryMultiplier3;
	MassGaussianNode* userHistoryMultiplier5;
	MassGaussianNode* userHistoryMultiplier10;
	MassGaussianNode* userHistoryMultiplier30;
	GaussianNode* allAdds;
	GammaNode* precision;
	DotProductGaussianAddNG* castFreeDPGA;

	// Variables for internal use
	bool debug;
	bool test;
	int startIteration;
	double oldBound;

	// Training dataset variables
	q** qsByUser;
	ushort* noViewingsPerUser;

	user** usersByMovie;
	rate** ratingsByMovie;
	ushort** whatIndexAMovieIsToAUser;
	ulong* noViewingsPerMovie;

public:
	// Constructor and Destructor
	InferenceEngineNG(bool isDebug, bool isTest, bool isOldSave,
			string &loadPath, int xiStartIteration);
	~InferenceEngineNG();

	// This is really the main function
	int go();

private:
	// Important functions for the engine
	void init();
	void joinMeUp();

	// Helper functions to get output file names.
	string getKString();
	string getSubmissionFileName(int iterations);
	string getStatsFileName(int iterations);
	string getSaveFileName(int iterations);

	// IO Helper functions
	void makeSubAndStatsFiles(int iterations);
	void doSave(int its);
	void save(string fileName);
	void loadThisModel(string fileName, bool loadPriors);
	void loadPreviousModel(string fileName, bool loadPriors);
	void IO(string fileName, bool isSave, bool loadPriors);
	void MakeStringStreamsSameLength(vector<stringstream*> xiSS);
	void AppendToAll(vector<stringstream*> xiSS, string xiContent);
	void OutputCurrentModelSummary();

	// Inference functions
	void update(int numIterations);
	double calculateBound(int itNumber);
	void getRMSEandBoundForObservedDataMatrix(double &bound, double &rmse,
			double &obsVariance);

	// Text output functions
	void nB(string output);
	void nb(string output);

	// Ratings data helper functions
	void readInViewings();
	void dataIntegrityCheck();

};

}

#endif /*INFERENCEENGINENG_H_*/
