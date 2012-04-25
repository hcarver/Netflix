#include "InferenceEngineNG.h"
#include <stdio.h>

// Main Function
// Process the inputted options
int main(int argc, char** args) {
	// Variables to store inputs
	bool argError = false;
	bool isDebug = false;
	bool isTest = false;
	bool isPreProc = false;
	bool isOldSave = false;
	int startIteration = -1;
	string loadPath = "";

	// arg 0 is the program name, ignore it
	// Process the rest
	for (int ii = 1; ii < argc; ii++) {
		if (args[ii][0] == '-') {
			switch (args[ii][1]) {
			case 'd':
			case 'D':
				isDebug = true;
				break;
			case 't':
			case 'T':
				isTest = true;
				break;
			case 'p':
			case 'P':
				isPreProc = true;
				break;
			case 'o':
			case 'O':
				isOldSave = true;
				break;
			default:
				// If it's a number, it's the start iteration number. If not, there must
				// be a problem.
				if (args[ii][1] >= '0' && args[ii][1] <= '9'
						&& startIteration < 0) {
					char * ptr = &(args[ii][1]);
					startIteration = strtol(ptr, &ptr, 10);
				} else {
					argError = true;
					cerr << args[ii]
							<< " didn't work properly - not recognised option\n";
				}
			}
		}
		// If we've already got a load path and there's another one in the input
		// there's a problem.
		else if (loadPath.length() > 0) {
			argError = true;
			cerr << args[ii]
					<< " didn't work properly - only one file path may be entered\n"
					<< "File Path " << loadPath << " already present\n";
		} else {
			loadPath = args[ii];
		}
	}

	// In case of argError, output a helpful string and return
	if (argError) {
		cerr
				<< "\nArgument error\n Only -d/D, -t/T, -p/P, -o/O -### flags recognised\nOnly one non flag may be entered (the path to load from)\n\n";
		return -1;
	}

	// Take action dependent on inputs
	if (isPreProc) {
		tool::PreProc pp(isDebug, isTest);
		pp.go();
	} else {
		infer::InferenceEngineNG eng(isDebug, isTest, isOldSave, loadPath,
				startIteration > 0 ? startIteration : 1);
		eng.go();
	}
}

namespace infer {

// Constructor and Destructor
InferenceEngineNG::InferenceEngineNG(bool isDebug, bool isTest, bool isOldSave,
		string & loadPath, int xiStartIteration) {
	// Take input and set debug and test as appropriate.
	startIteration = xiStartIteration;
	debug = isDebug;
	test = isTest;

	if (debug) {
		data.setDebug();
	}

	if (test) {
		data.setTest();
	}

	// Set member variable
	oldBound = -1E20;

	// Get data set mean and variance and set them to be correct if in test mode
	double mAndV[2];
	data.getDataMeanAndVariance(mAndV);

	if (test) {
		mAndV[0] = 3.0;
		mAndV[1] = 2.0;
	}

	// Set priors
	double allAddVariance = 0.0001;
	double movieMeanPriorVariance = 0.27675;
	double userMeanPriorVariance = 0.209642;
	double dpPriorVariance = pow(
			((mAndV[1] - (movieMeanPriorVariance + userMeanPriorVariance))
					/ ((double) K)), 0.5) / pow(2.0, 0.5);
	double wPriorVariance = dpPriorVariance;
	double userHistoryMultiplierVariance = 0.05;

	double obsVariance = mAndV[1] / 1.2;

	// Create the actual model nodes
	allAdds = new GaussianNode(mAndV[0], 1.0 / allAddVariance);
	precision = new GammaNode(5.0, 5.0 * obsVariance); // Variance of a U(-0.5, 0.5) distribution is 1/12 Precision new Constant(1.0/mAndV[1]);
	userAddition = new MassGaussianNode(U, 1, new GaussianNode(0, 20.0),
			new GammaNode(5.0, 5.0 * userMeanPriorVariance));
	movieAddition = new MassGaussianNode(M, 1, new GaussianNode(0, 20.0),
			new GammaNode(5.0, 5.0 * movieMeanPriorVariance));

	nB("\n\nGenerating dot-product Gaussian nodes... ");
	DPusers = new MassGaussianNode(U, K, new GaussianNode(0, 20.0),
			new GammaNode(5.0, 5.0 * dpPriorVariance));
	DPmovies = new MassGaussianNode(M, K, new GaussianNode(0, 20.0),
			new GammaNode(5.0, 5.0 * dpPriorVariance));
	W = new MassGaussianNode(M, K, new GaussianNode(0, 20.0),
			new GammaNode(5.0, 5.0 * wPriorVariance));
	userHistoryMultiplier1 = new MassGaussianNode(U, 5,
			new GaussianNode(0, 20.0),
			new GammaNode(5.0, 5.0 * userHistoryMultiplierVariance));
	userHistoryMultiplier3 = new MassGaussianNode(U, 5,
			new GaussianNode(0, 20.0),
			new GammaNode(5.0, 5.0 * userHistoryMultiplierVariance));
	userHistoryMultiplier5 = new MassGaussianNode(U, 5,
			new GaussianNode(0, 20.0),
			new GammaNode(5.0, 5.0 * userHistoryMultiplierVariance));
	userHistoryMultiplier10 = new MassGaussianNode(U, 5,
			new GaussianNode(0, 20.0),
			new GammaNode(5.0, 5.0 * userHistoryMultiplierVariance));
	userHistoryMultiplier30 = new MassGaussianNode(U, 5,
			new GaussianNode(0, 20.0),
			new GammaNode(5.0, 5.0 * userHistoryMultiplierVariance));

	// Create the DotProductGaussianAddNG object
	castFreeDPGA = new DotProductGaussianAddNG(debug);

	// If we have a loadPath, load the save according to whether it's from an old or current model
	if (loadPath.length() > 1) {
		if (isOldSave) {
			loadPreviousModel(loadPath, false);
		} else {
			loadThisModel(loadPath, false);
		}
	}

	// Read in viewings data
	nB("Reading viewings data...");
	readInViewings();
}

InferenceEngineNG::~InferenceEngineNG() {
	delete DPusers;
	delete DPmovies;
	delete userAddition;
	delete movieAddition;
	delete allAdds;
	delete precision;
	delete castFreeDPGA;
	delete W;
	delete userHistoryMultiplier1;
	delete userHistoryMultiplier3;
	delete userHistoryMultiplier5;
	delete userHistoryMultiplier10;
	delete userHistoryMultiplier30;

	for (user u = 0; u < U; u++) {
		delete[] qsByUser[u];
	}
	for (movie m = 0; m < M; m++) {
		delete[] usersByMovie[m];
		delete[] ratingsByMovie[m];
		delete[] whatIndexAMovieIsToAUser[m];
	}

	delete[] qsByUser;
	delete[] usersByMovie;
	delete[] ratingsByMovie;
	delete[] whatIndexAMovieIsToAUser;
	delete[] noViewingsPerMovie;
	delete[] noViewingsPerUser;
}

// This is really the main function
int InferenceEngineNG::go() {
	// Initialise the InferenceEngine
	init();
	// Update a lot
	update(100000);
	return 0;
}

// Important functions for the engine
void InferenceEngineNG::init() {
	nb("\nBeginning Initialisation... ");

	// Join up elements of model (i.e. give DPGA pointers to nodes)
	joinMeUp();

	// Initialise the DPGA
	castFreeDPGA->init(qsByUser, noViewingsPerUser);
	nB("Done. Engine initialised and ready to go.");
}

void InferenceEngineNG::joinMeUp() {
	castFreeDPGA->setUsers(DPusers);
	castFreeDPGA->setMovie(DPmovies);
	castFreeDPGA->setUserAdd(userAddition);
	castFreeDPGA->setMovieAdd(movieAddition);
	castFreeDPGA->setAllAdd(allAdds);
	castFreeDPGA->setPrecision(precision);
	castFreeDPGA->setW(W);
	// castFreeDPGA->setUserHistoryMultipliers(
	//   userHistoryMultiplier1,
	//   userHistoryMultiplier3,
	//   userHistoryMultiplier5,
	//   userHistoryMultiplier10,
	//   userHistoryMultiplier30);
}

// Helper functions to get output file names
string InferenceEngineNG::getKString() {
	stringstream s("");
	s << "K";
	s << K;
	return s.str();
}

string InferenceEngineNG::getSubmissionFileName(int iterations) {
	stringstream s("");
	s << data.getBaseSaveFileName() << "Sub_" << MODELSTRING << "_"
			<< getKString() << "_It" << iterations;
	return s.str();
}

string InferenceEngineNG::getStatsFileName(int iterations) {
	stringstream s("");
	s << data.getBaseSaveFileName() << "Stats_" << MODELSTRING << "_"
			<< getKString() << "_It" << iterations;
	return s.str();
}

string InferenceEngineNG::getSaveFileName(int iterations) {
	stringstream s("");
	s << data.getBaseSaveFileName() << "Save_" << MODELSTRING << "_"
			<< getKString() << "_It" << iterations;
	return s.str();
}

// IO Helper functions
void InferenceEngineNG::makeSubAndStatsFiles(int iterations) {
	data.writeQualifyingFile(getSubmissionFileName(iterations),
			getStatsFileName(iterations), castFreeDPGA);

	// If this isn't the first save since beginning the run, we delete the last save files.
	if (iterations - startIteration > saveAndBoundEvery) {
		data.deleteFileAsynch(
				getSubmissionFileName(iterations - saveAndBoundEvery));
		data.deleteFileAsynch(getStatsFileName(iterations - saveAndBoundEvery));
	}
}

void InferenceEngineNG::doSave(int iterations) {
	save(getSaveFileName(iterations));

	// If this isn't the first save since beginning the run, we delete the last save file.
	if (iterations - startIteration > saveAndBoundEvery) {
		data.deleteFileAsynch(getSaveFileName(iterations - saveAndBoundEvery));
	}
}

void InferenceEngineNG::save(string fileName) {
	cerr << "Save file: " << fileName << "\n";

	IO(fileName, true, true);
}

void InferenceEngineNG::loadThisModel(string fileName, bool loadPriors) {
	if (debug) {
		cerr << "Testing Load function\n";
	}

	IO(fileName, false, loadPriors);
}

void InferenceEngineNG::IO(string fileName, bool isSave, bool loadPriors) {
	if (isSave) {
		// Open output stream
		ofstream f;
		f.open(&fileName[0]);

		// Do IO for each model element in turn
		allAdds->IO("AllAdds", f, loadPriors, debug);
		userAddition->IO("UserAdds", f, loadPriors, debug);
		movieAddition->IO("MovieAdds", f, loadPriors, debug);
		DPusers->IO("DPusers", f, loadPriors, debug);
		DPmovies->IO("DPmovies", f, loadPriors, debug);
		W->IO("W", f, loadPriors, debug);
		userHistoryMultiplier1->IO("UserHistoryMultiplier1", f, loadPriors,
				debug);
		userHistoryMultiplier3->IO("UserHistoryMultiplier3", f, loadPriors,
				debug);
		userHistoryMultiplier5->IO("UserHistoryMultiplier5", f, loadPriors,
				debug);
		userHistoryMultiplier10->IO("UserHistoryMultiplier10", f, loadPriors,
				debug);
		userHistoryMultiplier30->IO("UserHistoryMultiplier30", f, loadPriors,
				debug);

		f.close();
	} else {
		// Open input stream
		ifstream f;
		f.open(&fileName[0]);

		// Do IO for each model element in turn
		allAdds->IO("AllAdds", f, loadPriors, debug);
		userAddition->IO("UserAdds", f, loadPriors, debug);
		movieAddition->IO("MovieAdds", f, loadPriors, debug);
		DPusers->IO("DPusers", f, loadPriors, debug);
		DPmovies->IO("DPmovies", f, loadPriors, debug);
		W->IO("W", f, loadPriors, debug);
		userHistoryMultiplier1->IO("UserHistoryMultiplier1", f, loadPriors,
				debug);
		userHistoryMultiplier3->IO("UserHistoryMultiplier3", f, loadPriors,
				debug);
		userHistoryMultiplier5->IO("UserHistoryMultiplier5", f, loadPriors,
				debug);
		userHistoryMultiplier10->IO("UserHistoryMultiplier10", f, loadPriors,
				debug);
		userHistoryMultiplier30->IO("UserHistoryMultiplier30", f, loadPriors,
				debug);

		f.close();
	}
}

void InferenceEngineNG::loadPreviousModel(string fileName, bool loadPriors) {
	cerr << "\n\n!!!!!loadPrevousModel not yet defined!!!!\n\n";
}

void InferenceEngineNG::MakeStringStreamsSameLength(
		vector<stringstream*> xiSS) {
	// NB. LEAVE THIS FUNCTION USING INTs NOT UINTs
	int maxSoFar = -1;

	for (int ii = 0; ii < xiSS.size(); ii++) {
		int thisLength = xiSS[ii]->str().length();

		if (thisLength > maxSoFar) {
			maxSoFar = thisLength;
		}
	}

	for (int ii = 0; ii < xiSS.size(); ii++) {
		while ((xiSS[ii]->str().length()) < maxSoFar) {
			(*xiSS[ii]) << " ";
		}
	}
}

void InferenceEngineNG::AppendToAll(vector<stringstream*> xiSS,
		string xiContent) {
	for (uint ii = 0; ii < xiSS.size(); ii++) {
		*xiSS[ii] << xiContent;
	}
}

void InferenceEngineNG::OutputCurrentModelSummary() {
	// HeaderRow, DPU, DPM, W, UA, MA, AA
	vector<stringstream*> lSS;
	for (uint ii = 0; ii < 7; ii++) {
		lSS.push_back(new stringstream(""));
	}

	*lSS[0] << "";
	*lSS[1] << "DPU";
	*lSS[2] << "DPM";
	*lSS[3] << "W";
	*lSS[4] << "UA";
	*lSS[5] << "MA";
	*lSS[6] << "AA";

	if (debug) {
		for (uint ii = 0; ii < lSS.size(); ii++) {
			cerr << lSS[ii]->str();
		}
	}

	MakeStringStreamsSameLength(lSS);

	if (debug) {
		for (uint ii = 0; ii < lSS.size(); ii++) {
			cerr << lSS[ii]->str();
		}
	}

	AppendToAll(lSS, "|");

	if (debug) {
		for (uint ii = 0; ii < lSS.size(); ii++) {
			cerr << lSS[ii]->str();
		}
	}

	// Now obtain all the statistics
	double prior11[6];
	double prior12[6];
	double prior21[6];
	double prior22[6];
	double meanMean[6];
	double varOfMean[6];
	double averageVariance[6];

	int eltIndex = 0;
	DPusers->getVitalStatistics(prior11[eltIndex], prior12[eltIndex],
			prior21[eltIndex], prior22[eltIndex], meanMean[eltIndex],
			varOfMean[eltIndex], averageVariance[eltIndex]);

	eltIndex = 1;
	DPmovies->getVitalStatistics(prior11[eltIndex], prior12[eltIndex],
			prior21[eltIndex], prior22[eltIndex], meanMean[eltIndex],
			varOfMean[eltIndex], averageVariance[eltIndex]);

	eltIndex = 2;
	W->getVitalStatistics(prior11[eltIndex], prior12[eltIndex],
			prior21[eltIndex], prior22[eltIndex], meanMean[eltIndex],
			varOfMean[eltIndex], averageVariance[eltIndex]);

	eltIndex = 3;
	userAddition->getVitalStatistics(prior11[eltIndex], prior12[eltIndex],
			prior21[eltIndex], prior22[eltIndex], meanMean[eltIndex],
			varOfMean[eltIndex], averageVariance[eltIndex]);

	eltIndex = 4;
	movieAddition->getVitalStatistics(prior11[eltIndex], prior12[eltIndex],
			prior21[eltIndex], prior22[eltIndex], meanMean[eltIndex],
			varOfMean[eltIndex], averageVariance[eltIndex]);

	eltIndex = 5;
	prior11[eltIndex] = allAdds->getMean();
	prior12[eltIndex] = 0.0;
	prior21[eltIndex] = allAdds->getPrecision();
	prior22[eltIndex] = log(allAdds->getPrecision());
	meanMean[eltIndex] = allAdds->getEX();
	varOfMean[eltIndex] = 0.0;
	averageVariance[eltIndex] = allAdds->getEX2()
			- allAdds->getEX() * allAdds->getEX();

	// Now output each statistic in turn
	// Prior11
	*lSS[0] << "E(E(X))";
	for (int ii = 1; ii < 7; ii++) {
		*lSS[ii] << prior11[ii - 1];
	}
	MakeStringStreamsSameLength(lSS);

	AppendToAll(lSS, "|");

	// Prior12
	*lSS[0] << "Var(E(X))";
	for (int ii = 1; ii < 7; ii++) {
		*lSS[ii] << prior12[ii - 1];
	}
	MakeStringStreamsSameLength(lSS);

	AppendToAll(lSS, "|");

	// Prior21
	*lSS[0] << "E(Var(X))";
	for (int ii = 1; ii < 7; ii++) {
		*lSS[ii] << 1.0 / prior21[ii - 1];
	}
	MakeStringStreamsSameLength(lSS);

	AppendToAll(lSS, "|");

	// Prior22
	*lSS[0] << "e^(E(Ln(Var(X))))";
	for (int ii = 1; ii < 7; ii++) {
		*lSS[ii] << exp(-prior22[ii - 1]);
	}
	MakeStringStreamsSameLength(lSS);

	AppendToAll(lSS, "|");

	// MeanMean
	*lSS[0] << "E(E(x))";
	for (int ii = 1; ii < 7; ii++) {
		*lSS[ii] << meanMean[ii - 1];
	}
	MakeStringStreamsSameLength(lSS);

	AppendToAll(lSS, "|");

	// varOfMean
	*lSS[0] << "Var(E(x))";
	for (int ii = 1; ii < 7; ii++) {
		*lSS[ii] << varOfMean[ii - 1];
	}
	MakeStringStreamsSameLength(lSS);

	AppendToAll(lSS, "|");

	// averageVariance
	*lSS[0] << "E(Var(x))";
	for (int ii = 1; ii < 7; ii++) {
		*lSS[ii] << averageVariance[ii - 1];
	}
	MakeStringStreamsSameLength(lSS);

	AppendToAll(lSS, "\n");

	// Output!
	cerr
			<< "\nMODEL Summary (which should be refactored, and should include observed / predicted sets)\n";
	for (int ii = 0; ii < 7; ii++) {
		cerr << lSS[ii]->str();
		delete lSS[ii];
	}
	cerr << "\n";
}

// Inference functions 
void InferenceEngineNG::update(int numIterations) {
	// Always begin with bound calculation
	OutputCurrentModelSummary();
	calculateBound(0);

	// If debugging on real data, test outputs too
	if (debug && !test) {
		makeSubAndStatsFiles(0);
		doSave(0);
	}

	// Iterate and update
	for (int i = startIteration; i <= numIterations + startIteration; i++) {
		cerr << "\nIteration number " << i << "...\n";
		castFreeDPGA->doAnIteration(qsByUser, noViewingsPerUser,
				noViewingsPerMovie, usersByMovie, whatIndexAMovieIsToAUser,
				false); //test || debug) || ((i - startIteration) > 3));

		// Always calculate bound if operating on small dataset
		if (test) {
			calculateBound(i);
		}

		// Otherwise, calculate bound, do outputs and make an entry
		else if ((i % saveAndBoundEvery == 0)) {
			calculateBound(i);
			OutputCurrentModelSummary();
			makeSubAndStatsFiles(i);
			data.makeEntryAsynch(getSubmissionFileName(i));
			doSave(i);
		}

		// Output model statistics in the first few iterations to see how new model's doing
		else if ((i - startIteration) < saveAndBoundEvery) {
			OutputCurrentModelSummary();
		}
	}
}

double InferenceEngineNG::calculateBound(int itNumber) {
	// Get the rmse and variance from the observed matrix
	double obs, obsRMSE, obsVariance;
	getRMSEandBoundForObservedDataMatrix(obs, obsRMSE, obsVariance);

	// Bounds for each node
	double aa = allAdds->getBound();
	double um = userAddition->getBound();
	double mm = movieAddition->getBound();
	double mdp = DPmovies->getBound();
	double udp = DPusers->getBound();
	double w = W->getBound();
	double uhm1 = userHistoryMultiplier1->getBound();
	double uhm3 = userHistoryMultiplier3->getBound();
	double uhm5 = userHistoryMultiplier5->getBound();
	double uhm10 = userHistoryMultiplier10->getBound();
	double uhm30 = userHistoryMultiplier30->getBound();
	double p = precision->getBound();

	// Bounds for all the means
	double meanBounds = userAddition->getMean()->getBound()
			+ movieAddition->getMean()->getBound()
			+ DPmovies->getMean()->getBound() + DPusers->getMean()->getBound()
			+ W->getMean()->getBound()
			+ userHistoryMultiplier1->getMean()->getBound()
			+ userHistoryMultiplier3->getMean()->getBound()
			+ userHistoryMultiplier5->getMean()->getBound()
			+ userHistoryMultiplier10->getMean()->getBound()
			+ userHistoryMultiplier30->getMean()->getBound();

	// Bounds for all the precisions
	double precisionBounds = userAddition->getPrecision()->getBound()
			+ movieAddition->getPrecision()->getBound()
			+ DPmovies->getPrecision()->getBound()
			+ DPusers->getPrecision()->getBound()
			+ W->getPrecision()->getBound()
			+ userHistoryMultiplier1->getPrecision()->getBound()
			+ userHistoryMultiplier3->getPrecision()->getBound()
			+ userHistoryMultiplier5->getPrecision()->getBound()
			+ userHistoryMultiplier10->getPrecision()->getBound()
			+ userHistoryMultiplier30->getPrecision()->getBound();

	// Sum the bounds
	double b = obs + aa + um + mm + mdp + udp + w + uhm1 + uhm3 + uhm5 + uhm10
			+ uhm30 + p + meanBounds + precisionBounds;

	cerr << "\nRMSE on Training Matrix: " << obsRMSE;

	cerr << "\n\nBOUNDS:\n";
	cerr
			<< "Obs    ; All Adds ; U Mean   ; M Mean ;   M DP   ;   U DP    ;    W    ;  UHMs  ; Precision ; Prior Means ; Prior Precs   ;\n";
	cerr << obs << " ; " << aa << " ; " << um << " ; " << mm << " ; " << mdp
			<< " ;" << udp << " ;" << w << ";"
			<< (uhm1 + uhm3 + uhm5 + uhm10 + uhm30) << ";" << p << " ;"
			<< meanBounds << " ; " << precisionBounds << ";\n";
	cerr << "Total Bound of " << b;

	// If the oldBound isn't a number (and the new bound is), or is smaller than the new bound, say so
	if ((isnan(oldBound) && !isnan(b)) || oldBound <= b) {
		cerr << " which is greater than the last bound (" << oldBound
				<< ") :D\n";
	} else {
		cerr << " which is LOWER than the last bound (" << oldBound
				<< ") :[(\n";
	}

	// Update old bound
	oldBound = b;

	return b;
}

void InferenceEngineNG::getRMSEandBoundForObservedDataMatrix(double &bound,
		double &rmse, double &obsVariance) {
	// Useful scalars
	double b = 0.0;

	// Scalars to track statistics
	double numberObservations = 0.0;
	double totalSize = 0.0;
	double totalSizeSquared = 0.0;
	double totalErrorSquared = 0.0;

	// Iterate over each user and rating and calculate stats
#pragma omp parallel for reduction(+:b) reduction(+:numberObservations) reduction(+:totalErrorSquared) reduction(+:totalSize) reduction(+:totalSizeSquared) schedule(guided,5)
	for (user u = 0; u < U; u++) {
		double prior[] = { 0.0, 0.0 };
		double stats[] = { 0.0, 0.0 };
		prior[1] = -0.5 * precision->getEX();

		for (int index = 0; index < noViewingsPerUser[u]; index++) {
			// Get rating and corresponding prediction
			double r = (double) (qsByUser[u][index] % RATE);
			castFreeDPGA->predictKnownIndex(stats, u, index, qsByUser);

			// Do obs matrix bound
			prior[0] = precision->getEX() * stats[0];

			b += (r * prior[0]) + (r * r * prior[1]);
			b += 0.5
					* (precision->getELnX() - precision->getEX() * stats[1]
							- LOG2PI);

			// Do stats bit
			numberObservations += 1.0;
			totalErrorSquared += pow((r - stats[0]), 2.0);
			totalSize += stats[0];
			totalSizeSquared += pow(stats[0], 2.0);
		}
	}
#pragma omp barrier

	// Set return values
	bound = b;
	rmse = pow((totalErrorSquared / numberObservations), 0.5);
	obsVariance = (totalSizeSquared / numberObservations)
			- pow(totalSize / numberObservations, 2.0);

	cerr << "\nOBSERVED data set summary\n";
	cerr << "totalSizeSquared: " << totalSizeSquared << ", totalSize: "
			<< totalSize << ", numberObservations: " << numberObservations
			<< "\n";
	cerr << "EX2 : " << (totalSizeSquared / numberObservations) << "\n";
	cerr << "EX: " << (totalSize / numberObservations) << "\n";
	cerr << "Variance: " << obsVariance << "\n";
}

// Text output functions
void InferenceEngineNG::nB(string output) {
	nb(output);
	nb("\n");
}

void InferenceEngineNG::nb(string output) {
	cerr << output;
}

// Ratings data helper functions
void InferenceEngineNG::readInViewings() {
	// Create first level of arrays
	qsByUser = new q*[U];
	noViewingsPerUser = new unsigned short[U];
	usersByMovie = new user*[M];
	noViewingsPerMovie = new unsigned long[M];
	whatIndexAMovieIsToAUser = new unsigned short*[M];
	ratingsByMovie = new rate*[M];

	data.readBigRatings(usersByMovie, noViewingsPerMovie, ratingsByMovie,
			qsByUser, noViewingsPerUser, whatIndexAMovieIsToAUser);

	// If in debug mode check both sides of data match
	if (debug) {
		dataIntegrityCheck();
	}
}

void InferenceEngineNG::dataIntegrityCheck() {
	cerr << "Checking data integrity...";

	// Create some monitoring variables
	int noViewings = 0;
	double totalRatings = 0;
	double totalSquared = 0;

	// Create some useful variable for the iteration
	rate r;
	movie m;
	user u;

	// Iterate over ratings on the user side. Calculate statistics, check validity of rating and movie values
#pragma omp parallel for private(r, m) reduction(+:noViewings) reduction(+:totalRatings) reduction(+:totalSquared) schedule(guided,5)
	for (u = 0; u < U; u++) {
		for (int index = 0; index < noViewingsPerUser[u]; index++) {
			noViewings++;
			r = qsByUser[u][index] % RATE;
			m = qsByUser[u][index] % MOV;
			totalRatings += r;
			totalSquared += r * r;

			if (r > 5 || r < 1) {
				cerr << "bad rating for user " << u << " movie " << m << " of "
						<< (int) r << "\n";
			}

			if (m >= M) {
				cerr << "bad movie for user " << u << " index " << index
						<< " thought to be " << m << "\n";
			}
		}
	}
#pragma omp barrier

	// Output statistics from the user side
	double mean = (totalRatings / ((double) noViewings));
	cerr << "\nUser side shows " << noViewings
			<< " total ratings and an average rating of " << mean
			<< " with variance "
			<< (totalSquared / ((double) noViewings) - mean * mean) << "\n";

	noViewings = 0.0;
	totalRatings = 0.0;
	totalSquared = 0.0;

	// Now on the movie side, go through all ratings, check statistics and validity of user and rating values
#pragma omp parallel for private(r, u) reduction(+:noViewings) reduction(+:totalRatings) reduction(+:totalSquared)
	for (m = 0; m < M; m++) {
		for (unsigned int index = 0; index < noViewingsPerMovie[m]; index++) {
			// Update statistics variables
			noViewings++;
			r = ratingsByMovie[m][index];
			u = usersByMovie[m][index];
			totalRatings = totalRatings + (double) r;
			totalSquared += pow((double) r, 2.0);

			// Check validity of values
			if (r > 5 || r < 1) {
				cerr << "bad rating for movie " << m << " user " << u << " of "
						<< (int) r << "\n";
			}
			if (u >= U || u < 0) {
				cerr << "bad user for movie " << m << " index " << index
						<< " thought to be " << u << "\n";
			}
		}
	}
#pragma omp barrier

	// Output statistics from the movie side
	mean = (totalRatings / ((double) noViewings));
	cerr << "Movie side shows " << noViewings
			<< " total ratings and an average rating of " << mean
			<< " with variance "
			<< (totalSquared / ((double) noViewings) - mean * mean) << "\n";

	// Now we check each user side rating has a corresponding movie side rating of the same value
#pragma omp parallel for schedule(guided,5)
	for (m = 0; m < M; m++) {
		for (uint index = 0; index < noViewingsPerMovie[m]; index++) {
			// First check each movie-side rating corresponds to a user-side rating
			if ((qsByUser[usersByMovie[m][index]][whatIndexAMovieIsToAUser[m][index]]
					% MOV) != m) {
				cerr << "movie " << m << " thought it was index "
						<< whatIndexAMovieIsToAUser[m][index] << " to user "
						<< usersByMovie[m][index] << " but it wasn't";
			}
			// Then check that the rating at each point had the same value
			if (ratingsByMovie[m][index]
					!= (qsByUser[usersByMovie[m][index]][whatIndexAMovieIsToAUser[m][index]]
							% RATE)) {
				cerr
						<< "ratings did not match in user and movie stores for movie "
						<< m << " index " << index << " user "
						<< usersByMovie[m][index] << "Pred1:"
						<< (int) ratingsByMovie[m][index] << "Pred2:"
						<< (int) (qsByUser[usersByMovie[m][index]][whatIndexAMovieIsToAUser[m][index]]
								% RATE);
			}
		}
	}
#pragma omp barrier
}

}
