#include "DotProductGaussianAddNG.h"
#include "omp.h"

namespace model {

// Constructor and Destructor
DotProductGaussianAddNG::DotProductGaussianAddNG(bool isDebug) {
	debug = isDebug;
}

DotProductGaussianAddNG::~DotProductGaussianAddNG() {
	// Need number of ratings per user to properly delete some of these arrays
	for (user u = 0; u < U; u++) {
		delete[] messages0[u];
		delete[] bigTempMessageArray[u];
		delete[] dpMomRecUrIx0[u];
	}

	delete[] messages0;
	delete[] bigTempMessageArray;
	delete[] dpMomRecUrIx0;

	for (int k = 0; k < K; k++) {
		for (user u = 0; u < U; u++) {
			delete[] wRecKUr[k][u];
		}
		delete[] wRecKUr[k];
	}
	delete[] wRecKUr;
}

// Initialisation function
void DotProductGaussianAddNG::init(q** &qsByUser,
		unsigned short* &ratingsPerUser) {
	// Start by initialising all the K-sized arrays
	wRecKUr = new double**[K];

#pragma omp parallel for  
	for (uint thisK = 0; thisK < K; thisK++) {
		wRecKUr[thisK] = new double *[U];
	}
#pragma omp barrier

	// Then initialise all the U-sized arrays
	messages0 = new double*[U];
	bigTempMessageArray = new double*[U];
	dpMomRecUrIx0 = new double*[U];

#pragma omp parallel for
	for (user u = 0; u < U; u++) {
		for (uint thisK = 0; thisK < K; thisK++) {
			wRecKUr[thisK][u] = new double[2];
		}

		messages0[u] = new double[ratingsPerUser[u]];
		dpMomRecUrIx0[u] = new double[ratingsPerUser[u]];
		bigTempMessageArray[u] = new double[2];

		for (int index = 0; index < ratingsPerUser[u]; index++) {
			messages0[u][index] = 0.0;
		}
	}
#pragma omp barrier

	// Then we fill the more important arrays with data about current estimates
	refreshWStore(qsByUser, ratingsPerUser);
	refreshDPStore(qsByUser, ratingsPerUser);
}

// Functions to set the nodes
void DotProductGaussianAddNG::setUsers(MassGaussianNode* gn) {
	DPusers = gn;
}

void DotProductGaussianAddNG::setMovie(MassGaussianNode* gn) {
	DPmovies = gn;
}

void DotProductGaussianAddNG::setUserAdd(MassGaussianNode* n) {
	userAdds = n;
}

void DotProductGaussianAddNG::setMovieAdd(MassGaussianNode* n) {
	movieAdds = n;
}

void DotProductGaussianAddNG::setW(MassGaussianNode* n) {
	W = n;
}

void DotProductGaussianAddNG::setUserHistoryMultipliers(MassGaussianNode* in1,
		MassGaussianNode* in3, MassGaussianNode* in5, MassGaussianNode* in10,
		MassGaussianNode* in30) {
	uhm1 = in1;
	uhm3 = in3;
	uhm5 = in5;
	uhm10 = in10;
	uhm30 = in30;
}

void DotProductGaussianAddNG::setAllAdd(EXandX2* n) {
	allAdds = n;
}

void DotProductGaussianAddNG::setPrecision(EXandLnX* n) {
	precision = n;
}

// Inference function
void DotProductGaussianAddNG::doAnIteration(q** &qsByUser,
		unsigned short* &noViewingsPerUser, unsigned long* &ratingsPerMovie,
		user ** &usersByMovie, unsigned short ** &whatIndexAMovieIsToAUser,
		bool updatePrecision) {
	if (updatePrecision) {
		// First update the overall precision (which should hopefully avoid the problems of
		// overfitting on a single model element)
		cerr << "Precision\n";
		toPrecision(qsByUser, noViewingsPerUser);
	}

	// Then the all adds
	cerr << "AllAdds\n";
	toAllAdds(qsByUser, noViewingsPerUser);

	// Then the user- and movie-specific additions
	cerr << "UserAdds\n";
	pastUserAdds(qsByUser, noViewingsPerUser);
	cerr << "MovieAdds\n";
	pastMovieAdds(qsByUser, noViewingsPerUser);

	// Then each element contributing to the dot-product node in turn
	cerr << "DP\n";
	pastDP(qsByUser, noViewingsPerUser, ratingsPerMovie, usersByMovie,
			whatIndexAMovieIsToAUser);

	// First the parents of the addition nodes
	userAdds->updateMeanParent();
	userAdds->updatePrecisionParent();
	movieAdds->updateMeanParent();
	movieAdds->updatePrecisionParent();

	// Then the parents of the dot-product nodes
	DPusers->updateMeanParent();
	DPusers->updatePrecisionParent();
	DPmovies->updateMeanParent();
	DPmovies->updatePrecisionParent();
	W->updateMeanParent();
	W->updatePrecisionParent();

	refreshDPStore(qsByUser, noViewingsPerUser);
}

// Functions to make predictions (either for training or test ratings)
void DotProductGaussianAddNG::predictKnownIndex(double* stats, user userID,
		uint index, q** &qsByUser) {
	predictUnknownIndex(stats, userID, (qsByUser[userID][index] % MOV));
}

void DotProductGaussianAddNG::predictUnknownIndex(double* stats, user userID,
		movie movieID) {
	getDPMomentsKnownMovie(userID, movieID, stats);
	AStats(stats, movieAdds->getMoments(movieID, 0));
	AStats(stats, userAdds->getMoments(userID, 0));
	AStats(stats, allAdds->getMoments());
}

// Parts of the inference algorithm

// Each of these follows a similar pattern:
// - Do any 'helper' calculations
// - Calculcate messages
// - If in debug, check messages
// - Update
// - Maintain current estimates
// - If in debug mode, check current estimates are correct

void DotProductGaussianAddNG::toPrecision(q **&qsByUser,
		unsigned short *&ratingsPerUser) {
	// Message to the precision node (only do when we do bound - it's just too expensive otherwise)
	double m0 = 0.0;
	double m1 = 0.0;

#pragma omp parallel for reduction(+:m0) reduction(+:m1) schedule(guided,5)
	for (user u = 0; u < U; u++) {
		double stats[] = { 0.0, 0.0 };

		for (int index = 0; index < ratingsPerUser[u]; index++) {
			double r = (double) (qsByUser[u][index] % RATE);
			predictKnownIndex(stats, u, index, qsByUser);
			m0 += r * stats[0] - 0.5 * r * r * stats[1];
		}
		m1 += 0.5 * ((double) ratingsPerUser[u]);
	}
#pragma omp barrier

	double message[] = { m0, m1 };

	precision->update(0, 0, 0, message);
}

void DotProductGaussianAddNG::toAllAdds(q** &qsByUser,
		unsigned short* &ratingsPerUser) {
	// Helper calculations
	double precEX = precision->getEX();
	double message1 = -0.5 * precEX;

	double m0 = 0.0;
	double m1 = 0.0;

	// Get the message to the allAdds node
#pragma omp parallel for reduction(+:m0) reduction(+:m1)
	for (user u = 0; u < U; u++) {
		int noViewings = ratingsPerUser[u];
		q* myQs = qsByUser[u];
		double* myDpMoms = dpMomRecUrIx0[u];
		double myUserAdd0 = userAdds->getMoments(u, 0)[0];

		m1 += message1 * ((double) noViewings);

		for (int i = 0; i < noViewings; i++) {
			q currentQ = myQs[i];
			m0 +=
					precEX * ((double) (currentQ % RATE))
							+ 2 * message1
									* (myUserAdd0
											+ movieAdds->getMoments(
													(currentQ % MOV), 0)[0]
											+ myDpMoms[i]);
		}
	}
#pragma omp barrier

	double message[] = { m0, m1 };

	// If in debug, check the validity of the message
	if (debug) {
		double correctMsg[2] = { 0, 0 };
		calcMessageToAllAdd(qsByUser, ratingsPerUser, correctMsg);

		double mod = fabs(correctMsg[0] - message[0])
				+ fabs(correctMsg[1] - message[1]);
		cerr << "Average absolute error in all adds messages is " << (mod / 2.0)
				<< "\n";
	}

	// Update the node
	allAdds->update(0, 0, 0, message);

	// Get message to pass on up
	double allAddsValueTimesMessage1Times2 = allAdds->getEX() * message1 * 2.0;

#pragma omp parallel for
	for (user u = 0; u < U; u++) {
		int noViewings = ratingsPerUser[u];
		q* myQs = qsByUser[u];

		double* myMessages = messages0[u];

		for (int i = 0; i < noViewings; i++) {
			myMessages[i] = precEX * ((double) (myQs[i] % RATE))
					+ allAddsValueTimesMessage1Times2;
		}
	}
#pragma omp barrier
}

void DotProductGaussianAddNG::pastUserAdds(q** &qsByUser,
		unsigned short* &ratingsPerUser) {
	// Helper variable - access precision only once
	double message1 = -0.5 * precision->getEX();

	// Can treat each user separately
#pragma omp parallel for
	for (user u = 0; u < U; u++) {
		// Make some helpful variables
		int noViewings = ratingsPerUser[u];
		q* myQs = qsByUser[u];
		double* myMessages0 = messages0[u];
		double* myDpMoms = dpMomRecUrIx0[u];

		double msg[2] = { 0.0, 0.0 };

		// Make the message for this user
		for (int i = 0; i < noViewings; i++) {
			msg[0] += myMessages0[i]
					+ 2 * message1
							* (movieAdds->getMoments((myQs[i] % MOV), 0)[0]
									+ myDpMoms[i]);
		}
		msg[1] = message1 * ((double) noViewings);

		// If debug (but not for every user), check the message
		if (debug && (u % 50 == 0)) {
			double correctMsg[2] = { 0, 0 };
			calcMessageToUserAdd(qsByUser, ratingsPerUser, u, correctMsg);

			double mod = fabs(correctMsg[0] - msg[0])
					+ fabs(correctMsg[1] - msg[1]);
			if (mod / 2.0 > 1e-9) {
				cerr << "Big error in user adds messages is " << (mod / 2.0)
						<< " for user " << u << "\n";
			}
		}

		// Perform the update
		userAdds->update(u, 0, 0, msg);

		// Update the current estimates
		double uNewValTimes2TimesMessage1 = userAdds->getMoments(u, 0)[0] * 2.0
				* message1;

		for (int i = 0; i < noViewings; i++) {
			myMessages0[i] = myMessages0[i] + uNewValTimes2TimesMessage1;
		}
	}
#pragma omp barrier
}

void DotProductGaussianAddNG::pastMovieAdds(q** &qsByUser,
		unsigned short* &ratingsPerUser) {
	// Helper variables - an array to store the messages for each thread, the value
	// of the precision and the value of the movieAdd elements
	double *** perThreadMessages = new double**[4];
	double message1 = -0.5 * precision->getEX();

	// Populate the per thread message array with 0s and populate the movie add local array
#pragma omp parallel for
	for (int thread = 0; thread < 4; thread++) {
		perThreadMessages[thread] = new double*[M];
		for (movie m = 0; m < M; m++) {
			perThreadMessages[thread][m] = new double[2];
			perThreadMessages[thread][m][0] = 0.0;
			perThreadMessages[thread][m][1] = 0.0;
		}
	}
#pragma omp barrier

	// Find the message for each movie, on each thread separately
#pragma omp parallel for
	for (user u = 0; u < U; u++) {
		int noViewings = ratingsPerUser[u];
		double* myMessages0 = messages0[u];
		double* myDpMoms = dpMomRecUrIx0[u];
		q* myQs = qsByUser[u];

		double** myMessages = perThreadMessages[omp_get_thread_num()];

		for (int i = 0; i < noViewings; i++) {
			movie m = myQs[i] % MOV;
			myMessages[m][0] += myMessages0[i] + 2 * message1 * myDpMoms[i];
			myMessages[m][1] += message1;
		}
	}
#pragma omp barrier

	// Add all messages from threads 1 to 3 on to thread 0.
#pragma omp parallel for
	for (movie m = 0; m < M; m++) {
		for (int thread = 1; thread < 4; thread++) {
			for (int moment = 0; moment < 2; moment++) {
				perThreadMessages[0][m][moment] +=
						perThreadMessages[thread][m][moment];
			}
		}
	}
#pragma omp barrier

	// If in debug mode, check the message for thread 0 are correct (but not for every movie)
	if (debug) {
		for (movie m = 0; m < M; m += 50) {
			double correctMsg[2];
			calcMessageToMovieAdd(qsByUser, ratingsPerUser, m, correctMsg);

			double mod = fabs(correctMsg[0] - perThreadMessages[0][m][0])
					+ fabs(correctMsg[1] - perThreadMessages[0][m][1]);
			if (mod / 2.0 > 1e-8) {
				cerr << "\nBig error in movie add messages of " << (mod / 2.0)
						<< " for movie " << m << "\n";
				cerr << "Got " << perThreadMessages[0][m][0] << ","
						<< perThreadMessages[0][m][1] << " should have been "
						<< correctMsg[0] << "," << correctMsg[1] << "\n";
			}
		}
	}

	// Update the movie adds with the messages and use a movie vector for access to the new movie addition values
	double movieAddTimes2TimesMessage1[M];

#pragma omp parallel for
	for (movie m = 0; m < M; m++) {
		movieAdds->update(m, 0, 0, perThreadMessages[0][m]);
		movieAddTimes2TimesMessage1[m] = movieAdds->getMoments(m, 0)[0] * 2.0
				* message1;
	}
#pragma omp barrier    

	// Update the messages0
#pragma omp parallel for
	for (user u = 0; u < U; u++) {
		int noViewings = ratingsPerUser[u];
		q* myQs = qsByUser[u];
		double* myMessages0 = messages0[u];

		for (int i = 0; i < noViewings; i++) {
			myMessages0[i] = myMessages0[i]
					+ movieAddTimes2TimesMessage1[myQs[i] % MOV];
		}
	}
#pragma omp barrier

	// Delete the per-thread message stores and the mov vector
	for (int thread = 0; thread < 4; thread++) {
		for (movie m = 0; m < M; m++) {
			delete[] perThreadMessages[thread][m];
		}
		delete[] perThreadMessages[thread];
	}
	delete[] perThreadMessages;
}

void DotProductGaussianAddNG::pastDP(q** &qsByUser,
		unsigned short* &noViewingsPerUser, unsigned long* &noViewingsPerMovie,
		user ** &usersByMovie, unsigned short ** &whatIndexAMovieIsToAUser) {
	// Update each element of the dot-product (0 to K-1) in order
	for (int i = 0; i < K; i++) {
		cerr << i << ",";
		pastDPUM(i, qsByUser, noViewingsPerUser, noViewingsPerMovie,
				usersByMovie, whatIndexAMovieIsToAUser);
	}
	cerr << '\n';
}

void DotProductGaussianAddNG::pastDPUM(int thisK, q** &qsByUser,
		unsigned short* &noViewingsPerUser, unsigned long* &noViewingsPerMovie,
		user ** &usersByMovie, unsigned short ** &whatIndexAMovieIsToAUser) {

	double message1times2 = 2 * (-0.5 * precision->getEX());
	double message1 = -0.5 * precision->getEX();

	// Since most iterations are over user, it helps to have an array access to the Movie DPs
	double* movDp[M];

//*****************************************************************************
// MOVIE DP    MOVIE DP    MOVIE DP    MOVIE DP    MOVIE DP    MOVIE DP
//*****************************************************************************

	double** perThreadMessage[4];

#pragma omp parallel for
	for (int i = 0; i < 4; i++) {
		perThreadMessage[i] = new double*[M];
		for (movie m = 0; m < M; m++) {
			perThreadMessage[i][m] = new double[2];
			perThreadMessage[i][m][0] = 0.0;
			perThreadMessage[i][m][1] = 0.0;

			movDp[m] = DPmovies->getMoments(m, thisK);
		}
	}
#pragma omp barrier

#pragma omp parallel for
	for (user u = 0; u < U; u++) {
		int noViewings = noViewingsPerUser[u];
		double** myMovieMessageArray = perThreadMessage[omp_get_thread_num()];
		q* myQs = qsByUser[u];
		double* myMessages = messages0[u];
		double* myDpMoms = dpMomRecUrIx0[u];

		double * us = DPusers->getMoments(u, thisK);
		double * w = wRecKUr[thisK][u];
		double uw0 = w[0] + us[0];
		double uw1 = us[1] + w[1] + 2 * us[0] * w[0];

		for (int index = 0; index < noViewings; index++) {
			movie m = myQs[index] % MOV;
			myDpMoms[index] -= (uw0) * movDp[m][0];
			double* movieMsgs = myMovieMessageArray[m];
			movieMsgs[0] += (myMessages[index]
					+ message1times2 * myDpMoms[index]) * uw0;
			movieMsgs[1] += uw1;
		}
	}
#pragma omp barrier

#pragma omp parallel for
	for (movie m = 0; m < M; m++) {
		for (int thread = 1; thread < 4; thread++) {
			perThreadMessage[0][m][0] += perThreadMessage[thread][m][0];
			perThreadMessage[0][m][1] += perThreadMessage[thread][m][1];
		}
		perThreadMessage[0][m][1] *= message1;
	}
#pragma omp barrier

	if (debug) {
		double correctMsg[2] = { 0, 0 };
		double mod = 0.0;
		for (movie m = 0; m < M; m += 50) {
			calcMessageToMovieDP(qsByUser, noViewingsPerUser, m, thisK,
					correctMsg);
			mod = fabs(correctMsg[0] - perThreadMessage[0][m][0]);
			mod += fabs(correctMsg[1] - perThreadMessage[0][m][1]);
			if (mod / 2.0 > 1e-10) {
				cerr << "\nBig error in message to MovieDP, got: "
						<< perThreadMessage[0][m][0] << ","
						<< perThreadMessage[0][m][1] << " instead of "
						<< correctMsg[0] << "," << correctMsg[1]
						<< " for movie " << m << '\n';
			}
		}
	}

#pragma omp parallel for
	for (movie m = 0; m < M; m++) {
		DPmovies->update(m, thisK, 0, perThreadMessage[0][m]);
	}
#pragma omp barrier

//*****************************************************************************
// USER DP    USER DP    USER DP    USER DP    USER DP    USER DP    USER DP    
//*****************************************************************************
#pragma omp parallel for
	for (user u = 0; u < U; u++) {
		bigTempMessageArray[u][0] = 0;
		bigTempMessageArray[u][1] = 0;
	}
#pragma omp barrier

#pragma omp parallel for 
	for (user u = 0; u < U; u++) {
		int myViewings = noViewingsPerUser[u];
		q* myQs = qsByUser[u];

		double* myMessages = messages0[u];
		double* thisUserMessage = bigTempMessageArray[u];
		double* myDpMoms = dpMomRecUrIx0[u];

		for (uint index = 0; index < myViewings; index++) {
			double* ms = movDp[myQs[index] % MOV];
			thisUserMessage[0] += (myMessages[index]
					+ message1times2 * myDpMoms[index]) * ms[0];
			thisUserMessage[1] += ms[1];
		}
		thisUserMessage[1] *= message1;

		// IMPORTANT POINT:
		// NOTE:
		// At this point the message is perfect to send towards the W matrix.
		thisUserMessage[0] += 2 * thisUserMessage[1] * wRecKUr[thisK][u][0];
	}
#pragma omp barrier

#pragma omp parallel for
	for (user u = 0; u < U; u++) {
		DPusers->update(u, thisK, 0, bigTempMessageArray[u]);
	}
#pragma omp barrier

	if (debug) {
		double correctMsg[2] = { 0, 0 };
		double mod = 0.0;
		for (user u = 0; u < U; u += 50) {
			calcMessageToUserDP(qsByUser, noViewingsPerUser, u, thisK,
					correctMsg);

			mod = fabs(correctMsg[0] - bigTempMessageArray[u][0]);
			mod += fabs(correctMsg[1] - bigTempMessageArray[u][1]);
			if (mod / 2.0 > 1e-10) {
				cerr << "\nBig error in message to UserDP, got: "
						<< bigTempMessageArray[u][0] << ","
						<< bigTempMessageArray[u][1] << " instead of "
						<< correctMsg[0] << "," << correctMsg[1] << " for user "
						<< u << '\n';
			}
		}
	}

#pragma omp parallel for
	for (user u = 0; u < U; u++) {
		// reconstruct the user-side message (see previous comment)
		bigTempMessageArray[u][0] -= 2.0 * bigTempMessageArray[u][1]
				* wRecKUr[thisK][u][0];
	}
#pragma omp barrier

//*****************************************************************************
// W   W   W   W   W   W   W   W   W   W   W   W   W   W   W   W   W   W   W   W    
//*****************************************************************************

	// Pass messages past DP users, and through the multiplier of 1/rootViewings
	double rootViewings[U];
#pragma omp parallel for
	for (user u = 0; u < U; u++) {
		double viewings = ((double) noViewingsPerUser[u]);
		double sqrtViewings = pow(viewings, 0.5);
		bigTempMessageArray[u][0] = (bigTempMessageArray[u][0]
				+ 2 * bigTempMessageArray[u][1]
						* DPusers->getMoments(u, thisK)[0]) / sqrtViewings;
		bigTempMessageArray[u][1] /= (viewings);
		rootViewings[u] = sqrtViewings;
	}
#pragma omp barrier

	double toSend[2] = { 0, 0 };

	for (movie wmovie = 0; wmovie < M; wmovie++) {
		double msg0 = 0.0;
		double msg1 = 0.0;
		int mViewings = noViewingsPerMovie[wmovie];
		double wMovie0 = W->getMoments(wmovie, thisK)[0];
		double** thisW = wRecKUr[thisK];
		user* myUsers = usersByMovie[wmovie];

#pragma omp parallel for reduction(+:msg0) reduction(+:msg1)
		for (int i = 0; i < mViewings; i++) {
			user u = myUsers[i];
			double sqrtViewings = rootViewings[u];
			double* thisMsgs = bigTempMessageArray[u];
			thisW[u][0] -= wMovie0 / sqrtViewings;
			msg0 += thisMsgs[0] + 2 * sqrtViewings * thisMsgs[1] * thisW[u][0]; // uViewings included because thisW includes the divide by number of viewings
			msg1 += thisMsgs[1]; // /(uViewings*uViewings);
		}
#pragma omp barrier

		toSend[0] = msg0;
		toSend[1] = msg1;
		W->update(wmovie, thisK, 0, toSend);

		if (debug && (wmovie % 50 == 0)) {
			double correctMsg[2] = { 0, 0 };
			calcMessageToW(wmovie, thisK, correctMsg, qsByUser,
					noViewingsPerUser, noViewingsPerMovie, usersByMovie,
					whatIndexAMovieIsToAUser);

			double mod0 = fabs(correctMsg[0] - toSend[0]);
			double mod1 = fabs(correctMsg[1] - toSend[1]);
			if ((mod0 > 1e-10 || mod1 > 1e-10)) {
				cerr << "\n Big error in messages to W for wmovie " << wmovie
						<< " K:" << thisK << " is (" << mod0 << "," << mod1
						<< ")\n";
				cerr << "(" << toSend[0] << "," << toSend[1] << ") should be: ("
						<< correctMsg[0] << "," << correctMsg[1] << ")\n";
				cerr << "This movie has " << noViewingsPerMovie[wmovie]
						<< " viewings.\n";
			}
		}

		wMovie0 = W->getMoments(wmovie, thisK)[0];

#pragma omp parallel for 
		for (int i = 0; i < mViewings; i++) {
			double* thisMsgs = bigTempMessageArray[myUsers[i]];
			thisMsgs[0] += 2.0 * thisMsgs[1] * wMovie0;
		}
#pragma omp barrier
	}

	refreshWStore(qsByUser, noViewingsPerUser, thisK);

#pragma omp parallel for 
	for (user u = 0; u < U; u++) {
		double* myMsgs = messages0[u];
		q* myQs = qsByUser[u];
		int myViewings = noViewingsPerUser[u];
		double uw0TimesMessage1Times2 = (DPusers->getMoments(u, thisK)[0]
				+ wRecKUr[thisK][u][0]) * message1times2;

		for (uint index = 0; index < myViewings; index++) {
			myMsgs[index] += movDp[myQs[index] % MOV][0]
					* uw0TimesMessage1Times2;
		}
	}
#pragma omp barrier

#pragma omp parallel for
	for (int i = 0; i < 4; i++) {
		for (movie m = 0; m < M; m++) {
			delete[] perThreadMessage[i][m];
		}
		delete[] perThreadMessage[i];
	}
#pragma omp barrier

}

// Functions for maintaining/accessing private statistics
void DotProductGaussianAddNG::refreshWStore(q** &qsByUser,
		ushort* &ratingsPerUser) {
	// For each user, for each K get the W contribution
#pragma omp parallel for 
	for (user u = 0; u < U; u++) {
		for (int thisK = 0; thisK < K; thisK++) {
			makeWContributionToUserSide(qsByUser, ratingsPerUser, u, thisK,
					wRecKUr[thisK][u]);
		}
	}
#pragma omp barrier
}

void DotProductGaussianAddNG::refreshWStore(q** &qsByUser,
		ushort* &ratingsPerUser, uint thisK) {
	// For each user, for this K, make the W contributions
#pragma omp parallel for
	for (user u = 0; u < U; u++) {
		makeWContributionToUserSide(qsByUser, ratingsPerUser, u, thisK,
				wRecKUr[thisK][u]);
	}
#pragma omp barrier
}

void DotProductGaussianAddNG::makeWContributionToUserSide(q** &qsByUser,
		ushort* &ratingsPerUser, user u, uint thisK, double inHere[]) {
	inHere[0] = 0.0;
	inHere[1] = 0.0;

	for (int index = 0; index < ratingsPerUser[u]; index++) {
		AStats(inHere, W->getMoments(qsByUser[u][index] % MOV, thisK));
	}

	MStats(inHere, 1.0 / pow((double) ratingsPerUser[u], 0.5));
}

void DotProductGaussianAddNG::addADpStatsBit(user userID, movie movieID,
		double statsSoFar[], uint kVal) {
	// Add one K's contribution to the DP moments for a specified (user, movie)
	double* dpu = DPusers->getMoments(userID, kVal);
	double * dpm = DPmovies->getMoments(movieID, kVal);
	double* wu = wRecKUr[kVal][userID];
	AStats(statsSoFar, (dpu[0] + wu[0]) * dpm[0],
			(dpu[1] + wu[1] + 2 * wu[0] * dpu[0]) * dpm[1]);
}

void DotProductGaussianAddNG::makeAllDPStats(user userID, movie movieID,
		double inHere[]) {
	// Get the DP stats for a given (user, movie)
	inHere[0] = 0;
	inHere[1] = 0;
	for (unsigned int i = 0; i < K; i++) {
		addADpStatsBit(userID, movieID, inHere, i);
	}
}

void DotProductGaussianAddNG::refreshDPStore(q** &qsByUser,
		unsigned short* &ratingsPerUser) {
#pragma omp parallel for
	for (user u = 0; u < U; u++) {
		int noViewings = ratingsPerUser[u];
		q* thisQs = qsByUser[u];
		double* thisDpMomRec = dpMomRecUrIx0[u];

		for (int i = 0; i < noViewings; i++) {
			thisDpMomRec[i] = 0.0;
			movie m = thisQs[i] % MOV;

			for (int k = 0; k < K; k++) {
				thisDpMomRec[i] += (DPusers->getMoments(u, k)[0]
						+ wRecKUr[k][u][0]) * DPmovies->getMoments(m, k)[0];
			}
		}
	}
#pragma omp barrier
}

void DotProductGaussianAddNG::getDPMomentsKnownMovie(user userID, movie movieID,
		double inHere[]) {
	makeAllDPStats(userID, movieID, inHere);
}

// Debugging functions to calculate (slowly) the correct message to each node
void DotProductGaussianAddNG::calcMessageToAllAdd(q** &qsByUser,
		unsigned short* &ratingsPerUser, double* message) {
	double m0 = 0;
	double m1 = 0;

#pragma omp parallel for reduction(+:m0) reduction(+:m1)
	for (user u = 0; u < U; u++) {
		for (int index = 0; index < ratingsPerUser[u]; index++) {
			double msg[2];
			q thisQ = qsByUser[u][index];
			getMsgFromG(thisQ % RATE, msg);
			passMsgThruA(msg, movieAdds->getMoments(thisQ % MOV, 0));
			passMsgThruA(msg, userAdds->getMoments(u, 0));

			for (int k = 0; k < K; k++) {
				double temp[2];
				makeWContributionToUserSide(qsByUser, ratingsPerUser, u, k,
						temp);
				AStats(temp, DPusers->getMoments(u, k));
				MStats(temp, DPmovies->getMoments(thisQ % MOV, k));
				passMsgThruA(msg, temp);
			}

			m0 += msg[0];
			m1 += msg[1];

		}
	}
#pragma omp barrier

	message[0] = m0;
	message[1] = m1;
}

void DotProductGaussianAddNG::calcMessageToUserAdd(q** &qsByUser,
		unsigned short* &ratingsPerUser, user u, double* message) {
	double *allAddsMoments = allAdds->getMoments();
	double m0 = 0.0;
	double m1 = 0.0;

#pragma omp parallel for reduction(+:m0) reduction (+:m1)
	for (int index = 0; index < ratingsPerUser[u]; index++) {
		double msgTemp[] = { 0.0, 0.0 };
		movie m = qsByUser[u][index] % MOV;
		getMsgFromG(qsByUser[u][index] % RATE, msgTemp);
		passMsgThruA(msgTemp, allAddsMoments);
		passMsgThruA(msgTemp, movieAdds->getMoments(m, 0));
		for (int k = 0; k < K; k++) {
			double temp[2];
			makeWContributionToUserSide(qsByUser, ratingsPerUser, u, k, temp);
			AStats(temp, DPusers->getMoments(u, k));
			MStats(temp, DPmovies->getMoments(m, k));
			passMsgThruA(msgTemp, temp);
		}
		m0 += msgTemp[0];
		m1 += msgTemp[1];
	}
#pragma omp barrier

	message[0] = m0;
	message[1] = m1;
}

void DotProductGaussianAddNG::calcMessageToMovieAdd(q** &qsByUser,
		unsigned short* &ratingsPerUser, movie m, double* message) {
	double* allAddsMoments = allAdds->getMoments();
	double m0 = 0.0;
	double m1 = 0.0;

#pragma omp parallel for reduction(+:m0) reduction(+:m1)
	for (user u = 0; u < U; u++) {
		for (int index = 0; index < ratingsPerUser[u]; index++) {
			if ((qsByUser[u][index] % MOV) == m) {
				double msgTemp[2] = { 0.0, 0.0 };

				getMsgFromG(qsByUser[u][index] % RATE, msgTemp);
				passMsgThruA(msgTemp, allAddsMoments);
				passMsgThruA(msgTemp, userAdds->getMoments(u, 0));
				for (int k = 0; k < K; k++) {
					double temp[2];
					makeWContributionToUserSide(qsByUser, ratingsPerUser, u, k,
							temp);
					AStats(temp, DPusers->getMoments(u, k));
					MStats(temp, DPmovies->getMoments(m, k));
					passMsgThruA(msgTemp, temp);
				}

				m0 += msgTemp[0];
				m1 += msgTemp[1];
			}
		}
	}
#pragma omp barrier

	message[0] = m0;
	message[1] = m1;
}

void DotProductGaussianAddNG::calcMessageToUserDP(q** &qsByUser,
		ushort* &ratingsPerUser, user u, uint thisK, double* message) {
	double* allAddsMoments = allAdds->getMoments();
	double* userAddsMoments = userAdds->getMoments(u, 0);

	double m0 = 0.0;
	double m1 = 0.0;

#pragma omp parallel for reduction(+:m0) reduction(+:m1)
	for (int index = 0; index < ratingsPerUser[u]; index++) {
		double msgTemp[2] = { 0.0, 0.0 };
		q thisQ = qsByUser[u][index];
		getMsgFromG(thisQ % RATE, msgTemp);
		passMsgThruA(msgTemp, allAddsMoments);
		passMsgThruA(msgTemp, movieAdds->getMoments(thisQ % MOV, 0));
		passMsgThruA(msgTemp, userAddsMoments);

		for (unsigned int kNum = 0; kNum < K; kNum++) {
			if (kNum == thisK) {
				continue;
			}

			double temp[2];
			makeWContributionToUserSide(qsByUser, ratingsPerUser, u, kNum,
					temp);
			AStats(temp, DPusers->getMoments(u, kNum));
			MStats(temp, DPmovies->getMoments(thisQ % MOV, kNum));
			passMsgThruA(msgTemp, temp);
		}

		passMsgThruM(msgTemp, DPmovies->getMoments(thisQ % MOV, thisK));

		m0 += msgTemp[0];
		m1 += msgTemp[1];
	}
#pragma omp barrier

	message[0] = m0;
	message[1] = m1;

	for (int index = 0; index < ratingsPerUser[u]; index++) {
		passMsgThruA(message,
				W->getMoments(qsByUser[u][index] % MOV, thisK)[0]
						/ pow((double) ratingsPerUser[u], 0.5));
	}
}

void DotProductGaussianAddNG::calcMessageToW(movie wmovie, uint thisK,
		double* message, q** &qsByUser, ushort* &ratingsPerUser,
		ulong* &ratingsPerMovie, user** &usersByMovie,
		ushort** &whatIndexAMovieIsToAUser) {
	double m0 = 0.0;
	double m1 = 0.0;

	double *allAddsMoments = allAdds->getMoments();

#pragma omp parallel for reduction(+:m0) reduction(+:m1)
	for (uint index = 0; index < ratingsPerMovie[wmovie]; index++) {
		user u = usersByMovie[wmovie][index];

		for (int ui = 0; ui < ratingsPerUser[u]; ui++) {
			double msgTemp[2] = { 0.0, 0.0 };
			movie m = qsByUser[u][ui] % MOV;
			getMsgFromG(qsByUser[u][ui] % RATE, msgTemp);
			passMsgThruA(msgTemp, allAddsMoments);
			passMsgThruA(msgTemp, movieAdds->getMoments(m, 0));
			passMsgThruA(msgTemp, userAdds->getMoments(u, 0));

			for (uint kNum = 0; kNum < K; kNum++) {
				if (kNum == thisK) {
					continue;
				}

				double temp[2];
				makeWContributionToUserSide(qsByUser, ratingsPerUser, u, kNum,
						temp);
				AStats(temp, DPusers->getMoments(u, kNum));
				MStats(temp, DPmovies->getMoments(m, kNum));
				passMsgThruA(msgTemp, temp);
			}

			passMsgThruM(msgTemp, DPmovies->getMoments(m, thisK));
			passMsgThruA(msgTemp, DPusers->getMoments(u, thisK));

			passMsgThruM(msgTemp, 1 / (pow((double) ratingsPerUser[u], 0.5)));

			for (uint ii = 0; ii < ratingsPerUser[u]; ii++) {
				if ((qsByUser[u][ii] % MOV) == wmovie) {
					continue;
				}
				passMsgThruA(msgTemp,
						W->getMoments(qsByUser[u][ii] % MOV, thisK)[0]);
			}

			m0 += msgTemp[0];
			m1 += msgTemp[1];
		}
	}
#pragma omp barrier

	message[0] = m0;
	message[1] = m1;
}

void DotProductGaussianAddNG::calcMessageToMovieDP(q** &qsByUser,
		unsigned short* &ratingsPerUser, movie m, unsigned int thisK,
		double* message) {
	double* allAddsMoments = allAdds->getMoments();
	double m0 = 0.0;
	double m1 = 0.0;

#pragma omp parallel for reduction(+:m0) reduction(+:m1)
	for (user u = 0; u < U; u++) {
		for (int index = 0; index < ratingsPerUser[u]; index++) {
			if ((qsByUser[u][index] % MOV) == m) {
				double msgTemp[2] = { 0.0, 0.0 };
				getMsgFromG(qsByUser[u][index] % RATE, msgTemp);
				passMsgThruA(msgTemp, allAddsMoments);
				passMsgThruA(msgTemp, userAdds->getMoments(u, 0));
				passMsgThruA(msgTemp, movieAdds->getMoments(m, 0));
				for (unsigned int kNum = 0; kNum < K; kNum++) {
					if (kNum == thisK) {
						continue;
					}

					double temp[2];
					makeWContributionToUserSide(qsByUser, ratingsPerUser, u,
							kNum, temp);
					AStats(temp, DPusers->getMoments(u, kNum));
					MStats(temp, DPmovies->getMoments(m, kNum));
					passMsgThruA(msgTemp, temp);
				}

				double temp[2];
				makeWContributionToUserSide(qsByUser, ratingsPerUser, u, thisK,
						temp);
				AStats(temp, DPusers->getMoments(u, thisK));
				passMsgThruM(msgTemp, temp);

				m0 += msgTemp[0];
				m1 += msgTemp[1];
			}
		}
	}
#pragma omp barrier

	message[0] = m0;
	message[1] = m1;
}

// Helper functions for calculating correct messages/ statistics
void DotProductGaussianAddNG::passMsgThruA(double*message,
		EXandX2* otherParent) {
	message[0] += 2 * message[1] * otherParent->getEX();
}

void DotProductGaussianAddNG::passMsgThruA(double* message, double*stats) {
	message[0] += 2 * message[1] * stats[0];
}

void DotProductGaussianAddNG::passMsgThruA(double* message, double stats0) {
	message[0] += 2 * message[1] * stats0;
}

void DotProductGaussianAddNG::passMsgThruM(double* message, double* stats) {
	message[0] *= stats[0];
	message[1] *= stats[1];
}

void DotProductGaussianAddNG::passMsgThruM(double* message, double stats0,
		double stats1) {
	message[0] *= stats0;
	message[1] *= stats1;
}

void DotProductGaussianAddNG::passMsgThruM(double* message, double constant0) {
	message[0] *= constant0;
	message[1] *= constant0 * constant0;
}

void DotProductGaussianAddNG::getMsgFromG(rate rating, double &message0) {
	message0 = precision->getEX() * (static_cast<double>(rating));
}

void DotProductGaussianAddNG::getMsgFromG(rate rating, double* message) {
	message[0] = precision->getEX() * ((double) rating);
	message[1] = (-0.5 * precision->getEX());
}

void DotProductGaussianAddNG::UnAddStats(double* stats, double other0,
		double other1) {
	stats[0] -= other0;
	stats[1] -= (other1 + 2 * other0 * stats[0]);
}

void DotProductGaussianAddNG::AStats(double* stats, EXandX2* otherAdd) {
	stats[1] += otherAdd->getEX2() + 2 * stats[0] * otherAdd->getEX();
	stats[0] += otherAdd->getEX();
}

void DotProductGaussianAddNG::AStats(double* stats, double* otherAdd) {
	stats[1] += otherAdd[1] + 2 * stats[0] * otherAdd[0];
	stats[0] += otherAdd[0];
}

void DotProductGaussianAddNG::AStats(double* stats, double other0,
		double other1) {
	stats[1] += other1 + 2 * stats[0] * other0;
	stats[0] += other0;
}

void DotProductGaussianAddNG::MStats(double*stats, double*otherStats) {
	stats[0] *= otherStats[0];
	stats[1] *= otherStats[1];
}

void DotProductGaussianAddNG::MStats(double*stats, double constant0) {
	stats[0] *= constant0;
	stats[1] *= constant0 * constant0;
}

void DotProductGaussianAddNG::MStats(double*stats, EXandX2* otherStats) {
	stats[0] *= otherStats->getEX();
	stats[1] *= otherStats->getEX2();
}

}

