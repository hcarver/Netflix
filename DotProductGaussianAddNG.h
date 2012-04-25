#ifndef DOTPRODUCTGAUSSIANADDNG_H_
#define DOTPRODUCTGAUSSIANADDNG_H_

#include "MassGaussianNode.h"
#include "EXandX2.h"
#include "EXandLnX.h"
#include "nfx.h"
#include "time.h"

namespace model {

class DotProductGaussianAddNG {

private:
	// Node variables
	MassGaussianNode* DPusers;
	MassGaussianNode* DPmovies;
	MassGaussianNode* userAdds;
	MassGaussianNode* movieAdds;
	MassGaussianNode* W;
	MassGaussianNode* uhm1;
	MassGaussianNode* uhm3;
	MassGaussianNode* uhm5;
	MassGaussianNode* uhm10;
	MassGaussianNode* uhm30;
	EXandX2* allAdds;
	EXandLnX* precision;

	// Private variable to set whether to run in debug mode
	bool debug;

	// Private variables to help calculate messages
	double*** wRecKUr; // W moments, by user and K, 2 moments
	double** messages0; // Current estimate of the rating, by user and index
	double** dpMomRecUrIx0; // Current estimate of the E(DP) for a certain user and index
	double** bigTempMessageArray; // Useful array for messages

public:
	// Constructor and Destructor
	DotProductGaussianAddNG(bool isDebug);
	~DotProductGaussianAddNG();

	// Initialisation function
	void init(q** &qsByUser, ushort* &ratingsPerUser);

	// Functions to set the nodes
	void setUsers(MassGaussianNode* gn);
	void setMovie(MassGaussianNode* gn);
	void setUserAdd(MassGaussianNode* n);
	void setMovieAdd(MassGaussianNode* n);
	void setW(MassGaussianNode* n);
	void setAllAdd(EXandX2* n);
	void setPrecision(EXandLnX* n);
	void setUserHistoryMultipliers(MassGaussianNode* uhm1,
			MassGaussianNode* uhm3, MassGaussianNode* uhm5,
			MassGaussianNode* uhm10, MassGaussianNode* uhm30);

	// Inference function
	void doAnIteration(q** &qsByUser, ushort* &noViewingsPerUser,
			ulong* &ratingsPerMovie, user ** &usersByMovie,
			ushort ** &whatIndexAMovieIsToAUser, bool updatePrecision);

	// Functions to make predictions (either for training or test ratings)
	void predictKnownIndex(double* stats, user userID, uint index,
			q** &moviesByUser);
	void predictUnknownIndex(double* stats, user userID, movie movieID);

private:
	// Parts of the inference algorithm
	void toPrecision(q** &qsByUser, ushort* &ratingsPerUser);
	void toAllAdds(q** &qsByUser, ushort* &ratingsPerUser);
	void pastUserAdds(q** &qsByUser, ushort* &ratingsPerUser);
	void pastMovieAdds(q** &qsByUser, ushort* &ratingsPerUser);
	void pastUHMs(q** &qsByUser, ushort* &ratingsPerUser);
	void pastDP(q** &qsByUser, ushort* &noViewingsPerUser,
			ulong* &ratingsPerMovie, user** &usersByMovie,
			ushort** &whatIndexAMovieIsToAUser);
	void pastDPUM(int thisK, q** &qsByUser, ushort* &ratingsPerUser,
			ulong* &ratingsPerMovie, user** &usersByMovie,
			ushort** &whatIndexAMovieIsToAUser);

	// Functions for maintaining/accessing private statistics
	void refreshWStore(q** &qsByUser, ushort* &ratingsPerUser);
	void refreshWStore(q** &qsByUser, ushort* &ratingsPerUser, uint thisK);
	void makeWContributionToUserSide(q** &qsByUser, ushort* &ratingsPerUser,
			user u, uint thisK, double inHere[]);
	void addADpStatsBit(user userID, movie movieID, double statsSoFar[],
			uint kVal);
	void makeAllDPStats(user userID, movie movieID, double *inHere);
	void refreshDPStore(q** &qsByUser, ushort* &ratingsPerUser);
	void getDPMomentsKnownMovie(user userID, movie movieID, double inHere[]);

	// Debugging functions to calculate (slowly) the correct message to each node
	void calcMessageToAllAdd(q** &qsByUser, ushort* &ratingsPerUser,
			double* message);
	void calcMessageToUserAdd(q** &qsByUser, ushort* &ratingsPerUser, user u,
			double* message);
	void calcMessageToMovieAdd(q** &qsByUser, ushort* &ratingsPerUser, movie m,
			double* message);
	void calcMessageToUHM1(q** &qsByUser, ushort* &ratingsPerUser, user u,
			double* message);
	void calcMessageToUHM3(q** &qsByUser, ushort* &ratingsPerUser, user u,
			double* message);
	void calcMessageToUHM5(q** &qsByUser, ushort* &ratingsPerUser, user u,
			double* message);
	void calcMessageToUHM10(q** &qsByUser, ushort* &ratingsPerUser, user u,
			double* message);
	void calcMessageToUHM30(q** &qsByUser, ushort* &ratingsPerUser, user u,
			double* message);
	void calcMessageToUserDP(q** &qsByUser, ushort* &ratingsPerUser, user u,
			uint thisK, double* message);
	void calcMessageToMovieDP(q** &qsByUser, ushort* &ratingsPerUser, movie m,
			uint thisK, double* message);
	void calcMessageToW(movie wmovie, uint thisK, double* message,
			q** &qsByUser, ushort* &ratingsPerUser, ulong* &ratingsPerMovie,
			user** &usersByMovie, ushort** &whatIndexAMovieIsToAUser);

	// Helper functions for calculating correct messages/ statistics
	void passMsgThruA(double* message, EXandX2* otherParent);
	void passMsgThruA(double* message, double* stats);
	void passMsgThruA(double* message, double stats0);
	void passMsgThruM(double* message, double* stats);
	void passMsgThruM(double* message, double stats0, double stats1);
	void passMsgThruM(double* message, double constant0);
	void getMsgFromG(rate rating, double &message0);
	void getMsgFromG(rate rating, double* stats);

	void UnAddStats(double* stats, double other0, double other1);
	void MStats(double* stats, double* otherStats);
	void MStats(double* stats, EXandX2* otherStats);
	void MStats(double* stats, double constant0);
	void AStats(double* stats, double* otherAdd);
	void AStats(double* stats, EXandX2* otherAdd);
	void AStats(double* stats, double other0, double other1);

};

}

#endif /*DOTPRODUCTGAUSSIANADDNG_H_*/

