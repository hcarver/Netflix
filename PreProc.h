#ifndef PREPROC_H_
#define PREPROC_H_

#include "nfx.h"
#include "NetflixDataNG.h"
using namespace utils;

namespace tool {

class PreProc {

private:
	NetflixDataNG data;

	// Variables for getting and storing ratings data.
	movie** moviesByUser;
	rate** ratingsByUser;
	unsigned short* noViewingsPerUser;
	unsigned int** datesByUser;

	unsigned short ** whatIndexAMovieIsToAUser;
	user** usersByMovie;
	unsigned long* noViewingsPerMovie;
	rate** ratesByMovie;
	unsigned int** datesByMovie;

	double** workingCopyOfRatingsByUser;

	q **bigRatings;
	user * bigPredictionUsers;
	q *bigPrediction;
	long int numPredictions;

	// Doubles for processing the ratings etc.
	double mAllAdd, mUserAddMean, mUserAddVar, mMovieAddMean, mMovieAddVar,
			mUserGivenMovieMean, mUserGivenMovieVariance, mWVar, mDPVar;

public:
	// Constructor and Destructor
	PreProc(bool isDebug, bool isTest);
	~PreProc();

	// Go function
	void go();

private:
	// Functions for processing & calculating the statistics of each model part
	void allAdds();
	void userSideAverages(bool IsWithMovieAveragesGiven);
	void movieSideAverages();
	void wAndDpFit();

	// Function for deciding what size variables to use for the larger model
	void memoryCapacityTest();

	// Function for working out the all-info rating, bit by bit
	void writeBigRatingFile();
	void testBigRatingFile();

	void writeBigPredictionFile();
	void testBigPredictionFile();

	q makeBig(user u, movie m, uint date, rate rating);

	void sortUserStuffByDate(movie* &thisMovies, rate* &thisRatings,
			uint* &thisDates, uint beginIndex, uint endIndex);

	uint partitionUserStuffByDate(movie* &thisMovies, rate* &thisRatings,
			uint* &thisDates, uint beginIndex, uint endIndex);

	void swap(movie* &thisMovies, rate* &thisRatings, uint* &thisDates,
			uint index1, uint index2);

	void sortMovieStuffByDate(user* &thisUsers, rate* &thisRatings,
			uint* &thisDates, uint beginIndex, uint endIndex);

	uint partitionMovieStuffByDate(user* &thisUsers, rate* &thisRatings,
			uint* &thisDates, uint beginIndex, uint endIndex);

	void swap(user* &thisUsers, rate* &thisRatings, uint* &thisDates,
			uint index1, uint index2);

	void checkUserDatesAreInOrder();
	void checkMovieDatesAreInOrder();

	uint getUserRatingHistoryFromDate(user u, uint targetDate,
			uint numOfDaysHistory);

	uint getUserEra(user u, uint ratingDate, uint upperLimit);
	uint getMovieEra(movie m, uint ratingDate, uint upperLimit);
};

}
#endif /*PREPROC_H_*/
