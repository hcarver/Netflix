#ifndef NETFLIXDATANG_H_
#define NETFLIXDATANG_H_

#include "nfx.h"
#include "DotProductGaussianAddNG.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>

using namespace std;

namespace utils {

class NetflixDataNG {

private:
	// String paths for accessing data
	string basePathToUse;
	string dataPath;
	string qualDataFile;

	// Array for conversion between internal and external user indices
	int userIDConversion[U];

	// Booleans for indicating whether running in test or debug mode
	bool test;
	bool debug;

	// Reusable stream reader
	ifstream individualReader;

	// Array for keeping calendar information
	short daysAfter30Sept1998[8][12];

public:
	// Constructor and Destructor
	NetflixDataNG();
	~NetflixDataNG();

	// File management public helper functions
	void deleteFileAsynch(string fileName);
	void zipFileAsynch(string fileName, string newFileName);
	void zipFileSynch(string fileName, string newFileName);
	void unzipFileSynch(string fileName, string newFileName);
	void makeEntryAsynch(string submissionFileName);

	// Data read functions
	void getDataMeanAndVariance(double mAndV[]);

	void allViewingsUsers(movie** &moviesByUser, rate** &ratesByUser,
			ushort* &noViewingsPerUser);

	void allViewingsMovies(user** &usersByMovie, ulong* &noViewingsPerMovie,
			rate** &ratesByMovie, rate** &ratesByUser, movie** &moviesByUser,
			ushort* &noViewingsPerUser, ushort ** &whatIndexAMovieIsToAUser);

	void allViewingsMovies(user** &usersByMovie, ulong* &noViewingsPerMovie,
			rate** &ratesByMovie, uint** &datesByMovie, rate** &ratesByUser,
			movie** &moviesByUser, ushort* &noViewingsPerUser,
			uint** &datesByUser, ushort** &whatIndexAMovieIsToAUser);

	void readBigRatings(ushort* &noViewingsPerUser, q** &bigRatingsByUser);

	void readBigRatings(user** &usersByMovie, ulong* &noViewingsPerMovie,
			rate** &ratesByMovie, q** &qsByUser, ushort* &noViewingsPerUser,
			ushort** &whatIndexAMovieIsToAUser);

	// Data utility functions
	string getBaseSaveFileName();

	void getAllBigPredictionData(vector<user> &users, vector<q> &preds);

	void getAllPredictionData(vector<user> &users, vector<movie> &movies,
			vector<uint> &dates);

	void writeQualifyingFile(string qualifyingFileName,
			string fullStatsFileName, model::DotProductGaussianAddNG* node);

	// Make big rating
	q GetBigRating(rate rating, movie m, uint userEra, uint movieEra,
			uint weekday, uint precisionNumber, uint useHistory1,
			uint useHistory3, uint useHistory5, uint useHistory10,
			uint useHistory30);

	void updateBigRating(q &currentBigRatingValue, q &currentBaseValue,
			q newValue, q newBaseValue);

	// Setting internal debug/test variables
	void setDebug();
	void setTest();

private:
	// Helper functions for IO
	void getQuartetFromUserData(movie& val, rate& r, string& fromThis,
			char* pEnd);
	void getQuartetFromMovieData(user& val, rate& r, string& fromThis,
			char* pEnd);
	void getInfoFromMovieData(user& val, rate& r,
			uint& daysAfter30September1998, string& fromThis, char* pEnd);
	uint getDaysAfter30September1998(uint year, uint month, uint day);

	uint parseForInt(string s, uint startAt, uint endAt);
	void sevenDigitFormat(stringstream &s, int number);

	// Helper functions
	void populateUserIdConverter(string path);
	string getTrainingDataPath(uint userMovieNumber, bool trueIfUser);

	// Functions for internal/external numbering
	user convertNxUserToHcUser(uint netflixUserID);
	uint convertHcUserToNxUser(user hcUserId);
	movie convertNxMovieToHcMovie(uint netflixMovieID);
	uint convertHcMovieToNxMovie(movie HcMovieID);

};

}

#endif /*NETFLIXDATANG_H_*/
