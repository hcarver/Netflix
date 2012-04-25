#include "PreProc.h"
#include <cstdlib>

namespace tool {

// Constructor and Destructor
PreProc::PreProc(bool isDebug, bool isTest) {
	if (isDebug) {
		data.setDebug();
	}
	if (isTest) {
		data.setTest();
	}
}

PreProc::~PreProc() {
	for (user u = 0; u < U; u++) {
		delete[] moviesByUser[u];
		delete[] ratingsByUser[u];
		delete[] workingCopyOfRatingsByUser[u];
	}

	delete[] moviesByUser;
	delete[] ratingsByUser;
	delete[] workingCopyOfRatingsByUser;

	delete[] noViewingsPerUser;
}

// Go function
void PreProc::go() {
	workingCopyOfRatingsByUser = new double*[U];

	moviesByUser = new movie*[U];
	ratingsByUser = new rate*[U];
	noViewingsPerUser = new ushort[U];
	datesByUser = new uint*[U];

	whatIndexAMovieIsToAUser = new ushort *[M];
	usersByMovie = new user*[M];
	noViewingsPerMovie = new ulong[M];
	ratesByMovie = new rate*[M];
	datesByMovie = new uint*[M];

	bigRatings = new q*[U];

	data.allViewingsMovies(usersByMovie, noViewingsPerMovie, ratesByMovie,
			datesByMovie, ratingsByUser, moviesByUser, noViewingsPerUser,
			datesByUser, whatIndexAMovieIsToAUser);

	for (user u = 0; u < U; u++) {
		sortUserStuffByDate(moviesByUser[u], ratingsByUser[u], datesByUser[u],
				0, noViewingsPerUser[u] - 1);
	}

	checkUserDatesAreInOrder();

	for (movie m = 0; m < M; m++) {
		sortMovieStuffByDate(usersByMovie[m], ratesByMovie[m], datesByMovie[m],
				0, noViewingsPerMovie[m] - 1);
	}

	checkMovieDatesAreInOrder();

	writeBigPredictionFile();
	testBigPredictionFile();

	writeBigRatingFile();
	testBigRatingFile();
}

void PreProc::testBigRatingFile() {
	unsigned short* testViewingsPerUser = new unsigned short[U];
	q** testBigRatings = new q*[U];

	cerr << "Testing big rating file\n";

	data.readBigRatings(testViewingsPerUser, testBigRatings);

	for (user u = 0; u < U; u++) {
		if (testViewingsPerUser[u] != noViewingsPerUser[u]) {
			cerr << "Num viewings wrong for user " << u << " "
					<< testViewingsPerUser[u] << " not " << noViewingsPerUser[u]
					<< "\n";
		}
		for (int i = 0; i < testViewingsPerUser[u]; i++) {
			if (testBigRatings[u][i] != bigRatings[u][i]) {
				cerr << "Big rating wrong for user i " << u << " " << i
						<< " was " << testBigRatings[u][i] << " not "
						<< bigRatings[u][i] << "\n";
			}
		}
	}

	cerr
			<< "Testing big rating file finished!!!\n\n ***If you haven't heard otherwise, everything was hunky-dory***\n";

}

void PreProc::memoryCapacityTest() {
	int NewK = 50;
	int UserTimes = 5;
	int MovieTimes = 30;

	// Ratings
	q** Viewings = new q*[U];
	for (user u = 0; u < U; u++) {
		Viewings[u] = new q[210];
	}

	cout << "Viewings done, size " << sizeof(Viewings) << "\n";
	int h;
	cin >> h;

	//Movie DP
	float **** movieDP = new float***[M];
	for (movie m = 0; m < M; m++) {
		movieDP[m] = new float**[NewK];
		for (int k = 0; k < NewK; k++) {
			movieDP[m][k] = new float*[MovieTimes];
			for (int time = 0; time < MovieTimes; time++) {
				movieDP[m][k][time] = new float[2];
			}
		}
	}

	cout << "movieDP done, size " << sizeof(movieDP) << "\n";
	cin >> h;

	//User DP
	double**** userDP = new double***[U];
	for (user u = 0; u < U; u++) {
		userDP[u] = new double**[NewK];
		for (int k = 0; k < NewK; k++) {
			userDP[u][k] = new double*[UserTimes];
			for (int time = 0; time < UserTimes; time++) {
				userDP[u][k][time] = new double[2];
			}
		}
	}

	cout << "userDP done, size " << sizeof(userDP) << "\n";
	cin >> h;

	//W
	float*** W = new float**[M];
	for (movie m = 0; m < M; m++) {
		W[m] = new float*[K];
		for (int k = 0; k < K; k++) {
			W[m][k] = new float[2];
		}
	}

	cout << "W done, size " << sizeof(W) << "\n";
	cin >> h;

	//Movie Adds
	float*** movieAdd = new float**[M];
	for (movie m = 0; m < M; m++) {
		movieAdd[m] = new float*[MovieTimes];
		for (int time = 0; time < MovieTimes; time++) {
			movieAdd[m][time] = new float[2];
		}
	}

	cout << "movieAdd done, size " << sizeof(movieAdd) << "\n";
	cin >> h;

	//User DP
	float *** userAdd = new float**[U];
	for (user u = 0; u < U; u++) {
		userAdd[u] = new float*[UserTimes];
		for (int time = 0; time < UserTimes; time++) {
			userAdd[u][time] = new float[2];
		}
	}

	cout << "userAdd done, size " << sizeof(userAdd) << "\n";
	cin >> h;
}

// Functions for processing & calculating the statistics of each model part
void PreProc::allAdds() {
	// First find the average value (i.e. numerical value of allAdd
	double allAdd = 0.0;
	int numRatings = 0;

	for (user u = 0; u < U; u++) {
		for (int index = 0; index < noViewingsPerUser[u]; index++) {
			allAdd += ratingsByUser[u][index];
			numRatings++;
		}
	}

	allAdd = allAdd / (double) (numRatings);
	cerr << "\n\nAll Add: " << allAdd << "\n";
	mAllAdd = allAdd;

	// Now remove the average value from all of the working copy
	for (user u = 0; u < U; u++) {
		for (uint index = 0; index < noViewingsPerUser[u]; index++) {
			workingCopyOfRatingsByUser[u][index] =
					((double) (ratingsByUser[u][index])) - allAdd;
		}
	}

	cerr << "Taken allAdds off the working copy\n";
}

void PreProc::userSideAverages(bool isWithMovieAveragesGiven) {
	// First calculate the average for each user
	double* userAverage = new double[U];

	double thisRates;
	int thisCount;

	for (user u = 0; u < U; u++) {
		thisRates = 0.0;
		thisCount = 0;
		for (int index = 0; index < noViewingsPerUser[u]; index++) {
			thisRates += workingCopyOfRatingsByUser[u][index];
			thisCount++;
		}
		userAverage[u] = thisRates / (double) (thisCount);
	}

	// Now calculate the statistics of the user averages
	thisRates = 0.0;
	double thisSquared = 0.0;

	for (user u = 0; u < U; u++) {
		thisRates += userAverage[u];
		thisSquared += userAverage[u] * userAverage[u];
	}

	double mean = thisRates / (double) (U);
	double biasedVar = (thisSquared / (double) (U)) - mean * mean;
	double unbiasedVar = biasedVar / (1.0 - 1.0 / (double) (U));

	// Output statistics, and save to member variables
	if (!isWithMovieAveragesGiven) {
		cerr << "User Add: " << mean << "," << unbiasedVar << "\n";
		mUserAddMean = mean;
		mUserAddVar = unbiasedVar;
	} else {
		cerr << "User Add with Movie Averages given: " << mean << ","
				<< unbiasedVar << "\n";
		mUserGivenMovieMean = mean;
		mUserGivenMovieVariance = unbiasedVar;

		// If this is the call where we've already removed movie averages, we remove the user averages too
		cerr
				<< "Removing the user averages (given movie averages) found from the working copy\n";
		for (user u = 0; u < U; u++) {
			for (int index = 0; index < noViewingsPerUser[u]; index++) {
				workingCopyOfRatingsByUser[u][index] -= userAverage[u];
			}
		}
	}
	delete[] userAverage;
}

void PreProc::movieSideAverages() {
	// First find total ratings for each movie and number of ratings
	double* movieAverage = new double[M];
	int* movieCount = new int[M];

	for (movie m = 0; m < M; m++) {
		movieAverage[m] = 0.0;
		movieCount[m] = 0;
	}

	for (user u = 0; u < U; u++) {
		for (int index = 0; index < noViewingsPerUser[u]; index++) {
			movieAverage[moviesByUser[u][index]] +=
					workingCopyOfRatingsByUser[u][index];
			movieCount[moviesByUser[u][index]]++;
		}
	}

	// Next find the average for each movie, and compute the statistics of the movie averages
	double total = 0.0;
	double totalSquared = 0.0;

	for (movie m = 0; m < M; m++) {
		movieAverage[m] /= (double) (movieCount[m]);
		total += movieAverage[m];
		totalSquared += movieAverage[m] * movieAverage[m];
	}

	double mean = total / (double) (M);
	double biasedVar = (totalSquared / (double) (M)) - mean * mean;
	double unbiasedVar = biasedVar / (1.0 - 1.0 / (double) (M));

	cerr << "Movie Add: " << mean << "," << unbiasedVar << "\n";
	mMovieAddMean = mean;
	mMovieAddVar = unbiasedVar;

	// Subtract the movie averages from the working copy of the ratings.
	cerr << "Taking movie adds off the working copy\n";

	for (user u = 0; u < U; u++) {
		for (int index = 0; index < noViewingsPerUser[u]; index++) {
			workingCopyOfRatingsByUser[u][index] -=
					movieAverage[moviesByUser[u][index]];
		}
	}

	// Next do the user side averages based on statistics with the movie averages removed.
	cerr << "The following are for user averages GIVEN the movie average\n";
	userSideAverages(true);

	delete[] movieAverage;
	delete[] movieCount;
}

void PreProc::wAndDpFit() {
	// First, find the statistics of what remains of the working copies
	double total = 0.0;
	double totalSquared = 0.0;
	int obs = 0;

	for (user u = 0; u < U; u++) {
		for (int index = 0; index < noViewingsPerUser[u]; index++) {
			total += workingCopyOfRatingsByUser[u][index];
			totalSquared += workingCopyOfRatingsByUser[u][index]
					* workingCopyOfRatingsByUser[u][index];
			obs++;
		}
	}

	double mean = total / (double) (obs);
	double biasedVar = (totalSquared / (double) (obs)) - mean * mean;
	double unbiasedVar = biasedVar / (1.0 - 1.0 / (double) (obs));

	cerr
			<< "Simple calculation gives SigmaDP^2 * (SigmaDP^2 + SigmaW ^ 2 / root(N) ) = "
			<< unbiasedVar << ". Incidentally, mean was " << mean << "\n";

	// Now compute estimates of Var(DP) and Var(W), based on fitting a line to the
	// variances of users ratings (y) and 1 over the square root of viewings made (x)
	double sigmaX = 0.0;
	double sigmaY = 0.0;
	double sigmaXY = 0.0;
	double sigmaXX = 0.0;
	double N = (double) (U);

	double* y = new double[U];
	double* x = new double[U];
	double tot = 0.0;

	for (user u = 0; u < U; u++) {
		y[u] = 0.0;
		x[u] = 1.0 / pow((double) (noViewingsPerUser[u]), 0.5);
		tot = 0.0;
		for (int index = 0; index < noViewingsPerUser[u]; index++) {
			y[u] += workingCopyOfRatingsByUser[u][index]
					* workingCopyOfRatingsByUser[u][index];
			tot += workingCopyOfRatingsByUser[u][index];
		}

		tot /= (double) (noViewingsPerUser[u]);
		y[u] /= (double) (noViewingsPerUser[u]);
		y[u] -= tot * tot;

		sigmaX += x[u];
		sigmaY += y[u];
		sigmaXY += x[u] * y[u];
		sigmaXX += x[u] * x[u];
	}

	// With these, we can find Maximum Likelihood values of sigmaDP and sigmaW
	double beta = (sigmaXY - sigmaX * sigmaY / N)
			/ (sigmaXX - sigmaX * sigmaX / N);
	double alpha = (sigmaY - beta * sigmaX) / N;

	cerr << "Alpha and beta found to be :" << alpha << " and " << beta << "\n";
	cerr << "As in Variances are modelled by alpha + beta / NumberOfViewings\n";
	cerr << "alpha = SigmaDP^2, beta = SigmaDP * SigmaW\n";

	double sigmaDP = pow(alpha, 0.5);
	double sigmaW = beta / sigmaDP;

	cerr << "SigmaDP = " << sigmaDP << ", sigmaW = " << sigmaW << ".\n\tBye!\n";

	delete[] x;
	delete[] y;
}

void PreProc::sortUserStuffByDate(movie* &thisMovies, rate* &thisRatings,
		uint* &thisDates, uint beginIndex, uint endIndex) {
	uint Piv_index;
	if (beginIndex < endIndex) {
		Piv_index = partitionUserStuffByDate(thisMovies, thisRatings, thisDates,
				beginIndex, endIndex);
		sortUserStuffByDate(thisMovies, thisRatings, thisDates, beginIndex,
				Piv_index - 1);
		sortUserStuffByDate(thisMovies, thisRatings, thisDates, Piv_index + 1,
				endIndex);
	}
}

uint PreProc::partitionUserStuffByDate(movie* &thisMovies, rate* &thisRatings,
		uint* &thisDates, uint beginIndex, uint endIndex) {
	uint pivotDate = thisDates[beginIndex];

	swap(thisMovies, thisRatings, thisDates, beginIndex, endIndex);

	uint storeIndex = beginIndex;

	for (uint i = beginIndex; i < endIndex; i++) {
		if (thisDates[i] <= pivotDate) {
			swap(thisMovies, thisRatings, thisDates, i, storeIndex);
			storeIndex++;
		}
	}

	swap(thisMovies, thisRatings, thisDates, storeIndex, endIndex);

	return storeIndex;
}

void PreProc::swap(movie* &thisMovies, rate* &thisRatings, uint* &thisDates,
		uint index1, uint index2) {
	movie tempM = thisMovies[index1];
	rate tempR = thisRatings[index1];
	uint tempD = thisDates[index1];

	thisMovies[index1] = thisMovies[index2];
	thisRatings[index1] = thisRatings[index2];
	thisDates[index1] = thisDates[index2];

	thisMovies[index2] = tempM;
	thisRatings[index2] = tempR;
	thisDates[index2] = tempD;
}

void PreProc::sortMovieStuffByDate(user* &thisUsers, rate* &thisRatings,
		uint* &thisDates, uint beginIndex, uint endIndex) {
	uint Piv_index;
	if (beginIndex < endIndex) {
		Piv_index = partitionMovieStuffByDate(thisUsers, thisRatings, thisDates,
				beginIndex, endIndex);
		sortMovieStuffByDate(thisUsers, thisRatings, thisDates, beginIndex,
				Piv_index - 1);
		sortMovieStuffByDate(thisUsers, thisRatings, thisDates, Piv_index + 1,
				endIndex);
	}
}

uint PreProc::partitionMovieStuffByDate(user* &thisUsers, rate* &thisRatings,
		uint* &thisDates, uint beginIndex, uint endIndex) {
	uint pivotDate = thisDates[beginIndex];
	swap(thisUsers, thisRatings, thisDates, beginIndex, endIndex);

	uint storeIndex = beginIndex;

	for (uint i = beginIndex; i < endIndex; i++) {
		if (thisDates[i] <= pivotDate) {
			swap(thisUsers, thisRatings, thisDates, i, storeIndex);
			storeIndex++;
		}
	}

	swap(thisUsers, thisRatings, thisDates, storeIndex, endIndex);

	return storeIndex;
}

void PreProc::swap(user* &thisUsers, rate* &thisRatings, uint* &thisDates,
		uint index1, uint index2) {
	user tempU = thisUsers[index1];
	rate tempR = thisRatings[index1];
	uint tempD = thisDates[index1];

	thisUsers[index1] = thisUsers[index2];
	thisRatings[index1] = thisRatings[index2];
	thisDates[index1] = thisDates[index2];

	thisUsers[index2] = tempU;
	thisRatings[index2] = tempR;
	thisDates[index2] = tempD;
}

void PreProc::checkUserDatesAreInOrder() {
	for (user u = 0; u < U; u++) {
		for (int i = 1; i < noViewingsPerUser[u]; i++) {
			if (datesByUser[u][i] < datesByUser[u][i - 1]
					|| datesByUser[u][i - 1] <= 0) {
				cerr << "Dates by user are NOT ordered correctly for user " << u
						<< "\n";
			}
		}
	}
}

void PreProc::checkMovieDatesAreInOrder() {
	for (movie m = 0; m < M; m++) {
		for (uint i = 1; i < noViewingsPerMovie[m]; i++) {
			if (datesByMovie[m][i] < datesByMovie[m][i - 1]
					|| datesByMovie[m][i - 1] <= 0) {
				cerr << "Dates by movie are NOT ordered correctly for movie "
						<< m << "\n";
			}
		}
	}
}

uint PreProc::getUserRatingHistoryFromDate(user u, uint targetDate,
		uint numOfDaysHistory) {
	uint total = 0;
	uint ratings = 0;

	for (uint index = 0; index < noViewingsPerUser[u]; index++) {
		if (datesByUser[u][index] > targetDate) {
			break;
		}

		if (datesByUser[u][index] + numOfDaysHistory > targetDate) {
			ratings++;
			total += ratingsByUser[u][index];
		}
	}

	if (total == 0) {
		return 0;
	}

	double exact = ((double) total) / ((double) ratings);

	// Rounded down if anything.
	uint aboutRight = ((exact - 0.98) / 0.04) / 1;

	double approxAboutRight = ((double) aboutRight) * 0.04 + 0.98;
	double approxAboutRightPlusOne = ((double) (aboutRight + 1)) * 0.04 + 0.98;

	if (fabs(exact - approxAboutRight)
			> fabs(exact - approxAboutRightPlusOne)) {
		aboutRight++;
	}

	if (fabs(((double) aboutRight) * 0.04 + 0.98 - exact) > 0.0205) {
		cerr << "The ratings history value was no good " << exact << " "
				<< ratings << " " << total << " " << aboutRight << "\n";
	}

	return min(max(aboutRight, (uint) 1), (uint) 100);
}

void PreProc::writeBigRatingFile() {
	cerr << "Big rating file started\n";

	ofstream f;
	f.open(HcTrainingFile);

	for (user u = 0; u < U; u++) {
		f << "!" << u << "!" << noViewingsPerUser[u] << "\n";
		bigRatings[u] = new q[noViewingsPerUser[u]];

		for (int i = 0; i < noViewingsPerUser[u]; i++) {
			if (i > 0) {
				f << ",";
			}

			bigRatings[u][i] = makeBig(u, moviesByUser[u][i], datesByUser[u][i],
					ratingsByUser[u][i]);
			f << bigRatings[u][i];
		}
		f << "\n";
	}

	f.close();

	cerr << "Big rating file finished\n";
}

void PreProc::writeBigPredictionFile() {
	cerr << "Big prediction file started\n";

	vector<user> users;
	vector<movie> movies;
	vector < uint > dates;

	data.getAllPredictionData(users, movies, dates);

	numPredictions = users.size();
	bigPredictionUsers = new user[numPredictions];
	bigPrediction = new q[numPredictions];

	ofstream f;
	f.open(HcPredictionFile);

	for (int i = 0; i < numPredictions; i++) {
		bigPredictionUsers[i] = users[i];
		bigPrediction[i] = makeBig(users[i], movies[i], dates[i], 0);
		f << bigPredictionUsers[i] << "," << bigPrediction[i] << "\n";
	}

	f.close();

	cerr << "Big prediction file finished after writing " << numPredictions
			<< " predictions\n";
}

void PreProc::testBigPredictionFile() {
	cerr << "Big prediction file test started\n";

	vector<user> users;
	vector<q> preds;

	data.getAllBigPredictionData(users, preds);

	if (users.size() != preds.size()
			|| users.size() != (unsigned int) numPredictions) {
		cerr << "Num predictions changed - " << users.size() << " "
				<< preds.size() << " " << numPredictions << "\n";
	}

	for (int i = 0; i < numPredictions; i++) {
		if (bigPredictionUsers[i] != users[i]) {
			cerr << "User changed for prediction " << i << " "
					<< bigPredictionUsers[i] << " " << users[i] << "\n";
		}

		if (bigPrediction[i] != preds[i]) {
			cerr << "Prediction value changed for prediction " << i << " "
					<< bigPrediction[i] << " " << preds[i] << "\n";
		}
	}

	cerr
			<< "Big prediction file test finished - if you didn't hear otherwise, everything was hunky-dory\n";
}

q PreProc::makeBig(user u, movie m, uint date, rate rating) {
	q bigRating = data.GetBigRating(rating, m, getUserEra(u, date, 5),
			getMovieEra(m, date, 50), date % 7, 0,
			getUserRatingHistoryFromDate(u, date, 1),
			getUserRatingHistoryFromDate(u, date, 3),
			getUserRatingHistoryFromDate(u, date, 5),
			getUserRatingHistoryFromDate(u, date, 10),
			getUserRatingHistoryFromDate(u, date, 30));

	// Check it's worked here.
	// Rating
	if (bigRating % RATE != rating) {
		cerr << "Rating wrong: " << (bigRating % RATE) << " " << rating << "\n";
	}

	// Movie
	if (bigRating % MOV != m) {
		cerr << "Movie wrong: " << (bigRating % MOV) << " " << m << "\n";
	}

	// User era
	if (bigRating % ERAU != getUserEra(u, date, 5)) {
		cerr << "User era wrong: " << (bigRating % ERAU) << " "
				<< getUserEra(u, date, 5) << "\n";
	}

	// Movie era
	if (bigRating % ERAM != getMovieEra(m, date, 50)) {
		cerr << "Movie era wrong: " << (bigRating % ERAM) << " "
				<< getMovieEra(m, date, 50) << "\n";
	}

	// Weekday
	if (bigRating % WKDY != (date % 7)) {
		cerr << "Weekday was wrong: " << (bigRating % WKDY) << " " << date
				<< "\n";
	}

	// Precision
	if (bigRating % PREC != 0) {
		cerr << "Precision was wrong: " << (bigRating % PREC) << "\n";
	}

	// History level 1
	if (bigRating % HIST1 != getUserRatingHistoryFromDate(u, date, 1)) {
		cerr << "History level 1 wrong: " << (bigRating % HIST1) << " "
				<< getUserRatingHistoryFromDate(u, date, 1) << "\n";
	}

	// History level 3
	if (bigRating % HIST3 != getUserRatingHistoryFromDate(u, date, 3)) {
		cerr << "History level 3 wrong: " << (bigRating % HIST3) << " "
				<< getUserRatingHistoryFromDate(u, date, 3) << "\n";
	}

	// History level 5
	if (bigRating % HIST5 != getUserRatingHistoryFromDate(u, date, 5)) {
		cerr << "History level 5 wrong: " << (bigRating % HIST5) << " "
				<< getUserRatingHistoryFromDate(u, date, 5) << "\n";
	}

	// History level 10
	if (bigRating % HIST10 != getUserRatingHistoryFromDate(u, date, 10)) {
		cerr << "History level 10 wrong: " << (bigRating % HIST10) << " "
				<< getUserRatingHistoryFromDate(u, date, 10) << "\n";
	}

	// History level 30
	if (bigRating % HIST30 != getUserRatingHistoryFromDate(u, date, 30)) {
		cerr << "History level 30 wrong: " << (bigRating % HIST30) << " "
				<< getUserRatingHistoryFromDate(u, date, 30) << "\n";
	}

	return bigRating;
}

unsigned int PreProc::getMovieEra(movie m, unsigned int ratingDate,
		unsigned int upperLimit) {
	// These cases are possible for predictions
	if (ratingDate <= datesByMovie[m][0]) {
		return 0;
	}

	if (ratingDate >= datesByMovie[m][noViewingsPerMovie[m] - 1]) {
		return upperLimit - 1;
	}

	double approx = (((double) (ratingDate - datesByMovie[m][0]))
			/ ((double) (datesByMovie[m][noViewingsPerMovie[m] - 1]
					- datesByMovie[m][0]))) * ((double) (upperLimit));
	unsigned int approxReturn = approx / 1;

	if (approxReturn < 0 || approxReturn > upperLimit) {
		cerr << "getMovieEra had bad approx return value of " << approxReturn
				<< " for movie " << m
				<< " (probably OK if it's in the prediction file)\n";
	}

	return min(approxReturn, upperLimit - 1);
}

unsigned int PreProc::getUserEra(user u, unsigned int ratingDate,
		unsigned int upperLimit) {
	if (ratingDate <= datesByUser[u][0]) {
		return 0;
	}

	if (ratingDate >= datesByUser[u][noViewingsPerUser[u] - 1]) {
		return upperLimit - 1;
	}

	double approx = (((double) (ratingDate - datesByUser[u][0]))
			/ ((double) (datesByUser[u][noViewingsPerUser[u] - 1]
					- datesByUser[u][0]))) * ((double) (upperLimit));
	unsigned int approxReturn = approx / 1;

	if (approxReturn < 0 || approxReturn > upperLimit) {
		cerr << "getUserEra had bad approx return value of " << approxReturn
				<< " for user " << u
				<< "(probably OK if it's in the prediction file)\n";
	}

	return min(approxReturn, upperLimit - 1);
}

}
