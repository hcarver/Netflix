#include "NetflixDataNG.h"

#include "MD5.h"

namespace utils {

// Constructor and Destructor
NetflixDataNG::NetflixDataNG() {
	debug = false;
	test = false;

	// Populate path variables
	basePathToUse = beastBasePath;
	dataPath = basePathToUse + pathToData;
	qualDataFile = basePathToUse + pathToQualifyingDataFile;

	// Populate the User id converter
	populateUserIdConverter(basePathToUse + pathToConversionFile);

	// Now need to populate the array of calendar information
	// First do 1998 itself. Populate with 0s up until November
	for (int month = 1; month <= 10; month++) {
		daysAfter30Sept1998[0][month - 1] = 0;
	}
	// November (October has 31)
	daysAfter30Sept1998[0][10] = 31;
	// December (November has 30)
	daysAfter30Sept1998[0][11] = 61;

	// Now iterate over years and populate forward depending on number of days in PREVIOUS month
	for (int year = 1999; year <= 2005; year++) {
		// Jan (Dec has 31, from previous year)
		daysAfter30Sept1998[year - 1998][0] =
				daysAfter30Sept1998[year - 1999][11] + 31;
		// Feb (Jan has 31)
		daysAfter30Sept1998[year - 1998][1] =
				daysAfter30Sept1998[year - 1998][0] + 31;
		// March (February has 28, or 29 in 2000 and 2004)
		if ((year % 4) == 0) {
			daysAfter30Sept1998[year - 1998][2] = daysAfter30Sept1998[year
					- 1998][1] + 29;
		} else {
			daysAfter30Sept1998[year - 1998][2] = daysAfter30Sept1998[year
					- 1998][1] + 28;
		}
		// April (March has 31)
		daysAfter30Sept1998[year - 1998][3] =
				daysAfter30Sept1998[year - 1998][2] + 31;
		// May (April has 30)
		daysAfter30Sept1998[year - 1998][4] =
				daysAfter30Sept1998[year - 1998][3] + 30;
		// June (May has 31)
		daysAfter30Sept1998[year - 1998][5] =
				daysAfter30Sept1998[year - 1998][4] + 31;
		// July (June has 30)
		daysAfter30Sept1998[year - 1998][6] =
				daysAfter30Sept1998[year - 1998][5] + 30;
		// August (July has 31)
		daysAfter30Sept1998[year - 1998][7] =
				daysAfter30Sept1998[year - 1998][6] + 31;
		// September (August has 31)
		daysAfter30Sept1998[year - 1998][8] =
				daysAfter30Sept1998[year - 1998][7] + 31;
		// October (September has 30)
		daysAfter30Sept1998[year - 1998][9] =
				daysAfter30Sept1998[year - 1998][8] + 30;
		// November (October has 31)
		daysAfter30Sept1998[year - 1998][10] =
				daysAfter30Sept1998[year - 1998][9] + 31;
		// December (November has 30)
		daysAfter30Sept1998[year - 1998][11] =
				daysAfter30Sept1998[year - 1998][10] + 30;
	}
}

NetflixDataNG::~NetflixDataNG() {
}

// File management public helper functions
void NetflixDataNG::deleteFileAsynch(string fileName) {
	stringstream s("");
	s << "rm " << fileName << " > /dev/null 2>&1 &"; // This extra bit silences it if the delete fails (after all, who cares?)
	string command = s.str();
	if (system(&command[0]) != 0) {
		// Don't care if file isn't deleted
	}
}

void NetflixDataNG::zipFileAsynch(string fileName, string newFileName) {
	stringstream s("");
	s << "gzip -f -q " << fileName << " " << newFileName << " &";
	string command = s.str();
	if (system(&command[0]) != 0) {
		cerr << "Couldn't zip file " << fileName << "\n";
	}
}

void NetflixDataNG::zipFileSynch(string fileName, string newFileName) {
	stringstream s("");
	s << "gzip -f -q " << fileName << " " << newFileName;
	string command = s.str();
	if (system(&command[0]) != 0) {
		cerr << "Couldn't synchronously zip file " << fileName << "\n";
	}
}

void NetflixDataNG::unzipFileSynch(string fileName, string newFileName) {
	stringstream s("");
	s << "gunzip " << fileName << " " << newFileName;
	string command = s.str();
	if (system(&command[0]) != 0) {
		cerr << "Couldn't synchronously unzip file " << fileName << "\n";
	}
}

void NetflixDataNG::makeEntryAsynch(string submissionFileName) {
	cerr << "Submitting latest prediction file\n";

	// Zip the file
	string zippedFileName = submissionFileName + ".gz";
	zipFileSynch(submissionFileName, zippedFileName);

	// Calculate the md5sum
	MD5 myWrapper;
	string hash = myWrapper.getHashFromFile(zippedFileName);

	// Command:
	// -s for silent, i.e. no progress or error messages.
	// -F for each input
	// URL to submit to
	// Send output html to null to avoid displaying.
	stringstream s("");
	s << "curl " << "-s " << "-F 'teamName=TeamOz' "
			<< "-F 'teamPassword=Selmer3' " << "-F 'MD5Hash=" << hash << "' "
			<< "-F 'file=@" << zippedFileName << "' "
			<< "http://www.netflixprize.com/submissions " << ">> null &";

	string command = s.str();

	if (system(&command[0]) != 0) {
		cerr << "Couldn't asynchronously submit file " << zippedFileName
				<< ", hash " << hash << "\n";
	}
}

// Data read functions
void NetflixDataNG::getDataMeanAndVariance(double mAndV[]) {
	try {
		// Open VitalStatistics file
		stringstream s("");
		s << basePathToUse;
		s << pathToVitalStatistics;
		string path = s.str();
		individualReader.open(&path[0]);

		// Read line which is "mean, variance"
		string line(100, '\0');
		individualReader.getline(&line[0], 99);
		char * pEnd;
		mAndV[0] = strtod(&line[0], &pEnd);
		pEnd = &pEnd[1];
		mAndV[1] = strtod(pEnd, NULL);
	} catch (...) {
		cerr << "Could not get the Data mean and variance";
	}
	individualReader.close();
}

//void NetflixDataNG::allViewingsUsers(
//  movie**                 &moviesByUser, 
//  rate**                    &ratesByUser,
//  unsigned short* &noViewingsPerUser)
//{
//    if(!test)
//    {
//        // Define variables for reading from files
//        string path;
//        ifstream reader;
//        string line(100, '\0');
//        char*pEnd;
//        
//        // Define variables for storing read ratings, movies and recording number of viewings etc.
//        movie mov;
//        rate rat;
//        vector<rate> rates;
//        vector<movie> moviesVec;
//        int views;
//        
//        // Parallelise making user-specific stuff private to each thread.
//#pragma omp parallel for private(pEnd, mov, rat, path, views, reader) firstprivate(rates, moviesVec, line) schedule(guided,5)
//
//        for (user u=0; u<U; u++)
//        {
//            if(u%48000==0)
//                cout<<u<<"\n";
//
//            rates.clear();
//            moviesVec.clear();
//            
//            // Open file and get to ratings data (past header line)
//            path = getTrainingDataPath(u, true);
//            reader.open(&path[0], ios_base::in);
//            reader.getline(&line[0], 99);
//            
//            while (true)
//            {
//                // Read lines until EOF
//                reader.getline(&line[0], 99);
//                if (line[0] == '\0')
//                    break;
//                    
//                // Obtain movie and rating
//                getQuartetFromUserData(mov, rat, line, pEnd);
//                
//                // Add these to the vector for this user.
//                moviesVec.push_back(mov);
//                rates.push_back(rat);
//            }
//            
//            // Close reader (ASAP)
//            reader.close();
//            
//            // Update data arrays with the number of viewings by this user
//            views = moviesVec.size();
//            noViewingsPerUser[u] = (short) views;
//            moviesByUser[u] = new movie[views];
//            ratesByUser[u] = new rate[views];
//            
//            // Update data arrays with the movies and ratings
//            for (int index=0; index<views; index++)
//            {
//                moviesByUser[u][index] = moviesVec[index];
//                ratesByUser[u][index] = rates[index];
//            }
//        }
//#pragma omp barrier
//    }
//    else
//    {
//        // If in test mode we populate the arrays with 3 random movies and ratings 
//        srand(time(0));
//        rand();
//        int views=3;
//        
//        // Iterate over users
//        for (user u=0; u<U; u++)
//        {
//            noViewingsPerUser[u] = views;
//            moviesByUser[u] = new movie[views];
//            ratesByUser[u] = new rate[views];
//
//            // First rating simply random
//            moviesByUser[u][0] = (rand() % M);
//            ratesByUser[u][0] = 1 + (rand() % 5);
//            
//            // Second rating only OK if movie isn't the same
//            do
//            {
//                moviesByUser[u][1] = (rand() % M);
//            }
//            while(moviesByUser[u][1] == moviesByUser[u][0]);
//            ratesByUser[u][1] = 1 + (rand() % 5);
//                        
//            // Third rating only OK if movie isn't the same as first 2
//            do
//            {
//                moviesByUser[u][2] = (rand() % M);
//            }
//            while(moviesByUser[u][2]==moviesByUser[u][1] || moviesByUser[u][2]==moviesByUser[u][0]);
//            ratesByUser[u][2] = 1 + (rand() % 5);
//        }
//    }
//}

void NetflixDataNG::allViewingsMovies(user** &usersByMovie,
		unsigned long* &noViewingsPerMovie, rate** &ratesByMovie,
		rate** &ratesByUser, movie** &moviesByUser,
		unsigned short* &noViewingsPerUser,
		unsigned short ** &whatIndexAMovieIsToAUser) {
	if (!test) {
		// Variables for each rating
		user use;
		rate rat;

		// Variables for each movie
		ifstream reader;
		string path;
		long views;
		char*pEnd;

#pragma omp parallel for private(use,rat,path,views, pEnd) private(reader) schedule(guided,5)
		for (movie m = 0; m < M; m++) {
			pEnd = 0;
			// Make vectors to store the ratings and user data
			vector<rate> rates;
			vector<user> usersVec;

			// IO stuff
			string line(100, '\0');
			path = getTrainingDataPath(m, false);
			reader.open(&path[0], ios_base::in);
			reader.getline(&line[0], 99); // Line of rubbish

			// Until EOF, read in line and get user and rating from it.
			while (true) {
				reader.getline(&line[0], 99);
				if (line[0] == '\0') {
					break;
				}
				getQuartetFromMovieData(use, rat, line, pEnd);
				usersVec.push_back(use);
				rates.push_back(rat);
			}
			// Close reader ASAP
			reader.close();

			// Use number of viewings to make arrays and populate noViewingsPerMovie
			views = (long) usersVec.size();

			noViewingsPerMovie[m] = views;
			usersByMovie[m] = new user[views];
			ratesByMovie[m] = new rate[views];
			whatIndexAMovieIsToAUser[m] = new unsigned short[views];

			// Populate usersByMovie & ratesByMovie
			for (int index = 0; index < views; index++) {
				usersByMovie[m][index] = usersVec[index];
				ratesByMovie[m][index] = rates[index];
			}
		}
#pragma omp barrier

		// This is new stuff for getting user-perspective ratings from movie-perspective ones.
		// Make vector arrays for storing the ratings user-side.
#define USERSATONCE 1000
		vector<rate> urates[USERSATONCE];
		vector<movie> umovies[USERSATONCE];

#pragma omp parallel for firstprivate(urates, umovies)
		for (uint i = 0; i <= U / USERSATONCE; i++) {
			// Limit to number of users dealt with by one iteration is either USERSATONCE or
			// on last iteration 189 ( = U - i*USERSATONCE)
			uint jLim = min((uint) USERSATONCE, U - i * USERSATONCE);

			// Make sure vectors are empty
			for (uint jj = 0; jj < jLim; jj++) {
				urates[jj].clear();
				umovies[jj].clear();
			}

			// Go through each movie and index & add rating and movie to the vectors
			for (movie m = 0; m < M; m++) {
				for (uint index = 0; index < noViewingsPerMovie[m]; index++) {
					user u = usersByMovie[m][index];
					if (u / USERSATONCE == i) {
						urates[u % USERSATONCE].push_back(
								ratesByMovie[m][index]);
						umovies[u % USERSATONCE].push_back(m);
						whatIndexAMovieIsToAUser[m][index] = urates[u
								% USERSATONCE].size() - 1;
					}
				}
			}

			// Go through each user in this range
			for (uint jj = 0; jj < jLim; jj++) {
				// Get user number & number of ratings
				uint views = urates[jj].size();
				user u = i * USERSATONCE + jj;

				// Fill in user info
				noViewingsPerUser[u] = views;
				moviesByUser[u] = new movie[views];
				ratesByUser[u] = new rate[views];
				for (uint index = 0; index < views; index++) {
					moviesByUser[u][index] = umovies[jj][index];
					ratesByUser[u][index] = urates[jj][index];
				}
			}
		}
#pragma omp barrier

	} else {
		// If in test mode we populate the user arrays with 3 random movies and ratings
		srand(time(0));
		rand();
		int views = 3;

		// Iterate over users
		for (user u = 0; u < U; u++) {
			noViewingsPerUser[u] = views;
			moviesByUser[u] = new movie[views];
			ratesByUser[u] = new rate[views];

			// First rating simply random
			moviesByUser[u][0] = (rand() % M);
			ratesByUser[u][0] = 1 + (rand() % 5);

			// Second rating only OK if movie isn't the same
			do {
				moviesByUser[u][1] = (rand() % M);
			} while (moviesByUser[u][1] == moviesByUser[u][0]);

			ratesByUser[u][1] = 1 + (rand() % 5);

			// Third rating only OK if movie isn't the same as first 2
			do {
				moviesByUser[u][2] = (rand() % M);
			} while (moviesByUser[u][2] == moviesByUser[u][1]
					|| moviesByUser[u][2] == moviesByUser[u][0]);

			ratesByUser[u][2] = 1 + (rand() % 5);
		}

		// Vectors for each movie (rate, user, index)
		vector<rate> rates[M];
		vector<user> usersVec[M];
		vector<rate> indexes[M];
		for (user u = 0; u < U; u++) {
			for (int index = 0; index < noViewingsPerUser[u]; index++) {
				rates[moviesByUser[u][index]].push_back(ratesByUser[u][index]);
				usersVec[moviesByUser[u][index]].push_back(u);
				indexes[moviesByUser[u][index]].push_back(index);
			}
		}

		int length;
		for (movie m = 0; m < M; m++) {
			// Use length to make arrays
			length = usersVec[m].size();
			noViewingsPerMovie[m] = length;
			usersByMovie[m] = new user[length];
			ratesByMovie[m] = new rate[length];
			whatIndexAMovieIsToAUser[m] = new short unsigned int[length];

			// Populate arrays with user, rate, index
			for (int i = 0; i < length; i++) {
				usersByMovie[m][i] = usersVec[m][i];
				ratesByMovie[m][i] = rates[m][i];
				whatIndexAMovieIsToAUser[m][i] = indexes[m][i];
			}
		}
	}
}

void NetflixDataNG::allViewingsMovies(user** &usersByMovie,
		ulong* &noViewingsPerMovie, rate** &ratesByMovie, uint** &datesByMovie,
		rate** &ratesByUser, movie** &moviesByUser, ushort* &noViewingsPerUser,
		uint** &datesByUser, unsigned short ** &whatIndexAMovieIsToAUser) {
	if (!test) {
		// Variables for each rating
		user use;
		rate rat;
		uint date;

		// Variables for each movie
		ifstream reader;
		string path;
		long views;
		char*pEnd;

#pragma omp parallel for private(use,rat, date, path,views, pEnd) private(reader) schedule(guided,5)
		for (movie m = 0; m < M; m++) {
			pEnd = 0;
			// Make vectors to store the ratings and user data
			vector<rate> rates;
			vector<user> usersVec;
			vector < uint > dateVec;

			// IO stuff
			string line(100, '\0');
			path = getTrainingDataPath(m, false);
			reader.open(&path[0], ios_base::in);
			reader.getline(&line[0], 99); // Line of rubbish

			// Until EOF, read in line and get user and rating from it.
			while (true) {
				reader.getline(&line[0], 99);
				if (line[0] == '\0') {
					break;
				}
				getInfoFromMovieData(use, rat, date, line, pEnd);
				usersVec.push_back(use);
				rates.push_back(rat);
				dateVec.push_back(date);
			}
			// Close reader ASAP
			reader.close();

			// Use number of viewings to make arrays and populate noViewingsPerMovie
			views = (long) usersVec.size();

			noViewingsPerMovie[m] = views;
			usersByMovie[m] = new user[views];
			ratesByMovie[m] = new rate[views];
			datesByMovie[m] = new uint[views];
			whatIndexAMovieIsToAUser[m] = new ushort[views];

			// Populate usersByMovie & ratesByMovie
			for (int index = 0; index < views; index++) {
				usersByMovie[m][index] = usersVec[index];
				ratesByMovie[m][index] = rates[index];
				datesByMovie[m][index] = dateVec[index];
			}
		}
#pragma omp barrier

		// This is new stuff for getting user-perspective ratings from movie-perspective ones.
		// Make vector arrays for storing the ratings user-side.
#define USERSATONCE 1000
		vector<rate> urates[USERSATONCE];
		vector<movie> umovies[USERSATONCE];
		vector<int> udates[USERSATONCE];

#pragma omp parallel for firstprivate(urates, umovies, udates)
		for (uint i = 0; i <= U / USERSATONCE; i++) {
			// Limit to number of users dealt with by one iteration is either USERSATONCE or
			// on last iteration 189 ( = U - i*USERSATONCE)
			uint jLim = min((uint) USERSATONCE, U - i * USERSATONCE);

			// Make sure vectors are empty
			for (uint jj = 0; jj < jLim; jj++) {
				urates[jj].clear();
				umovies[jj].clear();
				udates[jj].clear();
			}

			// Go through each movie and index & add rating and movie to the vectors
			for (movie m = 0; m < M; m++) {
				for (uint index = 0; index < noViewingsPerMovie[m]; index++) {
					user u = usersByMovie[m][index];
					if (u / USERSATONCE == i) {
						urates[u % USERSATONCE].push_back(
								ratesByMovie[m][index]);
						umovies[u % USERSATONCE].push_back(m);
						udates[u % USERSATONCE].push_back(
								datesByMovie[m][index]);
						whatIndexAMovieIsToAUser[m][index] = urates[u
								% USERSATONCE].size() - 1;
					}
				}
			}

			// Go through each user in this range
			for (uint jj = 0; jj < jLim; jj++) {
				// Get user number & number of ratings
				uint views = urates[jj].size();
				user u = i * USERSATONCE + jj;

				// Fill in user info
				noViewingsPerUser[u] = views;
				moviesByUser[u] = new movie[views];
				ratesByUser[u] = new rate[views];
				datesByUser[u] = new uint[views];
				for (uint index = 0; index < views; index++) {
					moviesByUser[u][index] = umovies[jj][index];
					ratesByUser[u][index] = urates[jj][index];
					datesByUser[u][index] = udates[jj][index];
				}
			}
		}
#pragma omp barrier

	} else {
		// If in test mode we populate the user arrays with 3 random movies and ratings
		srand(time(0));
		rand();
		int views = 3;

		// Iterate over users
		for (user u = 0; u < U; u++) {
			noViewingsPerUser[u] = views;
			moviesByUser[u] = new movie[views];
			ratesByUser[u] = new rate[views];
			datesByUser[u] = new uint[views];

			// First rating simply random
			moviesByUser[u][0] = (rand() % M);
			ratesByUser[u][0] = 1 + (rand() % 5);
			datesByUser[u][0] = 1 + (rand() % 2600);

			// Second rating only OK if movie isn't the same
			do {
				moviesByUser[u][1] = (rand() % M);
			} while (moviesByUser[u][1] == moviesByUser[u][0]);

			ratesByUser[u][1] = 1 + (rand() % 5);
			datesByUser[u][1] = 1 + (rand() % 2600);

			// Third rating only OK if movie isn't the same as first 2
			do {
				moviesByUser[u][2] = (rand() % M);
			} while (moviesByUser[u][2] == moviesByUser[u][1]
					|| moviesByUser[u][2] == moviesByUser[u][0]);

			ratesByUser[u][2] = 1 + (rand() % 5);
			datesByUser[u][2] = 1 + (rand() % 2600);
		}

		// Vectors for each movie (rate, user, index)
		vector<rate> rates[M];
		vector<user> usersVec[M];
		vector<rate> indexes[M];
		vector<int> dates[M];

		for (user u = 0; u < U; u++) {
			for (uint index = 0; index < noViewingsPerUser[u]; index++) {
				dates[moviesByUser[u][index]].push_back(datesByUser[u][index]);
				rates[moviesByUser[u][index]].push_back(ratesByUser[u][index]);
				usersVec[moviesByUser[u][index]].push_back(u);
				indexes[moviesByUser[u][index]].push_back(index);
			}
		}
		int length;
		for (movie m = 0; m < M; m++) {
			// Use length to make arrays
			length = usersVec[m].size();
			noViewingsPerMovie[m] = length;
			usersByMovie[m] = new user[length];
			ratesByMovie[m] = new rate[length];
			datesByMovie[m] = new uint[length];
			whatIndexAMovieIsToAUser[m] = new ushort[length];

			// Populate arrays with user, rate, index
			for (int i = 0; i < length; i++) {
				usersByMovie[m][i] = usersVec[m][i];
				ratesByMovie[m][i] = rates[m][i];
				datesByMovie[m][i] = dates[m][i];
				whatIndexAMovieIsToAUser[m][i] = indexes[m][i];
			}
		}
	}
}

void NetflixDataNG::readBigRatings(unsigned short* &noViewingsPerUser,
		q** &bigRatingsByUser) {
	if (!test) {
		ifstream f;
		f.open(HcTrainingFile);

		string readLine(400000, '\0');

		for (user u = 0; u < U; u++) {
			f.getline(&readLine[0], 399999); //    f << "!" << u << "!" << noViewingsPerUser[u] <<"\n";
			char* pos = &readLine[0];

			pos = &pos[1];
			int views = strtoul(pos, &pos, 10); // Actually this is the user number
			pos = &pos[1];
			views = strtoul(pos, &pos, 10);

			noViewingsPerUser[u] = views;

			bigRatingsByUser[u] = new q[views];

			f.getline(&readLine[0], 3999999); //Comma-separated viewings

			pos = &readLine[0];

			for (int i = 0; i < views; i++) {
				bigRatingsByUser[u][i] = strtoull(pos, &pos, 10);
				pos = &pos[1];
			}
		}
		f.close();
	} else {
		srand(time(0));
		rand();

		for (user u = 0; u < U; u++) {
			noViewingsPerUser[u] = 3;
			bigRatingsByUser[u] = new q[3];

			movie m0, m1, m2;

			m0 = (rand() % M);
			do {
				m1 = (rand() % M);
			} while (m1 == m0);
			do {
				m2 = (rand() % M);
			} while (m2 == m0 || m2 == m1);

			bigRatingsByUser[u][0] = GetBigRating(1 + (rand() % 5), m0,
					rand() % 5, rand() % 50, rand() % 7, 0, rand() % 101,
					rand() % 101, rand() % 101, rand() % 101, rand() % 101);
			bigRatingsByUser[u][1] = GetBigRating(1 + (rand() % 5), m1,
					rand() % 5, rand() % 50, rand() % 7, 0, rand() % 101,
					rand() % 101, rand() % 101, rand() % 101, rand() % 101);
			bigRatingsByUser[u][2] = GetBigRating(1 + (rand() % 5), m2,
					rand() % 5, rand() % 50, rand() % 7, 0, rand() % 101,
					rand() % 101, rand() % 101, rand() % 101, rand() % 101);
		}
	}
}

void NetflixDataNG::readBigRatings(user** &usersByMovie,
		ulong* &noViewingsPerMovie, rate** &ratesByMovie, q** &qsByUser,
		ushort* &noViewingsPerUser, ushort ** &whatIndexAMovieIsToAUser) {
	readBigRatings(noViewingsPerUser, qsByUser);

	// This is new stuff for getting movie-perspective ratings from use-perspective ones.
	// Make vector arrays for storing the ratings movie-side.
#define MOVIESATONCE 500

	vector<user> musers[MOVIESATONCE];
	vector<rate> mrates[MOVIESATONCE];
	vector < ushort > mindices[MOVIESATONCE];

#pragma omp parallel for firstprivate(musers, mrates, mindices)
	for (uint i = 0; i <= M / MOVIESATONCE; i++) {
		// Limit to number of movies dealt with by one iteration is either MOVIESATONCE or
		// on last iteration 270 ( = M - i*MOVIESATONCE)
		uint jLim = min((uint) MOVIESATONCE, M - i * MOVIESATONCE);

		// Make sure vectors are empty
		for (uint jj = 0; jj < jLim; jj++) {
			musers[jj].clear();
			mrates[jj].clear();
			mindices[jj].clear();
		}

		// Go through each user and index & add rating and user to the vectors
		for (user u = 0; u < U; u++) {
			for (uint index = 0; index < noViewingsPerUser[u]; index++) {
				movie m = qsByUser[u][index] % MOV;
				if (m / MOVIESATONCE == i) {
					musers[m % MOVIESATONCE].push_back(u);
					mrates[m % MOVIESATONCE].push_back(
							qsByUser[u][index] % RATE);
					mindices[m % MOVIESATONCE].push_back(index);
				}
			}
		}

		// Go through each movie in this range
		for (uint jj = 0; jj < jLim; jj++) {
			// Get movie number & number of ratings
			uint views = mrates[jj].size();
			movie m = i * MOVIESATONCE + jj;

			// Fill in user info
			noViewingsPerMovie[m] = views;
			usersByMovie[m] = new user[views];
			ratesByMovie[m] = new rate[views];
			whatIndexAMovieIsToAUser[m] = new ushort[views];

			for (uint index = 0; index < views; index++) {
				usersByMovie[m][index] = musers[jj][index];
				ratesByMovie[m][index] = mrates[jj][index];
				whatIndexAMovieIsToAUser[m][index] = mindices[jj][index];
			}
		}
	}
#pragma omp barrier
}

// Data utility functions
string NetflixDataNG::getBaseSaveFileName() {
	stringstream s("");
	s << basePathToUse;
	s << pathToResults;
	return s.str();
}

void NetflixDataNG::writeQualifyingFile(string qualifyingFileName,
		string fullStatsFileName, model::DotProductGaussianAddNG* node) {
	try {
		// Open qualifying datafile and prepare both output file writers
		individualReader.open(&qualDataFile[0]);
		ofstream qualWriter;
		ofstream fullWriter;
		qualWriter.open(&qualifyingFileName[0]);
		fullWriter.open(&fullStatsFileName[0]);

		cerr << "\nQualifying data file: " << qualDataFile << "\n";
		cerr << "Qualifying file: " << qualifyingFileName << "\n";
		cerr << "Full Stats file: " << fullStatsFileName << "\n";

		// Start reading, make variables for the prediction
		char line[100];
		individualReader.getline(line, 99);
		unsigned int position;
		char* pEnd;
		bool movieNoLine = false;

		movie movieNumber = 0;
		user userNumber = 0;
		double rating = 0.0;
		double stats[] = { 0.0, 0.0 };

		int lObservations = 0;
		double lTotalMean = 0.0;
		double lTotalMeanSquared = 0.0;
		double lTotalBoundedMean = 0.0;
		double lTotalBoundedMeanSquared = 0.0;
		double lTotalPredictionVariance = 0.0;

		do {
			pEnd = NULL;
			// Read line & check if it's a line with a movie number on it
			for (position = 0; position < 100; position++) {
				if (line[position] == '\0') {
					movieNoLine = false;
					break;
				} else if (line[position] == ':') {
					//If it is, record the movie number, and output the line to both output files
					pEnd = &line[position];
					movieNumber = convertNxMovieToHcMovie(
							strtoul(line, &pEnd, 10));
					movieNoLine = true;
					qualWriter << line << "\n";
					fullWriter << line << "\n";
					break;
				}
			}
			if (!movieNoLine) {
				// On normal lines, get user and predict stats as appropriate
				userNumber = convertNxUserToHcUser(strtoul(line, &pEnd, 10));
				node->predictUnknownIndex(stats, userNumber, movieNumber);

				// Output rating (in right range) and stats
				rating = stats[0];
				if (rating > 5.0) {
					rating = 5.0;
				} else if (rating < 1.0) {
					rating = 1.0;
				}

				qualWriter << rating << '\n';
				fullWriter << stats[0] << ',' << stats[1] << '\n';

				lObservations++;
				lTotalMean += stats[0];
				lTotalMeanSquared += stats[0] * stats[0];
				lTotalBoundedMean += rating;
				lTotalBoundedMeanSquared += rating * rating;
				lTotalPredictionVariance += stats[1] - stats[0] * stats[0];
			}
			individualReader.getline(&line[0], 99);
		} while (line[0] != '\0');

		// At EOF close all streams
		qualWriter.close();
		fullWriter.close();
		individualReader.close();

		// Output statistics
		double lObsDouble = (double) lObservations;
		lTotalMean /= lObsDouble;
		lTotalMeanSquared /= lObsDouble;
		lTotalBoundedMean /= lObsDouble;
		lTotalBoundedMeanSquared /= lObsDouble;
		lTotalPredictionVariance /= lObsDouble;

		cerr << "\n\n===Prediction Statistics===\n";
		cerr << "Observations      = " << lObservations << "\n";
		cerr << "E(predicted)      = " << lTotalMean << "\n";
		cerr << "Var(predicted)    = "
				<< (lTotalMeanSquared - lTotalMean * lTotalMean) << "\n";
		cerr << "E(Var(predicted)) = " << lTotalPredictionVariance << "\n";
		cerr << "E(|predicted|)    = " << lTotalBoundedMean << "\n";
		cerr << "Var(|predicted|)  = "
				<< (lTotalBoundedMeanSquared
						- lTotalBoundedMean * lTotalBoundedMean) << "\n";
		cerr << "===========================\n";
	} catch (...) {
		cerr << ("Couldn't read the prediction data-file");
	}
}

void NetflixDataNG::getAllBigPredictionData(vector<user> &users,
		vector<q> &preds) {
	ifstream f;
	f.open(HcPredictionFile);

	char line[100];
	f.getline(line, 99);
	char* pEnd;

	do {
		users.push_back(strtoul(line, &pEnd, 10));
		pEnd = &pEnd[1];
		preds.push_back(strtoull(pEnd, &pEnd, 10));
		f.getline(line, 99);
	} while (line[0] != '\0' && line[1] != '\0');

	f.close();
}

void NetflixDataNG::getAllPredictionData(vector<user> &users,
		vector<movie> &movies, vector<uint> &dates) {
	// Open qualifying datafile
	individualReader.open(&qualDataFile[0]);

	// Start reading, make variables for the prediction
	char line[100];
	individualReader.getline(line, 99);
	unsigned int position;
	char* pEnd;
	bool movieNoLine = false;

	movie movieNumber = 0;
	user userNumber = 0;

	do {
		pEnd = NULL;
		// Read line & check if it's a line with a movie number on it
		for (position = 0; position < 100; position++) {
			if (line[position] == '\0') {
				movieNoLine = false;
				break;
			} else if (line[position] == ':') {
				// If it is, record the movie number, and output the line to both output files
				pEnd = &line[position];
				movieNumber = convertNxMovieToHcMovie(strtoul(line, &pEnd, 10));
				movieNoLine = true;
				break;
			}
		}
		if (!movieNoLine) {
			// On normal lines, get user and predict stats as appropriate
			userNumber = convertNxUserToHcUser(strtoul(line, &pEnd, 10));
			pEnd = &pEnd[1];
			int year = strtoul(pEnd, &pEnd, 10);
			pEnd = &pEnd[1];
			int month = strtoul(pEnd, &pEnd, 10);
			pEnd = &pEnd[1];
			int day = strtoul(pEnd, &pEnd, 10);

			int date = daysAfter30Sept1998[year - 1998][month - 1] + day;

			users.push_back(userNumber);
			movies.push_back(movieNumber);
			dates.push_back(date);
		}

		individualReader.getline(&line[0], 99);
	} while (line[0] != '\0');

	individualReader.close();
}

q NetflixDataNG::GetBigRating(rate rating, movie m, uint userEra, uint movieEra,
		uint weekday, uint precisionNumber, uint useHistory1, uint useHistory3,
		uint useHistory5, uint useHistory10, uint useHistory30) {
	q bigRatingSoFar = 0;
	q baseValue = 1;

	updateBigRating(bigRatingSoFar, baseValue, rating, RATE);
	updateBigRating(bigRatingSoFar, baseValue, m, MOV);
	updateBigRating(bigRatingSoFar, baseValue, userEra, ERAU);
	updateBigRating(bigRatingSoFar, baseValue, movieEra, ERAM);
	updateBigRating(bigRatingSoFar, baseValue, weekday, WKDY);
	updateBigRating(bigRatingSoFar, baseValue, precisionNumber, PREC);
	updateBigRating(bigRatingSoFar, baseValue, useHistory1, HIST1);
	updateBigRating(bigRatingSoFar, baseValue, useHistory3, HIST3);
	updateBigRating(bigRatingSoFar, baseValue, useHistory5, HIST5);
	updateBigRating(bigRatingSoFar, baseValue, useHistory10, HIST10);
	updateBigRating(bigRatingSoFar, baseValue, useHistory30, HIST30);

	return bigRatingSoFar;
}

void NetflixDataNG::updateBigRating(q &currentBigRatingValue,
		q &currentBaseValue, q newValue, q newBaseValue) {
	for (q i = 0; i < newBaseValue; i++) {
		if ((currentBigRatingValue % newBaseValue) == newValue) {
			break;
		}
		currentBigRatingValue += currentBaseValue;
	}

	if ((currentBigRatingValue % newBaseValue) != newValue) {
		// We went round the whole loop without success
		throw "BUGGER, check the bases are all prime";
	}

	currentBaseValue *= newBaseValue;
}

// Setting internal debug/test variables  
void NetflixDataNG::setDebug() {
	debug = true;
}

void NetflixDataNG::setTest() {
	test = true;
}

// Helper functions for IO
void NetflixDataNG::getQuartetFromUserData(movie& val, rate& r,
		string& fromThis, char* pEnd) {
	val = convertNxMovieToHcMovie(strtoul(&fromThis[0], &pEnd, 10));
	pEnd = &pEnd[1];
	r = (rate) strtoul(pEnd, &pEnd, 10);
}

void NetflixDataNG::getQuartetFromMovieData(user& val, rate& r,
		string& fromThis, char* pEnd) {
	val = convertNxUserToHcUser(strtoul(&fromThis[0], &pEnd, 10));
	pEnd = &pEnd[1];
	r = (rate) strtoul(pEnd, &pEnd, 10);
}

void NetflixDataNG::getInfoFromMovieData(user& val, rate& r,
		uint& daysAfter30September1998, string& fromThis, char* pEnd) {
	val = convertNxUserToHcUser(strtoul(&fromThis[0], &pEnd, 10));
	pEnd = &pEnd[1];
	r = (rate) strtoul(pEnd, &pEnd, 10);
	pEnd = &pEnd[1];
	int y = strtoul(pEnd, &pEnd, 10);
	pEnd = &pEnd[1];
	int m = strtoul(pEnd, &pEnd, 10);
	pEnd = &pEnd[1];
	int d = strtoul(pEnd, &pEnd, 10);

	daysAfter30September1998 = daysAfter30Sept1998[y - 1998][m - 1] + d;
}

uint NetflixDataNG::parseForInt(string s, uint startAt, uint endAt) {
	int toReturn = -1;
	char c;
	for (int i = startAt; i <= endAt; i++) {
		c = s[i];
		if (c <= '9' && c >= '0') {
			if (toReturn == -1) {
				toReturn = c - '0';
			} else {
				toReturn = toReturn * 10 + c - '0';
			}
		} else if (toReturn == -1) {
			continue;
		} else {
			return toReturn;
		}
	}
	return toReturn;
}

void NetflixDataNG::sevenDigitFormat(stringstream &s, int number) {
	int noDigitsInNumber = 0;
	int tenPower = 1;
	while (true) {
		if (tenPower > number) {
			break;
		}
		noDigitsInNumber++;
		tenPower *= 10;
	}
	for (; noDigitsInNumber <= 6; noDigitsInNumber++) {
		s << '0';
	}
	s << number;
}

// Helper functions
void NetflixDataNG::populateUserIdConverter(string path) {
	individualReader.open(&path[0], ifstream::in);
	string nextLine(100, '\0');
	while (true) {
		individualReader.getline(&nextLine[0], 99);
		if (nextLine[0] == '\0') {
			break;
		}
		unsigned int commaPosition;
		for (commaPosition = 0; commaPosition < nextLine.length();
				commaPosition++) {
			if (nextLine[commaPosition] == ',') {
				break;
			}
		}
		userIDConversion[parseForInt(nextLine, 0, commaPosition - 1) - 1] =
				parseForInt(nextLine, commaPosition + 1, nextLine.length() - 1);
	}
	individualReader.close();
}

string NetflixDataNG::getTrainingDataPath(uint userMovieNumber,
		bool trueIfUser) {
	stringstream s("");
	s << dataPath;
	if (trueIfUser) {
		s << "ur_";
		sevenDigitFormat(s, convertHcUserToNxUser(userMovieNumber));
		s << ".txt";
	} else {
		s << "mv_";
		sevenDigitFormat(s, convertHcMovieToNxMovie(userMovieNumber));
		s << ".txt";
	}
	return s.str();
}

// Functions for internal/external numbering
user NetflixDataNG::convertNxUserToHcUser(unsigned int netflixUserID) {
	unsigned int leftPosition = 0;
	unsigned int rightPosition = U - 1;
	unsigned int midPosition;
	unsigned int current;

	while (true) {
		midPosition = (rightPosition + leftPosition) / 2;
		current = userIDConversion[midPosition];

		if (current == netflixUserID) {
			return (user) midPosition;
		} else if (current < netflixUserID) {
			leftPosition = midPosition + 1;
		} else {
			rightPosition = midPosition - 1;
		}
	}
}

uint NetflixDataNG::convertHcUserToNxUser(user gaplessUserId) {
	return userIDConversion[gaplessUserId];
}

movie NetflixDataNG::convertNxMovieToHcMovie(unsigned int netflixMovieID) {
	return (movie) (netflixMovieID - 1);
}

uint NetflixDataNG::convertHcMovieToNxMovie(movie HcMovieID) {
	return HcMovieID + 1;
}

}
