#ifndef NFX_H_
#define NFX_H_

#include <string>

// MODEL STRING
#define MODELSTRING "NgPmfAddAllStochastic"

// Define the types of user, movie and rating, and the type of iterators (unsigned int)
typedef unsigned long user;
typedef unsigned short movie;
typedef unsigned char rate;
typedef unsigned long long int q;

// These are a series of coprime numbers which act as keys, allowing extraction of
// a whole load of data about each rating from a single 64-bit unsigned long
#define RATE 7
#define MOV 17783
#define ERAU 5
#define ERAM 53
#define WKDY 11
#define PREC 3

#define HIST1  101
#define HIST3  103
#define HIST5  107
#define HIST10 109
#define HIST30 113

// Array indices for user history access
#define HIST1ACCESS  0
#define HIST3ACCESS  1
#define HIST5ACCESS  2
#define HIST10ACCESS 3
#define HIST30ACCESS 4

// Define some important values for the algorithm and some useful constants
#define K 60
#define U 480189
#define M 17770

#define OMP_NUM_THREADS 4
#define LOG2PI 1.837877066409345
#define PI 3.141592653589793238462643383279502884197169399375

// Frequency to save progress and calculate the bound
#define saveAndBoundEvery 6

// Define useful path strings
#define beastBasePath            "/home/hywel/workspace/"
#define HcTrainingFile           "/home/hywel/workspace/data/HcTrainingFile.txt"
#define HcPredictionFile         "/home/hywel/workspace/data/HcPredictionFile.txt"
#define pathToData               "training_set/training_set/"
#define pathToConversionFile     "data/conversionfile"
#define pathToQualifyingDataFile "download/qualifying.txt"
#define pathToStats              "data/"
#define pathToResults            "Results/"
#define pathToVitalStatistics    "data/meandAndMore"

#endif /*NFX_H_*/
