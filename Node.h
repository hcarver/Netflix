#ifndef NODE_H_
#define NODE_H_

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "nfx.h"

using namespace std;

namespace model {

class Node {
public:
	// Node stuff
	virtual void IO(string description, istream& file, bool loadPriors,
			bool debug) = 0;
	virtual void IO(string description, ostream& file, bool loadPriors,
			bool debug) = 0;
	virtual double getBound() = 0;
	virtual void update(int index1, int index2, int index3,
			double message[]) = 0;
};

}

#endif /*NODE_H_*/
