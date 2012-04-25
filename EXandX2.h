#ifndef EXANDX2_H_
#define EXANDX2_H_

#include "EX.h"
#include "Node.h"

namespace model {

class EXandX2: public EX, public Node {
public:
	virtual double getEX2() = 0;
	virtual double* getMoments() = 0;
};

}

#endif /*EXANDX2_H_*/
