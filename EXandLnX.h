#ifndef EXANDLNX_H_
#define EXANDLNX_H_

#include "EX.h"
#include "Node.h"

namespace model {

class EXandLnX: public EX, public Node {
public:
	virtual double getELnX() = 0;
};
}

#endif /*EXANDLNX_H_*/
