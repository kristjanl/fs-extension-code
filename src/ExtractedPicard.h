#ifndef EXTRACTEDPICARD_H_
#define EXTRACTEDPICARD_H_

#include "Continuous.h"
#include "MyLogger.h"

class ExtractedPicard: public Flowpipe {
  public:
    ExtractedPicard(const vector<Interval> & box, const Interval & I);
};



#endif /* EXTRACTEDPICARD_H_ */
