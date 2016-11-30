/*
 * RDF2Walk.h
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#ifndef RDF2WALK_H_
#define RDF2WALK_H_

#include "Snap.h"
#include "GraphWeigher.h"
#include "GraphWalker.h"

class WalkSink {

protected:
	virtual ~WalkSink() {

	}
public:
	virtual void consume(TVec<TStr> walk) = 0;
};

//removed because the particular way is different each time.
//void computeWalks(TStr filename, GraphWeigher & weighingStrategy, GraphWalker& walker, int amount, WalkSink & sink);

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > buildWalkableGraph(TStr filename, GraphWeigher & weighingStrategy);

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > buildWalkableGraphIgnoringLiterals(TStr filename, GraphWeigher & weighingStrategy);


class TextFileSink: public WalkSink {
private:
	FILE *fout;

public:
	TextFileSink(FILE *fout) :
			fout(fout) {
	}

	virtual void consume(TVec<TStr> walk);
};

#endif /* RDF2WALK_H_ */
