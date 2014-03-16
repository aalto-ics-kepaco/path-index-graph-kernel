#ifndef _TBWTBUILDER_H_
#define _TBWTBUILDER_H_

#include "SimpleForest.h"
#include "BlockArray.h"
#include "Tools.h"

#include "BitRank.h" // remove

#include <iostream>
#include <vector>

class TBWTBuilder
{
public:
    typedef SimpleForest::node_p node_p;
    typedef std::vector<node_p> namev_t; // FIXME Uses 32 bits per node

    explicit TBWTBuilder(SimpleForest *t)
        : nf(t)
    { }

    void sort(bool verbose = false, bool debug = false);

    void debugTest()
    {
        ulong r = 0;
        ulong *bv = getLeaf();
        for (ulong i = 0; i < nf->numberOfNodes(); ++i)
            if (Tools::GetField(bv, 1, i))
                ++r;
        std::cerr << "myrank = " << r << std::endl;
        BitRank *b = new BitRank(bv, nf->numberOfNodes(), true);
        std::cerr << "debug: rank = " << b->rank(nf->numberOfNodes()-1) << std::endl;
        if (nf->numberOfLeaves() != b->rank(nf->numberOfNodes()-1))
        {
            std::cerr << "debug: assert failed: number of leaves was " << nf->numberOfLeaves() << " but rank was " << b->rank(nf->numberOfNodes()-1) <<  std::endl;
        }
        delete b;
    }

    // call sort() first
    ulong * getLast() const; 
    ulong * getLeaf() const;
    uchar * getTBWTInternal() const;
    uchar * getTBWTLeaf() const;
    void getCF(unsigned *, unsigned *) const;
    BlockArray * getLeafEntry() const;
    BlockArray * getLastEntry() const;

private:
    void init();
    void debugPrint(unsigned, std::string msg = "");

    SimpleForest *nf;
    namev_t name;
};

#endif
