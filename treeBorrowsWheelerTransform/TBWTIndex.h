#ifndef _TBWTINDEX_H_
#define _TBWTINDEX_H_

#include "BitRank.h"
#include "HuffWT.h"
#include "BlockArray.h"
#include "Tools.h"

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <utility>
#include <cassert>

class TBWTIndex
{
public:
    // Pointer to one node
    typedef ulong node_p;
    // and a range of nodes
    typedef std::pair<node_p, node_p> node_r;

    explicit TBWTIndex(std::string &);
    ~TBWTIndex();

    bool isLeaf(node_p node) const
    {
        if (leaf->IsBitSet(node))
            return true; // Silly, but test.cpp may fail for returning the value from IsBitSet()
        return false;
    }

    bool isRoot(node_p node) const
    {
        if (node < trees)
            return true;
        return false;
    }

    node_p getRoot(ulong tree = 0) const
    { return tree; }
    uchar getC(node_p node) const
    { 
        if (leaf->IsBitSet(node))
        {
            node = leaf->rank(node) - 1;
            return wtleaf->access(node); 
        }
        
        node = leaf->rank0(node) - 1;
        return wt->access(node);
    }

    void getSubtrees(node_r const &range, std::set<uchar> &C) const
    {
        // FIXME Is std::set too slow?
        C.clear();
        node_p sp = 0;
        if (range.first > 0)
            sp  = leaf->rank0(range.first - 1);
        node_p ep = leaf->rank0(range.second);
        for (; sp < ep; ++sp)
            C.insert(wt->access(sp));
    }
    // FIXME clean up!
    node_r getSubtree(node_r const &range, uchar c) const
    {
        node_p sp = 0;
        if (range.first > 0)
            sp  = leaf->rank0(range.first - 1);
        node_p ep = leaf->rank0(range.second) - 1;

        if (sp > ep)
            return std::make_pair(sp,ep);

        /* replace with rank and select sp = wt->rank(c, sp-1);
           ep = wt->rank(c, ep); */
        // Linear FIXME
/*        node_p nsp = sp;
        node_p nep = ep;
        std::cerr << "nsp = " << nsp << ", nep = " << nep << std::endl;
        while (wt->access(nsp) != c && nsp < nep)
            ++nsp;
        while (wt->access(nep) != c && nep >= nsp)
        --nep;*/

        if (sp > 0)
            sp = wt->rank(c, sp-1);
        ep = wt->rank(c, ep);
        if (sp >= ep)
            return std::make_pair(sp+1,ep);
        sp = wt->select(c, sp+1);
        ep = wt->select(c, ep);

/*        if (sp != nsp || ep != nep)
        {
            std::cerr << "c = " << c << ", sp = " << sp << ", ep = " << ep << ", nsp = " << nsp << ", nep = " << nep <<  std::endl;
            for (unsigned i = 0; i < 10; ++i)
                std::cerr << wt->access(i);
            std::cerr << std::endl;
            std::abort();
            }
*/

        assert (sp <= ep);
        // Map the range back
        sp = leaf->select0(sp+1);
        ep = leaf->select0(ep+1);
        node_r r;
        r.first = getRankedChild(sp, 1);
        r.second = getRankedChild(ep, getDegree(ep));
        assert(r.first <= r.second);
        return r;
    }

    unsigned getSubtreeSize(node_r const &range) const
    {
        ulong r = 0;
        if (range.first > 0)
            r = last->rank(range.first - 1);
        return last->rank(range.second) - r;
    }

    void getLeaves(node_r const &range, std::set<uchar> &C) const
    {
        // FIXME Is std::set too slow?
        C.clear();
        node_p sp = 0;
        if (range.first > 0)
            sp  = leaf->rank(range.first - 1);
        node_p ep = leaf->rank(range.second);
        for (; sp < ep; ++sp)
            C.insert(wtleaf->access(sp));
    }
    std::map<unsigned,unsigned> getLeafFreq(node_r const &range, uchar c)
    {
        node_p sp = 0;
        if (range.first > 0)
            sp  = leaf->rank(range.first - 1);
        node_p ep = leaf->rank(range.second);
        
// FIXME Linear...
        std::map<unsigned, unsigned> freq;
        for (; sp < ep; ++sp)
            if (wtleaf->access(sp) == c)
                freq[(*leafEntry)[sp]] ++;
        
        return freq;
    }
    unsigned getLeafCount(node_r const &range, uchar c)
    {
        node_p sp = 0;
        if (range.first > 0)
            sp  = leaf->rank(range.first - 1);
        node_p ep = leaf->rank(range.second);
        if (sp > ep)
            return 0;

        if (sp)
            sp = wtleaf->rank(c, sp - 1);
        return wtleaf->rank(c, ep - 1) - sp;
    }
    unsigned getLeafCount(node_r const &range)
    {
        node_p sp = 0;
        if (range.first > 0)
            sp  = leaf->rank(range.first - 1);
        node_p ep = leaf->rank(range.second);

        return ep-sp;
    }


    std::map<unsigned,unsigned> getInternalFreq(node_r const &range)
    {
        node_p sp = 0;
        if (range.first > 0)
            sp  = last->rank(range.first - 1) - 1;
        node_p ep = last->rank(range.second) - 1;
        
// FIXME Linear...
        std::map<unsigned, unsigned> freq;
        for (; sp < ep; ++sp)
            freq[(*lastEntry)[sp]] ++;
        
        return freq;
    }
    unsigned getInternalCount(node_r const &range)
    {
        node_p sp = 0;
        if (range.first > 0)
            sp  = last->rank(range.first - 1);
        node_p ep = last->rank(range.second);

        return ep-sp;
    }


    node_r getChildren(node_p) const;
    node_p getRankedChild(node_p node, ulong k) const
    {
        node_r c = getChildren(node);
        if (k > c.second - c.first + 1)
        {
            std::cerr << "TBWTIndex::getRankedChild(): assert failed!" <<std::endl;
            std::abort();
        }
        return c.first + k - 1;
    }

    unsigned getDegree(node_p node) const
    {
        node_r range = getChildren(node);
        return range.second - range.first + 1;
    }

    ulong numberOfNodes() const
    { return t; }
    ulong numberOfLeaves() const
    { return leaf->rank(t-1); }
    unsigned numberOfEntries() const
    { return entries; }
    ulong numberOfTrees() const
    { return trees; }

private:
    unsigned entries;
    ulong trees;
    ulong t;
    unsigned F[256];
    BitRank *leaf;
    BitRank *last;
    HuffWT *wt,  // For internal nodes 
        *wtleaf; // For leaf nodes
    BlockArray *leafEntry;
    BlockArray *lastEntry;
    TBWTIndex();
};

#endif
