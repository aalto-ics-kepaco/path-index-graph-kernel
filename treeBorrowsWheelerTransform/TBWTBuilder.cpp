#include "TBWTBuilder.h"

#include <algorithm>

using namespace std;

class CompareLess {
    SimpleForest *nf;
public:
    CompareLess(SimpleForest *nf_) 
        : nf(nf_)
    { }

    inline unsigned parent_name(TBWTBuilder::node_p p) const
    {
        if (nf->isRoot(p))
            return 0;
        return nf->getName(nf->getParent(p));
    }
    
    bool operator()(TBWTBuilder::node_p i, TBWTBuilder::node_p j) const
    {
        TBWTBuilder::node_p namei = nf->getName(i);
        TBWTBuilder::node_p namej = nf->getName(j);
        if (namei < namej)
            return true;

        if (namei == namej && parent_name(i) < parent_name(j))
            return true;
        return false;
    }
};


void TBWTBuilder::debugPrint(unsigned iter, string msg)
{
    cerr << "-------------------------------" << endl;
    if (msg != "")
        cerr << msg << endl;
    else
        cerr << "iter = " << iter << endl;
    for (ulong i = 0; i < nf->numberOfNodes(); ++i)
    {
        cerr << "i = " << i << ", " << nf->isLast(name[i]) << ", " << nf->isLeaf(name[i]) << ", " << nf->getC(name[i]) << ", ";
        node_p p = name[i];
        while (!nf->isRoot(p))
        {
            p = nf->getParent(p);
            cerr << nf->getC(p);
        }        
        cerr << ", " << nf->getName(name[i]) << endl;
    }
}

void TBWTBuilder::sort(bool verbose, bool debug)
{
    ulong t = nf->numberOfNodes();
    if (t >= (1lu << (8*sizeof(unsigned)))-1)
    {
        cerr << "TBWTBuilder::sort(): number of nodes is higher than data type allows! Switch to 64-bit datatype." <<  endl;
        abort();
    }

    if (!name.empty())
    {
        cerr << "TBWTBuilder::sort(): assert failed: sort() was called multiple times?" <<  endl;
        abort();
    }
    name = namev_t(t);

    if (verbose) cerr << "TBWTBuilder: Initializing..." << endl;
    init();

    if (debug)
        debugPrint(0, "After init");

    ulong *incr = new ulong[t/W+1];

    for (unsigned iter = 0; (1u << iter) < nf->getHeight(); ++iter)
    {
        if (verbose) cerr << "TBWTBuilder: Sorting step " << iter << " out of " << Tools::CeilLog2(nf->getHeight()) << endl;
        stable_sort (name.begin(), name.end(), CompareLess(nf));

        if (verbose) cerr << "TBWTBuilder: Sorting done, renaming..." << endl;

        /**
         * Renaming
         */
        for (ulong i = 0; i < t/W+1; ++i)
            incr[i] = 0;
        CompareLess cl(nf);
        for (ulong j = 0; j < t; ++j)
            if (j < t-1 && cl(name[j], name[j+1]))
                Tools::SetField(incr, 1, j, 1);

        node_p i = 0;
        for (ulong j = 0; j < t; ++j)
        {
            nf->setName(name[j], i);
            if (Tools::GetField(incr, 1, j))
                ++i;
        }

        /**
         * Update parent pointers
         */
        if ((1u << (iter+1)) < nf->getHeight())
            for (ulong j = t; j > 0;)
            {
                --j;
                if (!nf->isRoot(j))
                {
                    node_p parent = nf->getParent(j);
                    if (nf->isRoot(parent))
                        nf->setParent(j, j);
                    else
                        nf->setParent(j, nf->getParent(parent));
                }
            }

        if (debug)
            debugPrint(iter);
    }
    if (debug)
        debugPrint(0, "After sort");
    
    delete [] incr;

    if (verbose) cerr << "TBWTBuilder: Sorting completed." << endl;
}

ulong * TBWTBuilder::getLast() const
{
    if (name.empty())
    {
        cerr << "TBWTBuilder::getTBWT(): assert failed: sort() was not called?" <<  endl;
        abort();
    }
    ulong t = nf->numberOfNodes();
    ulong *bv = new ulong[t/W+1];
    for (ulong i = 0; i < t/W+1; ++i)
        bv[i] = 0;
    // Roots have last == 0
    // except the last root.
    for (ulong i = nf->numberOfTrees()-1; i < t; ++i) 
        Tools::SetField(bv, 1, i, nf->isLast(name[i]));
    return bv;
}

ulong * TBWTBuilder::getLeaf() const
{
    if (name.empty())
    {
        cerr << "TBWTBuilder::getTBWT(): assert failed: sort() was not called?" <<  endl;
        abort();
    }
    ulong t = nf->numberOfNodes();
    ulong *bv = new ulong[t/W+1];
    for (ulong i = 0; i < t/W+1; ++i)
        bv[i] = 0;
    
    for (ulong i = 0; i < t; ++i)
        Tools::SetField(bv, 1, i, nf->isLeaf(name[i]));
    return bv;
}

uchar * TBWTBuilder::getTBWTInternal() const
{
    if (name.empty())
    {
        cerr << "TBWTBuilder::getTBWTInternal(): assert failed: sort() was not called?" <<  endl;
        abort();
    }
    ulong t = nf->numberOfNodes() - nf->numberOfLeaves();
    uchar *tbwt = new uchar[t];
    ulong j = 0;
    for (ulong i = 0; i < nf->numberOfNodes(); ++i)
        if (!nf->isLeaf(name[i]))
            tbwt[j++] = nf->getC(name[i]);

    if (j != t)
    {
        cerr << "BWTBuilder::getTBWTInternal(): assert failed: j != t" << endl;
        abort();
    }
    return tbwt;
}

uchar * TBWTBuilder::getTBWTLeaf() const
{
    if (name.empty())
    {
        cerr << "TBWTBuilder::getTBWTLeaf(): assert failed: sort() was not called?" <<  endl;
        abort();
    }
    ulong t = nf->numberOfLeaves();
    uchar *tbwt = new uchar[t];
    ulong j = 0;
    for (ulong i = 0; i < nf->numberOfNodes(); ++i)
        if (nf->isLeaf(name[i]))
            tbwt[j++] = nf->getC(name[i]);

    if (j != t)
    {
        cerr << "BWTBuilder::getTBWTLeaf(): assert failed: j != t" << endl;
        abort();
    }
    return tbwt;
}

void TBWTBuilder::getCF(unsigned *C, unsigned *F) const
{
    if (name.empty())
    {
        cerr << "TBWTBuilder::getCF(): assert failed: sort() was not called?" <<  endl;
        abort();
    }
    C[0] = nf->numberOfTrees();
    for (unsigned i = 1; i < 256; ++i)
        C[i] = 0;

    ulong t = nf->numberOfNodes();
    for (ulong i = 0; i < t; ++i)
        if (!nf->isLeaf(name[i]))
            C[nf->getC(name[i])] ++;

    // Init F
    F[0] = 0;
    for (unsigned i = 0; i < 255; ++i)
    {
        unsigned s = 0, j = F[i];
        while (s != C[i])
            if (nf->isLast(name[j++]))
                ++s;
        F[i+1] = j;
    }
}


// Parameter entry is a bitvector marking the first tree in each entry
BlockArray * TBWTBuilder::getLeafEntry() const
{
    ulong t = nf->numberOfLeaves();
    BlockArray *le = new BlockArray(t, Tools::CeilLog2(t));

    ulong j = 0;
    for (ulong i = 0; i < nf->numberOfNodes(); ++i)
        if (nf->isLeaf(name[i]))
            (*le)[j++] = nf->getEntry(name[i]);
    
    if (j != t)
    {
        cerr << "BWTBuilder::getLeafEntry(): assert failed: j != t" << endl;
        abort();
    }

    return le;
}

BlockArray * TBWTBuilder::getLastEntry() const
{
    ulong t = nf->numberOfNodes() - nf->numberOfLeaves(); // Number of internal nodes
    BlockArray *le = new BlockArray(t, Tools::CeilLog2(t));

    ulong j = 0;
    // Skipping all root nodes
    for (ulong i = nf->numberOfTrees(); i < nf->numberOfNodes(); ++i)
        if (nf->isLast(name[i]))
            (*le)[j++] = nf->getEntry(name[i]);

    if (j != t)
    {
        cerr << "BWTBuilder::getLastEntry(): assert failed: j != t" << endl;
        abort();
    }
    return le;
}


void TBWTBuilder::init()
{
    //unsigned tree = 0;
    //unsigned trees = nf->numberOfTrees();
    ulong t = nf->numberOfNodes();

    for (node_p i = 0; i < t; ++i)
    {
        name[i] = i;

        /**
         * Initialize names
         *  -  Root nodes have name 0               (orig. in [0, trees-1])
         *  -  Other nodes have names in [0,255]    (orig. in [trees, trees + 255])
         */
        if (nf->isRoot(i))
            nf->setName(i, 0); // was tree++);
        else
        {
            node_p parent = nf->getParent(i);
            nf->setName(i, (unsigned)nf->getC(parent)); // was trees + ...
        }
    }

    /* Old version:
       name[i] = node;
    nsigned n = 0;
    if (!nf->isRoot(node))
        n = node->getParent()->getC();
    name[i]->setName(n);
    ++i;
    if (!node->isLeaf())
        init(name, node->getChild(), i);
    if (!node->isLast())
        init(name, node->getSibling(), i);
    */
}

