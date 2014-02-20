#ifndef _SIMPLEFOREST_H_
#define _SIMPLEFOREST_H_

#include "BitRank.h"
#include "Tools.h"

#include <iostream>
#include <string>

// Like NaiveForest but uses about half the space,
// (and supports less operations).
class SimpleForest
{
public:
    typedef unsigned node_p; // Data type for accessing bit-packed nodes

    // Accessors
    inline node_p getParent(node_p node) const
    {
        return Tools::GetVariableField(data, psize, (ulong)node*bsize);
    }
    // Parent pointers are modified during construction
    inline void setParent(node_p node, node_p p)
    {
        Tools::SetVariableField(data, psize, (ulong)node*bsize, p);
    }
    inline node_p getName(node_p node) const
    {
        return Tools::GetVariableField(data, psize, (ulong)node*bsize+psize);
    }
    inline void setName(node_p node, node_p n)
    {
        Tools::SetVariableField(data, psize, (ulong)node*bsize+psize, n);
    }
    inline uchar getC(node_p node) const
    {
        return Tools::GetVariableField(data, csize, (ulong)node*bsize+2*psize);
    }
    inline unsigned getEntry(node_p node) const
    {
        return fentry->rank(node) - 1;
    }

    inline bool isRoot(node_p n) const
    { return (getParent(n) >= n ? true : false); }

    // Parents are updated during construction, test to check for root in original tree
    inline bool isOrigRoot(node_p n) const
    { return (getParent(n) == t ? true : false); }

    inline bool isLeaf(node_p node) const
    { 
        return Tools::GetVariableField(data, 1, (ulong)node*bsize+2*psize+csize);  
    }
    inline bool isLast(node_p node) const
    {
        return Tools::GetVariableField(data, 1, (ulong)node*bsize+2*psize+csize+1); 
    }

private:
    // More accessors
    inline void setC(node_p node, uchar c)
    {
        Tools::SetVariableField(data, csize, (ulong)node*bsize+2*psize, c);
    }
    inline void setLeaf(node_p node, bool p)
    {
        Tools::SetVariableField(data, 1, (ulong)node*bsize+2*psize+csize, p);
    }
    inline void setLast(node_p node, bool p)
    {
        Tools::SetVariableField(data, 1, (ulong)node*bsize+2*psize+csize+1, p);
    }
    
    void initNode(node_p node, node_p parent, uchar c)
    {        
        setParent(node, parent);
        setC(node, c);
    }

    node_p parse(unsigned &i, std::string const &bp, node_p parent, unsigned depth)
    {
        /**
         * Sanity checks
         */
        if (i+1 >= bp.size())
        {
            std::cerr << "SimpleForest::parse(): unable to parse tree structure, unexpected end of string at " 
                      << i << " of string:" << std::endl << bp << std::endl;
            std::abort();
        }
        if (bp[i] != '(')
        {
            std::cerr << "SimpleForest::parse(): unable to parse tree structure, expecting '(' at position " 
                      << i << " of string:" << std::endl << bp << std::endl;
            std::abort();
        }

        /**
         * Init node at freeBlock
         */
        unsigned node = freeBlock++;
        initNode(node, parent, bp[++i]);

        if (depth > h)
            h = depth;

        // Parse children if any
        if (bp[++i] == '(')
            parse(i, bp, node, depth + 1);
        else
        {
            setLeaf(node, 1);
            ++l; // incr. number of leaves
        }

        /**
         * Sanity checks
         */
        if (i >= bp.size())
        {
            std::cerr << "SimpleForest::parse(): unable to parse tree structure, unexpected end of string at " 
                      << i << " of string:" << std::endl << bp << std::endl;
            std::abort();
        }
        if (bp[i] != ')')
        {
            std::cerr << "SimpleForest::parse(): unable to parse tree structure, expecting ')' at position " 
                      << i << " of string:" << std::endl << bp << std::endl;
                std::abort();
        }

        // Parse siblings if any
        if (i+1 < bp.size() && bp[++i] == '(')
            parse(i, bp, parent, depth);
        else
            setLast(node, 1);

        return node;
    }

    void printNode(node_p &n, unsigned d) const
    {
        for (unsigned i = 0; i < d; ++i)
            std::cerr << " ";
        std::cerr << getC(n) << std::endl;

        bool last = isLast(n);
        if (!isLeaf(n))
            printNode(++n, d+1);
        if (!last)
            printNode(++n, d);
    }
        
public:

    SimpleForest(unsigned entries_, ulong trees_, ulong t_)
        : data(0), freeBlock(0), psize(0), csize(8), bsize(0), 
        entries(entries_), trees(trees_), t(t_), h(0), l(0), fentry(0)
    {
        /**
         * One node consists of 
         *
         *    Node *parent;     ~0lu if no parent
         *    Node *name;       Lexicographical rank
         *    uchar c;
         *    bit leaf;
         *    bit last;
         *
         * we allocate ceillog(t) bits for each Node pointer,
         * plus 8 bits for the character.
         */
        psize = Tools::CeilLog2(t + 1 > 255 ? t + 1 : 256);
        bsize = 2 * psize + csize + 2;

        // Debug printing
        std::cerr << "SimpleForest() debug: Using bsize = " << bsize 
                  << ", psize = " << psize << ", csize = " << csize << std::endl
                  << "Allocating " << (t * bsize/W + 1) * sizeof(ulong) << " bytes for tree data." << std::endl;

        // Reserve enough data blocks
        data = new ulong[t * bsize/W + 1];
        for (ulong i = 0; i < t * bsize/W + 1; ++i)
            data[i] = 0;
    }

    void add(std::string const &bp)
    {
        unsigned i = 0; // Iterator over bp
        // t is given as "parent" of the root
        parse(i, bp, t, 1);
        if (i != bp.size() - 1)
        {
            std::cerr << "SimpleForest::parse(): unable to parse tree structure, expecting end of line at position " 
                      << i << " of string of size " << bp.size() << ":" << std::endl << bp << std::endl;
            std::abort();
        }
    }

    void setFirstEntry(ulong *entry)
    {
        fentry = new BitRank(entry, t, true);
    }

    ~SimpleForest()
    {
        delete [] data;
        data = 0;
        delete fentry;
        fentry = 0;
    }

    node_p getRoot() const
    { return 0; }

    ulong numberOfNodes() const
    { return t; }
    ulong numberOfLeaves() const
    { return l; }
    inline unsigned getHeight() const
    { return h; }
    unsigned numberOfEntries() const
    { return entries; }
    ulong numberOfTrees() const
    { return trees; }

    void debugPrint() const
    {
        node_p root = 0;
        while (root < t)
        {
            printNode(root, 0);
            do 
                ++root;
            while (root < t && !isRoot(root));
        }
    }
private:
    ulong *data;
    ulong freeBlock;  // Pointer to the next free data block
    unsigned psize;   // Node pointer size, ceil(log(t)) bits.
    unsigned csize;   // Character size, log(alphabet size) bits
    unsigned bsize;   // Block size, 4 * psize + csize bits.

    unsigned entries; // Number of entries
    ulong trees;      // Number of trees
    ulong t;          // Number of nodes
    unsigned h;       // Height
    ulong l;          // Number of leaves
    BitRank *fentry;  // Bitvector marking the first root node in each entry
};

#endif
