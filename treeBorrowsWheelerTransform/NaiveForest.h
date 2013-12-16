#ifndef _NAIVEFOREST_H_
#define _NAIVEFOREST_H_

#include "Tools.h"

#include <iostream>
#include <string>

class NaiveForest
{
public:
    typedef unsigned node_p; // Data type for node pointers

    // Accessors
    inline node_p getParent(node_p node) const
    {
        return Tools::GetVariableField(data, psize, (ulong)node*bsize);
    }
    inline node_p getChild(node_p node) const
    {
        return Tools::GetVariableField(data, psize, (ulong)node*bsize+psize);
    }
    inline node_p getSibling(node_p node) const
    {
        return Tools::GetVariableField(data, psize, (ulong)node*bsize+2*psize);
    }

    inline node_p getName(node_p node) const
    {
        return Tools::GetVariableField(data, psize, (ulong)node*bsize+3*psize);
    }
    inline void setName(node_p node, node_p n)
    {
        Tools::SetVariableField(data, psize, (ulong)node*bsize+3*psize, n);
    }

    inline unsigned getEntry(node_p node) const
    {
        return Tools::GetVariableField(data, esize, (ulong)node*bsize+4*psize);
    }
    inline uchar getC(node_p node) const
    {
        return Tools::GetVariableField(data, csize, (ulong)node*bsize+4*psize+esize);
    }

    inline bool isRoot(node_p n) const
    { return (getParent(n) > n ? true : false); }
    inline bool isLeaf(node_p n) const
    { return (getChild(n) == 0 ? true : false); }
    inline bool isLast(node_p n) const
    { return (getSibling(n) == 0 ? true : false); }

private:
    // More accessors
    inline void setParent(node_p node, node_p p)
    {
        Tools::SetVariableField(data, psize, (ulong)node*bsize, p);
    }
    inline void setChild(node_p node, node_p c)
    {
        Tools::SetVariableField(data, psize, (ulong)node*bsize+psize, c);
    }
    inline void setSibling(node_p node, node_p s)
    {
        Tools::SetVariableField(data, psize, (ulong)node*bsize+2*psize, s);
    }
    inline void setEntry(node_p node, unsigned e)
    {
        Tools::SetVariableField(data, esize, (ulong)node*bsize+4*psize, e);
    }
    inline void setC(node_p node, uchar c)
    {
        Tools::SetVariableField(data, csize, (ulong)node*bsize+4*psize+esize, c);
    }

    void initNode(node_p node, node_p parent, uchar c, unsigned entry)
    {        
        setParent(node, parent);
        setC(node, c);
        setEntry(node, entry);
    }

    node_p parse(unsigned &i, std::string const &bp, node_p parent, unsigned depth, unsigned entry)
    {
        /**
         * Sanity checks
         */
        if (i+1 >= bp.size())
        {
            std::cerr << "NaiveForest::parse(): unable to parse tree structure, unexpected end of string at " 
                      << i << " of string:" << std::endl << bp << std::endl;
            std::abort();
        }
        if (bp[i] != '(')
        {
            std::cerr << "NaiveForest::parse(): unable to parse tree structure, expecting '(' at position " 
                      << i << " of string:" << std::endl << bp << std::endl;
            std::abort();
        }

        /**
         * Init node at freeBlock
         */
        unsigned node = freeBlock++;
        initNode(node, parent, bp[++i], entry);

        if (depth > h)
            h = depth;

        // Parse children if any
        if (bp[++i] == '(')
            setChild(node, parse(i, bp, node, depth + 1, entry));
        else
            ++l; // incr. number of leaves

        /**
         * Sanity checks
         */
        if (i >= bp.size())
        {
            std::cerr << "NaiveForest::parse(): unable to parse tree structure, unexpected end of string at " 
                      << i << " of string:" << std::endl << bp << std::endl;
            std::abort();
        }
        if (bp[i] != ')')
        {
            std::cerr << "NaiveForest::parse(): unable to parse tree structure, expecting ')' at position " 
                      << i << " of string:" << std::endl << bp << std::endl;
                std::abort();
        }

        // Parse siblings if any
        if (i+1 < bp.size() && bp[++i] == '(')
            setSibling(node, parse(i, bp, parent, depth, entry));
       
        return node;
    }

    void printNode(node_p n, unsigned d) const
    {
        for (unsigned i = 0; i < d; ++i)
            std::cerr << " ";
        std::cerr << getC(n) << std::endl;
        
        if (!isLeaf(n))
            printNode(getChild(n), d+1);
        if (!isLast(n))
            printNode(getSibling(n), d);
    }
        
public:

    NaiveForest(unsigned entries_, ulong trees_, ulong t_)
        : data(0), freeBlock(0), psize(0), csize(8), bsize(0), 
        entries(entries_), trees(trees_), t(t_), h(0), l(0)
    {
        /**
         * One node consists of 
         *
         *    Node *parent;     ~0lu if no parent
         *    Node *child;      Pointer to the leftmost child (0 if none)
         *    Node *sibling;    0 if no sibling
         *    Node *name;       Lexicographical rank
         *    Entry *entry;     Pointer to FASTA entry
         *    uchar c;
         *
         * we allocate ceillog(t) bits for each Node pointer,
         * ceillog(entries) bits for the Entry pointer
         * plus 8 bits for the character.
         */
        psize = Tools::CeilLog2(t > 255 ? t : 256); // Was (t > trees+255 ? t : trees+255);
        esize = Tools::CeilLog2(entries);
        bsize = 4 * psize + esize + csize;

        // Debug printing
        std::cerr << "NaiveForest() debug: Using bsize = " << bsize << ", esize = " << esize
                  << ", psize = " << psize << ", csize = " << csize << std::endl
                  << "Allocating " << (t * bsize/W + 1) * sizeof(ulong) << " bytes for tree data." << std::endl;

        // Reserve enough data blocks
        data = new ulong[t * bsize/W + 1];
        for (ulong i = 0; i < t * bsize/W + 1; ++i)
            data[i] = 0;
    }

    void add(std::string const &bp, unsigned entry)
    {
        unsigned i = 0; // Iterator over bp
        // t-1 is given as "parent" of the root
        parse(i, bp, t-1, 1, entry);
        if (i != bp.size() - 1)
        {
            std::cerr << "NaiveForest::parse(): unable to parse tree structure, expecting end of line at position " 
                      << i << " of string of size " << bp.size() << ":" << std::endl << bp << std::endl;
            std::abort();
        }
    }

    ~NaiveForest()
    {
        delete [] data;
        data = 0;
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
    unsigned esize;   // Entry pointer size, ceil(log(entries)) bits
    unsigned csize;   // Character size, log(alphabet size) bits
    unsigned bsize;   // Block size, 4 * psize + csize bits.

    unsigned entries; // Number of entries
    ulong trees;      // Number of trees
    ulong t;          // Number of nodes
    unsigned h;       // Height
    ulong l;          // Number of leaves
};

#endif
