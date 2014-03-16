#ifndef _NAIVETREE_H_
#define _NAIVETREE_H_

#include "Tools.h"

#include <iostream>
#include <string>

class NaiveTree
{
public:
    class Node
    {
    public:
        Node(Node *p, Node *chld, Node *sblng, uchar chr)
            : parent(p), child(chld), sibling(sblng), c(chr)
        { }
        ~Node()
        {
            delete child; child = 0;
            delete sibling; sibling = 0;
            parent = 0; // Does not delete parent
        }

        inline Node * getParent() const
        { return parent; }
        inline Node * getChild() const
        { return child; }
        inline void setChild(Node *chld)
        { child = chld; }
        inline Node * getSibling() const
        { return sibling; }
        inline void setSibling(Node *s)
        { sibling = s; }
        inline uchar getC() const
        { return c; }
        inline bool isRoot() const
        { return (parent == 0 ? 1 : 0); }
        inline bool isLeaf() const
        { return (child == 0 ? 1 : 0); }
        inline bool isLast() const
        { return (sibling == 0 ? 1 : 0); }

        void printNode(unsigned d)
        {
            for (unsigned i = 0; i < d; ++i)
                std::cerr << " ";
            std::cerr << c << std::endl;
            
            if (!isLeaf())
                child->printNode(d+1);
            if (!isLast())
                sibling->printNode(d);
        }
        

        inline void setName(unsigned n)
        { name = n; }
        unsigned getName() const
        { return name; }

    private:
        Node();

        Node *parent;
        Node *child; // leftmost child
        Node *sibling;
        uchar c;

        unsigned name;
    };

private:
    Node * parse(ulong &i, std::string const &bp, Node * parent, unsigned depth)
    {
        if (i+1 >= bp.size())
        {
            std::cerr << "NaiveTree::parse(): unable to parse tree structure, unexpected end of string at " 
                      << i << " of string:" << std::endl << bp << std::endl;
            std::abort();
        }

        if (bp[i] != '(')
        {
            std::cerr << "NaiveTree::parse(): unable to parse tree structure, expecting '(' at position " 
                      << i << " of string:" << std::endl << bp << std::endl;
            std::abort();
        }

        Node *node = new Node(parent, 0, 0, bp[++i]);
        ++t;
        if (depth > h)
            h = depth;

        if (bp[++i] == '(')
            node->setChild(parse(i, bp, node, depth + 1));
        else
            ++l; // incr. number of leaves

        if (i >= bp.size())
        {
            std::cerr << "NaiveTree::parse(): unable to parse tree structure, unexpected end of string at " 
                      << i << " of string:" << std::endl << bp << std::endl;
            std::abort();
        }

        if (bp[i] != ')')
        {
            std::cerr << "NaiveTree::parse(): unable to parse tree structure, expecting ')' at position " 
                      << i << " of string:" << std::endl << bp << std::endl;
                std::abort();
        }

        if (i+1 < bp.size() && bp[++i] == '(')
            node->setSibling(parse(i, bp, parent, depth));
       
        return node;
    }

public:

    NaiveTree(std::string const &bp)
        : root(0), t(0), h(0), l(0)
    {
        ulong i = 0;
        root = parse(i, bp, 0, 1);
        if (i != bp.size() - 1)
        {
            std::cerr << "NaiveTree::parse(): unable to parse tree structure, expecting end of line at position " 
                      << i << " of string of size " << bp.size() << ":" << std::endl << bp << std::endl;
            std::abort();
        }
    }

    ~NaiveTree()
    {
        delete root;
        root = 0;
    }

    Node * getRoot() const
    { return root; }

    ulong numberOfNodes() const
    { return t; }
    ulong numberOfLeaves() const
    { return l; }
    inline unsigned getHeight() const
    { return h; }

    void debugPrint() const
    {
        root->printNode(0);
    }
private:
    

    Node *root;
    ulong t; // Number of nodes
    unsigned h; // Height
    ulong l; // Number of leaves
};

#endif
