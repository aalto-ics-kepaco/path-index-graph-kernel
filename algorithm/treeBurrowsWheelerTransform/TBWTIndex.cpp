#include "TBWTIndex.h"

#include <stdexcept>
#include <cstring>
using namespace std;

TBWTIndex::TBWTIndex(std::string &inputfile)
    : entries(0), trees(0), t(0), leaf(0), last(0), wt(0), wtleaf(0), leafEntry(0), lastEntry(0)
{
    FILE *input = fopen((inputfile+".tbwt").c_str(), "rb");
    if (!input)
    {
        input = fopen(inputfile.c_str(), "rb");
        if (!input)
            throw std::runtime_error("TBWTIndex::TBWTIndex(): could not read the index file");
    }

    // TBWT_SAVEFILE_MSG is defined in Tools.h
    char msg[10];
    if (10 <= strlen(TBWT_SAVEFILE_MSG))
        throw std::runtime_error("TBWTIndex::TBWTIndex(): assert failed, message is longer than temporary buffer!");

    if (fread(msg, sizeof(char), strlen(TBWT_SAVEFILE_MSG), input) 
        != strlen(TBWT_SAVEFILE_MSG))
        throw std::runtime_error("file read error (msg).");
    msg[strlen(TBWT_SAVEFILE_MSG)] = 0;
    if (string(msg) != string(TBWT_SAVEFILE_MSG))
        throw std::runtime_error("TBWTIndex::TBWTIndex(): detected incorrect index file or old index version; please, check the filename given and/or rebuild the index!");

    if (fread(&entries, sizeof(unsigned), 1, input) != 1)
        throw std::runtime_error("file read error (entries).");
    if (fread(&trees, sizeof(ulong), 1, input) != 1)
        throw std::runtime_error("file read error (trees).");
    if (fread(&t, sizeof(ulong), 1, input) != 1)
        throw std::runtime_error("file read error (t).");

    leaf = new BitRank(input);
    last = new BitRank(input);
    wt = HuffWT::load(input);
    wtleaf = HuffWT::load(input);

    if (fread(F, sizeof(unsigned), 256, input) != 256)
        throw std::runtime_error("file read error (F).");

    leafEntry = new BlockArray(input);
    lastEntry = new BlockArray(input);
    fclose(input);
}

TBWTIndex::~TBWTIndex()
{
    delete leaf; leaf = 0;
    delete last; last = 0;
    HuffWT::deleteHuffWT(wt); wt = 0;
    HuffWT::deleteHuffWT(wtleaf); wtleaf = 0;
}

TBWTIndex::node_r TBWTIndex::getChildren(node_p node) const
{
    if (isLeaf(node))
        return make_pair((node_p)1, (node_p)0);
    ulong r = 0;
    node = leaf->rank0(node) - 1; // Map to internal nodes
    uchar c = wt->access(node, r); // Returns rank value in &r
    ulong y = F[c];
    ulong z = last->rank(y-1);
    node_r range;
    range.first = last->select(z+r-1)+1;
    range.second = last->select(z + r);
    return range;
}
