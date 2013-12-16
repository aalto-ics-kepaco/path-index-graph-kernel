
#include "NaiveForest.h"
#include "TBWTIndex.h"

#include <sstream>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <ctime>
#include <cstring>
#include <cassert>
#include <getopt.h>

using namespace std;

/**
 * Flags set based on command line parameters
 */
bool verbose = false;
bool debug = false;

void traverse(TBWTIndex * tbwti, NaiveForest *nf, TBWTIndex::node_p &nodei, NaiveForest::node_p noden, unsigned d)
{
    if (debug)
    {
        for (unsigned i = 0; i < d; ++i)
            cout << " ";
        cout << tbwti->getC(nodei) << endl;
    }
    assert(tbwti->isRoot(nodei) == nf->isRoot(noden));
    assert(tbwti->isLeaf(nodei) == nf->isLeaf(noden));
    assert(tbwti->getC(nodei) == nf->getC(noden));
    
    if (!tbwti->isLeaf(nodei))
    {
        // Process children
        TBWTIndex::node_r range = tbwti->getChildren(nodei);
        //cout << "i = " << nodei << ", chilren = " << range.first << ", " << range.second << endl;
        NaiveForest::node_p childn = nf->getChild(noden);
        for (ulong i = range.first; i <= range.second; ++i)
        {
            traverse(tbwti, nf, i, childn, d+1);
            childn = nf->getSibling(childn);
        }
        assert(childn == 0);
    }
}

void compare(TBWTIndex * tbwti, NaiveForest *nf)
{
    NaiveForest::node_p noden = 0;

    for (ulong j = 0; j < tbwti->numberOfTrees(); ++j)
    {
        if (verbose)
            cout << "----------------------------------------------------------------------------" << endl
                 << "Tree number " << j << endl;
        while (!nf->isRoot(noden))
            ++noden;
        TBWTIndex::node_p nodei = tbwti->getRoot(j);
        traverse(tbwti, nf, nodei, noden, 0);
        ++noden;
    }
}

void input_size(istream *input, unsigned &entries, ulong &trees, ulong &nodes)
{
    entries = 0;
    trees = 0;
    nodes = 0;
    string row;
    while (getline(*input, row).good()) 
    {
        if (row[0] == '>')
            ++entries;
        else
        {
            ++trees;
            if (row.size() % 3 != 0)
            {
                cerr << "error: invalid input row length at row:" << endl
                     << row << endl;
                abort();
            }
            nodes += row.size()/3;
        }
    }
}

void parse_entries(istream *input, NaiveForest *nf, ulong estimatedLength, time_t wctime)
{
    unsigned entry = 0;
    ulong j = 0;
    string row;
    string name = "undef";
    while (getline(*input, row).good()) 
    {
        j += row.size();

        if (row[0] == '>')
        {
            name = row.substr(1);
            ++entry;
            if (verbose && entry % 1000 == 0 )
            {
                cerr << "Inserting entry n:o " << entry << ", name: " << name << " (";
                if (estimatedLength)
                    cerr << (100*j/estimatedLength) << "%, ";
                cerr << "elapsed " << std::difftime(time(NULL), wctime) << " s, " 
                     << std::difftime(time(NULL), wctime) / 3600 << " hours)" << endl;
            }
        }
        else
            nf->add(row, entry-1);
    }
    row.clear();
}

void print_usage(char const *name)
{
    cerr << "usage: " << name << " [options] <fasta> <index>" << endl
         << "Check README or `" << name << " --help' for more information." << endl;
}

void print_help(char const *name)
{
    cerr << "usage: " << name << " [options] <fasta> <index>" << endl << endl
         << "Give two filenames as input. "
         << "The first file must be in FASTA format." << endl
         << "The second file must be index built from the given FASTA." << endl
         << "The index is compared against the original tree to check validity." << endl << endl
         << "Options:" << endl
         << " -h, --help                    Display command line options." << endl
         << " -v, --verbose                 Print progress information." << endl;
}

int main(int argc, char **argv) 
{
    // Sanity check
#ifdef NDEBUG
    cerr << "error: compile test code without the -DNDEBUG flag!" << endl;
    abort();
#endif

    /**
     * Parse command line parameters
     */
    if (argc == 1)
    {
        print_usage(argv[0]);
        return 1;        
    }

    static struct option long_options[] =
        {
            {"help",        no_argument,       0, 'h'},
            {"verbose",     no_argument,       0, 'v'},
            {0, 0, 0, 0}
        };
    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "hv",
                            long_options, &option_index)) != -1) 
    {
        switch(c) 
        {
        case 'h':
            print_help(argv[0]);
            return 0;
        case 'v':
            verbose = true; break;
        case '?':
            print_usage(argv[0]);
            return 1;
        default:
            print_usage(argv[0]);
            std::abort();
        }
    }
 
    if (argc - optind < 2)
    {
        cerr << argv[0] << ": expecting two input files!" << endl;
        print_usage(argv[0]);
        return 1;
    }
        
    if (argc - optind > 2)
        cerr << "Warning: too many filenames given! Ignoring all but first two." << endl;

    string inputfile = string(argv[optind++]);
    string indexfile = string(argv[optind++]);

    istream *fp;
    /*if (inputfile == "-") 
        fp = &std::cin; // Not supported here!
    else*/
    fp = new std::ifstream(inputfile.c_str());

    if (!fp->good())
    {
        cerr << "test: unable to read input file " << inputfile << endl;
        exit(1); 
    }

    cerr << std::fixed;
    cerr.precision(2);
    time_t wctime = time(NULL);

    // estimate the total input sequence length
    fp->seekg(0, ios::end);
    long estLength = fp->tellg();
    if (estLength == -1)
    {
        cerr << "error: unable to seek input file" << endl;
        estLength = 0;
    }
    fp->seekg(0);
    if (verbose)
        cerr << "Input file size is " << estLength/1024 << " kb." << endl;

    /**
     * Count the input sizes
     */
    unsigned entries = 0;
    ulong trees = 0, nodes = 0;
    input_size(fp, entries, trees, nodes);
    fp->clear(); // forget previous EOF
    fp->seekg(0);

    if (verbose)
        cerr << "Number of FASTA entries: " << entries << endl 
             << "Total number of trees: " << trees << endl
             << "Total number of nodes: " << nodes << endl
             << "Parsing the input..." << endl;
    
    /**
     * Parse all trees
     */
    NaiveForest *nf = new NaiveForest(entries, trees, nodes);
    parse_entries(fp, nf, estLength, wctime);
    if (debug)
        nf->debugPrint();

    ulong t = nf->numberOfNodes();
    if (verbose)
        cerr << "Number of nodes: " << t << endl
             << "Number of leaves: " << nf->numberOfLeaves() << endl
             << "Tree height: " << nf->getHeight() << endl; 

    /**
     * Load the index
     */
    if (verbose) cerr << "------------------------------------------------------------------------" << endl 
                      << "Loading the index..." << endl;
    TBWTIndex * tbwti = new TBWTIndex(indexfile);


    /**
     * Comparing..
     */
    if (verbose) cerr << "------------------------------------------------------------------------" << endl 
                      << "Comparing..." << endl;
    assert(trees == tbwti->numberOfTrees());
    assert(entries == tbwti->numberOfEntries());
    assert(t == tbwti->numberOfNodes());
    assert(nf->numberOfLeaves() == tbwti->numberOfLeaves());

    compare(tbwti, nf);

    /**
     * Clean up
     */
    delete tbwti;
    tbwti = 0;
    delete nf;
    nf = 0;
    
    if (fp != &std::cin)
        delete fp;
    fp = 0;

    std::cerr << "Test complete. "
              << "(total wall-clock time " << std::difftime(time(NULL), wctime) << " s, " 
              << std::difftime(time(NULL), wctime) / 3600 << " hours)" << endl;
    return 0;
}
 
