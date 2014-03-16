
#include "TBWTIndex.h"

#include <sstream>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <map>
#include <set>
#include <deque>
#include <string>
#include <ctime>
#include <cstring>
#include <cassert>
#include <getopt.h>


using namespace std;

/**
 * Definitions for parsing command line options
 */
#define DEFAULT_MIN_LENGTH 1
enum parameter_t { long_opt_max = 256, long_opt_debug };

/**
 * Flags set based on command line parameters
 */
bool verbose = false;
bool debug = false;
unsigned minLength = DEFAULT_MIN_LENGTH,
    maxLength = ~0u;

// Variables used during the recursion
TBWTIndex * tbwti = 0; 
deque<uchar> path;
time_t wctime = 0;
ulong traversed = 0;
ulong totalOccs = 0;
typedef map<unsigned, unsigned> freq_map;

inline void sumFreqs(freq_map &to, freq_map const &from)
{
    for (freq_map::const_iterator it = from.begin(); it != from.end(); ++it)
        to[it->first] += it->second;
}

void output(freq_map const &freq)
{
    for (deque<uchar>::iterator it = path.begin(); it != path.end(); ++it)
        printf("%c", (char)*it);
    for (freq_map::const_iterator it = freq.begin(); it != freq.end(); ++it)
    {
        printf(" %u:%u", it->first, it->second);
        totalOccs += it->second;
    }
    printf("\n");
}

// FIXME std::map is too slow?
// FIXME Use simpler traverse when path.size is small
void traverseSubtree(TBWTIndex::node_r range, freq_map const &leafFreq = freq_map())
{
    if (debug)
    {
        cerr << "range = " << range.first << ", " << range.second << ": ";
        for (deque<uchar>::iterator it = path.begin(); it != path.end(); ++it)
            cerr << *it;
        cerr << endl;
    }

    // Number of traversed internal nodes
    traversed += tbwti->getSubtreeSize(range);
    if (verbose && traversed % 10000 == 0)
        cerr << "Traversed " << traversed << " internal nodes of " << (tbwti->numberOfNodes() - tbwti->numberOfLeaves() + 1) << " (" 
             << 100*traversed/(tbwti->numberOfNodes() - tbwti->numberOfLeaves()) << "%, "
              << "elapsed time " << std::difftime(time(NULL), wctime) << " s, " 
              << std::difftime(time(NULL), wctime) / 3600 << " hours)" << endl;
    
    set<uchar> C, leaves;
    set<uchar>::iterator Cit;

    // Traverse subtrees and collect frequencies
    tbwti->getSubtrees(range, C);
    tbwti->getLeaves(range, leaves);
    for (Cit = C.begin(); Cit != C.end(); ++Cit)
    {
        path.push_back(*Cit);
        set<uchar>::iterator leaf = leaves.find(*Cit);
        if (leaf != leaves.end())
        {
            // There is a leaf/leaves branching with the same symbol,
            // we have to include their frequencies in the subtree
            leaves.erase(leaf); // Remove since it is taken care of here!
            traverseSubtree(tbwti->getSubtree(range, *Cit), tbwti->getLeafFreq(range, *Cit));
        }
        else // No leaves branching with symbol *Cit
            traverseSubtree(tbwti->getSubtree(range, *Cit));
        path.pop_back();
    }

    // Process the remaining leaves (those that have no subtree branching with the same symbol)
    if (path.size() + 1 >= minLength && path.size() + 1 <= maxLength)
        for (Cit = leaves.begin(); Cit != leaves.end(); ++Cit)
        {
            path.push_back(*Cit);
            output(tbwti->getLeafFreq(range, *Cit));
            path.pop_back();
        }

    freq_map freq;
    if (range.first != 0) // FIXME Skip root earlier
        freq = tbwti->getInternalFreq(range);
    sumFreqs(freq, leafFreq);

    // Output current node
    if (path.size() >= minLength && path.size() <= maxLength)
        output(freq);
}

int atoi_min(char const *value, int min, char const *parameter, char const *name)
{
    std::istringstream iss(value);
    int i;
    char c;
    if (!(iss >> i) || iss.get(c))
    {
        cerr << name << ": argument of " << parameter << " must be of type <int>, and greater than or equal to " << min << endl
             << "Check README or `" << name << " --help' for more information." << endl;
        std::exit(1);
    }

    if (i < min)
    {
        cerr << name << ": argument of " << parameter << " must be greater than or equal to " << min << endl
             << "Check README or `" << name << " --help' for more information." << endl;
        std::exit(1);
    }
    return i;
}


void print_usage(char const *name)
{
    cerr << "usage: " << name << " [options] <index>" << endl
         << "Check README or `" << name << " --help' for more information." << endl;
}

void print_help(char const *name)
{
    cerr << "usage: " << name << " [options] <index>" << endl << endl
         << "The index file must be an index built with 'builder'." << endl
         << "Result is outputted to the standard output." << endl << endl
         << "Options:" << endl
         << " -h, --help                    Display command line options." << endl
         << " -v, --verbose                 Print progress information." << endl
         << " -m, --min                     Minimum path length in the" << endl
         << "                               result set. (default: no limit)" << endl
         << " --max                         Maximum path length in result."
         << "                               (default: no limit)" << endl;
/* not yet implemented, default is to search all paths starting from the root(s):
         << " --path                        Traverse only the given path." << endl
         << "                               (default: all paths)" << endl;*/
}

int main(int argc, char **argv) 
{
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
            {"min",         required_argument, 0, 'm'},
            {"max",         required_argument, 0, long_opt_max},
            {"debug",       no_argument,       0, long_opt_debug},
            {0, 0, 0, 0}
        };
    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "hvm:",
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
        case 'm':
            minLength = atoi_min(optarg, 1, "-m, --min", argv[0]);
            break;
        case long_opt_max:
            maxLength = atoi_min(optarg, 1, "--max", argv[0]);
            break;
        case long_opt_debug:
            debug = true; break;
        default:
            print_usage(argv[0]);
            std::abort();
        }
    }
 
    if (argc - optind < 1)
    {
        cerr << argv[0] << ": expecting one input file!" << endl;
        print_usage(argv[0]);
        return 1;
    }
         
    if (argc - optind > 1)
        cerr << "Warning: too many filenames given! Ignoring all but the first one." << endl;

    string indexfile = string(argv[optind++]);
    
    /**
     * Load the index
     */
    if (verbose) 
        cerr << "Loading the index..." << endl;
    tbwti = new TBWTIndex(indexfile);
    if (verbose)
        cerr << "Number of FASTA entries: " << tbwti->numberOfEntries() << endl 
             << "Total number of trees: " << tbwti->numberOfTrees() << endl
             << "Total number of nodes: " << tbwti->numberOfNodes() << endl
             << "Traversing the tree..." << endl;

    cerr << std::fixed;
    cerr.precision(2);
    wctime = time(NULL);

    /**
     * Traversing..
     *
     * Start from the root node of all trees:
     */
    TBWTIndex::node_r root = make_pair(0, tbwti->numberOfTrees() - 1);
    path.clear();
    traverseSubtree(root);

    /**
     * Clean up
     */
    delete tbwti;   

    if (verbose)
        cerr << "Traverse complete, " << traversed << " internal nodes checked, " 
             << "in total " << totalOccs << " occurrences found." << endl
              << "Total time " << std::difftime(time(NULL), wctime) << " s, " 
              << std::difftime(time(NULL), wctime) / 3600 << " hours." << endl;
    return 0;
}
 
