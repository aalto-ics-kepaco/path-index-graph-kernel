/**
 * Testing random queries
 *
 * subpath count
 * subpath freq
 * subpath subtree traverse
 */
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
typedef TBWTIndex::node_p node_p;

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
// FIXME Use simpler traverse when path.size is smaller than minLength
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

// Freq of root-originating path
// Note: assumes that path exists
void subpath_freq(string const &path, freq_map &result)
{
    // Start from root
    TBWTIndex::node_r range = make_pair(0, tbwti->numberOfTrees() - 1);
    for (unsigned i = 0; i < path.size() - 1; ++i)
        range = tbwti->getSubtree(range, path[i]);
    
    result = tbwti->getLeafFreq(range, path[path.size()-1]);
    range = tbwti->getSubtree(range, path[path.size()-1]);
    if (range.first <= range.second)
        sumFreqs(result, tbwti->getInternalFreq(range));
}

// Count of root-originating path
// Note: assumes that path exists
unsigned subpath_count(string const &path)
{
    // Start from root
    TBWTIndex::node_r range = make_pair(0, tbwti->numberOfTrees() - 1);
    for (unsigned i = 0; i < path.size() - 1; ++i)
        range = tbwti->getSubtree(range, path[i]);
    
    unsigned lc = tbwti->getLeafCount(range, path[path.size()-1]);
    range = tbwti->getSubtree(range, path[path.size()-1]);
    return tbwti->getInternalCount(range) + lc;
}

// Returns the number of leaves in a subtree
unsigned subpath_subtree_leaves(TBWTIndex::node_r range)
{
    set<uchar> C;
    set<uchar>::iterator Cit;

    // Traverse subtrees and collect frequencies
    tbwti->getSubtrees(range, C);
    unsigned leaves = 0;
    for (Cit = C.begin(); Cit != C.end(); ++Cit)
        leaves += subpath_subtree_leaves(tbwti->getSubtree(range, *Cit));

    return leaves + tbwti->getLeafCount(range);
}


// Count leaves of root-originating path's subtree
// Note: assumes that the given path exists
unsigned subpath_subtree(string const &path)
{
    // Start from root
    TBWTIndex::node_r range = make_pair(0, tbwti->numberOfTrees() - 1);
    for (unsigned i = 0; i < path.size() - 1; ++i)
        range = tbwti->getSubtree(range, path[i]);
    
    unsigned lc = tbwti->getLeafCount(range, path[path.size()-1]);
    range = tbwti->getSubtree(range, path[path.size()-1]);
    if (range.first <= range.second)
        return subpath_subtree_leaves(range) + lc;
    else
        return lc;
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
    cerr << "usage: " << name << " [options] <index> <result>" << endl
         << "Check README or `" << name << " --help' for more information." << endl;
}

void print_help(char const *name)
{
    cerr << "usage: " << name << " [options] <index> <result>" << endl << endl
         << "The index file must be an index built with 'builder'." << endl
         << "Result file must be generated by 'traverse'." << endl << endl
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
    unsigned nqueries = ~0u;
    unsigned patlen = 0;

    while ((c = getopt_long(argc, argv, "hvm:n:k:",
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
        case 'n':
            nqueries = atoi_min(optarg, 1, "-n, --nqueries", argv[0]);
            break;
        case 'k':
            patlen =  atoi_min(optarg, 1, "-k, --patlen", argv[0]);
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
 
    if (argc - optind < 2)
    {
        cerr << argv[0] << ": expecting two input files!" << endl;
        print_usage(argv[0]);
        return 1;
    }
         
    if (argc - optind > 2)
        cerr << "Warning: too many filenames given! Ignoring all but the first two." << endl;

    string indexfile = string(argv[optind++]);
    string resultfile = string(argv[optind++]);

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
     * Testing...
     */
    istream *fp;
    if (resultfile == "-") 
        fp = &std::cin;
    else
        fp = new std::ifstream(resultfile.c_str());

    if (!fp->good())
    {
        cerr << argv[0] << ": unable to read result file " << resultfile << endl;
        exit(1); 
    }
    if (verbose) cerr << "------------------------------------------------------------------------" << endl 
                      << "Testing queries..." << endl;

    ulong j = 0;
    string row;
    string pattern = "undef";
/*    while (getline(*fp, row).good() && j < nqueries) 
    {
        if (verbose && j && j % 100000 == 0)
            cerr << j << " queries..." << endl;
            
        // Parse row
        pattern = row.substr(0, row.find(' '));
        unsigned count = subpath_count(pattern);
        ++j;
        totalOccs += count;
    }
    cerr << "Test complete (subpath_count), " << j << " queries done, " 
             << "in total " << totalOccs << " occurrences found." << endl
              << "Total time " << std::difftime(time(NULL), wctime) << " s, " 
              << std::difftime(time(NULL), wctime) / 3600 << " hours." << endl
         << "queries per s: " << (double)j/std::difftime(time(NULL), wctime) << endl;

    fp->clear(); // forget previous EOF
    fp->seekg(0);
    wctime = time(NULL);
    j = 0;
    totalOccs = 0;
    while (getline(*fp, row).good() && j < nqueries) 
    {
        if (verbose && j && j % 100000 == 0)
            cerr << j << " results checked... OK" << endl;
            
        // Parse row
        pattern = row.substr(0, row.find(' '));
        freq_map result;
        subpath_freq(pattern, result);
        unsigned tc = 0;
        for (freq_map::const_iterator it = result.begin(); it != result.end(); ++it)
            tc += it->second;

        ++j;
        totalOccs += tc;
    }

    cerr << "Test complete (subpath_freq), " << j << " queries done, " 
             << "in total " << totalOccs << " occurrences found." << endl
              << "Total time " << std::difftime(time(NULL), wctime) << " s, " 
              << std::difftime(time(NULL), wctime) / 3600 << " hours." << endl
         << "queries per s: " << (double)j/std::difftime(time(NULL), wctime) << endl;
  
*/
    cerr << "subtree queries for k = " << patlen << endl;
    fp->clear(); // forget previous EOF
    fp->seekg(0);
    wctime = time(NULL);
    j = 0;
    totalOccs = 0;
    while (getline(*fp, row).good() && j < nqueries) 
    {
        // Parse row
        pattern = row.substr(0, row.find(' '));
        if (pattern.size() != patlen) continue;

        if (verbose && j && j % 10000 == 0)
            cerr << j << " queries... (subtree)" << endl;

        unsigned leaves = subpath_subtree(pattern);

        ++j;
        totalOccs += leaves;
    }

    cerr << "Test complete (subpath_subtree_leaves), " << j << " queries done, " 
             << "in total " << totalOccs << " occurrences found." << endl
              << "Total time " << std::difftime(time(NULL), wctime) << " s, " 
              << std::difftime(time(NULL), wctime) / 3600 << " hours." << endl
         << "queries per s: " << (double)j/std::difftime(time(NULL), wctime) << endl;


    /**
     * Clean up
     */
    delete tbwti;   
    return 0;
}
 
