
#include "SimpleForest.h"
#include "TBWTBuilder.h"
#include "BitRank.h"
#include "HuffWT.h"

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
 * Definitions for parsing command line options
 */
enum parameter_t { long_opt_rlcsa = 256, long_opt_debug };

/**
 * Flags set based on command line parameters
 */
bool verbose = false;
bool debug = false; 


int atoi_min(char const *value, int min, char const *parameter, char const *name)
{
    std::istringstream iss(value);
    int i;
    char c;
    if (!(iss >> i) || iss.get(c))
    {
        cerr << "builder: argument of " << parameter << " must be of type <int>, and greater than or equal to " << min << endl
             << "Check README or `" << name << " --help' for more information." << endl;
        std::exit(1);
    }

    if (i < min)
    {
        cerr << "builder: argument of " << parameter << " must be greater than or equal to " << min << endl
             << "Check README or `" << name << " --help' for more information." << endl;
        std::exit(1);
    }
    return i;
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

void parse_entries(istream *input, SimpleForest *nf, ulong *entry, ulong estimatedLength, time_t wctime)
{
    unsigned entries = 0;
    unsigned node = 0;
    ulong j = 0;
    string row;
    string name = "undef";
    while (getline(*input, row).good()) 
    {
        j += row.size();

        if (row[0] == '>')
        {
            // Mark first tree in this entry
            Tools::SetField(entry, 1, node, 1);

            name = row.substr(1);
            ++entries;
            if (verbose && entries % 100 == 0 )
            {
                cerr << "Inserting entry n:o " << entries << ", name: " << name << " (";
                if (estimatedLength)
                    cerr << (100*j/estimatedLength) << "%, ";
                cerr << "elapsed " << std::difftime(time(NULL), wctime) << " s, " 
                     << std::difftime(time(NULL), wctime) / 3600 << " hours)" << endl;
            }
        }
        else
        {
            node += row.size() / 3;
            nf->add(row);
        }
    }
    row.clear();
}

void print_usage(char const *name)
{
    cerr << "usage: " << name << " [options] <input> [output]" << endl
         << "Check README or `" << name << " --help' for more information." << endl;
}

void print_help(char const *name)
{
    cerr << "usage: " << name << " [options] <input> [output]" << endl << endl
         << "<input> is the input filename. "
         << "The input must be in FASTA format." << endl
         << "If no output filename is given, the index is stored as <input>.fmi" << endl
         << "or as <input>.rlcsa." << endl << endl
         << "Options:" << endl
         << " -h, --help                    Display command line options." << endl
         << " -v, --verbose                 Print progress information." << endl;
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
            {"debug",       no_argument,       0, long_opt_debug},
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
        case long_opt_debug:
            debug = true; break;
        case '?':
            print_usage(argv[0]);
            return 1;
        default:
            print_usage(argv[0]);
            std::abort();
        }
    }
 
    if (argc - optind < 1)
    {
        cerr << "builder: no input filename given!" << endl;
        print_usage(argv[0]);
        return 1;
    }
        
    if (argc - optind > 2)
        cerr << "Warning: too many filenames given! Ignoring all but first two." << endl;

    string inputfile = string(argv[optind++]);
    string outputfile = inputfile; 
    if (optind != argc)
        outputfile = string(argv[optind++]);
    
    outputfile += ".tbwt";

    istream *fp;
    /*if (inputfile == "-") 
        fp = &std::cin; // Not supported here!
    else*/
    fp = new std::ifstream(inputfile.c_str());

    if (!fp->good())
    {
        cerr << "builder: unable to read input file " << inputfile << endl;
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
    // Bitvector that marks the first root node in each entry
    ulong *entry = new ulong[nodes/W + 1];
    for (unsigned i = 0; i < nodes/W + 1; ++i)
        entry[i] = 0;

    SimpleForest *nf = new SimpleForest(entries, trees, nodes);

    parse_entries(fp, nf, entry, estLength, wctime);
    if (debug)
        nf->debugPrint();

    ulong t = nf->numberOfNodes();
    if (verbose)
        cerr << "Number of nodes: " << t << endl
             << "Number of leaves: " << nf->numberOfLeaves() << endl
             << "Tree height: " << nf->getHeight() << endl; 

    /**
     * Build TBWT
     */
    if (verbose) 
        cerr << "Sorting nodes..." << endl;
    TBWTBuilder tbwtb(nf);
    tbwtb.sort(verbose, debug);
    if (verbose) 
        cerr << "Sorting done! Saving TBWT to disk..." << endl;
    
    /**
     * Save to disk
     */
    {
        FILE *output = fopen(outputfile.c_str(), "wb");
        if (!output)
            throw std::runtime_error("builder: could not write file");

        // TBWT_SAVEFILE_MSG is defined in Tools.h
        if (fwrite(TBWT_SAVEFILE_MSG, sizeof(char), strlen(TBWT_SAVEFILE_MSG), output) 
            != strlen(TBWT_SAVEFILE_MSG))
            throw std::runtime_error("builder: file write error (msg).");

        if (fwrite(&entries, sizeof(unsigned), 1, output) != 1)
            throw std::runtime_error("builder: file write error (entries).");
        if (fwrite(&trees, sizeof(ulong), 1, output) != 1)
            throw std::runtime_error("builder: file write error (trees).");
        if (fwrite(&t, sizeof(ulong), 1, output) != 1)
            throw std::runtime_error("builder: file write error (t).");

        BitRank *b = new BitRank(tbwtb.getLeaf(), t, true);
        if (nf->numberOfLeaves() != b->rank(t-1))
        {
            cerr << "TBWTBuilder::getTBWT(): assert failed: number of leaves was " << nf->numberOfLeaves() << " but rank was " << b->rank(t-1) <<  endl;
            abort();
        }
        b->save(output);
        delete b; b = 0;
        
        b = new BitRank(tbwtb.getLast(), t, true);
        if (t - nf->numberOfLeaves() != b->rank(t-1) - 1)
        {
            cerr << "TBWTBuilder::getTBWT(): assert failed: number of internal nodes was " << t-nf->numberOfLeaves() << " but rank was " << b->rank(t-1) <<  endl;
            abort();
        }
        b->save(output);
        delete b; b = 0;
        
        uchar *tbwt = tbwtb.getTBWTInternal();
        HuffWT *hwt = HuffWT::makeHuffWT(tbwt, t-nf->numberOfLeaves());
        HuffWT::save(hwt, output);
        HuffWT::deleteHuffWT(hwt); 
        delete [] tbwt; tbwt = 0;
        
        tbwt = tbwtb.getTBWTLeaf();
        hwt = HuffWT::makeHuffWT(tbwt, nf->numberOfLeaves());
        HuffWT::save(hwt, output);
        HuffWT::deleteHuffWT(hwt); 
        delete [] tbwt; tbwt = 0;

        unsigned C[256];
        unsigned F[256];
        tbwtb.getCF(C, F);
        /**
         * C is not required
         * if (fwrite(C, sizeof(unsigned), 256, output) != 256)
         *     throw std::runtime_error("builder: file write error (C).");
         */

        if (fwrite(F, sizeof(unsigned), 256, output) != 256)
            throw std::runtime_error("builder: file write error (F).");

        if (debug)
            for (unsigned i = 0; i < 256; ++i)
                if (C[i])
                    cerr << "C[" << i << "] = " << C[i] << ", F[" << i << "] = " << F[i] << endl;

        nf->setFirstEntry(entry);

        BlockArray *leafEntry = tbwtb.getLeafEntry();
        leafEntry->Save(output);
        delete leafEntry;
        leafEntry = 0;

        BlockArray *lastEntry = tbwtb.getLastEntry();
        lastEntry->Save(output);
        delete lastEntry;
        lastEntry = 0;

        fflush(output);
        fclose(output);
    }


    /**
     * Clean up
     */
    // delete [] entry; already taken care of by BitRank in SimpleForest class.
    delete nf;
    nf = 0;
    
    if (fp != &std::cin)
        delete fp;
    fp = 0;

    if (verbose) 
        std::cerr << "Save complete. "
                  << "(total wall-clock time " << std::difftime(time(NULL), wctime) << " s, " 
                  << std::difftime(time(NULL), wctime) / 3600 << " hours)" << endl;
    return 0;
}
 
