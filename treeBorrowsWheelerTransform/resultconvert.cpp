// Convert result into the original alphabet
// FIXME Uses 32 bit variables
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <cassert>
#include <cstdlib> // abort()

using namespace std;

typedef unsigned char uchar;


unsigned str2unsigned(string const &str) 
{
    stringstream ss(str);
    unsigned n;
    ss >> n;
    return n;
}

vector<string> splitFreqs(string const &str) 
{
    // Separator is ' '
    stringstream ss(str);
    vector<string> tokens;
    while (ss.good ()) 
    {
        string token;
        ss >> token;
        tokens.push_back(token);
    }
    return tokens;
}

pair<unsigned,unsigned> splitFreq(string const &str)
{
    // Separator is ':'
    size_t i = str.find(':');
    if (i == string::npos)
    {
        cerr << "error: unable to tokenize the string \"" << str << "\"" << endl;
        abort();
    }

    unsigned d = str2unsigned(str.substr(0, i));
    unsigned o = str2unsigned(str.substr(i + 1));
    return make_pair(d,o);
}

int main (int argc, char **argv)
{
    if (argc != 3)
    {
        cerr << "usage: " << argv[0] << " <fasta> <mapping>  < result  > output" << endl
             << "where <fasta> is the original FASTA input, and " << endl
             << "      <mapping> is a file containing the std. error output from 'tconvert'." << endl;
        return 1;
    }

    /**
     * Parse mapping from log file
     */
    vector<string> otoi(256);    
    {
        ifstream mappingf(argv[2]);
        if (!mappingf.good())
        {
            cerr << argv[0] << ": unable to read mapping file " << argv[2] << endl;
            exit(1); 
        }    
        string row;
        unsigned mappings = 0;
        while (getline(mappingf, row).good()) 
        {
            /**
             * The mapping file should contain rows having this
             * kind of a format:
             *
             *       otoi[']'] = "{Mg";
             *
             * I.e. mapping from uchar to a string.
             */
            if (row.substr(0,4) == "otoi")
            {
                ++mappings;
                uchar c = row[6];
                assert (otoi[c].empty());
            
                string out = row.substr(row.find('\"'));
                assert (out.size() >= 5);
                out = out.substr(1, out.find('\"', 1) - 1);
                cerr << "otoi[" << c << "] = " << out << endl;
                otoi[c] = out;
            }
        }
        if (mappings == 0)
        {
            cerr << argv[0] << ": unable to parse the given mapping file!" << cerr;
            abort();
        }
        cerr << "Found " << mappings << " mappings." << endl
             << "Reading FASTA titles..." << endl;
    }
    /**
     * Parse FASTA titles
     */
    vector<string> title;
    {
        ifstream fastaf(argv[1]);
        if (!fastaf.good())
        {
            cerr << argv[0] << ": unable to read FASTA file " << argv[1] << endl;
            exit(1); 
        }
        string row;
        unsigned entries = 0;
        while (getline(fastaf, row).good()) 
        {
            if (row[0] == '>')
            {
                ++entries;
                string t = row.substr(row.find_first_not_of("> \t"));
                if (t.find(".seqs") == t.size() - 5)
                    t = t.substr(0, t.size() - 5); // Trim the suffix

                if (t.find(':') != string::npos)
                    cerr << "Warning: FASTA title \"" << t << "\" contains the delimiter character ':' !" << endl;
                title.push_back(t);
            }
        }
        if (entries == 0)
        {
            cerr << argv[0] << ": unable to parse the given FASTA file; no fasta title rows found!" << cerr;
            abort();
        }
        cerr << "Found " << entries << " FASTA entries." << endl
             << "Converting results..." << endl;
    }

    /**
     * Convert result to original alphabet
     */
    string row;
    unsigned rows = 0;
    while (getline(cin, row).good()) 
    {
        ++rows;
        if (rows % 100000 == 0)
            cerr << rows << " rows converted..." << endl;

        string pattern = row.substr(0, row.find(' '));
        string occs = row.substr(row.find(' ') + 1);

        // Convert pattern
        for (string::iterator it = pattern.begin(); it != pattern.end(); ++it)
            cout << otoi[*it];

        // Convert frequencies (insert fasta titles)
        vector<string> tokens = splitFreqs(occs);
        for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); ++it)
        {
            pair<unsigned,unsigned> token = splitFreq(*it);
            cout << ' ' << title[token.first] << ':' << token.second;
        }        
        cout << endl;
    }

    cerr << "Done after " << rows << endl;
    return 0;
}
