// Convert into single symbol alphabet
#include <iostream>
#include <string>
#include <set>
#include <map>
#include <cassert>

using namespace std;

typedef unsigned char uchar;

int main (int argc, char **argv)
{
    if (argc != 1)
    {
        cerr << "usage: " << argv[0] << " < input > output" << endl;
        return 1;        
    }

    // Output alphabet
    set<uchar> freec;
    cerr << "output:";
    for (unsigned i = 42; i < 256; ++i)
    {
        cerr << " " << (uchar)i;
        freec.insert(i);
    }
    cerr << endl;

    // Mapping from input to output
    map<string, uchar> itoo;

    string row;
    unsigned entries = 0;

    while (getline(cin, row).good()) 
    {
        if (row[0] == '>')
        {
            ++entries;
            if (entries % 1000 == 0)
                cerr << "after " << entries << " (" << row << ")" << endl;
            cout << row << endl;
            continue;
        }
        
        while (!row.empty())
        {
            if (row.find_first_of("})]") == 0)
            {
                do
                {
                    row = row.substr(1);
                    cout << ")";
                } while (row.find_first_of("})]") == 0);
                continue;
            }

            string in = row.substr(0, row.find_first_of("{([])}", 1));
            if (itoo.count(in) == 0)
            {
                // pick symbol
                assert( !freec.empty() );
                set<uchar>::iterator it = freec.begin();
                uchar c = *it;
                freec.erase(it);
             cerr << "using itoo[\"" << in << "\"] = '" << c << "';" << endl;
                itoo[in] = c;
            }

            cout << "(" << itoo[in];

            row = row.substr(in.size());
        }
        cout << endl;
    }

    // Output mapping
    for (map<string, uchar>::iterator it = itoo.begin(); it != itoo.end(); ++it)
        cerr << "itoo[\"" << it->first << "\"] = '" << it->second << "';" << endl;

    for (map<string, uchar>::iterator it = itoo.begin(); it != itoo.end(); ++it)
        cerr << "otoi['" << it->second << "'] = \"" << it->first << "\";" << endl;

    cerr << "total number of entries " << entries << endl;
    return 0;
}
