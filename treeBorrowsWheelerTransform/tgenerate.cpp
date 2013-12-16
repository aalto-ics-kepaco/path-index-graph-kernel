// Generates "random" trees into stdout
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>

using namespace std;


void output_node(int &t)
{
    if (t <= 0)
        return;
    --t;

    cout << "("
         << (char)(rand() % (73-65) + 65); //(char)(rand() % (90-65) + 65);

    if ((rand() % 100) < 50)
        output_node(t);  // Generate children

    cout << ")";

    if ((rand() % 100) < 50)
        output_node(t); // Generate siblings
}

int main (int argc, char **argv)
{
    if (argc != 3)
    {
        cerr << "usage: " << argv[0] << " [size] [entries]" << endl;
        return 1;        
    }
    int t = atoi(argv[1]); // Size
    int entries = atoi(argv[2]);

    srand ( time(NULL) );

    for (int i = 0; i < entries; ++i)
    {
        // Fasta header
        cout << ">" << i << endl;
  
        // Generate t trees of size t^2
        for (int j = 0; j < t; ++j)
        {
            int tt = t*t-1; //t-1;
            // Root
            cout << "("
                 << (char)(rand() % (90-65) + 65);
            
            while (tt > 0)
                output_node(tt);

            // root ends
            cout << ")" << endl;
        }
    }
    return 0;
}
