// c++ code for bron-kerbosch

#include <vector>
#include <string>
#include <fstream>
#include <ios>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <map>
#include <exception>
#include <list>
#include <set>
#include <utility>
#include <stack>


using std::istream;
using std::ostream;
using std::ifstream;
using std::fstream;
using std::string;
using std::ios;
using std::cout;
using std::endl;
using std::exception;
using std::map;
using std::vector;
using std::list;
using std::set;
using std::pair;
using std::make_pair;
using std::stack;


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



// for cands and x a hashtable would be best
// hold an array where a fixed slot for each node
// another array where for each node (slot), there's next/prev nodes
/*

*/

// class for hashtable with fixed set of possible elements
// contains three arrays, example:
//
// bits = [ 0, 1, 0, 0, 1, 0, 1]
// ptrs = [-1, 2,-1,-1, 1,-1, 0]
// list = [ 6, 4, 1,-1,-1,-1,-1]
//
// bits = [ 0, 1, 0, 0, 0, 0, 1]
// ptrs = [-1, 1,-1,-1,-1,-1, 0]
// list = [ 6, 1,-1,-1,-1,-1,-1]
//
// bits is a bitvector for O(1) search
// ptrs is a middle array for indicating where corresponding elements are on "list"
// list is a basic array of elements
// 'size' tells the position of first new elemenent
class bithashtable
{
	int size;
	int maxsize;
	vector<bool> bits;
	int[] ptrs;
	int[] list;
	
	bithashtable(int size)
	{
		this->maxsize = size;
		this->size = 0;
		bits = vector<bool>(size);
		ptrs = int[size];
		list = int[size];
		
		for (int i=0;i<size;i++)
		{
			bits[i] = 0;
			ptrs[i] = -1;
			list[i] = -1;
		}
	}
	
	void fill()
	{
		for (int i = 0; i<maxsize; i++)
		{
			bits[i] = 1;
			ptrs[i] = i;
			list[i] = i;
		}
		
		size = maxsize;
	}
	
	void add(int id)
	{
		if (bits[id] == 0)
		{
			bits[id] = 1;
			list[size] = id;
			ptrs[id] = size;
			size++;
		}
	}
	
	void remove(int id)
	{
		if (bits[id] == 1)
		{
			bits[id] = 0;
			
			// swap removable with last
			list[ptrs[id]] = list[size-1];
			list[size-1] = -1;
			ptrs[list[ptrs[id]]] = ptrs[id]]
			ptrs[id] = -1;
			
			size--;
		}
	}
};



// global
int[] cliq;
int cliqptr;
//bithashtable cands;
//vector<bool> x;
//int xsize;

// graph
int[] vertices;
int vertexcount;
int[][] N; // sorted list of neighs
int[][] Nbits; // bitvector of neighs
int[] Nsizes;

// on each round for cands:
//  1x length()
//  1x iteration over
//  2x remove a set
//  1x remove a single elem
// on each round for x
//  1x length()
//  1x remove a set
//  1x remove a single elem

// assume neighs are sorted

void bk(vector<int> & cands, vector<int> & x)
{
	if (cands.empty() && x.empty())
		cout << "clique! size:" << cliqptr-1 << endl; // clique is cliq[0...clipqtr]
	
	int u = choose_pivot(cands, x);
	int[] uNeighs = N[u];
	int uNeighsize = Nsizes[u];
	
	while (!cands.empty())
	{
		int v = -1;
		// get first candidate not in uNeighs
		for (int i=0; i<cands.size(); i++)
		{
			if (Nbits[u][cands[i]] == 0)
			{
				v = cands[i];
				cands.erase(i);
				break;
			}
		}
		
		if (v == -1)
			break;
		
		int[] vNeighs = N[v];
		int vNeighsize = Nsizes[v];
		
		// getting ready to call bk(cliq+v, cands .intersect. N(v), x .intersect. N(v))
		// clique
		cliq[cliqptr++] = v;
		
		// newcands contains everything form cands that's also in neighs
		vector<int> newcands;
		for (int i=0; i<cands.size(); i++)
			if (Nbits[u][cands[i]] == 1)
				newcands.push_back(cands[i]);
		
		// x
		vector<int> newx;
		for (int i=0; i<x.size(); i++)
			if (Nbits[u][x[i]] == 1)
				newx.push_back(x[i]);
		
		// calling bk(newcands, newx)
		bk(newcands, newx);
		
		// returning
		cliqptr--;
		
		// add to exclusion
		x.push_back(v);
	}
}


void bronkerbosch(int nodecount, vector<int> nodes, vector< vector<int> > neighbors)
{
	vertexcount = nodecount;
	vertices = int[vertexcount];
	
	N = vertex< int[] >;
	
	for (int i=0; i<vertexcount; i++)
	{
		vertices[i] = nodes[i];
		N[i] = int[ neighbors[i].size() ];
		
		for 
	}
	vertices = vertices;
	N = neighbors;
	
	Nbits
	
	
	// initialize graph data structures
	nodes = nodes;
	nodesize = nodesize;
	neighs = neighs;
	neighsizes = neighsizes;
	
	// initializing bk data structures
	cands = bithashtable(nodesize);
	cands.fill()
	
	x = vector<bool>();
	cliq = int[nodesize];
	
	// initialize cands to contain all nodes, x and cliq nothing
	for (int i=0;i < nodesize; i++)
	{
		x[i] = 0;
		cliq[i] = 0;
	}
	
	xsize = 0;
	cliqptr = 0;
	
	bk();
}

// reads the graph and places result in global variables
void read_graph(string fn)
{
	// FILEFORMAT:
	// 
	// N1 label
	// N2 label
	// ...
	// N1 N2 label
	// ...
	
	
	fstream f(fn.c_str(), ios::read);
	vector<string> lines;
	
	string temp;
	while (std::getline(f, temp).eof() == false)
	{
		lines.push_back( string(temp.c_str() );
	}
	
	for (int i=0; i<lines.size(); i++)
	{
		// 2 words in the row
		if (lines[i].find("_") == lines[i].rfind("_"))
		{
			nodes[atoi(lines[i].substr(0,lines[i].find("_")))]
		}
		
	}
	
	
}



int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		cout << "arg error" << endl;
		return EXIT_FAILURE;
	}
	
	string fn = string(argv[1]);
	read_graph(fn);
	bron_kerbosch();
	
	return EXIT_SUCCESS;
}



