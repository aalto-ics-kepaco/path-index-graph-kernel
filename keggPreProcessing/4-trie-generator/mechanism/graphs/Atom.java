package mechanism.graphs;

import java.util.*;


public class Atom extends Node implements Cloneable
{
	public Atom(Graph parent, int id, String c)
	{
		this.parent = parent;
		this.id = id;
		this.symbol = c;
	}

	public Node clone()
	{
		Node x = new Atom(parent, id, symbol);
		x.coredist = coredist;
		x.nodeneighs = new HashSet<Node>();
		x.edgeneighs = new HashMap<Node,Edge>();
		return x;
	}
	
	// 4:P
	public String toString()
	{
		return id + ":" + symbol;
	}
}

class AtomSymbolComparator implements Comparator<Atom>
{
	public int compare(Atom x, Atom y)
	{
		return x.getSymbol().compareTo(y.getSymbol());
	}
}


