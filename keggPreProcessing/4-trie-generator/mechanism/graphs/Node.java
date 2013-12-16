package mechanism.graphs;

import java.util.*;


public abstract class Node
{
	protected Graph parent = null;
	protected String symbol = "";
	protected int id = -1;
	protected int coredist = Integer.MAX_VALUE;
	protected Set<Node> nodeneighs = new HashSet<Node>();
	protected Map<Node,Edge> edgeneighs= new HashMap<Node,Edge>();
	

	public int getId()
	{
		return id;
	}

	public String getSymbol()
	{
		return this.symbol;
	}

	public int getCoreDist()
	{
		return coredist;
	}
	
	public void setCoreDist(int value)
	{
		coredist = value;
	}
	
	public Collection<? extends Node> getNodeNeighbors()
	{
		return nodeneighs;
	}
	
	public Collection<? extends Edge> getEdgeNeighbors()
	{
		return edgeneighs.values();
	}
	
	public boolean isNeighbor(Node other)
	{
		return getEdge(other) != null;
	}
	
	public Edge getEdge(Node other)
	{
		Edge val = edgeneighs.get(other);
		boolean res = edgeneighs.containsKey(other);
		return edgeneighs.get(other);
	}

	public void addNeighbor(Edge e, Node n)
	{
		nodeneighs.add(n);
		edgeneighs.put(n, e);
	}
	
	public Graph getParent()
	{
		return parent;
	}

	public String toString()
	{
		return symbol + "(" + id + "/" + parent.getLigand() + ")";
	}
	
	// requires id's to be unique!
	public boolean equals(Node at)
	{
		return id == at.getId() && this.getParent() == at.getParent();
	}

	public int hashCode()
	{
		return id;
	}
}

class NodeSymbolComparator implements Comparator<Node>
{
	public int compare(Node x, Node y)
	{
		return x.getSymbol().compareTo(y.getSymbol());
	}
}