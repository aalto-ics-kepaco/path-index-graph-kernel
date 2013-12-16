package mechanism.graphs;

import java.util.*;


public class PGNode extends Node
{
	public Node a1, a2; // pair
	protected Set<PGNode> nodeneighs = new HashSet<PGNode>(); // overwrite superclasses containers, more specific here
	protected Map<PGNode,PGEdge> edgeneighs = new HashMap<PGNode,PGEdge>();
	
	public PGNode(Graph parent, Node a1, Node a2)
	{
		this.parent = parent;
		this.a1 = a1;
		this.a2 = a2;
	}

	public String toString()
	{
		return a1 + "-" + a2;
	}

	public Node getNode1()
	{
		return a1;
	}
	
	public Node getNode2()
	{
		return a2;
	}
	
	public Set<PGNode> getNodeNeighbors()
	{
		return nodeneighs;
	}

	public Collection<PGEdge> getEdgeNeighbors()
	{
		return edgeneighs.values();
	}
	
	public boolean isNeighbor(PGNode other)
	{
		return getEdge(other) != null;
	}
	
	public PGEdge getEdge(PGNode other)
	{
		return edgeneighs.get(other);
	}

	public void addNeighbor(PGEdge e, PGNode n)
	{
		nodeneighs.add(n);
		edgeneighs.put(n,e);
	}
}
