package mechanism.graphs;



public class PGEdge extends Edge
{
	protected PGNode source, target;
	
	public PGEdge(Graph parent, PGNode n1, PGNode n2)
	{
		this.parent = parent;
		this.source = n1;
		this.target = n2;
		
		source.addNeighbor(this, target);
		target.addNeighbor(this, source);
	}
	
	public PGNode getSource()
	{
		return source;
	}

	public PGNode getTarget()
	{
		return target;
	}	
	
	public PGNode getOther(PGNode n)
	{
		if (source == n)
			return target;
		if (target == n)
			return source;
		
		return null;
	}
}