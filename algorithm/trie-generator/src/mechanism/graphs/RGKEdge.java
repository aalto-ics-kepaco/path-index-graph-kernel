package mechanism.graphs;

public class RGKEdge extends Edge
{
	public String type;
	
	public RGKEdge(RGKGraph parent, String type)
	{
		this.parent = parent;
		this.type = type;
		this.source = null; // no neighbors
		this.target = null;
	}

	public RGKEdge(RGKGraph parent, String type, RGKNode n1, RGKNode n2)
	{
		this.parent = parent;
		this.type = type;
		this.source = n1;
		this.target = n2;
		
		source.addNeighbor(this, target);
		target.addNeighbor(this, source);
	}

	public String toString()
	{
		return id + ":" + type + ":" + source + "<->" + target;
	}

	public int getType()
	{
		if (type.equals("main"))
			return 1;
		if (type.equals("leave"))
			return 2;
		if (type.equals("cofac"))
			return 3;
		if (type.equals("trans"))
			return 4;
		if (type.equals("ligase"))
			return 5;
		
		return 0;
	}
}