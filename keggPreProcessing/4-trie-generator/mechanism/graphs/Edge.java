package mechanism.graphs;


public abstract class Edge
{
	protected Node source;
	protected Node target;
	protected int id;
	protected int type;
	protected int change;  // the change pattern, 0 (no change), +1 (new bond), -1 (cleaved bond)
	protected int oldtype; // the bond's old type, e.g. 2 for double bond
	protected int newtype; // the bond's new type, e.g. 0 for no bond (when cleavage happens)
	protected Graph parent;
	
	public int getId()
	{
		return id;
	}
	
	public Node getSource()
	{
		return source;
	}

	public Node getTarget()
	{
		return target;
	}

	public int getChangetype()
	{
		return change;
	}

	public int getOldtype()
	{
		return oldtype;
	}

	public int getNewtype()
	{
		return newtype;
	}
	
	public int getType()
	{
		return type;
	}
	
	public Graph getParent()
	{
		return parent;
	}
	
	public String toString()
	{
		return source.getId() + "<->" + target.getId() + "(" + change + ":" + oldtype + "->" + newtype + ")";
	}

	public Node getOther(Node n)
	{
		if (source == n)
			return target;
		if (target == n)
			return source;
		
		return null;
	}
}
