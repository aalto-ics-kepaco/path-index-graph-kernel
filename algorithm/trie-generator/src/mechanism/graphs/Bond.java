package mechanism.graphs;


public class Bond extends Edge implements Cloneable
{
	public Bond(int id, Atom source, Atom target, int type, int change, int oldtype, int newtype, Graph parent)
	{
		this.id = id;
		this.source = source;
		this.target = target;
		this.change = change;
		this.type = type;
		this.oldtype = oldtype;
		this.newtype = newtype;
		this.parent = parent;
		
		source.addNeighbor(this, target);
		target.addNeighbor(this, source);
	}
	
	public Edge clone()
	{
		return new Bond(id, (Atom)source, (Atom)target, type, change, oldtype, newtype, parent);
	}
	
	// 4:+1:4P-2O
	public String toString()
	{
		return id + ":" + source + "<->" + target;
	}
}






