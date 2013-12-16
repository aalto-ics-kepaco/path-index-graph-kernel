package mechanism.graphs;

public class RGKNode extends Node
{
	public MoleculeGraph molecule;
	
	public RGKNode(RGKGraph parent)
	{
		this.parent = parent;
		id = 0;
	}

	public int getMolIndex()
	{
		return molecule.getIndex();
	}
	
	public String getSymbol()
	{
		return molecule.getLigand();
	}

	public String toString()
	{
		return id + ":" + molecule.getLigand();
	}
}