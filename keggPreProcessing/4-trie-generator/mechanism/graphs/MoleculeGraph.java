package mechanism.graphs;


import java.util.*;
import java.io.*;


public class MoleculeGraph extends Graph
{
	public MoleculeGraph(String filename)
	{
		formula = new HashMap<String,Integer>();
		ligand = filename.substring(filename.lastIndexOf("/")+1, filename.lastIndexOf("."));
		
		try
		{
			read(filename);
		} catch (IOException e)
		{
			System.out.println("error " + e.getMessage());
		}
	}
	
	private MoleculeGraph(Graph parent, BitSet nodebits)
	{
		super(parent, nodebits);
	}
	
	// id:ligand:size:formula
	public String toString()
	{
		StringBuffer sbf = new StringBuffer();
		sbf.append(id + ":" + this.ligand + ":" + this.getSize() + ":");
		for (String s : formula.keySet())
		{	if (formula.get(s) > 1)
				sbf.append(s + formula.get(s));
			else
				sbf.append(s);
		}		
		return sbf.toString();
	}

	public void read(String filename) throws IOException
	{
		List<String> lines = new ArrayList<String>();
		
		Scanner sc = new Scanner(new File(filename));
		while (sc.hasNextLine())
			lines.add(sc.nextLine());
		sc.close();

		
		int id = 0; // start molecule's internal numbering from zero

		int atomcount = Integer.parseInt(lines.get(3).substring(0, 3).trim());
		int bondcount = Integer.parseInt(lines.get(3).substring(3, 6).trim());
		
		Node[] tmpatoms = new Node[atomcount];
		Edge[] tmpbonds = new Edge[bondcount];
		int aid = 0;
		int bid = 0;
		
		// contains information whether atom is hydrogen or not
		int[] validatoms = new int[atomcount];
		
		// atom block
		for (int i = 4; i < 4 + atomcount; i++)
		{
			String line = lines.get(i);
			String[] words = line.split("\\s+");
			String symbol = words[4].trim();
			
			// don't take hydrogens
			if (!symbol.equals("H") && !symbol.equals("H+"))
			{
				Atom a = new Atom(this, id++, symbol);
				tmpatoms[aid++] = a;
				validatoms[i-4] = id-1;
				
				if (formula.containsKey(symbol))
					formula.put(symbol, formula.get(symbol) + 1);
				else
					formula.put(symbol, 1);
			}
			else
				// all hydrogens into rawlist
				validatoms[i-4] = -1;
		}
		
		// bond block, ignore hydrogen bonds
		for (int i = 4 + atomcount; i < 4 + atomcount + bondcount; i++)
		{
			String line = lines.get(i);
			int source = Integer.parseInt(line.substring(0, 3).trim())-1; // numbering correction
			int target = Integer.parseInt(line.substring(3, 6).trim())-1; // mol-files start from 1, we from 0
			int type   = Integer.parseInt(line.substring(6, 9).trim());
			
			// hydrogen bonds are not created
			if (validatoms[source] == -1 || validatoms[target] == -1)
				continue;
			
			Node src = tmpatoms[validatoms[source]];
			Node tgt = tmpatoms[validatoms[target]];
			
			Bond b = new Bond(bid, (Atom)src, (Atom)tgt, type, 0, type, type, this);
			tmpbonds[bid++] = b;

			src.addNeighbor(b, tgt);
			tgt.addNeighbor(b, src);
		}
		
		// create the actual atom and bond arrays
		nodes = new Node[aid];
		edges = new Edge[bid];
		
		// ...copy non-hydrogen stuff to them
		System.arraycopy(tmpatoms, 0, nodes, 0, aid);
		System.arraycopy(tmpbonds, 0, edges, 0, bid);
		
		// node id's are compressed (consecutive)
	}
	
	public void computeDistances()
	{
		for (Node n : nodes)
			n.setCoreDist(0);
	}
	
	public List<Graph> getMolecules()
	{
		List<Graph> x = new ArrayList<Graph>();
		x.add(this);
		return x;
	}

	public Graph createSubgraph(BitSet nodebits)
	{
		return new MoleculeGraph(this, nodebits);
	}
}	


