package mechanism.graphs;

import java.io.File;
import java.io.IOException;
import java.util.*;

/*
 * ReactionGraph, read from mol-files
 * 
 */
public class ReactionGraph extends Graph
{
	private int mapnum;
	private List<String> ec_codes = new ArrayList<String>();
	private List<String> reactant_ligands = new ArrayList<String>();
	private List<String> product_ligands = new ArrayList<String>();
	private List<MoleculeGraph> reactants = new ArrayList<MoleculeGraph>();
	private List<MoleculeGraph> products = new ArrayList<MoleculeGraph>();
	

	private ReactionGraph(Graph parent, BitSet nodebits)
	{
		super(parent,nodebits);
	}
		
	public ReactionGraph(String filename)
	{
		try
		{
			read(filename);
		} catch (IOException e)
		{
			System.out.println("error " + e.getMessage());
		}
	}
	
	
	public List<String> getReactantLigands()
	{
		return reactant_ligands;
	}

	public List<String> getProductLigands()
	{
		return product_ligands;
	}
	
	public List<MoleculeGraph> getReactants()
	{
		return reactants;
	}
	public List<MoleculeGraph> getProducts()
	{
		return products;
	}
	
	public void addReactant(MoleculeGraph mg)
	{
		reactants.add(mg);
	}
	public void addProduct(MoleculeGraph mg)
	{
		products.add(mg);
	}
	
	public int getMapNum()
	{
		return mapnum;
	}
	
	// id:ligand:size
	public String toString()
	{
		return id + ":" + ligand + ":" + nodes.length + "/" + edges.length;
	}
	

	public void read(String filename) throws IOException
	{
		String[] fn_splitted = filename.split("/");		
		if (fn_splitted[fn_splitted.length-1].contains("_f"))
			direction = 1;
		else if (fn_splitted[fn_splitted.length-1].contains("_b"))
			direction = -1;
		
		String[] parts = fn_splitted[fn_splitted.length-1].split("_");
		if (parts.length == 3)
			mapnum = Integer.parseInt(parts[1]);
		
		List<String> lines = new ArrayList<String>();
		
		Scanner sc = new Scanner(new File(filename));
		while (sc.hasNextLine())
			lines.add(sc.nextLine());
		sc.close();
		
		// first line contains ligand and ec-codes
		if (!lines.get(0).trim().isEmpty())
		{
			ligand = lines.get(0).substring(0,6).trim();
			String[] ecwords = lines.get(0).split("[ \t]");
			for (int i = 1; i < ecwords.length; i++)
				ec_codes.add(ecwords[i]);
		}
		
		// second line contains the reaction equation as C1 C2 => C3 C4
		if (!lines.get(1).trim().isEmpty())
		{
			String[] substs = lines.get(1).trim().split("[ ]");
			int j;
			for (j = 0; j < substs.length; j++)
			{
				if (!substs[j].contains("=>"))
					reactant_ligands.add(substs[j]);
				else
					break;
			}
			for (j = j + 1; j < substs.length; j++)
				product_ligands.add(substs[j]);
		}		
		
		int aid = 0;
		int bid = 0;
		int atomcount = Integer.parseInt(lines.get(3).substring(0, 3).trim());
		int bondcount = Integer.parseInt(lines.get(3).substring(3, 6).trim());
		
		this.formula = new HashMap<String,Integer>();

		
		nodes = new Node[atomcount];
		edges = new Edge[bondcount];
		
		for (int i = 4; i < 4 + atomcount; i++)
		{
			String line = lines.get(i);
			String[] words = line.split("\\s+");
			String symbol = words[4].trim();
			
			Node a = new Atom(this, aid, symbol);
			
			nodes[aid++] = a;
			
			if (formula.containsKey(symbol))
				formula.put(symbol, formula.get(symbol)+1);
			else
				formula.put(symbol, 1);
		}
		
		for (int i = 4 + atomcount; i < 4 + atomcount + bondcount; i++)
		{
			String line = lines.get(i);
			int source = Integer.parseInt(line.substring(0, 3).trim());
			int target = Integer.parseInt(line.substring(3, 6).trim());
			int event = Integer.parseInt(line.substring(6, 9).trim());
			int oldtype = Integer.parseInt(line.substring(9, 12).trim());
			int newtype = Integer.parseInt(line.substring(12, 15).trim());
			
			Node src = nodes[source-1];
			Node tgt = nodes[target-1];
			
			Edge b = new Bond(bid, (Atom)src, (Atom)tgt, 0, event, oldtype, newtype, this);
			edges[bid++] = b;
			
			src.addNeighbor(b, tgt);
			tgt.addNeighbor(b, src);
		}
		
		atomcount = nodes.length;
		bondcount = edges.length;
		
		// compress nodes so that id's are consecutive
		
		
//		// Sort atoms by their symbols into order
//		Arrays.sort(atoms, new AtomSymbolComparator());
		
//		String dir = direction == 1 ? "+1" : "-1";
	}

	public Graph createSubgraph(BitSet atombits)
	{
		return new ReactionGraph(this, atombits);
	}

	public List<MoleculeGraph> getSubstrates()
	{
		List<MoleculeGraph> substrates = new ArrayList<MoleculeGraph>(reactants);
		substrates.addAll(products);
		return substrates;
	}
}






