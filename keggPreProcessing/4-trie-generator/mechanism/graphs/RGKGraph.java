package mechanism.graphs;

import java.io.*;
import java.util.*;

import mechanism.Kernel;


// represents the tsuda's reaction graph,
// ie a standard undirected graph where nodes are 'MoleculeGraph's and edges Strings (types)

public class RGKGraph extends Graph
{
	private RGKNode[] subnodes, prodnodes, nodes;
	private RGKEdge[] edges;
	private List<RGKEdge> edgelist;
	private List<String> reactant_ligands;
	private List<String> product_ligands;
	
	
	public RGKGraph(String filename)
	{
		reactant_ligands = new ArrayList<String>();
		product_ligands = new ArrayList<String>();
		edgelist = new ArrayList<RGKEdge>();
		
		// construct the graph
		try
		{
			read_reaction(filename);
			read_rpair();
		} catch (IOException e)
		{
			System.out.println(e.getMessage());
		}
		
		edges = new RGKEdge[edgelist.size()];
		int i = 0;
		for (RGKEdge e : edgelist)
			edges[i++] = e;
	}
	
	
	private void read_reaction(String filename) throws IOException
	{
		List<String> lines = new ArrayList<String>();
		
		Scanner sc = new Scanner(new File(filename));
		while (sc.hasNextLine())
			lines.add(sc.nextLine());
		sc.close();
		
		// first line contains ligand and ec-codes
		ligand = lines.get(0).substring(0,6).trim();
		
		// second line contains the reaction equation as C1 C2 => C3 C4
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
		
		nodes = new RGKNode[reactant_ligands.size() + product_ligands.size()];
		subnodes = new RGKNode[reactant_ligands.size()];
		prodnodes = new RGKNode[product_ligands.size()];
		int i = 0;
		
		// construct nodes
		for (String s : reactant_ligands)
		{
			RGKNode n = new RGKNode(this);
			n.molecule = new MoleculeGraph(Kernel.MOL_FOLDER + s + ".mol");
			subnodes[i] = n;
			nodes[i++] = n;
		}
		
		i = 0;
		for (String s : product_ligands)
		{
			RGKNode n = new RGKNode(this);
			n.molecule = new MoleculeGraph(Kernel.MOL_FOLDER + s + ".mol");
			prodnodes[i] = n;
			nodes[subnodes.length + i++] = n;
		}
	}
	
	
	private void read_rpair() throws IOException
	{
		List<List<String>> rpairs = new ArrayList<List<String>>();
		
		boolean region = false;
		Scanner sc = new Scanner(new File("/group/home/icomic/data/kegg/ligand-2010.07.01-results/rpairdefs.txt"));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			
			if (line.startsWith(ligand))
			{
				region = true;
				String[] words = line.split(" ");
				
				String sub = words[1];
				String prod = words[2];
				String type = words[3];
				
				List<String> rpair = new ArrayList<String>();
				rpair.add(sub);
				rpair.add(prod);
				rpair.add(type);
				
				rpairs.add(rpair);
			}
			else if (region) // break when first non-ligand line is read with region TRUE
				break;
		}
		sc.close();
		
		// go through all node pairs
		// -> add node if they are on both sides
		//    or if they can be found from rpair
		
		for (RGKNode n1 : nodes)
		{
			for (RGKNode n2 : nodes)
			{
				for (List<String> rpair : rpairs)
				{
					// if exists in rpair
					if ((rpair.get(0).equals(n1.molecule.getLigand()) && rpair.get(1).equals(n2.molecule.getLigand())) || (rpair.get(1).equals(n1.molecule.getLigand()) && rpair.get(0).equals(n2.molecule.getLigand())))
					{
						if (n1.isNeighbor(n2))
							continue;
						
						RGKEdge e = new RGKEdge(this, rpair.get(2), n1, n2);
						edgelist.add(e);
						break;
					}
				}
			}
		}
		
		// clique both sides unless already a edge
		for (RGKNode n1 : subnodes)
			for (RGKNode n2 : subnodes)
				if (!n1.isNeighbor(n2) && n1 != n2)
					edgelist.add(new RGKEdge(this, "group", n1, n2));
		
		for (RGKNode n1 : prodnodes)
			for (RGKNode n2 : prodnodes)
				if (!n1.isNeighbor(n2) && n1 != n2)
					edgelist.add(new RGKEdge(this, "group", n1, n2));
	}
	
	// mask superclass's getters to specify these 
	public RGKNode[] getNodes()
	{
		return nodes;
	}

	public RGKEdge[] getEdges()
	{
		return edges;
	}

	public int getNodeCount()
	{
		return nodes.length;
	}

	public int getEdgeCount()
	{
		return edges.length;
	}
	
	public int getSize()
	{
		return nodes.length;
	}
	
	public void computeDistances()
	{
		// set coredists as zero
		for (Node n : nodes)
			n.setCoreDist(0);
	}

	public Graph createSubgraph(BitSet nodebits)
	{
		return null; // do nothing
	}
	
	public String toString()
	{
		return id + ":" + ligand;
	}
}




