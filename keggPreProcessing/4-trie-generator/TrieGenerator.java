import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;

import mechanism.*;
import mechanism.graphs.Atom;
import mechanism.graphs.Graph;
import mechanism.graphs.MoleculeGraph;
import mechanism.graphs.Node;
import mechanism.graphs.ReactionGraph;


public class TrieGenerator
{
	/**
	 * Generates a Trie out of the paths of the reaction graphs
	 */
	public static void main(String[] args)
	{
		Graph[] graphs;
		
		// Read reaction file names
		List<String> files = new ArrayList<String>();
		for (int i = 1; i < args.length; i++)
			if (new File(args[i]).isFile())
				files.add(args[i]);
		
		int maxdepth = Integer.parseInt(args[0]);
		
		int gc = files.size();
		
		// always use sorted indices
		Collections.sort(files);
		
		System.out.println(gc + " reaction graphs");
		
		graphs = new MoleculeGraph[gc];
		for (int i = 0; i < gc; i++)
		{
			graphs[i] = new MoleculeGraph(files.get(i));
			graphs[i].setIndex(i);
		}

		TG tg = new TG(graphs, maxdepth);
		tg.Generate();
	}
}



class TG
{
	private BitSet currbits = null;
	private Trie T = null;
	private long sz = 0;
	private long oldsz = 0;
	private int maxdepth = 0;
	private int depth = 0;
	private List<List<Character>> seqs;
	
	private ArrayList<Character> seqstr;
	
//	private char[] seqstr;
//	private int seqptr;
	private Graph graphs[];
	
	public TG(Graph[] graphs2, int maxdepth)
	{
		this.graphs = graphs2;
		this.maxdepth = maxdepth;
	}
	
	public void Generate()
	{
		seqstr = new ArrayList<Character>();
		seqs = new ArrayList<List<Character>>();
		
		for (Graph g : graphs)
		{
			System.out.print("starting rg " + g + " ");

			T = new Trie();
			TrieNode root = T.getRoot();
			
			currbits = new BitSet(g.getSize());
			seqs.clear();
			sz = 0;
			
			for (int i = 0; i < g.getSize(); i++)
			{
				depth = 0;
				seqstr.clear();
				DFS(g.getNodes()[i], root, 0);
				
				seqs.add((ArrayList<Character>)seqstr.clone());
				
//				seqs[i] = Arrays.copyOf(seqstr, seqptr);
			}
			
			System.out.println(" size " + sz);
			
			String fn;
			if (g.getDirection() == 1)
				fn = g.getLigand() + "_f_" + g.getMapNum() + ".seqs";
			else if (g.getDirection() == -1)
				fn = g.getLigand() + "_b_" + g.getMapNum() + ".seqs";
			else
				fn = g.getLigand() + ".seqs";
			
			write(fn);
		}
	}
	
	private void write(String filename)
	{
		BufferedWriter out;
		
		try
		{
			out = new BufferedWriter(new FileWriter(filename));
			
			out.write(">" + filename + "\n");
			
			for (int i = 0; i < seqs.size(); i++)
			{
				for (int j = 0; j < seqs.get(i).size(); j++)
				{
					out.write(seqs.get(i).get(j));
				}
				
				out.newLine();
			}
			
			out.close();
		}
		catch (IOException e)
		{
			System.out.println(e.getMessage());
		}
	}
	
	
	public void DFS(Node v, TrieNode currtn, int type)
	{
		if (depth+1 > maxdepth)
			return;
		
		depth++;

		if (type == 0)
			seqstr.add('(');
		else if (type == -1)
			seqstr.add('{');
		else if (type == 1)
			seqstr.add('[');
		
//		currtn = currtn.AddChild(v);
//		T.addNode(currtn);
		currbits.set(v.getId(), true);

		char[] symb = v.getSymbol().toCharArray();
		
		for (int i = 0; i < symb.length; i++)
			seqstr.add(symb[i]);
		
		sz++;
		
//		System.out.println(" seq " + currbits);
		
		// go through neighbors, push to stack if not already in curr
		for (Node ne : v.getNodeNeighbors())
			if (currbits.get(ne.getId()) == false)
				DFS(ne, currtn, ne.getEdge(v).getChangetype());
		
		currbits.set(v.getId(), false);
		
		depth--;
		
		if (type == 0)
			seqstr.add(')');
		else if (type == -1)
			seqstr.add('}');
		else if (type == 1)
			seqstr.add(']');
	}
	
	
}





