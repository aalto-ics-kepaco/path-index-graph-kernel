package mechanism.kernels;

import mechanism.*;
import mechanism.graphs.Graph;
import mechanism.graphs.Node;
import mechanism.graphs.MoleculeGraph;
import mechanism.graphs.ReactionGraph;

import java.text.DecimalFormat;
import java.util.*;

public class SubgraphsKernel extends Kernel
{
//	private List<List<BitSet>> sgbitsets = null;
	
	private Subgraphs current_sg = null;
	private Map<Graph, Integer> subs = null;
	
	private int[] counts1;
	private int[] counts2;
	
	public SubgraphsKernel(Graph[] graphs, KernelParams params)
	{
		super(graphs, params);
	}

	public void compute()
	{
		// optimally we would compute subgraphbits for all graphs
		// and then do N^2 kernel computation using that
		// resulting in a single sweep
		//
		// that takes too much memory, tho
		// instead
		//
		// start...end  (i.e. 1000..1030)
		// do:
		//  compute start...end subgraphbits into memory  (1000...1030)
		//  
		//  for j = 1...end:               #    (1...1030)
		//   for i = start...end:          # (1000...1030)
		//    kernel(i,j) = compute(i,j)   # using 'i's in memory
		//
		// now this results in a single sweep over graphs
		//
		
//		// compute bitsets for our range
//		sgbitsets = new ArrayList<List<BitSet>>(params.end - params.start);
//		compute_bitsets();
		
		// compute only the lower left triangle
		for (int i = params.start-1; i < params.end; i++)
		{
			pgtime = 0;
			walktime = 0;
			
			// compute i only once
			current_sg = new Subgraphs(graphs[i]);
			subs = current_sg.enumerateDistinct(params.maxlen);
			
			for (int j = 0; j <= i; j++)
			{
				matrix[i-params.start+1][j] = (float)compute(graphs[i], graphs[j]);
			}
			
			totalpgtime += pgtime;
			totalwalktime += walktime;
			
			String dir = "";
			if (graphs[i].getDirection() == 1)
				dir = "_+1";
			else if (graphs[i].getDirection() == -1)
				dir = "_-1";
			
			if (graphs[i] instanceof ReactionGraph)
				System.out.println(graphs[i].getLigand() + "_" + ((ReactionGraph)graphs[i]).getMapNum() + dir + ": row " + (i+1) + " computed");
			else if (graphs[i] instanceof MoleculeGraph)
				System.out.println(graphs[i].getLigand() + ": row " + (i+1) + " computed");
		}
	}

	
	public double compute(Graph g1, Graph g2)
	{
		// two ways to compute:
		// (1) product graph and enumerate subgraphs:
		//     ends up with all permutations of a single subgraph
		// (2) enumerate subgraphs separately, check for common ones
		//     needs isomorphism
		
		Subgraphs sg2 = new Subgraphs(g2);
		Map<Graph,Integer> res2 = sg2.enumerateDistinct(params.maxlen);
		
		if (params.op == KernelOperationType.MinNormalized)
		{
			// compute the number of walks of different lenghts to use in normalization
			counts1 = new int[params.maxlen+1];
			counts2 = new int[params.maxlen+1];
			Arrays.fill(counts1, 0);
			Arrays.fill(counts2, 0);
			
			for (Graph g : subs.keySet())
				counts1[g.getSize()] += subs.get(g);
			for (Graph g : res2.keySet())
				counts2[g.getSize()] += res2.get(g);
		}
		
		
		double kvalue = 0.0;
		for (Graph g : subs.keySet())
		{
			if (res2.containsKey(g))
			{
				int val1 = subs.get(g);
				int val2 = res2.get(g);
//				int sz = g.getSize();
				
				// min comparison
				if (params.op == KernelOperationType.Min)
					kvalue += Math.min(val1, val2); 
				// min comparison + normalization
				else if (params.op == KernelOperationType.MinNormalized)
					kvalue += Math.min(val1, val2) / (Math.sqrt(counts1[g.getSize()]) * Math.sqrt(counts2[g.getSize()]));
				// indicator comparison
				else if (params.op == KernelOperationType.Indicator)
					kvalue += 1;
				// dot product
				else if (params.op == KernelOperationType.DotProduct)
					kvalue += val1 * val2; 
			}
		}
		
		return kvalue;
	}
}

class Subgraphs
{
	private Graph g;
	private BitSet fragAtoms;
	private List<BitSet> fragments;
	private Map<Graph,Integer> fragdict;
	private int[] forbidden;
	private int atomptr; // Indeksi viimeisimp�n� lis�ttyyn kaareen
	private Node[] lastAtom; // Taulukko fragmenttiin lis�tyille solmuille
	private int nodes;
	private int limit;
	private int resultcount;
	private boolean countOnly;
	
	public Subgraphs(Graph g)
	{
		this.g = g;
		nodes = g.getNodes().length;
		
		fragments = new ArrayList<BitSet>();
		fragdict = null;
		fragAtoms = new BitSet(nodes);
		forbidden = new int[nodes];
		limit = Integer.MAX_VALUE;
		countOnly = false;
	}
	
	// count just the number of subgraphs
	public int countSubgraphs(int limit)
	{
		this.limit = limit;
		return countSubgraphs();
	}
	
	// count just the number of subgraphs
	public int countSubgraphs()
	{
		countOnly = true;
		if (limit > 0)
			DFS();
		return resultcount;
	}
	
	// enumerate all node-subgraphs
	public List<BitSet> enumerateSubgraphs(int limit)
	{
		this.limit = limit;
		return enumerateSubgraphs();
	}
	
	public List<BitSet> enumerateSubgraphs()
	{
		if (limit > 0)
			DFS();
		
		return fragments;
	}
	
	// enumerate the isomorphic subgraphs with counts
	public Map<Graph,Integer> enumerateDistinct(int limit)
	{
		this.limit = limit;
		return enumerateDistinct();
	}

	// enumerate the isomorphic subgraphs with counts
	public Map<Graph,Integer> enumerateDistinct()
	{
		if (limit > 0)
		{
			DFS();
			separateDistinct();
		}
		
		return fragdict;
	}
	
	private void DFS()
	{
		atomptr = -1; // Solmujen lukumäärä
		lastAtom = new Node[nodes];

		// for each node, we have forbidden-status
		Arrays.fill(forbidden, 0);
		fragments.clear();
		resultcount = 0;

		// start enumerating subgraphs from each node at a time
		// after each round, prevous start-node should be in forbidden
		for (Node a : g.getNodes())
		{
			Add(a);
			while (Forward() || Backward()) 
			{ 
				;
			}
		}
	}

	private void separateDistinct()
	{
		// construct molecular graphs out of the BitSets, test for isomorphism on the fly
		// and aggregate resulting molecular graphs into a dict with counts
		
		// Tree data structure
		// root -> size1,size2,...
		// sizek -> formula1,formula2,...
		// formulak -> abdist1,abdist2,...
		// 
		
		fragdict = new HashMap<Graph,Integer>(fragments.size()/2);
		for (int i = 0; i < fragments.size(); i++)
		{
			// take fragment bits
			BitSet fbits = fragments.get(i);
			
			// construct moleculegraph out of it
			Graph f = g.createSubgraph(fbits);

			// test its isomorphism against all previous fragments
			if (fragdict.containsKey(f))
				fragdict.put(f, fragdict.get(f)+1);
			else
				fragdict.put(f, 1);
		}
	}
	
	private boolean Forward()
	{
		if (limit <= atomptr+1)
			return false;
		
		// loop through lastNode[nodeptr]
		int i = atomptr;
		while (i >= 0)
		{
			Node last = lastAtom[i];
			for (Node ne : last.getNodeNeighbors())
			{
				if (fragAtoms.get(ne.getId()) == false && forbidden[ne.getId()] == 0) // forbidden.defined(current[an]) == false)
				{
					Add(ne);
					return true;
				}
			}
			
			i--;
		}

		return false;
	}

	private void Add(Node a)
	{
		// add node/atom
		fragAtoms.set(a.getId());
		lastAtom[++atomptr] = a;
		if (!countOnly)
			fragments.add((BitSet)fragAtoms.clone());
		resultcount++;
	}

	private boolean Backward()
	{
		if (atomptr == -1) // Alkioita nodeptr+1 eli nolla
			return false;

		forbidden[lastAtom[atomptr].getId()] = atomptr + 1;
//		forbidden_size++;

		Remove(lastAtom[atomptr]);

		// Liian isot forbidit karsitaan
		for (int i = 0; i < nodes; i++)
		{
			if (forbidden[i] > atomptr + 2)
			{
				forbidden[i] = 0;
//				forbidden_size--;
			}
		}

		return true;
	}

	private void Remove(Node a)
	{
		fragAtoms.set(a.getId(), false);
		atomptr--;
	}
}
