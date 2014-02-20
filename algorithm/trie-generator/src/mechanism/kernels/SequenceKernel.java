package mechanism.kernels;

import java.util.*;

import mechanism.*;
import mechanism.graphs.*;

/*
 *  Superclass of sequence kernels (i.e. walks/paths)
 */
public abstract class SequenceKernel extends Kernel
{
	public SequenceKernel(Graph[] graphs, KernelParams params)
	{
		super(graphs, params);
		
		if (params.kw == KernelWeight.Diffusion)
		{
			diffs = new HashMap<Graph,Diffusion>(); // from 'Kernel' superclass
		}
	}

	public double compute(Graph g1, Graph g2)
	{
		// if diffusion model used, precompute them
		if (params.kw == KernelWeight.Diffusion)
		{
			// compute diffusion matrices
			diffs.put(g1, new Diffusion(g1, Math.abs(params.beta))); // from 'Kernel' superclass
			diffs.put(g2, new Diffusion(g2, Math.abs(params.beta)));
		}
		
		// construct product graph
		long pgtime = System.currentTimeMillis();
		
		ProductGraph pg;
		if (params.reduced)
			pg = new ReducedProductGraph(g1, g2, params);
		else
			pg = new ProductGraph(g1, g2, params);
		
		this.pgtime += System.currentTimeMillis() - pgtime;

		// use dynamic programming
		long walktime = System.currentTimeMillis();
		
		double value;
		if (params.paths)
			value = paths(pg);
		else if (params.nontottering)
			value = nontottering_walks(pg);
		else
			value = walks(pg);
		
		this.walktime += System.currentTimeMillis() - walktime;
		
		counted++;
		
		return value;
	}
	
	// standard weighting function gives singular weight
	protected double weight(Node[] atoms)
	{
		return 1.0;
	}
	
	abstract protected double walks(ProductGraph pg);
	abstract protected double nontottering_walks(ProductGraph pg);
	
	
	/*
	 * Double Rucker's algorithm
	 * 
	 * INPUT: (G,G')
	 * OUTPUT: sum of products of path counts in G and G'
	 *  
	 * R = TreeMap<String,int>             # results
	 * C = []                              # current path
	 * F = set()                           # forbidden nodes
	 * last = [] # table
	 * for (s,s') pairs in G \times G':    # start from all pairs 
	 *   while |C| > 0: 
	 *   
	 *     for (v,v') \in border(C) - F - C:
	 *       add pair to C
	 *       increment R[C.half]
	 *     else:
	 *       C.remove(last)
	 *       add last pair to F with value len(C)
	 *       remove all forbids that are larger than last forbid
	 * 
	 * all references to pairs can be precomputed into a pair table and index suffices
	 * 
	 */
	protected double paths(ProductGraph pg)
	{
		// doesn't use dynamic programming, instead we use rucker to compute the 
		// paths into a table and sum it to get dot product
		Rucker r = new Rucker(pg.getG1(), pg.getG2(), params.reduced);
		r.DFS();
		
		return r.kvalue;
	}
	
	
	
	
	
	
	
	
	
	class Rucker
	{
		class AtomPair
		{
			public Node a1, a2;
			public AtomPair(Node a1, Node a2)
			{
				this.a1=a1;
				this.a2=a2;
			}
		}
		
		private Graph g1, g2;
//		private boolean countOnly;
		private int limit;
		private AtomPair[] atompairs;
		private AtomPair[] lastPair;
		private BitSet seq1atoms;
		private BitSet seq2atoms;
		private int pairptr;
		private int resultcount;
		private int[] forbidden1;
		private int[] forbidden2;
		private int nc1, nc2;
		private double kvalue = 0;
		private boolean reduced;
		
		public Rucker(Graph g1, Graph g2, boolean reduced)
		{
			this.g1 = g1;
			this.g2 = g2;
			nc1 = g1.getNodes().length;
			nc2 = g2.getNodes().length;
			this.reduced = reduced;
		}
			
		public void DFS()
		{
			// generate pairs, all against all (also duplicates)
			int i = 0;
			for (Node a1 : g1.getNodes())
				for (Node a2 : g2.getNodes())
					if (a1.getSymbol().equals(a2.getSymbol()) && reducedvalid(a1,a2))
						i++;
			atompairs = new AtomPair[i];
			
			i = 0;
			for (Node a1 : g1.getNodes())
				for (Node a2 : g2.getNodes())
					if (a1.getSymbol().equals(a2.getSymbol()) && reducedvalid(a1,a2))
						atompairs[i++] = new AtomPair((Atom)a1,(Atom)a2);
			
			limit = Math.min(g1.getSize(), g2.getSize());
			pairptr = -1;
			lastPair = new AtomPair[limit];
			forbidden1 = new int[nc1];
			forbidden2 = new int[nc2];
			resultcount = 0;
			seq1atoms = new BitSet(nc1);
			seq2atoms = new BitSet(nc2);
			
			// start enumerating paths from each node at a time
			// after each round, previous start-node should be in forbidden
			for (AtomPair ap : atompairs)
			{
				// reset forbids
				Arrays.fill(forbidden1, 0);
				Arrays.fill(forbidden2, 0);
				
				Add(ap);
				while (Forward() || Backward()) 
				{
					;
				}
			}
		}
		
		private boolean reducedvalid(Node n1, Node n2)
		{
			if (!reduced)
				return true;
			
			if (n1.getCoreDist() != n2.getCoreDist() && (n1.getCoreDist() == 0 || n2.getCoreDist() == 0))
				return false;
			
			return true;
		}
		
		private boolean Forward()
		{
			// stop if full
			if (pairptr+1 == limit)
				return false;
			if (pairptr == -1)
				return false;		
			
			// find next pair to add
			AtomPair last = lastPair[pairptr];
			for (Node ne1 : last.a1.getNodeNeighbors())
			{
				if (seq1atoms.get(ne1.getId()) || (forbidden1[ne1.getId()] > 0 && forbidden1[ne1.getId()] <= pairptr+2)) // already in seq, or forbidden
					continue;
				
				for (Node ne2 : last.a2.getNodeNeighbors())
				{
					if (seq2atoms.get(ne2.getId()) || (forbidden2[ne2.getId()] > 0 && forbidden2[ne2.getId()] <= pairptr+2)) // already in seq, or forbidden
						continue;
					
					// another way is to hold pairs instead of single atoms
					if (ne1.getSymbol().equals(ne2.getSymbol()) && reducedvalid(ne1,ne2))
					{
						Add(new AtomPair(ne1,ne2));
						return true;
					}
				}
			}
			
			return false;
		}
	
		private void Add(AtomPair ap)
		{
			Node a1 = ap.a1;
			Node a2 = ap.a2;
			
			// add pair to lastAtom
			lastPair[++pairptr] = ap;
			seq1atoms.set(a1.getId());
			seq2atoms.set(a2.getId());
			resultcount++;
			
			// construct the atom sequences to weight the paths
			Node[] atoms1 = new Atom[pairptr+1];
			Node[] atoms2 = new Atom[pairptr+1];
			
			for (int i = 0; i < pairptr + 1; i++)
			{
				atoms1[i] = lastPair[i].a1;
				atoms2[i] = lastPair[i].a2;
			}
			
			double w1 = weight(atoms1);
			double w2 = weight(atoms2);
			
			kvalue += w1*w2;
		}
	
		private boolean Backward()
		{
			if (pairptr == -1) // Alkioita nolla
				return false;
	
			AtomPair last = lastPair[pairptr];
			
			forbidden1[last.a1.getId()] = pairptr + 1;
			forbidden2[last.a2.getId()] = pairptr + 1;
	
			Remove(last);
			
			return true;
		}
	
		private void Remove(AtomPair ap)
		{
			seq1atoms.set(ap.a1.getId(), false);
			seq2atoms.set(ap.a2.getId(), false);
			pairptr--;
		}
	}
}

