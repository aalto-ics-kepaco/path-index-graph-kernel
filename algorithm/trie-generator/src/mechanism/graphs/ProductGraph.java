package mechanism.graphs;

import java.io.*;
import java.util.*;

import mechanism.*;


public class ProductGraph extends Graph
{
	protected Graph g1,g2;
	
	// overrides similarly named guys from superclass
	protected PGNode[] nodes;
	protected PGEdge[] edges;
	
	protected KernelParams params;
	
	// 3d array giving the set of neighboring nodes of 'x' given 'y' at
	// neighnodes[x][y] = [ nodes with (*|x) or (x|*) or (x|×) ]
	public Map<PGNode, Map<PGNode,List<PGNode>> > neighnodes; // 3d-array giving the set of neighboring nodes
	
	// node/edge counts
	protected int nc;
	protected int ec; 
	
	
	
	protected ProductGraph()
	{
		// empty
	}
	
	
	/*
	 * forms a direct product graph by matching all nodes against each other
	 * and only matching
	 */
	public ProductGraph(Graph g1, Graph g2, KernelParams params)
	{
		this(g1,g2,params,true);
	}
	
	
	/*
	 * Forms a direct product graph (g1 X g2)
	 * First, create nodes for all pairs of nodes (N_1 X N_2) with matching labels
	 * Then, create edges between the pair-nodes if they are neighbors in the orinigal graphs
	 * 
	 * Creation of product graphs is a major bottleneck in the performance, as they can become quite large
	 * Optimizations:
	 *  1) count the number of PGNodes by multiplying matching label atom counts (fast)
	 *  2) form the PGNodes by pre-sorting the atoms on both sides by label, and forming PGNodes
	 *     only inside these blocks (diagonal blocks on the (N_1 X N_2) matrix
	 *  3) its impossible to count the number of PGEdges quickly, so do it by enumerating
	 *  4) PGEdges are computed by checking neighborhoods of Atoms from both sides
	 *  
	 *  
	 */
	public ProductGraph(Graph g1, Graph g2, KernelParams params, boolean nodematch)
	{
		this.g1 = g1;
		this.g2 = g2;
		this.params = params;
		
		createNodes(nodematch);
		createEdges();
		

//		int count = 0;
//		for (PGNode v : nodes)
//		{
//			for (PGNode u : nodes)
//			{
//				// normal C-C pairing
//				if (v.a1.getSymbol().equals("C") && 
//					u.a1.getSymbol().equals("C") && 
//					v.a1.isNeighbor(u.a1) && 
//					v.a2.isNeighbor(u.a2) && 
//					v.a1.getEdge(u.a1).getChangetype() == 0 && 
//					v.a1.getId() == 2 &&
//					u.a1.getId() == 3)
//				{
//					count++;
//				}
//			}
//		}
//		
//		int z = 3;
	}
	
	protected void createNodes(boolean nodematch)
	{
		nc = 0;
		ec = 0;
		
		// count first the number of PG-nodes and edges, and then create correct size
		// arrays, probably faster than adding right away to ArrayList
		
		// optimized case for symbolmatching pg-nodes
		if (nodematch)
		{
			// count the number of nodes by multiplying symbol-counts
			for (String s : g1.getFormula().keySet())
			{
				if (g2.getFormula().containsKey(s))
					nc += g1.getFormula().get(s) * g2.getFormula().get(s);
			}
			nodes = new PGNode[nc];
			
			// new try: v2
			// 'eat' from front labels until both arrays contain same fronts
			//  then loop them through and go into 'eat' mode again
			int iptr = 0;
			int jptr = 0;
			int nextjptr = 0;
			String front = "";
			int index = 0;
			
			Node[] g1nodes = g1.getNodes().clone();
			Node[] g2nodes = g2.getNodes().clone();
			Arrays.sort(g1nodes, new NodeSymbolComparator());
			Arrays.sort(g2nodes, new NodeSymbolComparator());
	
			while (iptr < g1nodes.length && jptr < g2nodes.length && nextjptr < g2nodes.length)
			{
				Node a1 = g1nodes[iptr];
				Node a2 = g2nodes[jptr];
				
				// if fronts out of sync -> chomp 'em until match
				if (a1.getSymbol().compareTo(a2.getSymbol()) < 0)
					iptr++;
				else if (a1.getSymbol().compareTo(a2.getSymbol()) > 0)
					jptr++;
				else
				{
					front = a1.getSymbol();
					int i;
					// fronts match -> enumerate combinations
					for (i = iptr; i < g1nodes.length; i++)
					{
						a1 = g1nodes[i];
						
						// quit if symbol changes
						if (!a1.getSymbol().equals(front))
						{
							iptr = i;
							jptr = nextjptr;
							break;
						}
						
						for (int j = jptr; j < g2nodes.length; j++)
						{
							a2 = g2nodes[j];
							
							if (!a2.getSymbol().equals(front))
							{
								nextjptr = j;
								break;
							}
							
							if (match(a1,a2))
							{
								PGNode next = new PGNode(this, a1,a2);
								next.id = index;
								nodes[index++] = next;
							}
						}
					}
					iptr = i;
				}
			}
		}
		else // directly create the pgnodes by matching everything
		{
			nc = g1.getSize() * g2.getSize();
			nodes = new PGNode[nc];
			int index = 0;
			
			for (Node v1 : g1.getNodes())
			{
				for (Node v2 : g2.getNodes())
				{
					if (match(v1,v2))
					{
						PGNode next = new PGNode(this,v1,v2);
						next.id = index;
						nodes[index++] = next;						
					}
				}
			}
		}
	}
	
	protected void createEdges()
	{
		
//		int z = 0;
		
		// count the number of necessary edges
		for (int i = 0; i < nodes.length-1; i++)
		{
			for (int j = i+1; j < nodes.length; j++)
			{
				PGNode n1 = nodes[i];
				PGNode n2 = nodes[j];
				
//				if (n1.a1.isNeighbor(n2.a1) && n1.a2.isNeighbor(n2.a2) && match(n1.a1.getEdge(n2.a1), n1.a2.getEdge(n2.a2)) && n1.a1.getSymbol().equals("C") && n1.a1.getEdge(n2.a1).getChangetype() == 0)
//				{
//					z++;
//				}
				
				if (n1.a1.isNeighbor(n2.a1) && n1.a2.isNeighbor(n2.a2) && match(n1.a1.getEdge(n2.a1), n1.a2.getEdge(n2.a2)))
					ec++;
			}
		}
		
		// construct the edge array (suitable size)
		edges = new PGEdge[ec];
		
		// create the edges
		id = 0;
		for (int i = 0; i < nodes.length-1; i++)
		{	
			for (int j = i+1; j < nodes.length; j++)
			{ 
				PGNode n1 = nodes[i];
				PGNode n2 = nodes[j];
				
				if (n1.a1.isNeighbor(n2.a1) && n1.a2.isNeighbor(n2.a2) && match(n1.a1.getEdge(n2.a1), n1.a2.getEdge(n2.a2)))
				{
					PGEdge next = new PGEdge(this, n1,n2);
					next.id = id;
					edges[id++] = next;
				}
			}
		}
		
		// neighbors array
		neighnodes = new HashMap<PGNode, Map<PGNode,List<PGNode>>>();
		for (PGNode n1 : nodes)
		{
			neighnodes.put(n1, new HashMap<PGNode,List<PGNode>>());
			
			for (PGNode v2 : n1.nodeneighs)
			{
				neighnodes.get(n1).put(v2, new ArrayList<PGNode>());
				
				for (PGNode v : n1.nodeneighs)
				{
					// contains same on either side
					if (v.a1 == v2.a1 || v.a1 == v2.a2 || v.a2 == v2.a1 || v.a2 == v2.a2)
					{
						neighnodes.get(n1).get(v2).add(v);
					}
				}
			}
		}
	}
	
	public Graph getG1()
	{
		return g1;
	}
	
	public Graph getG2()
	{
		return g2;
	}
	
	// overwritten getters to specify the PGNode instead of Node
	public PGNode[] getNodes()
	{
		return nodes;
	}
	
	public PGEdge[] getEdges()
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
	
	protected boolean match(Edge b1, Edge b2)
	{
		if (params.edgematch)
			return b1.getType() == b2.getType() && b1.getChangetype() == b2.getChangetype() ? true : false;
		
		return b1.getChangetype() == b2.getChangetype() ? true : false;
	}
	
	protected boolean match(Node a1, Node a2)
	{
		if (params.nodematch)
			return a1.getSymbol().equals(a2.getSymbol()) ? true : false;
		
		return true;
	}
	
	public String toString()
	{
		return g1 + "<->" + g2;
	}
	
	/*
	 * The walk/path computation routines
	 * all routines expect the dynamic programming matrix to have
	 * nodes as rows and columns as counts of 'k'-length walk/paths
	 * thus D_ij corresponds to walks of length 'j' ending at node 'i'
	 */
	
	
	// compute number common walks upto 'maxlength'
	// the counts are inherently infinite, thus maxlength should be 
	// relatively small (say the count of nodes)
	public long computeWalkCounts(int maxlength)
	{
		// use dynamic programming with n * k matrix
		long[][] D = new long[nc][maxlength+1];
		long sum = 0;
		
		// initialize
		for (int i = 0; i < nc; i++)
		{
			D[i][0] = 0;
			D[i][1] = 1;
			sum += 1;
		}
		
		for (int l = 2; l < maxlength+1; l++)
		{
			for (int i = 0; i < nc; i++)
			{
				long val = 0;
				for (PGNode ne : nodes[i].nodeneighs)
				{
					val += D[ne.id][l-1];
				}
				
				D[i][l] = val;
				sum += val;
			}
		}
		
		return sum;
	}
	
	public long[] computeWalkCountArray(int maxlength)
	{
		// use dynamic programming with n * k matrix
		long[][] D = new long[nc][maxlength+1];
		long[] counts = new long[maxlength+1];
		
		// initialize
		for (int i = 0; i < nc; i++)
		{
			D[i][0] = 0;
			D[i][1] = 1;
			counts[1] += 1;
		}
		
		for (int l = 2; l < maxlength+1; l++)
		{
			for (int i = 0; i < nc; i++)
			{
				long val = 0;
				for (PGNode ne : nodes[i].nodeneighs)
				{
					val += D[ne.id][l-1];
				}
				
				D[i][l] = val;
				counts[l] += val;
			}
		}
		
		return counts;
	}
	
	
	// compute the weighted sum of common walks upto 'maxlength'
	// or until convergence if lambda is sufficiently small
	// all walks are weighted such that each node (a,b) contributes 
	//  lambda^d(a)*lambda^d(b) to the score
	public double computeWalkScore(int maxlength)
	{
		// use dynamic programming with n * k matrix
		double[][] D = new double[nc][maxlength+1];
		double sum = 0.0;
		
		// initialize
		for (int i = 0; i < nc; i++)
		{
			PGNode v = nodes[i];
			
			D[i][0] = 0.0;
			D[i][1] = logistic(v);
			sum += D[i][1];
		}
		
		double oldsum = 0.0;
		for (int l = 2; l < maxlength + 1; l++)
		{
			oldsum = sum;
			for (int i = 0; i < nc; i++)
			{
				PGNode v = nodes[i];
				
				double val = 0.0;
				for (PGNode ne : nodes[i].nodeneighs)
				{
					// ne.match computes the bond-match
					val += D[ne.id][l-1];
				}
				
				D[i][l] = logistic(v) * (1.0 - params.lambda) * (1.0 - params.lambda) * val;
				sum += D[i][l];
			}
			
			double increase = (sum - oldsum) / oldsum; // in percentages / 100
			
			// stop if change in score is negligible -> convergence occurred
			if (increase < 0.00001)
				break;
		}
		
		return sum;
	}
	
	// compute the number of paths in the product graph,
	// these don't correspond to paths in the original graphs
	public double computePathCount(int maxlength)
	{
		long sum = 0;
		
		if (maxlength > nc)
			maxlength = nc;
		
		// use dynamic programming with n * k matrix
		long[][] D = new long[nc][maxlength+1];  // +1 for zero length paths
		
		// initialize
		for (int i = 0; i < nc; i++)
		{
			D[i][0] = 0; // zero-length paths, none of them
			D[i][1] = 1; //  one-length paths,  one of them
			sum += 1;
		}
		
		boolean stop;
		for (int l = 2; l <= maxlength; l++)
		{
			stop = true;
			for (int i = 0; i < nc; i++)
			{
//				PGNode v = nodes[i];
				
				long val = 0;
				for (PGNode ne : nodes[i].nodeneighs)
				{
					for (int k = 1; k < l; k++)
					{
						if (k % 2 != 0) // odd
							val += D[ne.id][l-k];
						else // even
							val -= D[i][l-k];
					}
				}
				D[i][l] = val;
				sum += D[i][l];
				
				if (val != 0)
					stop = false;
			}
			
			if (stop)
				break;
		}
			
		return sum;
	}
	
	
	// computes weighted sum of common paths in the product graph
	// these paths don't correspond to paths in the original graphs
	public double computePathScore(int maxlength)
	{
		double sum = 0.0;
		
		if (maxlength > nc)
			maxlength = nc;
		
		// use dynamic programming with n * k matrix
		double[][] D = new double[nc][maxlength+1];  // +1 for zero length paths
		
		// initialize
		for (int i = 0; i < nc; i++)
		{
			PGNode v = nodes[i];
			
			D[i][0] = 0; // zero-length paths, none of them
			D[i][1] = lambda(v) * 1; // one-length paths, one of them
			sum += D[i][1];
		}
		
		double oldsum = 0.0;
		for (int l = 2; l <= maxlength; l++)
		{
			oldsum = sum;
			for (int i = 0; i < nc; i++)
			{
				PGNode v = nodes[i];
				
				double val = 0.0;
				for (PGEdge ne : nodes[i].edgeneighs.values())
				{
					for (int k = 1; k < l; k++)
					{
						if (k % 2 != 0) // odd
							val += D[ne.getOther(v).id][l-k];
						else // even
							val -= D[i][l-k];
					}
					
					// multiply the edge weights in
//					val = ne.match * val;
				}
				
				D[i][l] = lambda(v) * val;
				sum += D[i][l];
			}	

			double increase = (sum - oldsum) / oldsum; // in percentages / 100
				
			// stop if change in score is negligible -> convergence occurred
			if (increase < 0.00001)
				break;
		}
			
		return sum;
	}
	
	// computes the corrected path count
	public long computeCorrectedPathCount(int maxlength)
	{
		long sum = 0;
		
		if (maxlength > Math.max(g1.getNodes().length, g2.getNodes().length))
			maxlength = Math.max(g1.getNodes().length, g2.getNodes().length);
		
		// use dynamic programming with n * k matrix
		long[][] D = new long[nc][maxlength+1];  // +1 for zero length paths
		
		// initialize
		for (int i = 0; i < nc; i++)
		{
			D[i][0] = 0; // zero-length paths, none of them
			D[i][1] = 1; //  one-length paths,  one of them
			sum += 1;
		}
		
		boolean stop;
		for (int l = 2; l <= maxlength; l++)
		{
			stop = true;
			for (int i = 0; i < nc; i++)
			{
				PGNode v = nodes[i];
				
				long val = 0;
				for (PGNode ne : nodes[i].nodeneighs)
				{
					// add the first term of each neighbor
					val += D[ne.id][l-1];
					
//					int s = neighnodes.get(ne).get(v).size();
					
					// add the sum rest of the positive terms multiplied by neighcounts
//					for (int k = 3; k < l; k+=2)
//					{
//						val += neighnodes.get(ne).get(v).size() * D[ne.id][l-k];
//					}
					
					// minus the negative term serie for each similar neighbor
					for (PGNode ine : neighnodes.get(ne).get(v)) // neighbors of 'ne' which contain some 'v'
					{
						for (int k = 2; k < l; k++)
						{
							if (k % 2 == 0)  // even
								val -= D[ine.id][l-k];
							else  // odd
							{
								for (PGNode ine2 : neighnodes.get(v).get(ne))
								{
									val += D[ine2.id][l-k];
								}
							}
						}
					}
				}
				D[v.id][l] = val;
				sum += val;
				
				if (val != 0)
					stop = false;
			}
			
			if (stop)
				break;
		}
			
		return sum;
	}
	
	// computes the corrected path count
	public double computeCorrectedPathScore(int maxlength)
	{
		double sum = 0.0;
		
		if (maxlength > Math.max(g1.getNodes().length, g2.getNodes().length))
			maxlength = Math.max(g1.getNodes().length, g2.getNodes().length);
		
		// use dynamic programming with n * k matrix
		double[][] D = new double[nc][maxlength+1];  // +1 for zero length paths
		
		// initialize
		for (int i = 0; i < nc; i++)
		{
			PGNode v = nodes[i];
			
			D[i][0] = 0.0; // zero-length paths, none of them
			D[i][1] = lambda(v) * 1.0; // one-length paths, one of them
			sum += D[i][1];
		}
		
		double oldsum = 0.0;
		for (int l = 2; l <= maxlength; l++)
		{
			oldsum = sum;
			for (int i = 0; i < nc; i++)
			{
				PGNode v = nodes[i];
				
				double val = 0.0;
				for (PGEdge ne : nodes[i].edgeneighs.values())
				{
					// add the first term of each neighbor
					val += D[ne.getOther(v).id][l-1];
					
					// add the sum rest of the positive terms multiplied by neighcounts
					for (int k = 3; k < l; k+=2)
					{
						val += neighnodes.get(ne.getOther(v)).get(v).size() * D[ne.getOther(v).id][l-k];
					}
					
					// minus the negative term serie for each similar neighbor
					for (PGNode ine : neighnodes.get(ne.getOther(v)).get(v)) // neighbors of 'ne' which contain some 'v'
					{
						for (int k = 2; k < l; k+=2)
						{
							val -= D[ine.id][l-k];
						}
					}
					
//					val = val;
				}
				
				D[i][l] = lambda(v) * val;
				sum += D[i][l];
			}

			double increase = (sum - oldsum) / oldsum; // in percentages / 100
				
			// stop if change in score is negligible -> convergence occurred
			if (increase < 0.00001)
				break;
		}
			
		return sum;
	}


	public void writeDotFile(String filename)
	{
		FileWriter fw;
		try
		{
			fw = new FileWriter(filename);
		
			fw.write("digraph G {\n");
			
			for (PGNode n : nodes)
			{
				fw.write(n.id + " [title=" + n.a1.getSymbol() + n.a1.getId() + "|" + n.a2.getSymbol() + n.a2.getId() + "];\n");
			}
			
			for (PGEdge e : edges)
			{
				fw.write(e.source.id + " -> " + e.target.id + ";\n");
			}
				
			fw.close();
		} catch (IOException e1)
		{
			;
		}
	}
	
	private double logistic(double x)
	{
		return 1.0 / (1.0 + Math.exp(params.beta * (x - 0.5))); // use logistic function instead of gaussian
	}
	
	private double logistic(PGNode v)
	{
		return logistic(v.a1.getCoreDist()) * logistic(v.a2.getCoreDist());
	}
	
	private double lambda(PGNode v)
	{
		return (1.0 - params.lambda) * (1.0 - params.lambda);
	}

	// do nothing
	public Graph createSubgraph(BitSet nodebits)
	{
		return null;
	}
	
}





