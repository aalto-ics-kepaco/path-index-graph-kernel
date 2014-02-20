package mechanism.kernels;


import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import mechanism.*;
import mechanism.graphs.*;

// computes shortest paths kernel
//
// K(G,G') = sum_{v1,v2} sum_{v1',v2'} K_sp(v1<->v2, v1'<->v2'),
// where K_sp is the kernel for two individual shortest paths,
// namely: K_sp(e,e') = prod_i L(e_i) = L(e'_i),
// basically we only match path pairs with identical length, and identical atoms
//
//

public class ShortestPathsKernel extends Kernel
{
	public ShortestPathsKernel(Graph[] graphs, KernelParams params)
	{
		super(graphs, params);
	}

	public double compute(Graph g1, Graph g2)
	{
		// first compute list of shortest paths in both input graphs
		// i.e. take all atom pairs in G1 and G2 and compute the atom sequence
		// -> we get size(g1)*size(g1)-1/2 pairs 
		
		Map<String,Integer> g1sp = FloydWarshall.floydwarshall(g1);
		Map<String,Integer> g2sp = FloydWarshall.floydwarshall(g2);
		
		double k = 0.0;
		for (String path : g1sp.keySet())
		{
			if (g2sp.containsKey(path))
			{
				// add weight decay
				if (params.lambda != 1.0) // dangerous comparison..
				{
					// count the number of edges
					int count = 1;
					for (int i = 0; i < path.length(); i++)
						if (path.charAt(i) == '1' || path.charAt(i) == '0')
							count++;
					
					k += Math.pow(params.lambda, count) * (g1sp.get(path) * g2sp.get(path));
				}
				else // don't use weight decay
					k += (g1sp.get(path) * g2sp.get(path));
			}
		}
		
		return k;
	}
}

class FloydWarshall
{
	private static int[][] cost;
	private static int[][] next;
	private static Graph g;
	
	public static Map<String,Integer> floydwarshall(Graph g)
	{
		FloydWarshall.g = g;
		
		// initializes to zeros
		int N = g.getSize();
		cost = new int[N][N];
		next = new int[N][N];
		
		for (int i = 0; i < N; i++)
			Arrays.fill(next[i], -1);
		
		// initialize costs to 0 (identity), 1 (neighbors) or infinity (rest)
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				Node v = g.getNodes()[i];
				Node u = g.getNodes()[j];
				
				if (v.isNeighbor(u) && i != j)
					cost[i][j] = 1;
				else
				{
					if (i == j)
						cost[i][j] = 0;
					else
						cost[i][j] = Integer.MAX_VALUE/2;
				}
			}
		}
		
		// compute in dynamic programming fashion the shortest distance between
		// all nodes with 'k' steps starting with k=1,...
		//
		// if we can shorten the distance between two nodes by including an additional
		// node, we go for it
		for (int k = 0; k < N; k++)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					if (cost[i][k] + cost[k][j] < cost[i][j])
					{
						int oldcost = cost[i][j];
						int newcost = cost[i][k] + cost[k][j];
						
						cost[i][j] = cost[i][k] + cost[k][j];
						next[i][j] = k;
					}
				}
			}
		}
		
		// get all paths
		Map<String,Integer> spcounts = new HashMap<String,Integer>();
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				if (i==j)
					continue;
				
				String path = getpath(i,j);
				
				// don't put nulls into map
				if (path == null)
					continue;
				
				if (spcounts.containsKey(path))
					spcounts.put(path, spcounts.get(path) + 1);
				else
					spcounts.put(path, 1);
			}
		}
		
		return spcounts;
	}
	
	private static String getpath(int i, int j)
	{
		// no path
		if (cost[i][j] == Integer.MAX_VALUE/2)
			return null;
		
		int mid = next[i][j];
		
		// neighbors
		if (mid == -1)
			return g.getNodes()[i].getSymbol() + g.getNodes()[i].getEdge(g.getNodes()[j]).getChangetype() + g.getNodes()[j].getSymbol();
		else
			return getpath(i,mid) + getpath(mid,j);
	}
}

