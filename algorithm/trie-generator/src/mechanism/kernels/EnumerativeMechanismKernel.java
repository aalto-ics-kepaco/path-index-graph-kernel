package mechanism.kernels;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;

import mechanism.*;
import mechanism.graphs.*;

/*
 * Enumerative non-probabilistic walk kernel
 * two variants: 
 *    (1) non-decomposable gaussian weighting function
 *        f(w) = gauss(avg_d)
 *    (2) decomposable gaussian weighting function
 *        f(w) = \prod gauss(d(v))
 * 
 *    The first variant requires a DP approach where a histogram of distancesums are kept
 *    The second one is trivial: F_k(i) = g(d_v) * \sum_j F_(k-1)(j)
 * 
 *    The second one requires finetuning of the sigma to quite large value
 *    
 */

public class EnumerativeMechanismKernel extends SequenceKernel
{
	private double[][] Z; // Z[i][k] equals Z-constant of reaction 'i' at level 'k'
	private int maxdistsum = 100;
	private int maxdv = 0;

	
	public EnumerativeMechanismKernel(Graph[] graphs, KernelParams params)
	{
		super(graphs, params);
		
		Z = new double[graphs.length][params.maxlen+1];
	}

	public double compute(Graph g1, Graph g2)
	{
		return super.compute(g1,g2);
	}
	
	private void compute_Zconstants(Graph g)
	{
		// count distance-sums of all walks ending at 'i' of length 'k'
		// also, to compute distance sums, we need two columns of them
		// thus -> DSprev, DScurr
		//
		// we have: DScurr[i][k  ][distsum]
		//          DSprev[i][k-1][distsum]
		//
		
		int[][] DSprev = new int[g.getSize()][maxdistsum+maxdv];
		int[][] DScurr = new int[g.getSize()][maxdistsum+maxdv];
		double zsum = 0.0;
		
		// initialize
		// DSprev is now k=1
		// DScurr is now k=2
		for (int i = 0; i < g.getSize(); i++)
		{
			Node a = g.getNodes()[i];
			DSprev[a.getId()][a.getCoreDist()] = 1;
			zsum += logistic(a);
		}
		
		// kernel value for length k=1 walks 
		Z[g.getIndex()][1] = zsum;
		
		for (int l = 2; l <= params.maxlen; l++)
		{
			zsum = 0.0;
			
			for (int i = 0; i < g.getSize(); i++)
			{
				Node v = g.getNodes()[i];
				int dv = v.getCoreDist();
				
				int[] newds = new int[maxdistsum+maxdv]; // assume zero-initialized
				
				for (Node ne : v.getNodeNeighbors())
				{
					// shift DS by d(v_i), add to newds'
					
					for (int ds = 0; ds < maxdistsum; ds++)
					{
						newds[ds+dv] += DSprev[ne.getId()][ds];
					}
				}
				
				DScurr[i] = newds;
				
				// compute Z's
				zsum += computeZ(newds, l);
			}
			
			// f(w)'s are counted twice (both directions)
			zsum = zsum/2;
			
			// swap currents into prevs
			DSprev = DScurr;
			
			// add zsum
			Z[g.getIndex()][l] = zsum;
		}
	}
	
	protected double weight(Node[] nodes)
	{
		// apply the logistic function
		double val = 1.0;
		
		for (Node x : nodes)
		{
			val *= lambda(x) * weight(x);
		}
		
		return val;
	}
	
	protected double nondecomposable_weight(Node[] nodes)
	{
		// TODO
		return 1.0;
	}
	
	protected double walks(ProductGraph pg)
	{
		double[][] F = new double[pg.getNodeCount()][params.maxlen+1];
		double sum = 0.0;
		
		// initialize
		for (PGNode v : pg.getNodes())
		{
			F[v.getId()][0] = 0.0;
			F[v.getId()][1] = lambda(v) * weight(v);
			sum += F[v.getId()][1];
		}
		
//		for (Edge e : pg.getG1().getEdges())
//			if (e.getSource().getSymbol().equals("C") && e.getTarget().getSymbol().equals("C") && e.getChangetype() == 0)
//				System.out.println(e.getSource().getId() + " - " + e.getTarget().getId());
//		System.out.println();
		
		
		double oldsum = 0.0;
		for (int l = 2; l <= params.maxlen; l++)
		{
//			int z = 0;
			
			
			oldsum = sum;
			for (PGNode v : pg.getNodes())
			{
				
				double val = 0.0;
				for (PGNode u : v.getNodeNeighbors())
				{
					val += F[u.getId()][l-1];
					
//					if (v.a1.getSymbol().equals("C") && u.a1.getSymbol().equals("C") && v.a1.getEdge(u.a1).getChangetype() == 0 &&
//						v.a1.getId()==2 && u.a1.getId()==6)
//					{
//						z++;
//						System.out.println(v.a2.getId() + " - " + u.a2.getId());
//					}
				}
				
				// add gaussian function value at the end
				F[v.getId()][l] = lambda(v) * weight(v) * val;
				sum += F[v.getId()][l];
			}
			
			double increase = (sum - oldsum) / oldsum; // in percentages / 100
			
			// stop if change in score is negligible -> convergence occurred
			if (increase < params.epsilon)
				break;
		}
		
		return sum;
	}
	
	protected double nontottering_walks(ProductGraph pg)
	{
		// uses message passing scheme to compute non-tottering walks
		//
		// keeps a M-matrix for messages across edges, these are walk counts of everything
		// 'behind' that node.
		//
		
		int EDGES = pg.getEdgeCount();
		
		double[][] F = new double[pg.getNodeCount()][params.maxlen+1];
		double[][] M = new double[pg.getEdgeCount()*2][params.maxlen+1]; // both directions
		
		double sum = 0.0;
		double oldsum = 0.0;
		
		// initialize cases 0 and 1 of probability sums (F)
		for (PGNode v : pg.getNodes())
		{
			F[v.getId()][0] = 0.0;
			F[v.getId()][1] = lambda(v) * weight(v);
			sum += F[v.getId()][1];
		}
		
		// initialize trivial cases for messages (M)
		for (PGEdge e : pg.getEdges())
		{
			M[e.getId()][0] = 0.0;
			M[e.getId()][1] = weight(e.getSource()) * weight(e.getSource());
			M[e.getId()+EDGES][0] = 0.0;
			M[e.getId()+EDGES][1] = weight(e.getTarget()) * weight(e.getTarget());
		}
		
		for (int l = 2; l <= params.maxlen; l++)
		{
			oldsum = sum;

			for (PGNode v : pg.getNodes())
			{
				double val = 0.0;
				
				for (PGNode u : v.getNodeNeighbors())
				{
					PGEdge e = u.getEdge(v);
					
					int eid = e.getId();
					if (e.getSource() == v)
						eid += EDGES;
					
					val += M[eid][l-1];
				}
				
				F[v.getId()][l] = lambda(v) * weight(v) * val;
				sum += F[v.getId()][l];
			}
			
			// define new messages (M)
			for (PGEdge e : pg.getEdges())
			{
				PGNode src = e.getSource();
				PGNode tgt = e.getTarget();
				
				// forward direction
				double val = 0.0;
				for (PGNode u : src.getNodeNeighbors())
				{
					if (u.a1 == tgt.a1 || u.a2 == tgt.a2)
						continue;
					
					PGEdge e2 = u.getEdge(src);
					
					int eid = e2.getId();
					if (e2.getSource() == src)
						eid += EDGES;

					val += M[eid][l-1];
				}
				
				M[e.getId()][l] = lambda(src) * weight(src) * val;
				
				// backward direction
				val = 0.0;
				for (PGNode u : tgt.getNodeNeighbors())
				{
					if (u.a1 == src.a1 || u.a2 == src.a2)
						continue;
					
					PGEdge e2 = u.getEdge(tgt);
					
					int eid = e2.getId();
					if (e2.getSource() == tgt)
						eid += EDGES;

					val += M[eid][l-1];
				}
				
				M[e.getId()+EDGES][l] = lambda(tgt) * weight(tgt) * val;
			}
			
			double increase = (sum - oldsum) / oldsum; // in percentages / 100
			
			// stop if change in score is negligible -> convergence occurred
			if (increase < params.epsilon)
				break;
		}
		
		return sum;
	}
	
	
	
	protected double nondecomposable_walks(ProductGraph pg)
	{
		// count distance-sums of all walks ending at 'i' of length 'k'
		// distance-sums for both w and w' needed (DS1, DS2)
		// also, to compute distance sums, we need two columns of them
		// thus -> DS1prev, DS1curr, DS2prev, DS2curr
		//
		// we have: DS1curr[i][k  ][distsum]
		//          DS2curr[i][k  ][distsum']
		//          DS1prev[i][k-1][distsum]
		//          DS2prev[i][k-1][distsum']
		
		int[][] DS1prev = new int[pg.getNodeCount()][maxdistsum+maxdv];
		int[][] DS1curr = new int[pg.getNodeCount()][maxdistsum+maxdv];
		int[][] DS2prev = new int[pg.getNodeCount()][maxdistsum+maxdv];
		int[][] DS2curr = new int[pg.getNodeCount()][maxdistsum+maxdv];
		double zsum = 0.0;
		double kvalue = 0.0;
		
		// initialize
		// DSprev is now k=1
		// DScurr is now k=2
		for (int i = 0; i < pg.getNodeCount(); i++)
		{
			PGNode v = pg.getNodes()[i];
			DS1prev[v.getId()][v.a1.getCoreDist()] = 1;
			DS2prev[v.getId()][v.a2.getCoreDist()] = 1;
			zsum += logistic(v);
		}
		
		// kernel value for length k=1 walks 
		kvalue += params.lambda * zsum / (Z[pg.getG1().getIndex()][1] * Z[pg.getG2().getIndex()][1]);
		
		for (int l = 2; l <= params.maxlen; l++)
		{
			zsum = 0.0;
			
			for (int i = 0; i < pg.getNodeCount(); i++)
			{
				PGNode v = pg.getNodes()[i];
				
				int dv1 = v.a1.getCoreDist();
				int dv2 = v.a2.getCoreDist();
				
				int[] newds1 = new int[maxdistsum+maxdv]; // assume zero-initialized
				int[] newds2 = new int[maxdistsum+maxdv];
				
				for (PGNode ne : v.getNodeNeighbors())
				{
					// shift DS by d(v_i), add to newds'
					
					array_add(newds1, DS1prev[ne.getId()], dv1);
					array_add(newds2, DS2prev[ne.getId()], dv2);
				}
				
				DS1curr[i] = newds1;
				DS2curr[i] = newds2;
				
				// compute Z's
				zsum += computeZx(newds1, newds2, l);
			}
			
			// f(w)'s are counted twice (both directions)
			zsum = zsum/2;
			
			// swap currents into prevs
			DS1prev = DS1curr;
			DS2prev = DS2curr;
			
			// add zsum
			kvalue += Math.pow(params.lambda,l) * zsum / (Z[pg.getG1().getIndex()][l] * Z[pg.getG2().getIndex()][l]);
		}
		
		return kvalue;
	}
	
	protected double nondecomposable_nontottering_walks(ProductGraph pg)
	{
		// count distance-sums of all walks ending at 'i' of length 'k'
		// distance-sums for both w and w' needed (DS1, DS2)
		// also, to compute distance sums, we need two columns of them
		// thus -> DS1prev, DS1curr, DS2prev, DS2curr
		//
		// we have: DS1curr[i][k  ][distsum]
		//          DS2curr[i][k  ][distsum']
		//          DS1prev[i][k-1][distsum]
		//          DS2prev[i][k-1][distsum']
		
		
		pg.computeCorrectedPathCount(10);
		
		
		int[][] DS1prevprev = new int[pg.getNodeCount()][maxdistsum+maxdv];
		int[][] DS2prevprev = new int[pg.getNodeCount()][maxdistsum+maxdv];
		int[][] DS1prev = new int[pg.getNodeCount()][maxdistsum+maxdv];
		int[][] DS2prev = new int[pg.getNodeCount()][maxdistsum+maxdv];
		int[][] DS1curr = new int[pg.getNodeCount()][maxdistsum+maxdv];
		int[][] DS2curr = new int[pg.getNodeCount()][maxdistsum+maxdv];
		double zsum = 0.0;
		double kvalue = 0.0;
		
		// initialize
		// DSprevprev is now k=1
		// DSprev is now k=2
		// DScurr is now k=3
		// compute only prev and prevprev, curr is applied at the next stage
		for (PGNode v : pg.getNodes())
		{
			DS1prevprev[v.getId()][v.a1.getCoreDist()] = 1;
			DS2prevprev[v.getId()][v.a2.getCoreDist()] = 1;
			zsum += logistic(v);
		}
		
		for (PGNode v : pg.getNodes())
		{
			for (PGNode u : v.getNodeNeighbors())
			{
				array_add(DS1prev[v.getId()], DS1prevprev[u.getId()], v.a1.getCoreDist());
				array_add(DS2prev[v.getId()], DS2prevprev[u.getId()], v.a2.getCoreDist());
			}
			
			zsum += computeZx(DS1prev[v.getId()], DS2prev[v.getId()], 2);
		}
		
		// kernel value for length k=1 walks 
		kvalue += params.lambda * zsum / (Z[pg.getG1().getIndex()][1] * Z[pg.getG2().getIndex()][1]);
		
		for (int l = 3; l <= params.maxlen; l++)
		{
			zsum = 0.0;
			
			for (PGNode v : pg.getNodes())
			{
				int dv1 = v.a1.getCoreDist();
				int dv2 = v.a2.getCoreDist();
				
				int[] newds1 = new int[maxdistsum+maxdv]; // assume zero-initialized
				int[] newds2 = new int[maxdistsum+maxdv];
				
				for (PGNode u : v.getNodeNeighbors())
				{
					// shift DS by d(v_i), add to newds'
					
					// first add the neighboring dist-sums with offset
					array_add(newds1, DS1prev[u.getId()], dv1);
					array_add(newds2, DS2prev[u.getId()], dv2);
					
					for (PGNode s : pg.neighnodes.get(u).get(v))
					{
						if (s.a1 == v.a1)
							array_minus(newds1, DS1prevprev[s.getId()], u.a1.getCoreDist()+dv1);
						if (s.a2 == v.a2)
							array_minus(newds2, DS2prevprev[s.getId()], u.a2.getCoreDist()+dv2);
					}
				}
				
				DS1curr[v.getId()] = newds1;
				DS2curr[v.getId()] = newds2;
				
				// compute Z's
				zsum += computeZx(newds1, newds2, l);
			}
			
			// f(w)'s are counted twice (both directions)
			zsum = zsum/2;
			
			// swap currents into prevs
			DS1prevprev = DS1prev;
			DS2prevprev = DS2prev;
			
			DS1prev = DS1curr;
			DS2prev = DS2curr;
			
			DS1curr = new int[pg.getNodeCount()][maxdistsum+maxdv];
			DS2curr = new int[pg.getNodeCount()][maxdistsum+maxdv];
			
			// add zsum
			kvalue += Math.pow(params.lambda,l) * zsum / (Z[pg.getG1().getIndex()][l] * Z[pg.getG2().getIndex()][l]);
		}
		
		return kvalue;
	}
	
	private void array_add(int[] newds, int[] array, int offset)
	{
		for (int ds = 0; ds < maxdistsum; ds++)
		{
			newds[ds+offset] += array[ds];
		}
	}
	
	private void array_minus(int[] newds, int[] array, int offset)
	{
		for (int ds = 0; ds < maxdistsum; ds++)
		{
			newds[ds+offset] -= array[ds];
		}
	}
	
	
	protected double weight(PGNode v)
	{
		if (params.kw == KernelWeight.Exponential)
			return exp(v);
		else if (params.kw == KernelWeight.Logistic)
			return logistic(v);
		else if (params.kw == KernelWeight.Diffusion)
			return diffusion(v);
		
		return 1.0;
	}
	
	protected double weight(Node v)
	{
		if (params.kw == KernelWeight.Exponential)
			return exp(v);
		else if (params.kw == KernelWeight.Logistic)
			return logistic(v);
		else if (params.kw == KernelWeight.Diffusion)
			return diffusion(v);
		
		return 1.0;
	}
	

	private double computeZ(int[] DS, int k)
	{
		// in: array DScurr[i] and k
		// out: sum DS[i]*f(i/k)
		
		double sum = 0.0;
		for (int i = 0; i < maxdistsum+maxdv; i++)
		{
			sum += (DS[i] * logistic(1.0*i/k));
		}
		
		return sum;
	}
	
	private double computeZx(int[] DS1, int[] DS2, int k)
	{
		return computeZ(DS1, k) * computeZ(DS2, k);
	}
}

