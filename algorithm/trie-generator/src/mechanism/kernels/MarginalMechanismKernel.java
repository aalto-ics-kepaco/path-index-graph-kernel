package mechanism.kernels;

import mechanism.*;
import mechanism.graphs.Atom;
import mechanism.graphs.Graph;
import mechanism.graphs.Node;
import mechanism.graphs.PGEdge;
import mechanism.graphs.PGNode;
import mechanism.graphs.ProductGraph;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;


public class MarginalMechanismKernel extends SequenceKernel
{
	/*
	 * From first principles:
	 * 
	 * Assume that both original graphs g1 and g2 have their own proper markov probability models
	 * Construct a product graph and use random walks on that with
	 *  
	 *  Ps(v,v') = Ps(v) * Ps(v')
	 *  Pt(v,v'|u,u') = Pt(v|u) * Pt(v'|u')
	 *  Pe(v,v') = Pe(v) * Pe(v')
	 * 
	 *  -> P(w) = Ps(v_0)*Ps(v_0') \prod Pt(v_i|v_{i-1})*Pt(v_i'|v_{i-1'}) Pe(v_k)*Pe(v_k')
	 * 
	 * Normalization required afterwards
	 * 
	 * Two variants: standard 1st order random walks and 2nd order non-tottering walks
	 * The non-tottering variant produces walks which don't repeat in the original graphs (not in the pg)
	 * 
	 */
	
	
	// probability structures
	protected double ps;  // P_s(v) * P_s(v')
	protected double pe;  // P_e(v) * P_e(v) 
	protected List<Double> ps1;  // Ps(v)*Ps(v)
	protected List<Double> pe1;  // Pe(v)*Pe(v')
	protected List<Map<Integer,Double>> pe2;  // Pe(v|u)*Pe(v'|u')
	protected List<Map<Integer,Double>> pt1;  // Pt(v|u)*Pt(v'|u')
	protected List<Map<Integer,Map<Integer,Double>>>  pt2; // Pt(v|u,s)*Pt(v'|u',s')
	
	// Z contains K_k(i,i) kernel values as Z[i][k]
	protected double[][] Z;
	
	public MarginalMechanismKernel(Graph[] graphs, KernelParams params)
	{
		super(graphs, params);
	}

	public void compute()
	{
		if (params.partialnorm)
			computeZ();
		
		if (params.kw == KernelWeight.Exponential && params.alpha <= 0.0)
			return;
		
		super.compute();
	}

	
	// (1) unnormalized probability functions

	protected double start_value(Node a)
	{
		if (params.kw == KernelWeight.Exponential)
			return exp(a);
		else if (params.kw == KernelWeight.Logistic)
			return logistic(a);
		else if (params.kw == KernelWeight.Diffusion)
			return diffusion(a);
		
		return 1.0;
	}
	
	protected double trans_value(Node curr, Node next)
	{
		if (params.kw == KernelWeight.Exponential)
			return exp(next);
		else if (params.kw == KernelWeight.Logistic)
			return logistic(next);
		else if (params.kw == KernelWeight.Diffusion)
			return diffusion(next);
		
		return 1.0;
	}
	
	protected double trans_value(Node prev, Node curr, Node next)
	{
		if (prev == next || prev == curr || curr == next)
			return 0.0;
		
		return trans_value(curr,next);
	}
	
	protected double end_value(Node a)
	{
		return 1.0 - params.lambda;
	}

		
	// (2) normalized probability functions
	
	protected double Pstart(Node a)
	{
//		// no prob if not connected to anything
//		if (a.getCoreDist() == Integer.MAX_VALUE)
//			return 0.0;
		
		double sum = 0.0;
		for (Node x : a.getParent().getNodes())
		{
			sum += start_value(x);
		}
		
		if (sum == 0.0)
			return 0.0;
		
		return start_value(a) / sum;
	}
	
	protected double Ptransition(Node from, Node to)
	{
		if (from == to || !from.isNeighbor(to))
			return 0.0;
		
		// no weighting if not connected
		if (from.getCoreDist() == Integer.MAX_VALUE)
			return 1.0;
		
		double sum = 0.0;
		for (Node a : from.getNodeNeighbors())
		{
			sum += trans_value(from,a);
		}
		
		if (sum == 0.0)
			return 0.0;
		
		return (1.0 - end_value(from)) * (trans_value(from,to) / sum);
	}	
	
	protected double Ptransition(Node prev, Node curr, Node next)
	{
		if (prev == next || curr == next || !curr.isNeighbor(next) || !prev.isNeighbor(curr))
			return 0.0;
		
		// no weighting if not connected
		if (curr.getCoreDist() == Integer.MAX_VALUE)
			return 1.0;
		
		double sum = 0.0;
		for (Node a : curr.getNodeNeighbors())
		{
			sum += trans_value(prev,curr,a);
		}
		
		if (sum == 0.0)
			return 0.0;
		
		return (1.0 - end_value(curr)) * (trans_value(prev,curr,next) / sum);
	}
	
	// n-th order markov transition probability, 
	protected double Ptransition(Node[] atoms, int nextid)
	{
		Node curr = atoms[nextid-1];
		Node next = atoms[nextid];
		
		// no weighting if not connected
		if (curr.getCoreDist() == Integer.MAX_VALUE)
			return 1.0;
		
		// check that next is not visited
		for (int i = 0; i < nextid; i++)
			if (atoms[i] == next)
				return 0.0;
		
		// otherwise next is not visited yet
		double sum = 0.0;
		for (Node ne : curr.getNodeNeighbors())
		{
			// check that 'ne' is not in atoms
			boolean neigh = false;
			for (int i = 0; i < nextid-1; i++)
			{
				if (ne == atoms[i])
				{
					neigh = true;
					break;
				}
			}
			if (!neigh)
				sum += trans_value(curr, ne);
		}
		
		if (sum == 0.0)
			return 0.0;
		
		assert (1.0 - end_value(curr)) * (trans_value(curr, next) / sum) <= 1.0;
		
		return (1.0 - end_value(curr)) * (trans_value(curr, next) / sum);
	}
	
	protected double Pend(Node curr)
	{
		// p_e = 1 - \sum p_t(out)
		double end = 1.0;
		for (Node ne : curr.getNodeNeighbors())
			end -= Ptransition(curr,ne);
		
		return end;
	}
	
	protected double Pend(Node prev, Node curr)
	{
		// p_e = 1 - \sum p_t(out)
		double end = 1.0;
		for (Node ne : curr.getNodeNeighbors())
			end -= Ptransition(prev,curr,ne);
		
		return end;
	}
	
	protected double Pend(Node[] atoms, int currid)
	{
		// p_e = 1 - \sum p_t(out)
		Node curr = atoms[currid];
		double end = 1.0;
		for (Node ne : curr.getNodeNeighbors())
		{
			Node[] newatoms = Arrays.copyOf(atoms, atoms.length+1);
			newatoms[currid+1] = ne;
			end -= Ptransition(newatoms, currid+1);
		}
		
		return end;
	}
	
	
	// (3) cached probability functions
	
	protected double Pstart(PGNode v)
	{
		return ps1.get(v.getId());
	}
	
	protected double Ptransition(PGNode curr, PGNode next)
	{
		return pt1.get(curr.getId()).get(next.getId());
	}
	
	protected double Ptransition(PGNode prev, PGNode curr, PGNode next)
	{
		return pt2.get(prev.getId()).get(curr.getId()).get(next.getId());
	}
	
	protected double Pend(PGNode v)
	{
		return pe1.get(v.getId());
	}
	
	protected double Pend(PGNode prev, PGNode curr)
	{
		return pe2.get(prev.getId()).get(curr.getId());
	}
	
	// weights the whole walk, assume path-model (full order markov)
	protected double weight(Node[] walk)
	{
		double result = 1.0;
		
		result *= Pstart(walk[0]);
		
		// 1st order
		if (!params.nontottering && !params.paths)
		{
			for (int i = 1; i < walk.length; i++)
				result *= Ptransition(walk[i-1], walk[i]);
			result *= Pend(walk[walk.length-1]);
		}
		// 2nd order
		else if (params.nontottering && !params.paths)
		{
			result *= Ptransition(walk[0], walk[1]);
			for (int i = 2; i < walk.length; i++)
				result *= Ptransition(walk[i-2], walk[i-1], walk[i]);
			result *= Pend(walk[walk.length-2], walk[walk.length-1]);
		}
		// full order
		else if (params.paths)
		{
			for (int i = 1; i < walk.length; i++)
				result *= Ptransition(walk, i);
			result *= Pend(walk, walk.length-1);
		}
		
		return result;
	}

	
	// initialize probabilites
	protected void initialize_probabilities(ProductGraph pg)
	{
		double svalue, evalue, tvalue;
		
		// standard ps/pe values
		ps = 1.0 / (pg.getG1().getSize() * pg.getG2().getSize());
		pe = (1.0 - params.lambda) * (1.0 - params.lambda);
		
		// fill with zeros
		ps1 = new ArrayList<Double>(pg.getNodeCount());
		for (int i = 0; i < pg.getNodeCount(); i++) // prefill
			ps1.add(ps);
		
		pe1 = new ArrayList<Double>(pg.getNodeCount());
		for (int i = 0; i < pg.getNodeCount(); i++) // prefill
			pe1.add(pe);

		pe2 = new ArrayList<Map<Integer,Double>>(pg.getNodeCount());
		
		pt1 = new ArrayList<Map<Integer,Double>>(pg.getNodeCount());
		pt2 = new ArrayList<Map<Integer,Map<Integer,Double>>>(pg.getNodeCount());
		
		for (PGNode v1 : pg.getNodes())
		{
			// 1st order start prob
			svalue = Pstart(v1.a1) * Pstart(v1.a2);
			ps1.set(v1.getId(), svalue);
			
			// 1st order end prob
			evalue = Pend(v1.a1) * Pend(v1.a2);
			pe1.set(v1.getId(), evalue);
			
			// make room
			pe2.add(new HashMap<Integer,Double>(v1.getNodeNeighbors().size()));
			pt1.add(new HashMap<Integer,Double>(v1.getNodeNeighbors().size()));
			pt2.add(new HashMap<Integer, Map<Integer,Double>>(v1.getNodeNeighbors().size()));
			
			for (PGNode v2 : v1.getNodeNeighbors())
			{
				// 1st order transition prob
				tvalue = Ptransition(v1.a1, v2.a1) * Ptransition(v1.a2, v2.a2);
				pt1.get(v1.getId()).put(v2.getId(), tvalue);
				
				// 2nd order end prob
				evalue = Pend(v1.a1, v2.a1) * Pend(v1.a2, v2.a2);
				pe2.get(v1.getId()).put(v2.getId(), evalue);
				
				// make room
				pt2.get(v1.getId()).put(v2.getId(), new HashMap<Integer,Double>(v2.getNodeNeighbors().size()));
				
				for (PGNode v3 : v2.getNodeNeighbors())
				{
					// 2nd order transition prob
					tvalue = Ptransition(v1.a1, v2.a1, v3.a1) * Ptransition(v1.a2, v2.a2, v3.a2);
					pt2.get(v1.getId()).get(v2.getId()).put(v3.getId(), tvalue);
				}
			}
		}
	}
	
	protected void computeZ()
	{
		// generates K_k(i,i) kernel values for normalization
		//  ie. Z[i][k] = K_k(i,i) value
		
		Z = new double[graphs.length][params.maxlen+1];
		
		for (Graph g : graphs)
		{
			if (params.kw == KernelWeight.Diffusion)
				diffs.put(g, new Diffusion(g, Math.abs(params.beta)));
			
			ProductGraph pg = new ProductGraph(g, g, params);
			if (params.nontottering)
				Z[g.getIndex()] = nontottering_walk_counts(pg);
			else
				Z[g.getIndex()] = walk_counts(pg);
		}
	}
	
	// tests that probs sum to one and are valid
	protected boolean test_probabilities(ProductGraph pg)
	{
		// start probs
		double sum = 0.0;
		for (Node v : pg.getG1().getNodes())
		{
			sum += Pstart(v);
		}
		assert sum < 1.01 && sum > 0.99;
		
		sum = 0.0;
		for (Node v : pg.getG2().getNodes())
		{
			sum += Pstart(v);
		}
		assert sum < 1.01 && sum > 0.99;

		// transition probs
		for (Node v1 : pg.getG1().getNodes())
		{
			sum = 0.0;
			for (Node v2 : pg.getG1().getNodes())
			{
				sum += Ptransition(v1, v2);
				
				if (Ptransition(v1, v2) > 0.0)
					assert v1.isNeighbor(v2);
			}
			sum += Pend(v1);
			assert sum < 1.01 && sum > 0.99;
		}
		
		for (Node v1 : pg.getG2().getNodes())
		{
			sum = 0.0;
			for (Node v2 : pg.getG2().getNodes())
			{
				sum += Ptransition(v1, v2);
				
				if (Ptransition(v1, v2) > 0.0)
					assert (v1).isNeighbor(v2);
			}
			sum += Pend(v1);
			assert sum < 1.01 && sum > 0.99;
		}
		
		// end probs
		for (Node v : pg.getG1().getNodes())
		{	
			if (Pend(v) == 1.0)
			{	
				sum = 0.0;
				for (Node v2 : v.getNodeNeighbors())
				{
					sum += Ptransition(v,v2);
				}

				assert v.getNodeNeighbors().size() == 0 || sum == 0.0;
			}
		}
		
		for (Node v : pg.getG2().getNodes())
		{
			if (Pend(v) == 1.0)
			{
				sum = 0.0;
				for (Node v2 : v.getNodeNeighbors())
				{
					sum += Ptransition(v,v2);
				}

				assert v.getNodeNeighbors().size() == 0 || sum == 0.0;
			}
		}
		
		// second order stuff
		// transition probs
		for (Node v3 : pg.getG1().getNodes())
		{
			for (Node v1 : pg.getG1().getNodes())
			{
				sum = 0.0;
				for (Node v2 : pg.getG1().getNodes())
				{
					sum += Ptransition(v3, v1, v2);
				
					if (Ptransition(v3, v1, v2) > 0.0)
					{
						assert (v1).isNeighbor(v2);
						assert v3 != v1;
					}
				}
				sum += Pend(v3, v1);
				assert (sum < 1.01 && sum > 0.99) || sum == Pend(v3, v1);
			}
		}
		
		for (Node v3 : pg.getG2().getNodes())
		{
			for (Node v1 : pg.getG2().getNodes())
			{
				sum = 0.0;
				for (Node v2 : pg.getG2().getNodes())
				{
					sum += Ptransition(v3, v1, v2);
				
					if (Ptransition(v3, v1, v2) > 0.0)
					{
						assert (v1).isNeighbor(v2);
						assert v3 != v1;
					}
				}
				sum += Pend(v3, v1);
				assert (sum < 1.01 && sum > 0.99) || sum == Pend(v3, v1);
			}
		}
		
		// end probs
		for (Node v1 : pg.getG1().getNodes())
		{
			for (Node v2 : pg.getG1().getNodes())
			{
				if (Pend( v1, v2) == 1.0)
				{
					sum = 0.0;
					for (Node v3 : v2.getNodeNeighbors())
					{
						sum += Ptransition(v1, v2, v3);
					}

					assert v2.getNodeNeighbors().size() <= 1 || sum == 0.0;
				}
			}
		}
		
		for (Node v1 : pg.getG2().getNodes())
		{
			for (Node v2 : pg.getG2().getNodes())
			{
				if (Pend( v1, v2) == 1.0)
				{
					sum = 0.0;
					for (Node v3 : v2.getNodeNeighbors())
					{
						sum += Ptransition(v1, v2, v3);
					}
	
					assert v2.getNodeNeighbors().size() <= 1 || sum == 0.0;
				}
			}
		}
		
		// product graph stuff

//		for (PGNode v : pg.getNodes())
//		{
//			if (Pend(v) == 1.0)
//				assert v.a1.getNodeNeighbors().size() == 0 && v.a2.getNodeNeighbors().size() == 0;
//			else if (Pend(v) == 1 - params.lambda)
//				assert v.a1.getNodeNeighbors().size() == 0 || v.a2.getNodeNeighbors().size() == 0;
//			else if (Pend(v) == (1-params.lambda)*(1-params.lambda))
//				assert v.a1.getNodeNeighbors().size() > 0 && v.a2.getNodeNeighbors().size() > 0;
//		}
		
		for (PGNode v1 : pg.getNodes())
			for (PGNode v2 : pg.getNodes())
				if (pt1.get(v1.getId()).containsKey(v2) && Ptransition(v1,v2) > 0)
					assert Ptransition(v1,v2) == Ptransition(v1.a1,v2.a1) * Ptransition(v1.a2,v2.a2);
		
		for (PGNode v3 : pg.getNodes())
		{
			for (PGNode v1 : pg.getNodes())
			{
				if (!pt2.get(v3.getId()).containsKey(v1.getId()))
					continue;
				
				for (PGNode v2 : pg.getNodes())
				{
					if (pt2.get(v3.getId()).get(v1.getId()).containsKey(v2.getId()) && Ptransition(v3,v1,v2) > 0)
					{
						// no common stuff
						assert v3.a1 != v2.a1 && v3.a2 != v2.a2;
						// decomposition
						assert Ptransition(v3,v1,v2) == Ptransition(v3.a1,v1.a1,v2.a1) * Ptransition(v3.a2,v1.a2,v2.a2);
					}
				}
			}
		}
		
		return true;
	}
	
	
	protected double walks(ProductGraph pg)
	{
		double sum = 0.0;
		for (Double d : walk_counts(pg))
			sum += d;
		return sum;
	}
	
	protected double[] walk_counts(ProductGraph pg)
	{
		// change the definition of the kernel to:
		// K(G,G') = sum beta_i K_i(G,G'),
		//  where K_i() is kernel for i-len walks
		// we normalize each subkernel by dividing with 1-norm
		// -> K_i(G,G') / sum_i phi(w_i)
		// 
		// we need function double[] levelsums = pg.getG1().computeWalksArray();
		// or more generally we need a matrix Z with Z(i,j) = ||phi_i(G_j)||_1, i.e.
		//   the 'i'-len walk feature sum of reaction 'j'
		
		DecimalFormat df = new DecimalFormat("0.00000", new DecimalFormatSymbols(Locale.US));
		
		initialize_probabilities(pg);

		assert test_probabilities(pg);
		
		// use dynamic programming with n * k matrix
		double[][] F = new double[pg.getNodeCount()][params.maxlen+1];
		double[] res = new double[params.maxlen+1];
		double levelsum = 0.0;
		double sum = 0.0;
		
		// initialize
		for (PGNode v : pg.getNodes())
		{
			F[v.getId()][0] = 0.0;
			F[v.getId()][1] = Pstart(v) * Pend(v);
			levelsum += F[v.getId()][1];
		}
		
		// divide by partial normalization
		if (params.partialnorm)
			levelsum /= Math.sqrt(Z[pg.getG1().getIndex()][1] * Z[pg.getG2().getIndex()][1]);
		
		res[1] = levelsum;
		sum += res[1];
		
		for (int l = 2; l <= params.maxlen; l++)
		{
			levelsum = 0.0;
			for (PGNode v : pg.getNodes())
			{
				double val = 0.0;
				for (PGNode u : v.getNodeNeighbors())
				{
					val += F[u.getId()][l-1] * Ptransition(u,v);
				}
				
				F[v.getId()][l] = val;
				levelsum += F[v.getId()][l];
			}

			if (params.partialnorm)
				levelsum /= Math.sqrt(Z[pg.getG1().getIndex()][l] * Z[pg.getG2().getIndex()][l]);
			
			res[l] = levelsum;
			
			double increase = res[l] / sum; // in percentages / 100
			
			sum += res[l];
			
			// stop if change in score is negligible -> convergence occurred
			if (increase < params.epsilon)
				break;
		}
		
		return res;
	}
	
	protected double nontottering_walks(ProductGraph pg)
	{
		double sum = 0.0;
		for (Double d : nontottering_walk_counts(pg))
			sum += d;
		return sum;
	}
	
	protected double[] nontottering_walk_counts(ProductGraph pg)
	{
		// lets use messages with second order markov model
		// 
		// each edge has a message (i -> j)_k which means the k-length walks
		//  ending at 'i' without anything coming from 'j'
		// 
		// The equations are:
		// 
		// F_k(i) = \sum_j M_{k-1}(j -> i) Pt(j -> i)
		// F_1(i) = PsPe
		// 
		// M_{1}(i -> j) = Ps(i) Pt(i -> j)
		// M_{k}(i -> j) = \sum_s M_{k-1}(s -> i) Pt(s -> i -> j)
		//
		
		DecimalFormat df = new DecimalFormat("0.000", new DecimalFormatSymbols(Locale.US));
		
		initialize_probabilities(pg);
		
		assert test_probabilities(pg);
		
		int EDGES = pg.getEdgeCount();
		
		double[][] F = new double[pg.getNodeCount()][params.maxlen+1];
		double[][] M = new double[pg.getEdgeCount()*2][params.maxlen+1]; // both directions
		
		double totalmass, remainmass, levelmass, contmass, ratio;
		double[] res = new double[params.maxlen+1];
		double sum = 0.0;
		double levelsum = 0.0;
		
		// initialize cases 0 and 1 of probability sums (F)
		for (PGNode v : pg.getNodes())
		{
			F[v.getId()][0] = 0.0;
			F[v.getId()][1] = Pstart(v) * Pend(v);
			
			levelsum += F[v.getId()][1];
		}
		
		if (params.partialnorm)
			levelsum /= Math.sqrt(Z[pg.getG1().getIndex()][1] * Z[pg.getG2().getIndex()][1]);
		
		res[1] = levelsum;
		sum += levelsum;
				
		// initialize trivial cases for messages (M)
		for (PGEdge e : pg.getEdges())
		{
			PGNode src = e.getSource();
			PGNode tgt = e.getTarget();
			
			M[e.getId()][0] = 0.0;
			M[e.getId()][1] = Pstart(src) * Ptransition(src, tgt);
			M[e.getId()+EDGES][0] = 0.0;
			M[e.getId()+EDGES][1] = Pstart(tgt) * Ptransition(tgt, src);
		}
		
		levelmass = levelsum;
		totalmass = sum;
		remainmass = 1.0 - sum;
		ratio = levelmass / remainmass;
//		System.out.println("sum for level 1 is " + df.format(levelmass) + " of total " + df.format(totalmass) + " (" + df.format(ratio) + "% of remaining mass " + df.format(remainmass));

		for (int l = 2; l <= params.maxlen; l++)
		{
			levelsum = 0.0;
			for (PGNode v : pg.getNodes())
			{
				double val = 0.0;
				
				for (PGNode u : v.getNodeNeighbors())
				{
					PGEdge e = u.getEdge(v);
					
					int eid = e.getId();
					if (e.getSource() == v)
						eid += EDGES;
					
					val += M[eid][l-1] * Pend(u,v);
					
//					int dt = (int)Math.round(0.90/Ptransition(u,v));
					int dt = count;
					int de = (int)Math.round(1.00/Pend(u,v));
					
//					System.out.println("" + (u.a1.getId())+(u.a2.getId()) + " -> " + v.a1.getId()+v.a2.getId() + ": Pt=1/" + dt + " Pe=1/" + de + " M=" + M[eid][l-1]);
				}
				
//				System.out.println("node " + v + " gets " + val + " prev:" + F[v.getId()][l-1]);
				
				F[v.getId()][l] = val;
				levelsum += val;
			}
			
			if (params.partialnorm)
				levelsum /= Math.sqrt(Z[pg.getG1().getIndex()][1] * Z[pg.getG2().getIndex()][1]);

			res[l] = levelsum;
			
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

					val += M[eid][l-1] * Ptransition(u,src,tgt);
				}
				
				M[e.getId()][l] = val;
				
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

					val += M[eid][l-1] * Ptransition(u,tgt,src);
				}
				
				M[e.getId()+EDGES][l] = val;
			}
			
			levelmass = levelsum;
			totalmass = sum+levelsum;
			remainmass = 1.0 - sum;
			ratio = levelmass / remainmass;
//			System.out.println("sum for level " + l + " is " + df.format(levelmass) + " of total " + df.format(totalmass) + " (" + df.format(ratio) + "% of remaining mass " + df.format(remainmass));
			
			double increase = res[l] / sum; // in percentages / 100
			
			sum += res[l];
			
			// stop if change in score is negligible -> convergence occurred
			if (increase < params.epsilon)
				break;
		}
		
		return res;
	}
}


