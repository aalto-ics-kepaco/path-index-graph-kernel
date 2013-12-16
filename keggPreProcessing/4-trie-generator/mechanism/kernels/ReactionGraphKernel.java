package mechanism.kernels;

import java.util.*;

import mechanism.*;
import mechanism.graphs.*;



public class ReactionGraphKernel extends MarginalMechanismKernel
{
	/*
	 * inner kernel is a molecule walk kernel where atom AND bond symbols
	 * have to match
	 * 
	 * outer kernel works on the reaction graph as a walk kernel with weight 
	 *  K_z(v,v') = K_wk(v,v')
	 *  
	 *  i.e. each node is a molecule and kernel between nodes is the inner walk kernel
	 * 
	 * in outer kernel the edges have to match
	 * 
	 * edges are 'main', 'leave', 'cofac', 'trans', 'ligase', 'group' (from RPAIR)
	 * 
	 */
	
	// inner level walk kernel for molecules
	private Kernel moleculekernel;
	
	
	public ReactionGraphKernel(Graph[] graphs, KernelParams params)
	{
		super(graphs, params);
		
		// fix parameters to match Tsuda's implementation
		// only free parameter is the lambda
		params.edgematch = true; // match edges both in inner and outer kernel
		params.nodematch = false;
		params.maxlen = 50; // run until convergence, 50 should be enough
	}

	public void compute()
	{
		// precompute the inner walk kernels
		
		// count the number of molecules used
		Set<String> ligands = new HashSet<String>();
		for (int i = 0; i < graphs.length; i++)
		{
			RGKGraph g = (RGKGraph)graphs[i];
			
			for (int j = 0; j < g.getNodes().length; j++)
				ligands.add(g.getNodes()[j].getSymbol());
		}
		
		String[] ligs = new String[ligands.size()];
		ligs = ligands.toArray(ligs);
		Arrays.sort(ligs);
		
		MoleculeGraph[] mols = new MoleculeGraph[ligands.size()];
		
		for (int i = 0; i < ligs.length; i++)
		{
			mols[i] = new MoleculeGraph(Kernel.MOL_FOLDER + ligs[i] + ".mol");
			mols[i].setIndex(i);
			mols[i].computeDistances();
		}
		
		// set MG indices for all reactions to match the 'mols'
		for (Graph g : graphs)
		{
			for (Node n : g.getNodes())
			{
				int id = Arrays.binarySearch(ligs, ((RGKNode)n).molecule.getLigand());
				((RGKNode)n).molecule.setIndex(id);
			}
		}
		
		System.out.println("Precomputing inner molecular walk kernel with " + mols.length + " mols");
		
		// now compute the kernel: uniform marginal with k=50
		KernelParams wkparams = params.clone();
		wkparams.start = 1;
		wkparams.end = mols.length;
		wkparams.walks = true;
		wkparams.nontottering = false;
		wkparams.paths = false;
//		wkparams.lambda = 0.90;
		wkparams.alpha = 1.0;
		wkparams.kw = KernelWeight.Exponential;
		wkparams.normalize = true;
		wkparams.maxlen = 50;
		wkparams.reduced = false;
		wkparams.edgematch = true; // match also bonds
		wkparams.nodematch = true;
		moleculekernel = new MarginalMechanismKernel(mols, wkparams);
		moleculekernel.compute();
		moleculekernel.normalize();
		
		// compute normally from now on
		System.out.println("Computing outer reaction walk kernel..");
		super.compute();
	}

	// parameters are now 'RGKGraph's
	public double compute(Graph g1, Graph g2)
	{
		// use standard walk kernel on top of these two graphs with
		// inner walk kernels as matchings
		
		// construct product graph
		long pgtime = System.currentTimeMillis();
		ProductGraph pg = new ProductGraph(g1, g2, params, false);
		this.pgtime += System.currentTimeMillis() - pgtime;

		// use dynamic programming
		long walktime = System.currentTimeMillis();
		double value = walks(pg);
		this.walktime += System.currentTimeMillis() - walktime;
		
		counted++;
		
		return value;	
	}

	protected double walks(ProductGraph pg)
	{
		initialize_probabilities(pg);
		
		// use dynamic programming with n * k matrix
		double[][] F = new double[pg.getNodeCount()][params.maxlen+1];
		double sum = 0.0;
		double oldsum = 0.0;
		
		// initialize
		for (PGNode v : pg.getNodes())
		{
			F[v.getId()][0] = 0.0;
			F[v.getId()][1] = Pstart(v) * Pend(v) * moleculekernel.getValue(((RGKNode)v.a1).getMolIndex(), ((RGKNode)v.a2).getMolIndex());
			sum += F[v.getId()][1];
		}
		
		for (int l = 2; l <= params.maxlen; l++)
		{
			oldsum = sum;
			for (PGNode v : pg.getNodes())
			{
				double val = 0.0;
				for (PGNode u : v.getNodeNeighbors())
				{
					val += F[u.getId()][l-1] * Ptransition(u,v);
				}
				
				F[v.getId()][l] = val * moleculekernel.getValue(((RGKNode)v.a1).getMolIndex(), ((RGKNode)v.a2).getMolIndex());
				sum += F[v.getId()][l];
			}

			double increase = (sum - oldsum) / oldsum; // in percentages / 100
			
			// stop if change in score is negligible -> convergence occurred
			if (increase < params.epsilon)
				break;
		}
		
		return sum;
	}
	
}
