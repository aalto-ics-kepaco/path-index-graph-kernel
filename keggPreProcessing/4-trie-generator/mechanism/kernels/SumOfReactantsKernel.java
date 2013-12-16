package mechanism.kernels;

import java.text.DecimalFormat;

import mechanism.*;
import mechanism.graphs.Graph;
import mechanism.graphs.MoleculeGraph;
import mechanism.graphs.ReactionGraph;

public class SumOfReactantsKernel extends Kernel
{
	// sum-of-reactants kernel
	// basically the sum of all substrate-pair kernels [s X s'] -> for 3 sub -> 9
	// the molecule-kernel used is subgraph kernel up to k=10
	
	// molecule-kernel matrix
	private Kernel moleculekernel;
	
	public SumOfReactantsKernel(Graph[] graphs, KernelParams params)
	{
		super(graphs, params);
	}

	public void compute()
	{
//		System.out.println("Computing molecule kernel first");
		MoleculeGraph[] substrates = parseSubstrates();
//		System.out.println(substrates.length + " substrates");
		
		// compute molecular kernel
		KernelParams sgparams = params.clone();
		sgparams.start = 1;
		sgparams.end = substrates.length;
		moleculekernel = new SubgraphsKernel(substrates, sgparams);
		moleculekernel.compute();
		moleculekernel.normalize();

		super.compute();
	}

	public double compute(Graph g1, Graph g2)
	{
		// cast to reactiongraphs first
		ReactionGraph rg1 = (ReactionGraph)g1;
		ReactionGraph rg2 = (ReactionGraph)g2;
		
		double result = 0.0;
		
		// all substrates against all
		for (MoleculeGraph m1 : rg1.getSubstrates())
			for (MoleculeGraph m2 : rg2.getSubstrates())
				result += moleculekernel.getValue(m1.getIndex(),m2.getIndex());		
		
		return result;
	}
}




