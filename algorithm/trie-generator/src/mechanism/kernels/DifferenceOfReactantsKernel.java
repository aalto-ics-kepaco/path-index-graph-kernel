package mechanism.kernels;

import mechanism.*;
import mechanism.graphs.*;

public class DifferenceOfReactantsKernel extends Kernel
{
	// molecule-kernel matrix
	private Kernel moleculekernel;
	
	public DifferenceOfReactantsKernel(Graph[] graphs, KernelParams params)
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

	// use uni-directional version, take positive sum over S-S' and P-P' pairs
	// add negative sums over S-P' and P-S' pairs
	public double compute(Graph g1, Graph g2)
	{
		// cast to reactiongraphs first
		ReactionGraph rg1 = (ReactionGraph)g1;
		ReactionGraph rg2 = (ReactionGraph)g2;
		
		double kpos = 0.0;
		double kneg = 0.0;
		
		// all reactants against all multiplied by all products against all
		for (MoleculeGraph m1 : rg1.getReactants())
			for (MoleculeGraph m2 : rg2.getReactants())
				kpos += moleculekernel.getValue(m1.getIndex(),m2.getIndex());
		for (MoleculeGraph m1 : rg1.getProducts())
			for (MoleculeGraph m2 : rg2.getProducts())
				kpos += moleculekernel.getValue(m1.getIndex(),m2.getIndex());
		for (MoleculeGraph m1 : rg1.getReactants())
			for (MoleculeGraph m2 : rg2.getProducts())
				kneg += moleculekernel.getValue(m1.getIndex(),m2.getIndex());
		for (MoleculeGraph m1 : rg1.getProducts())
			for (MoleculeGraph m2 : rg2.getReactants())
				kneg += moleculekernel.getValue(m1.getIndex(),m2.getIndex());
		
		return kpos - kneg;
	}
}
