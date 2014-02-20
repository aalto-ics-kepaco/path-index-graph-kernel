package mechanism.kernels;

import java.text.DecimalFormat;

import mechanism.*;
import mechanism.graphs.Graph;
import mechanism.graphs.MoleculeGraph;
import mechanism.graphs.ReactionGraph;

public class ReactantMatchingKernel extends Kernel
{
	// molecule-kernel matrix
	private Kernel moleculekernel;
	
	public ReactantMatchingKernel(Graph[] graphs, KernelParams params)
	{
		super(graphs, params);
	}

	public void compute()
	{
		System.out.println("Computing molecule sg kernel first");
		MoleculeGraph[] substrates = parseSubstrates();
		System.out.println(substrates.length + " substrates");
		
		// compute molecular kernel
		KernelParams sgparams = params.clone();
		sgparams.start = 1;
		sgparams.end = substrates.length;
//		sgparams.maxlen = 10;
		sgparams.op = KernelOperationType.DotProduct;
		moleculekernel = new SubgraphsKernel(substrates, sgparams);
		moleculekernel.compute();
		moleculekernel.normalize();
		
		System.out.println("Computing RM kernel");
		
		super.compute();
	}

	// use uni-directional version, where we sum the substrates and products separately
	public double compute(Graph g1, Graph g2)
	{
		// cast to reactiongraphs first
		ReactionGraph rg1 = (ReactionGraph)g1;
		ReactionGraph rg2 = (ReactionGraph)g2;
		
		double klhs = 0.0;
		double krhs = 0.0;
		
		// all reactants against all multiplied by all products against all
		for (MoleculeGraph m1 : rg1.getReactants())
			for (MoleculeGraph m2 : rg2.getReactants())
				klhs += moleculekernel.getValue(m1.getIndex(),m2.getIndex());		
		for (MoleculeGraph m1 : rg1.getProducts())
			for (MoleculeGraph m2 : rg2.getProducts())
				krhs += moleculekernel.getValue(m1.getIndex(),m2.getIndex());		
		
		return klhs*krhs;
	}
	
}
