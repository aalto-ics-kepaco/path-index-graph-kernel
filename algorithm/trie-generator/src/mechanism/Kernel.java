package mechanism;


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.*;

import mechanism.graphs.*;


/*
 * Kernels are symmetric, thus only lower triangle of the matrix is necessary to compute and store
 * Also, for memory-efficiency, floats are used instead of doubles
 */


public abstract class Kernel
{
	public static String MOL_FOLDER = "/group/home/icomic/data/kegg/ligand/LATEST/mol/";
	
	protected double matrix[][] = null;
	protected Graph[] graphs;
	protected int count;
	protected KernelParams params;
	protected int ec;
	
	protected int counted = 0;
	protected long totalpgtime = 0;
	protected long totalwalktime = 0;
	protected long pgtime = 0;
	protected long walktime = 0;
	
	protected Map<Graph, Diffusion> diffs = null;;
	
		
	public Kernel(Graph[] graphs, KernelParams params)
	{
		this.graphs = graphs;
		this.count = graphs.length;
		this.params = params;
		
		ec = params.end - params.start + 1;
		
		// jagged (lower left) triangle matrix
		matrix = new double[ec][];
		for (int i = 0; i < ec; i++)
			matrix[i] = new double[params.start + i];
	}
	
	public void compute()
	{
		
//		Map<String,Integer> spectrum1 = new HashMap<String,Integer>();
//		for (Edge e : graphs[0].getEdges())
//		{
//			String estr = "";
//			if (e.getChangetype() == 1)
//				estr += "+";
//			else if (e.getChangetype() == -1)
//				estr += "-";
//			else if (e.getChangetype() == 0)
//				estr += " ";
//			
//			if (e.getSource().getSymbol().compareTo(e.getTarget().getSymbol()) < 0)
//				estr += e.getSource().getSymbol() + e.getTarget().getSymbol();
//			else
//				estr += e.getTarget().getSymbol() + e.getSource().getSymbol();
//			
//			if (spectrum1.containsKey(estr))
//				spectrum1.put(estr, spectrum1.get(estr) + 1);
//			else
//				spectrum1.put(estr, 1);
//		}
//
//		Map<String,Integer> spectrum2 = new HashMap<String,Integer>();
//		for (Edge e : graphs[1].getEdges())
//		{
//			String estr = "";
//			if (e.getChangetype() == 1)
//				estr += "+";
//			else if (e.getChangetype() == -1)
//				estr += "-";
//			else if (e.getChangetype() == 0)
//				estr += " ";
//			
//			if (e.getSource().getSymbol().compareTo(e.getTarget().getSymbol()) < 0)
//				estr += e.getSource().getSymbol() + e.getTarget().getSymbol();
//			else
//				estr += e.getTarget().getSymbol() + e.getSource().getSymbol();
//			
//			if (spectrum2.containsKey(estr))
//				spectrum2.put(estr, spectrum2.get(estr) + 1);
//			else
//				spectrum2.put(estr, 1);
//		}

		
		
		
		
		long count = 0;
		long tocompute = ((params.end-params.start+1)*params.end)/2;
		
		long starttime = System.currentTimeMillis();
		
		// compute only the lower left triangle and for normalization purposes
		// also [i,i] diagonal
		for (int i = params.start-1; i < params.end; i++)
		{
//			pgtime = 0;
//			walktime = 0;
			
			for (int j = 0; j <= i; j++)
			{
				matrix[i-params.start+1][j] = compute(graphs[i], graphs[j]);
				
				// nan's not allowed at this stage, infinities are ok
				// after normalization nans are ok
				assert !Double.isNaN(matrix[i-params.start+1][j]);
				
				if (Double.isNaN(matrix[i-params.start+1][j]))
					System.out.println("Error: nan at " + i + " " + j);
				
				count++;
			}
			
//			totalpgtime += pgtime;
//			totalwalktime += walktime;
			
			if (graphs[i] instanceof ReactionGraph)
			{
				String dir = "";
				if (graphs[i].getDirection() == 1)
					dir = "_+1";
				else if (graphs[i].getDirection() == -1)
					dir = "_-1";
				
				// also show percentage
				// i.e. [34.3%] 
				System.out.println("row " + (i+1) + "\t" + graphs[i].getLigand() + "_" + ((ReactionGraph)graphs[i]).getMapNum() + dir + ": computed " + getUpdateStr(starttime, 1.0*count/tocompute));
			}
			else
				System.out.println("row " + (i+1) + "\t" + graphs[i].getLigand() + ": computed " + getUpdateStr(starttime, 1.0*count/tocompute));
		}
	}
	
	public abstract double compute(Graph g1, Graph g2);
	
	public double getValue(int i, int j)
	{
		// lower triangle matrix, if the (i,j) pair hits upper right triangle, return (j,i) cell
		if (j < matrix[i].length)
			return matrix[i][j];
		return matrix[j][i];
	}
	
	public void normalize()
	{
		double[][] normalized_matrix = new double[count][];
		for (int i = 0; i < count; i++)
			normalized_matrix[i] = new double[i+1];
		
		for (int i = 0; i < count; i++)
			for (int j = 0; j <= i; j++)
				if (Math.sqrt(matrix[i][i]*matrix[j][j]) != 0.0)
					normalized_matrix[i][j] = matrix[i][j] / Math.sqrt( matrix[i][i] * matrix[j][j] );
		
		matrix = normalized_matrix;
	}
	
	public MoleculeGraph[] parseSubstrates()
	{
		// i have {reac: subs}, i need {sub: reacs}

		// get the set of substrates in all reactions
		Set<String> ligandset = new HashSet<String>();
		for (Graph g : graphs)
		{
			ReactionGraph rg = (ReactionGraph)g;
			ligandset.addAll(rg.getReactantLigands());
			ligandset.addAll(rg.getProductLigands());
		}
		
		// produce a sorted list of substrates
		List<String> ligandlist = new ArrayList<String>(ligandset);
		Collections.sort(ligandlist);
		
		// construct a list of reacs for each substrate
		Map<String,List<ReactionGraph>> reactant_reacs = new HashMap<String,List<ReactionGraph>>(ligandlist.size());
		Map<String,List<ReactionGraph>> product_reacs = new HashMap<String,List<ReactionGraph>>(ligandlist.size());
		for (int i = 0; i < ligandlist.size(); i++)
		{
			reactant_reacs.put(ligandlist.get(i), new ArrayList<ReactionGraph>());
			product_reacs.put(ligandlist.get(i), new ArrayList<ReactionGraph>());
		}
		
		// attach moleculegraphs to reactiongraphs
		for (Graph g : graphs)
		{
			ReactionGraph rg = (ReactionGraph)g;
			for (String s : rg.getReactantLigands())
				reactant_reacs.get(s).add(rg);
			for (String s : rg.getProductLigands())
				product_reacs.get(s).add(rg);
			
		}
		
		// parse all moleculegraphs used in our reaction set
		MoleculeGraph[] substrates = new MoleculeGraph[ligandlist.size()];
		for (int i = 0; i < ligandlist.size(); i++)
		{
			MoleculeGraph mg = new MoleculeGraph(MOL_FOLDER + ligandlist.get(i) + ".mol");
			mg.setIndex(i);
			substrates[i] = mg;
			
			for (ReactionGraph rg : reactant_reacs.get(ligandlist.get(i)))
				rg.addReactant(mg);
			for (ReactionGraph rg : product_reacs.get(ligandlist.get(i)))
				rg.addProduct(mg);
		}
		
		return substrates;
	}
	
	
	public void writeToFile(String filename)
	{
		BufferedWriter out;
		DecimalFormat df = new DecimalFormat("0.#####E0", new DecimalFormatSymbols(Locale.US));
		
		try
		{
			out = new BufferedWriter(new FileWriter(filename));
			
			for (int i = 0; i < count; i++)
			{
				for (int j = 0; j <= i; j++)
				{
					if (Double.isInfinite(matrix[i][j]))
						out.write("Inf");
					else if (Double.isNaN(matrix[i][j]))
						out.write("NaN");
					else
						out.write( df.format(matrix[i][j]));
					
					out.write("\t");
				}
				
				out.newLine();
			}
			
			out.close();
		}
		catch (IOException e)
		{
			System.out.println(e.getMessage());
		}
		
		System.out.println("Kernel written to file " + filename);
	}

	public void writePartialToFile(String filename)
	{
		BufferedWriter out;
		DecimalFormat df = new DecimalFormat("0.#####E0", new DecimalFormatSymbols(Locale.US));
		
		try
		{
			out = new BufferedWriter(new FileWriter(filename));
			
			for (int i = params.start-1; i < params.end; i++)
			{
				for (int j = 0; j <= i; j++)
				{
					if (Double.isInfinite(matrix[i-params.start+1][j]))
						out.write("Inf");
					else if (Double.isNaN(matrix[i-params.start+1][j]))
						out.write("NaN");
					else
						out.write(df.format(matrix[i-params.start+1][j]));
					
					out.write("\t");
				}
				
				out.newLine();
			}
			
			out.close();
		}
		catch (IOException e)
		{
			System.out.println(e.getMessage());
		}
		
		System.out.println("Partial kernel (" + params.start + ".." + params.end + ") written to file " + filename);
	}
	
	public String toString()
	{
		return Arrays.deepToString(matrix);
	}
	
	public void print(int start, int end)
	{
		System.out.print("              ");
		for (int i = start-1; i < end; i++)
			System.out.print(graphs[i].getLigand() + " ");
		System.out.println();
		
		// output directly
		for (int i = start-1; i < end; i++)
		{
			System.out.print(graphs[i].getLigand() + ": [ ");
			for (int j = start-1; j <= i; j++)
				System.out.format(Locale.ENGLISH, "%.3f       ", matrix[i][j]);
			System.out.println("");
		}
	}

    // return phi(x) = standard Gaussian pdf
    public static double phi(double x) {
        return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
    }

    // return phi(x, mu, signma) = Gaussian pdf with mean mu and stddev sigma
    public static double pdf(double x, double mu, double sigma) {
        return phi((x - mu) / sigma) / sigma;
    }

    
	// Exponential model
	protected double exp(double x)
	{
		if (params.alpha > 0.0)
			return Math.pow(params.alpha, -x);
		
		return 0.0;
	}
	
	protected double exp(Node v)
	{
		return exp(v.getCoreDist());
	}
	
	protected double exp(PGNode v)
	{
		return exp(v.a1) * exp(v.a2);
	}
	
	// Logistic model
	protected double logistic(double x)
	{
		return 1.0 / (1.0 + Math.exp(params.beta * (x - 0.5))); // use logistic function instead of gaussian
	}
	
	protected double logistic(Node v)
	{
		return logistic(v.getCoreDist());
	}
	
	protected double logistic(PGNode v)
	{
		return logistic(v.a1) * logistic(v.a2);
	}

	// Diffusion model
	protected double diffusion(PGNode v)
	{
		return diffusion(v.a1) * diffusion(v.a2);
	}
	
	protected double diffusion(Node v)
	{
		if (params.beta >= 0.0) // positive
			return maxdiffvalue(v);
		else // negative
			return 1.0 - maxdiffvalue(v);
	}
	
	// find max diffusion value against reaction core
	protected double maxdiffvalue(Node a)
	{
		Diffusion diff = diffs.get(a.getParent());
		
		double maxval = 0.0;
		for (Node x : a.getParent().getNodes())
		{
			if (x.getCoreDist() == 0 && diff.get(x.getId(), a.getId()) > maxval)
				maxval = diff.get(x.getId(), a.getId());
		}
		
		return maxval;
	}
	
	
	// Lambda model
	protected double lambda(PGNode v)
	{
		return params.lambda;
	}
	
	protected double lambda(Node v)
	{
		return params.lambda;
	}

    
	public static String getUpdateStr(long starttime, double doneratio)
	{
		DecimalFormat df = new DecimalFormat("00.00%", new DecimalFormatSymbols(Locale.US));
		
		long runtime = (long)((System.currentTimeMillis()-starttime)/1000);
		long yettime = (long)((double)runtime / doneratio - runtime);
		
		String runstr = "";
		if (runtime >= (60*60*24))
		{
			runstr += runtime / (60*60*24) + "d ";
			runtime -= (runtime / (60*60*24)) * (60*60*24);
		}
		if (runtime >= (60*60))
		{
			runstr += runtime / (60*60) + "h ";
			runtime -= (runtime / (60*60)) * (60*60);
		}
		if (runtime >= 60)
		{
			runstr += runtime / 60 + "m ";
			runtime -= (runtime / 60) * 60;
		}
		runstr += runtime + "s";
		
		String yetstr = "";
		if (yettime >= (60*60*24))
		{
			yetstr += yettime / (60*60*24) + "d ";
			yettime -= (yettime / (60*60*24)) * (60*60*24);
		}
		if (yettime >= (60*60))
		{
			yetstr += yettime / (60*60) + "h ";
			yettime -= (yettime / (60*60)) * (60*60);
		}
		if (yettime >= 60)
		{
			yetstr += yettime / 60 + "m ";
			yettime -= (yettime / 60) * 60;
		}
		yetstr += yettime + "s";
		
		
		String res = "[" + df.format(doneratio) + " elapsed:" + runstr + " ETA:" + yetstr + "]";
		
		return res;
	}
	
	
}



