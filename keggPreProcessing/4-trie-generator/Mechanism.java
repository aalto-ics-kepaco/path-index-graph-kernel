//package mechanism;

/*
 * Mechanism.java - a mechanism kernel computation
 * 
 * Mechanism kernel is a reaction kernel (based on KEGG molecule/reaction definitions),
 * where common walks are counted and weighted by distance decay
 * 
 * 
 * The program takes as input a set of reactions defined as MOL-files, a label-matching matrix,
 * some parameters and computer the kernel
 * 
 */


import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;

import mechanism.*;
import mechanism.graphs.Graph;
import mechanism.graphs.MoleculeGraph;
import mechanism.graphs.RGKGraph;
import mechanism.graphs.ReactionGraph;
import mechanism.kernels.*;

public class Mechanism
{
	// Starting point of the program
	public static void main(String[] args) throws IOException
	{
		System.out.println("Parsing arguments..");
		
		// construct argument string
		StringBuilder argsb = new StringBuilder("");
		for (String s : args)
			argsb.append(s + " ");

		String argstr = argsb.toString();
		
		// variables with default values
		double lambda = 0.90;
		double alpha = 1.0;
		double beta = Double.MAX_VALUE;
		double epsilon = 0.0001;
		boolean normalize = false;
		boolean partnorm = false;
		boolean batch = false;
		boolean nodematch = true;
		boolean edgematch = false;
//		String substmatrix = null;
		boolean walks = true;
		boolean paths = false;
		boolean nontottering = false;
		boolean reduced = false;
		boolean moleculegraph = false;
		int maxlen = 20;
		int start = 1;
		int end = Integer.MAX_VALUE;
		KernelType type = KernelType.MMECH;
		KernelParams params = new KernelParams();
		KernelOperationType op = KernelOperationType.DotProduct;
		KernelWeight kw = KernelWeight.Exponential;
		String outputdir = "./";

		// want to compute different types of kernels on reactions:
		// - sucky kernels (reactant-matching, tsuda's kernel, etc
		// - simple walk kernel
		// - exponential kernel
		// - (diffusion kernel)
		// - mechanism kernel (marginal)
		// - mechanism kernel (enumerative)
		
		// help
		if (argstr.indexOf("-h") >= 0 || argstr.indexOf("--help") >= 0)
		{
			// print help
			System.out.println("Mechanism kernel computer\n"
							+ "java Mechanism [-opts] REACTIONGRAPH-MOLFILES\n\n"
							+ "Arguments to the program are:\n"
							+ " Kernel selection\n"
							+ " -t KERNEL        - kernel type [params]\n"
							+ "    WK     gartner  exponential walk kernel [lambda, k]\n"
							+ "    RWK    kashima  marginal walk kernel [lambda, k]\n"
							+ "    SG              subgraphs kernel [k]\n"
							+ "    SP     borgward shortest-paths kernel\n"
							+ "    DOR    rousu's  difference-of-reactants [k]\n"
							+ "    SOR    rousu's  sum-of-reactants [k]\n"
							+ "    RM     rousu's  reactant-matching [k]\n"
							+ "    RGK    tsuda's  reaction graph kernel\n"
							+ "    MMECH  heino's  marginal mechanism kernel [lambda, k]\n"
							+ "    EMECH  heino's  enumerative mechanism kernel [lambda, k]\n"
							+ "    MG              molecule kernel [k]\n"
							+ " Dot product type (only used for SG)\n"
							+ " -y OPERATION     - kernel feature operation type\n"
							+ "    DOT             standard dot product [default]\n"
							+ "    IND             boolean 0/1 indicator\n"
							+ "    MIN             minimum (count of shared features)\n"
							+ "    MINNORM         length-normalized minimum\n"
							+ " Weighting function\n"
							+ " -f FUNCTION      - weight function [default=exp]\n"
							+ "    EXP             exponential attraction [alpha]\n"
							+ "    LOGISTIC        logistic attraction [beta]\n"
							+ "    DIFFUSION       diffusion attraction [beta]\n"
							+ " Model parameters\n"
							+ " -l --lambda      - lambda parameter for length decay [default=" + lambda + "]\n"
							+ " -a --alpha       - alpha parameter for random walk core attraction [default=" + alpha + "]\n"
							+ " -b --beta        - beta parameter for diffusion core attraction [default=" + beta + "]\n"
							+ " -e --epsilon     - epsilon convergence condition [default=" + epsilon + "]\n"
							+ " -k --maxlen      - max walk/sg size [default=" + maxlen + "]\n"
							+ " Other parameters\n"
							+ " -w --walks       - walks [default]\n"
							+ " -g               - non-tottering walks\n"
							+ " -p --paths       - use paths instead of walks\n"
							+ " -n --normalize   - normalize kernel matrix\n"
							+ " -m --mol         - molecular inputs\n"
							+ " -pn              - partial normalization\n"
							+ " -r --reduced     - use reduced product graphs\n"
							+ "    --start       - first index\n"
							+ "    --end         - last index\n"
							+ " -o dir           - output dir\n"
							+ " -h --help        - this help");

			System.exit(0);
		}
		
		// arguments
		
		if (argstr.indexOf("-t") >= 0)
		{
			String m = getStrParam(argstr, "-t");
			
			if (m.equals("RM"))
				type = KernelType.RM;
			else if (m.equals("SOR"))
				type = KernelType.SOR;
			else if (m.equals("DOR"))
				type = KernelType.DOR;
			else if (m.equals("SG"))
				type = KernelType.SG;
			else if (m.equals("RGK"))
				type = KernelType.RGK;
			else if (m.equals("WK"))
				type = KernelType.WK;
			else if (m.equals("RWK"))
				type = KernelType.RWK;
			else if (m.equals("MMECH"))
				type = KernelType.MMECH;
			else if (m.equals("EMECH"))
				type = KernelType.EMECH;
			else if (m.equals("MG"))
				type = KernelType.MG;
			else if (m.equals("SP"))
				type = KernelType.SP;
		}
		
		if (argstr.indexOf("-y ") >= 0)
		{
			String m = getStrParam(argstr, "-y");
			
			if (m.equals("DOT"))
				op = KernelOperationType.DotProduct;
			else if (m.equals("IND"))
				op = KernelOperationType.Indicator;
			else if (m.equals("MIN"))
				op = KernelOperationType.Min;
			else if (m.equals("MINNORM"))
				op = KernelOperationType.MinNormalized;
		}
		
		if (argstr.indexOf("-f ") >= 0)
		{
			String m = getStrParam(argstr, "-f");
			
			if (m.equals("EXP"))
				kw = KernelWeight.Exponential;
			else if (m.equals("LOGISTIC"))
				kw = KernelWeight.Logistic;
			else if (m.equals("DIFFUSION"))
				kw = KernelWeight.Diffusion;
		}		
		
		if (argstr.indexOf("-l") >= 0)
			lambda = getDoubleParam(argstr, "-l");
		else if (argstr.indexOf("--lambda") >= 0)
			lambda = getDoubleParam(argstr, "--lambda");
		
		if (argstr.indexOf("-b ") >= 0)
			beta = getDoubleParam(argstr, "-b");
		else if (argstr.indexOf("--beta") >= 0)
			beta = getDoubleParam(argstr, "--beta");

		if (argstr.indexOf("-a") >= 0)
			alpha = getDoubleParam(argstr, "-a");
		else if (argstr.indexOf("--alpha") >= 0)
			alpha = getDoubleParam(argstr, "--alpha");
		
		if (argstr.indexOf("-e ") >= 0)
			epsilon = getDoubleParam(argstr, "-e");
		else if (argstr.indexOf("--epsilon") >= 0)
			epsilon = getDoubleParam(argstr, "--epsilon");
		
		if (argstr.indexOf("-k") >= 0)
			maxlen = getIntParam(argstr, "-k");
		
		normalize = argstr.indexOf("-n ") >= 0 || argstr.indexOf("--normalize") >= 0;
		partnorm = argstr.indexOf("-pn ") >= 0;
		paths = argstr.indexOf("-p ") >= 0 || argstr.indexOf("--paths") >= 0;
		walks = argstr.indexOf("-w ") >= 0 || argstr.indexOf("--walks") >= 0;
		reduced = argstr.indexOf("-r ") >= 0 || argstr.indexOf("--reduced") >= 0;
		nontottering = argstr.indexOf("-g ") >= 0;
		moleculegraph = argstr.indexOf("-m ") >= 0 || argstr.indexOf("--mol") >= 0;

		
		
		if (argstr.indexOf("-o") >= 0)
			outputdir = getStrParam(argstr, "-o");
		
		if (argstr.indexOf("--start") >= 0)
			start = getIntParam(argstr, "--start");
		if (argstr.indexOf("--end") >= 0)
			end = getIntParam(argstr, "--end");

		
		// normalization requires square matrix
		if (normalize)
			start = 1;
		
		////////////////////////////////////////////////////////////////////////////
		//
		// Params read
		//
		
		params.alpha = alpha;
		params.beta = beta;
		params.lambda = lambda;
		params.epsilon = epsilon;
		params.maxlen = maxlen;
		params.walks = walks;
		params.nontottering = nontottering;
		params.paths = paths;
		params.start = start;
		params.reduced = reduced;
		params.normalize = normalize;
		params.partialnorm = partnorm;
		params.op = op;
		params.kw = kw;
		
		Graph[] graphs;
		
		// Read reaction file names
		List<String> files = new ArrayList<String>(args.length);
		for (String s : args)
			if (new File(s).isFile())
				files.add(s);
		
		// always use sorted indices
		Collections.sort(files);
		
//		System.out.println(files);
				
			
		// update 'end'
		end = Math.min(end, files.size());
		params.end = end;

		if (type == KernelType.MG)
			params.end = files.size(); // magic constant
			
		System.out.println("Parameter string code: " + params.toString());
		System.out.println("Reading reaction graphs 1.." + params.end + " of " + files.size() + " reaction graphs.."); 
		
		graphs = new Graph[params.end];
		for (int i = 0; i < params.end; i++)
		{
			if (type == KernelType.RGK)
				graphs[i] = new RGKGraph(files.get(i));
			else if (!moleculegraph)
				graphs[i] = new ReactionGraph(files.get(i));
			else
				graphs[i] = new MoleculeGraph(files.get(i));
			
			graphs[i].setIndex(i);
		}

		
		System.out.println("Precomputing core distances...");
		// Precompute the distances
		for (Graph g : graphs)
			g.computeDistances();
		
		System.out.println("Computing rows " + start+".."+end+" from [" + graphs.length + " x " + graphs.length + "] kernel matrix (lower triangle only, approx " + ((end-start+1)*end)/2 + " cells)");
		
		Kernel k = null;
		if (type == KernelType.MMECH)     // heinonen's marginal and enumerative
			k = new MarginalMechanismKernel(graphs, params);
		else if (type == KernelType.EMECH)
			k = new EnumerativeMechanismKernel(graphs, params);
		else if (type == KernelType.WK)  // basic exp (enumerative)
		{
			// unweighted enumerative mechanism kernel
			params.kw = KernelWeight.Exponential;
			params.alpha = 1.0;
			k = new EnumerativeMechanismKernel(graphs, params);
		}
		else if (type == KernelType.RWK)  // basic random walk kernel (marginal)
		{
			// unweighted marginal mechanism kernel
			params.kw = KernelWeight.Exponential;
			params.alpha = 1.0;
			k = new MarginalMechanismKernel(graphs, params);
		}
		else if (type == KernelType.SOR)  // rousu's SoR, DoR and RM
			k = new SumOfReactantsKernel(graphs, params);
		else if (type == KernelType.DOR)
			k = new DifferenceOfReactantsKernel(graphs, params);
		else if (type == KernelType.RM)
			k = new ReactantMatchingKernel(graphs, params);
		else if (type == KernelType.RGK)  // tsuda's reaction kernel
			k = new ReactionGraphKernel(graphs, params);
		else if (type == KernelType.MG)   // molecule kernel
		{
			k = new SubgraphsKernel(graphs, params);
			MoleculeGraph[] mgraphs = k.parseSubstrates();
			
			params.start = start;
			params.end = end;
			k = new SubgraphsKernel(mgraphs, params);
		}
		else if (type == KernelType.SG)   // subgraph kernel
			k = new SubgraphsKernel(graphs, params);
		else if (type == KernelType.SP)
			k = new ShortestPathsKernel(graphs, params);
		
		
		// batch mode
		if (batch)
		{
			batch(params, graphs, files);
			return;
		}
		
		// normal mode of operation
		// Compute kernels
		k.compute();
		System.gc();

		if (normalize)
		{
			System.out.println("Normalizing kernel matrix...");
			k.normalize();
		}
		
		if (end == files.size() && start == 1)
		{
			System.out.println("Writing matrix...");
			k.writeToFile(outputdir + "kernel-" + type + "-" + params + "-full.txt");
		}
		else
		{
			System.out.println("Writing partial matrix...");
			
			k.writePartialToFile(outputdir + "kernel-" + type + "-" + params + "-" + start + "-" + end + ".txt");
		}
		
		System.out.println("Done");
	}
	
	private static void batch(KernelParams params, Graph[] graphs, List<String> files)
	{
		DecimalFormat df = new DecimalFormat("0.000", new DecimalFormatSymbols(Locale.US));
		DecimalFormat df2 = new DecimalFormat("0.###", new DecimalFormatSymbols(Locale.US));
		Kernel k = null;
		
		// for a marginal kerneltype, compute following
		// - {exp,logistic,diff} x {pg, reduced} x {walks,nontot,paths}
		//
		// -> gets 3x2x3 = 18 datasets with each having a set of parameters
		//
		// ENUMERATIVE
		// - {pg,reduced} x {walks,nontot,paths}
		// -> 2x3=6 datasets with parameter sigma
		// RGK
		//- 1 dataset with lambda 0.9
		// DOR/SOR/RM
		// SG
		// WK (basically enumerative mech with sigma=1.0)
		//
		//
		//
		// thus -> 6 exp marginal [k=15, lambda=0.9, reduced|w|nt|p, alphas=...]
		//         6 log marginal [k=15, lambda=0.9, reduced|w|nt|p, betas=...]
		//         6 diff marginal [k=15, lambda=0.9, reduced|w|nt|p, betas=...]
		//         6 exp enumerative [k=15, lambda=0.9, reduced|w|nt|p, alphas=...]
		//         6 log enumerative [k=15, lambda=0.9, reduced|w|nt|p, betas=...]
		//         6 diff enumerative [k=15, lambda=0.9, reduced|w|nt|p, betas=...]
		//         4 SOR [k=10, dp|min|minnorm|ind] (no params)
		//         4 DOR [k=10, dp|min|minnorm|ind] (no params)
		//         4 RM [k=10, dp|min|minnorm|ind] (no params)
		//         4 SG [k=10, dp|min|minnorm|ind] (no params)
		//         1 RGK [lambda=0.9] (no params)
		//         1 WK  [lambda=0.9, em] (no params)
		//         1 RWK [lambda=0.9, mm] (no params)
		

		int cc = 36; // configcount
		int SORi = cc;
		int DORi = cc+4;
		int RMi = cc+8;
		int SGi = cc+12;
		int RGKi = cc+16;

		
		int zeroindex;
		
		// {reduced,walks,nontot,paths} quintuples
		boolean[][] fields = {{false,true,false,false},
	                          {false,false,true,false},
	                          {false,false,false,true},
	                          {true,true,false,false},
	                          {true,false,true,false},
	                          {true,false,false,true}};
		
		// generate values
		List<Double> values = new ArrayList<Double>();
		values.add(-1000.0);
		values.add(-100.0);
		values.add(-50.0);
		values.add(-20.0);
		values.add(-10.0);
		for (double d = -5.0; d < -0.01; d += 0.10)
			values.add(d);
		values.add(-0.01);
		values.add(0.0);
		values.add(0.01);
		for (double d = 0.1; d < 5.0; d += 0.10)
			values.add(d);
		values.add(10.0);
		values.add(20.0);
		values.add(50.0);
		values.add(100.0);
		values.add(1000.0);
		
		zeroindex = values.size() / 2;
		
		double[][] results = new double[values.size()][cc+19];
		
		int i = 0;
		for (i = 0; i < cc; i++)
		{
			params.maxlen = 15;
			params.lambda = 0.9;
			params.reduced = fields[i%6][0];
			params.walks = fields[i%6][1];
			params.nontottering = fields[i%6][2];
			params.paths = fields[i%6][3];
			
			for (int j = 0; j < values.size(); j++)
			{
				// set both, only use either one
				params.beta = values.get(j);
				params.alpha = values.get(j);
				
				if (i/6 == 0)
				{
					params.kw = KernelWeight.Exponential;
					k = new MarginalMechanismKernel(graphs, params);
				}
				if (i/6 == 1)
				{
					params.kw = KernelWeight.Logistic;
					k = new MarginalMechanismKernel(graphs, params);
				}
				if (i/6 == 2)
				{
					params.kw = KernelWeight.Diffusion;
					k = new MarginalMechanismKernel(graphs, params);
				}
				if (i/6 == 3)
				{
					params.kw = KernelWeight.Exponential;
					k = new EnumerativeMechanismKernel(graphs, params);
				}
				if (i/6 == 4)
				{
					params.kw = KernelWeight.Logistic;
					k = new EnumerativeMechanismKernel(graphs, params);
				}
				if (i/6 == 5)
				{
					params.kw = KernelWeight.Diffusion;
					k = new EnumerativeMechanismKernel(graphs, params);
				}
				
				k.compute();
				k.normalize();

				results[j][i] = k.getValue(0,1);
			}
			System.out.println(i + k.getClass().toString());
		}

		
		params.maxlen = 10;
		
		// SOR
		for (; i < SORi+4; i++)
		{
			if (i == SORi)
				params.op = KernelOperationType.DotProduct;
			if (i == SORi+1)
				params.op = KernelOperationType.Min;
			if (i == SORi+2)
				params.op = KernelOperationType.MinNormalized;
			if (i == SORi+3)
				params.op = KernelOperationType.Indicator;
			
			k = new SumOfReactantsKernel(graphs, params);
			k.compute();
			k.normalize();
			results[zeroindex][i] = k.getValue(0,1);
			System.out.println(i + k.getClass().toString());
		}
		
		// DOR
		for (; i < DORi+4; i++)
		{
			if (i == DORi)
				params.op = KernelOperationType.DotProduct;
			if (i == DORi+1)
				params.op = KernelOperationType.Min;
			if (i == DORi+2)
				params.op = KernelOperationType.MinNormalized;
			if (i == DORi+3)
				params.op = KernelOperationType.Indicator;
			
			k = new DifferenceOfReactantsKernel(graphs, params);
			k.compute();
			k.normalize();
			results[zeroindex][i] = k.getValue(0,1);
			System.out.println(i + k.getClass().toString());
		}
		
		// RM
		for (; i < RMi+4; i++)
		{
			if (i == RMi)
				params.op = KernelOperationType.DotProduct;
			if (i == RMi+1)
				params.op = KernelOperationType.Min;
			if (i == RMi+2)
				params.op = KernelOperationType.MinNormalized;
			if (i == RMi+3)
				params.op = KernelOperationType.Indicator;
			
			k = new ReactantMatchingKernel(graphs, params);
			k.compute();
			k.normalize();
			results[zeroindex][i] = k.getValue(0,1);
			System.out.println(i + k.getClass().toString());
		}
		
		// SG
		for (; i < SGi+4; i++)
		{
			if (i == SGi)
				params.op = KernelOperationType.DotProduct;
			if (i == SGi+1)
				params.op = KernelOperationType.Min;
			if (i == SGi+2)
				params.op = KernelOperationType.MinNormalized;
			if (i == SGi+3)
				params.op = KernelOperationType.Indicator;
			
			k = new SubgraphsKernel(graphs, params);
			k.compute();
			k.normalize();
			results[zeroindex][i] = k.getValue(0,1);
			System.out.println(i + k.getClass().toString());
		}
		
		// WK
		params.maxlen = 15;
		params.lambda = 0.9;
		params.alpha = 1.0;
		params.kw = KernelWeight.Exponential;
		params.reduced = false;
		params.walks = true;
		params.nontottering = false;
		params.paths = false;
		params.nodematch = true;
		params.edgematch = false;
		k = new EnumerativeMechanismKernel(graphs, params);
		k.compute();
		k.normalize();
		results[zeroindex][i] = k.getValue(0,1);
		System.out.println(i + k.getClass().toString());

		i++;
		
		// RWK
		k = new MarginalMechanismKernel(graphs, params);
		k.compute();
		k.normalize();
		results[zeroindex][i] = k.getValue(0,1);
		System.out.println(i + k.getClass().toString());
		
		i++;
		
		// RGK
		for (int j = 0; j < files.size(); j++)
		{
			graphs[j] = new RGKGraph(files.get(j));
		}
		params.maxlen = 50;
		params.lambda = 0.9;
		k = new ReactionGraphKernel(graphs, params);
		k.compute();
		k.normalize();
		results[zeroindex][i] = k.getValue(0,1);
		System.out.println(i + k.getClass().toString());
	
		
		
		
		// write into gnuplottable file
		BufferedWriter out;
		try 
		{
			out = new BufferedWriter(new FileWriter("marginals.data"));
			
			out.write("par");
			out.write("\tem:w\tem:nt\tem:p\tem:wr\tem:ntr\tem:pr");
			out.write("\tlm:w\tlm:nt\tlm:p\tlm:wr\tlm:ntr\tlm:pr");
			out.write("\tdm:w\tdm:nt\tdm:p\tdm:wr\tdm:ntr\tdm:pr");
			out.write("\tee:w\tee:nt\tee:p\tee:wr\tee:ntr\tee:pr");
			out.write("\tle:w\tle:nt\tle:p\tle:wr\tle:ntr\tle:pr");
			out.write("\tmle:w\tme:nt\tme:p\tme:wr\tme:ntr\tme:pr");
			out.write("\tsor:dp\tsor:m\tsor:mn\tsor:i");
			out.write("\tdor:dp\tdor:m\tdor:mn\tdor:i");
			out.write("\trm:dp\trm:m\trm:mn\trm:i");
			out.write("\tsg:dp\tsg:m\tsg:mn\tsg:i");
			out.write("\twk\trwk\trgk");
			out.newLine();
			
			for (i = 0; i < values.size(); i++)
			{
				String s = "";
				
				s += df2.format(values.get(i)) + "\t";
				for (int j = 0; j < cc+19; j++)
				{
					if (results[i][j] == 0.0)
						s += "\t";
					else
						s += df.format(results[i][j]) + "\t";
				}
				
				s = s.trim();
				
				out.write(s);
				out.newLine();
			}
			
			out.close();
		} 
		catch (IOException e)
		{
			System.out.println(e.getMessage());
		}
	}
	
	
	public static double getDoubleParam(String arg, String param)
	{
		double st = -1;
		try
		{
			int begin = arg.indexOf(param) + param.length() + 1;
			int stop = arg.indexOf(" ", begin);
			st = Double.parseDouble(arg.substring(begin, stop).trim());
		} catch (NumberFormatException e)
		{
			System.out.println("Awkward argument for param " + arg + ", " + e);
			System.exit(0);
		}
		return st;
	}

	public static int getIntParam(String arg, String param)
	{
		int st = 0;
		try
		{
			int begin = arg.indexOf(param) + param.length() + 1;
			int stop = arg.indexOf(" ", begin);
			st = Integer.parseInt(arg.substring(begin, stop).trim());
		} catch (NumberFormatException e)
		{
			System.out.println("Awkward argument for param " + arg + ", " + e);
			System.exit(0);
		}
		return st;
	}

	public static long getLongParam(String arg, String param)
	{
		long st = 0;
		try
		{
			int begin = arg.indexOf(param) + param.length() + 1;
			int stop = arg.indexOf(" ", begin);
			st = Long.parseLong(arg.substring(begin, stop).trim());
		} catch (NumberFormatException e)
		{
			System.out.println("Awkward argument for param " + arg + ", " + e);
			System.exit(0);
		}
		return st;
	}

	public static String getStrParam(String arg, String param)
	{
		String st = null;
		try
		{
			int begin = arg.indexOf(param) + param.length() + 1;
			int stop = arg.indexOf(" ", begin);
			st = arg.substring(begin, stop).trim();
		} catch (NumberFormatException e)
		{
			System.out.println("Awkward argument for param " + arg + ", " + e);
			System.exit(0);
		}
		return st;
	}
	


}
