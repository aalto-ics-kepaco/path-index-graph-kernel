package mechanism;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

public class KernelParams
{
	public double alpha, beta, epsilon, lambda;
	public int maxlen, start, end;
	public boolean walks, nontottering, paths, reduced, normalize, partialnorm, nodematch, edgematch;
	public KernelOperationType op;
	public KernelWeight kw;
	
	
	public String toString()
	{
		// format for resulting files is:
		// kernel-TYPE-K-LAMBDA-ALPHA-BETA-SIGMA-NORMALIZE-TOTTERING-REDUCED-OPERATION-START-END.txt
		//
		// where NORMALIZE = "normalized" or "raw"
		//       TOTTERING = "notottering" or "tottering"
		//       REDUCED = "reducedpg" or "plainpg"
		//       OPERATION = dot | indicator | min | minnorm

		DecimalFormat df = new DecimalFormat("0.00", new DecimalFormatSymbols(Locale.US));
		
		String str = "";
		str += maxlen + "-" + df.format(lambda) + "-" + df.format(alpha) + "-";
		
		if (beta == Double.MAX_VALUE)
			str +="inf-";
		else
			str += df.format(beta) + "-";
		if (normalize)
			str += "normalized-";
		else
			str += "raw-";
		
		if (nontottering)
			str += "nontottering-";
		else if (paths)
			str += "paths-";
		else
			str += "walks-";
		
		if (reduced)
			str += "reducedpg-";
		else
			str += "plainpg-";
		
		if (op == KernelOperationType.DotProduct)
			str += "dot-";
		else if (op == KernelOperationType.Indicator)
			str += "ind-";
		else if (op == KernelOperationType.Min)
			str += "min-";
		else if (op == KernelOperationType.MinNormalized)
			str += "minnorm-";
		
		if (kw == KernelWeight.Exponential)
			str += "exp";
		else if (kw == KernelWeight.Logistic)
			str += "logistic";
		if (kw == KernelWeight.Diffusion)
			str += "diff";
		
		return str;
	}
	
	public KernelParams clone()
	{
		KernelParams x = new KernelParams();
		x.alpha = alpha;
		x.beta = beta;
		x.lambda = lambda;
		x.epsilon = epsilon;
		x.maxlen = maxlen;
		x.start = start;
		x.end = end;
		x.nontottering = nontottering;
		x.reduced = reduced;
		x.paths = paths;
		x.normalize = normalize;
		x.partialnorm = partialnorm;
		x.op = op;
		x.kw = kw;
		return x;
	}	
}
