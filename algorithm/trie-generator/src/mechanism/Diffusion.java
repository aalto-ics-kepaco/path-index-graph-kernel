package mechanism;

import mechanism.graphs.*;
import Jama.*;

public class Diffusion
{
	private Graph g;
	private int n;
	private double[][] H;
	private double[][] E;
	private double beta;

	public Diffusion(Graph g, double beta)
	{
		this.g = g;
		this.beta = beta;

		constructH();
		expm(); // compute expm
	}

	// computes the matrix: beta*H = beta*(A - D)
	// uses this.g to compute this
	private void constructH()
	{
		n = g.getSize();
		H = new double[n][n];

		for (Node n : g.getNodes())
		{
			// fill in the diagonal value (-degree)
			H[n.getId()][n.getId()] = - beta * n.getNodeNeighbors().size();

			// fill in the adjacency matrix
			for (Node n2 : n.getNodeNeighbors())
			{
				H[n.getId()][n2.getId()] = beta*1.0;
			}
		}
	}

	// compute e^H
	public void expm()
	{
		E = expm_matlab();
	}

	public double get(int i, int j)
	{
		return E[i][j];
	}

	private double[][] matrixmult(double[][] A, double[][] B)
	{
		double[][] X = new double[n][n];

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				for (int k = 0; k < n; k++)
					X[i][j] += A[i][k] * B[k][j];
		return X;
	}

    /*
     * code translated to java from
     *  '/matlab/toolbox/matlab/matfun/expm.m'
     *
     * some optimization and cleanup done
     *
     */
    private double[][] expm_matlab()
    {
    	// mvals
    	int[] mvals = {3, 5, 7, 9, 13};
    	int mvals_last = mvals[mvals.length-1];

    	// theta
    	double[] theta = {//3.650024139523051e-008,
                //5.317232856892575e-004,
                1.495585217958292e-002,  // m_vals = 3
                //8.536352760102745e-002,
                2.539398330063230e-001,  // m_vals = 5
                //5.414660951208968e-001,
                9.504178996162932e-001,  // m_vals = 7
                //1.473163964234804e+000,
                2.097847961257068e+000,  // m_vals = 9
                //2.811644121620263e+000,
                //3.602330066265032e+000,
                //4.458935413036850e+000,
                5.371920351148152e+000}; // m_vals = 13
    	double theta_last = theta[theta.length-1];


    	// transform H into object of type 'Matrix A'
    	Matrix A = new Matrix(n,n);
    	for (int i = 0; i < n; i++)
    		for (int j = 0; j < n; j++)
    			A.set(i,j,H[i][j]);

    	// 1-norm
    	double normA = A.norm1();

    	if (normA <= theta_last)
    	{
    		// no scaling and squaring is required
    	    for (int i = 0; i < mvals.length; i++)
    	    {
    	    	// use first large enough theta
    	    	if (normA <= theta[i])
    	    	{
    	    		return PadeApproximantOfDegree(A, mvals[i]);
    	    	}
    	    }
    	}

		// scaling required
	    MantissaExp ts = frexp(normA/theta_last);
	    double t = ts.mantissa;
	    int s = ts.exp;

   	    // adjust s if normA/theta(end) is a power of 2
	    if (t == 0.5) // accuracy problems?
	    	s -= 1;
	    A = A.times(1.0 / Math.pow(2,s)); // Scaling
	    double[][] F = PadeApproximantOfDegree(A, mvals_last);
	    for (int i = 0; i < s; i++)
	    {
	        F = matrixmult(F, F); // Squaring
	    }

    	return F;
    }

    private double[][] PadeApproximantOfDegree(Matrix A, int m)
    {
    	// coefficients
    	long[][] C = {{},
			         {},
			         {},
			         {120,60,12,1}, // m == 3
			         {},
			         {30240,15120,3360,420,30,1}, // m == 5
			         {},
			         {17297280,8648640,1995840,277200,25200,1512,56,1}, // m == 7
			         {},
			         {17643225600L,8821612800L,2075673600,302702400,30270240,2162160,110880,3960,90,1}, // m == 9
			         {},
			         {},
			         {},
			         {64764752532480000L,32382376266240000L,7771770303897600L,1187353796428800L,129060195264000L,10559470521600L,670442572800L,33522128640L,1323241920L,40840800L,960960,16380,182,1}}; // m == 13

    	// correct coefficient by degree 'm'
    	long[] c = C[m];

    	Matrix F = new Matrix(n,n);
    	Matrix I = Matrix.identity(n,n);

    	if (m <= 9)
    	{
    		// create a set of matrices
    		int ceil = (int)Math.ceil((m+1)/2);
    		Matrix[] Ap = new Matrix[ceil];

    		Ap[0] = I;
    		Ap[1] = A.times(A);

    		for (int j = 2; j < ceil; j++)
    			Ap[j] = Ap[j-1].times(Ap[1]);

    		Matrix U = new Matrix(n,n); // zeros
    		Matrix V = new Matrix(n,n); // zeros

    		for (int j = m+1; j >= 2; j -= 2) // even: [m+1, m-1, m-3, ... , 2]
    			U = U.plus( Ap[j/2-1].times(c[j-1]) );

    		U = A.times(U);

    		for (int j = m; j >= 1; j -= 2) // odd: [m,m-2,m-4, ... , 1]
    			V = V.plus( Ap[(j+1)/2-1].times(c[j-1]) );

    		F = V.minus(U).solve(V.plus(U));
    	}
        // For optimal evaluation need different formula for m >= 12.
    	else if (m >= 12)
    	{
	        Matrix A2 = A.times(A);
	        Matrix A4 = A2.times(A2);
	        Matrix A6 = A2.times(A4);

	        Matrix U = A.times( A6.times( A6.times(c[13]).plus(A4.times(c[11])).plus(A2.times(c[9])) ).plus(A6.times(c[7])).plus(A4.times(c[5])).plus(A2.times(c[3])).plus(I.times(c[1])) );
	        Matrix V = A6.times( A6.times(c[12]).plus(A4.times(c[10])).plus(A2.times(c[8])) ).plus(A6.times(c[6])).plus(A4.times(c[4])).plus(A2.times(c[2])).plus(I.times(c[0]));

	        F = V.minus(U).solve(V.plus(U));
    	}

    	// convert back to double[][] array
        double[][] res = new double[n][n];
        for (int i = 0; i < n; i++)
        	for (int j = 0; j < n; j++)
        		res[i][j] = F.get(i,j);

        return res;
    }

    private MantissaExp frexp(double x)
    {
    	MantissaExp me = new MantissaExp();
    	if (x == 0)
    	{
    	    me.mantissa = 0;
    	    me.exp = 0;
    	    return me;
    	}

    	long bits = Double.doubleToLongBits(x);
    	me.mantissa = Double.longBitsToDouble((0x800fffffffffffffL & bits) | 0x3fe0000000000000L);
    	me.exp = (int)((0x7ff0000000000000L & bits) >> 52) - 1022;
    	return me;
    }

    class MantissaExp
    {
    	public double mantissa;
    	public int exp;
    }
}
