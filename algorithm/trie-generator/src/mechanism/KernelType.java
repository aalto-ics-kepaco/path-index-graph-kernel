package mechanism;

public enum KernelType
{
	WK, // standard gartner's exponential walk kernel (equals unweighted enumerative mechanism)
	RWK, // standard kashima's marginal walk kernel (equals unweighted marginal mechanism)
	RGK, // Tsuda's reaction graph kernel (unweighted RWK inside)
	SG, // standard subgraph kernel
	RM, // three rousu's subgraph kernels
	SOR,
	DOR,
	MMECH, // marginal mechanism
	EMECH, // enumerative mechanism
	MG, // molecule kernel
	SP
}
