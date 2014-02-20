package mechanism;

public enum KernelOperationType
{
	Indicator, // boolean indicator <a,b> == 1 if a > 0 and b > 0
	DotProduct, // basic dot product
	Min, // Histogram intersection = minimum of <a,b>
	MinNormalized
}