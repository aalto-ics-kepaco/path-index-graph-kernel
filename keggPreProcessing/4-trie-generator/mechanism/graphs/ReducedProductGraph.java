package mechanism.graphs;

import mechanism.KernelParams;

public class ReducedProductGraph extends ProductGraph
{
	public ReducedProductGraph(Graph g1, Graph g2, KernelParams params)
	{
		this(g1,g2,params,true);
	}
	
	public ReducedProductGraph(Graph g1, Graph g2, KernelParams params, boolean nodematch)
	{
		this.g1 = g1;
		this.g2 = g2;
		this.params = params;
		
		createNodes(nodematch);
		removeNC();
		createEdges();
	}
	
	private void removeNC()
	{
		// remove nodes with NC pairings (non-core/core)
		int toremove = 0;
		for (PGNode v : nodes)
		{
			if ((v.a1.getCoreDist() == 0 || v.a2.getCoreDist() == 0) && (v.a1.getCoreDist() != v.a2.getCoreDist()))
			{
				nodes[v.id] = null;
				toremove++;
			}
		}
		
		// create a new nodearray
		PGNode[] nodes2 = new PGNode[nodes.length - toremove];
		// copy stuff to new array and compress id's
		int k = 0;
		for (PGNode v : nodes)
		{
			if (v != null)
			{
				v.id = k;
				nodes2[k++] = v;
			}
		}
		
		// change pointer
		nodes = nodes2;
		nc = nodes.length;		
	}
}
