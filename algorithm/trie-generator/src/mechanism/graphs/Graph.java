package mechanism.graphs;

import java.util.*;

import mechanism.Isomorphism;

/*
 * General abstract graph class which provides roughly the stuff needed for
 * molecular graphs (ligand,direction,formula,..)
 */
public abstract class Graph
{
	protected int id = -1;
	protected String ligand = null;
	protected Map<String,Integer> formula = null;
	protected Map<String,Integer> nodespectrum = null;
	protected Map<String,Integer> edgespectrum = null;
	protected int direction = 0;
	protected int hash = 0;

	protected Node[] nodes = null;
	protected Edge[] edges = null;
	
	// empty constructor
	public Graph()
	{
		
	}
	
	// abstract constructor
	// create a copy of 'parent' using only 'nodebits' nodes
	protected Graph(Graph parent, BitSet nodebits)
	{
		id = -1;
		ligand = parent.ligand;
		direction = parent.direction;
		
		formula = new HashMap<String,Integer>();
		
		// if parent has spectra, also put here
		if (parent.getNodeSpectrum() != null)
		{
			nodespectrum = new HashMap<String,Integer>();
			edgespectrum = new HashMap<String,Integer>();
		}
		
		nodes = new Node[nodebits.cardinality()];
		
		int[] nodemap = new int[parent.getNodeCount()];
		
		// place atoms into subgraph
		int aid = 0;
		for (int i = nodebits.nextSetBit(0); i >= 0; i = nodebits.nextSetBit(i+1))
		{
			Node a = ((Atom)parent.nodes[i]).clone(); //new Atom(aid, parent.nodes[i].getSymbol());
//			a.setParent(parent);
			
			nodemap[i] = aid;
			
			nodes[aid++] = a;
			
			if (formula.containsKey(a.getSymbol()))
				formula.put(a.getSymbol(), formula.get(a.getSymbol())+1);
			else
				formula.put(a.getSymbol(), 1);
		}
		
		// count number of bonds
		int edgecount = 0;
		for (Edge b : parent.edges)
			if (nodebits.get(b.getSource().getId()) && nodebits.get(b.getTarget().getId()))
				edgecount++;
		
		edges = new Edge[edgecount];
		
		// place bonds into subgraph
		int bid = 0;
		for (Edge b : parent.edges)
		{
			if (nodebits.get(b.getSource().getId()) && nodebits.get(b.getTarget().getId()))
			{
//				Bond newb = new Bond(bid, nodes[ nodemap[b.getSource().getId()] ], nodes[nodemap[b.getTarget().getId()]], b.getType(), b.getChangetype(), b.getOldtype(), b.getNewtype(), this);
				Edge newe = ((Bond)b).clone();
				newe.getSource().addNeighbor(newe, newe.getTarget());
				newe.getTarget().addNeighbor(newe, newe.getSource());
				edges[bid++] = newe;
				
				if (edgespectrum != null)
				{
					String bstr;
					if (newe.getSource().getSymbol().compareTo(newe.getTarget().getSymbol()) < 0)
						bstr = newe.getSource().getSymbol() + b.getChangetype() + newe.getTarget().getSymbol();
					else
						bstr = newe.getTarget().getSymbol() + b.getChangetype() + newe.getSource().getSymbol();
					
					if (edgespectrum.containsKey(bstr))
						edgespectrum.put(bstr, edgespectrum.get(bstr)+1);
					else
						edgespectrum.put(bstr, 1);
				}
			}
		}
		
		if (nodespectrum != null)
		{
			// collect atom spectrum 
			for (Node a : nodes)
			{
				String astr = "";
				for (Node ne_a : a.getNodeNeighbors())
				{
					astr += ne_a.getSymbol();
				}
				
				// sort internally
				char[] chars = astr.toCharArray();
				Arrays.sort(chars);
				astr = new String(chars);
				astr = a.getSymbol() + astr;
				
				
				if (nodespectrum.containsKey(astr))
					nodespectrum.put(astr, nodespectrum.get(astr)+1);
				else
					nodespectrum.put(astr, 1);
			}
		}
	}
	
	
	public int getIndex()
	{
		return id;
	}
	
	public void setIndex(int i)
	{
		id = i;
	}
	
	public String getLigand()
	{
		return ligand;
	}
		
	public Map<String,Integer> getFormula()
	{
		return formula;
	}
	
	public int getDirection()
	{
		return direction;
	}
	public Node[] getNodes()
	{
		return nodes;
	}
	public Edge[] getEdges()
	{
		return edges;
	}
	
	public int getSize()
	{
		return nodes.length;
	}
	public int getNodeCount()
	{
		return nodes.length;
	}
	public int getEdgeCount()
	{
		return edges.length;
	}

	public void computeDistances()
	{
		List<Node> corenodes = new ArrayList<Node>();
		
		// find core atoms
		for (Node a : nodes)
		{
			for (Edge b : a.getEdgeNeighbors())
			{
				if (b.getChangetype() != 0 && !corenodes.contains(a))
				{
					a.setCoreDist(0);
					corenodes.add(a);
				}
			}
		}
		
		// compute min distances to core atoms using multiple dijkstra
		// start from core atoms
		for (Node src : corenodes)
		{
			int[] dist = new int[nodes.length];
			Node[] previous = new Node[nodes.length];
			List<Node> Q = new ArrayList<Node>();
			
			for (Node a : nodes)
			{
				dist[a.getId()] = Integer.MAX_VALUE;
				previous[a.getId()] = null;
			}
			dist[src.getId()] = 0;
			
			for (Node a : nodes)
				Q.add(a);
			while (!Q.isEmpty())
			{
				// find one with lowest dist in Q
				Node u = Q.get(0);
				int mindist = dist[u.getId()];
				for (Node cand : Q)
				{
					if (dist[cand.getId()] < mindist)
					{
						mindist = dist[cand.getId()];
						u = cand;
					}
				}
				
				if (dist[u.getId()] == Integer.MAX_VALUE)
					break;
				
				Q.remove(u);
				for (Node ne : u.getNodeNeighbors())
				{
					int alt = dist[u.getId()] + 1;
					if (alt < dist[ne.getId()])
					{
						dist[ne.getId()] = alt;
						previous[ne.getId()] = u;
					}
				}
				
				// node 'u' processed
				if (dist[u.getId()] < u.getCoreDist())
					u.setCoreDist(dist[u.getId()]);
			}
		}
		
		// no core -> set all nodes at core
		// this way all weightings become uniform
		if (corenodes.size() == 0)
		{
			for (Node a : nodes)
				a.setCoreDist(0);
		}

		// remove unconnected areas, i.e. if coredist == MAX_INT
		Node[] nodes2 = new Node[nodes.length];
		Edge[] edges2 = new Edge[edges.length];
		formula.clear();
		
		int i = 0;
		int j = 0;
		
		for (Node n : nodes)
		{
			if (n.getCoreDist() != Integer.MAX_VALUE) // reachable
			{
				nodes2[i++] = n;
				
				if (formula.containsKey(n.getSymbol()))
					formula.put(n.getSymbol(), formula.get(n.getSymbol()) + 1);
				else
					formula.put(n.getSymbol(), 1);
				
				for (Edge e : n.getEdgeNeighbors())
					if (e.source == n)
						edges2[j++] = e;
			}
		}
		
		// compress
		nodes = Arrays.copyOf(nodes2, i);
		edges = Arrays.copyOf(edges2, j);
		
		// compress id's
		int id = 0;
		for (Node n : nodes)
			n.id = id++;
		
		id = 0;
		for (Edge e : edges)
			e.id = id++;
		
		// rehash .nodeneighs and .edgeneighs
		for (Node n : nodes)
		{
			// make a copy of neighbor-edges
			List<Edge> nedges = new ArrayList<Edge>(n.edgeneighs.values());
			
			n.nodeneighs.clear(); // clear previous ones
			n.edgeneighs.clear();
			
			for (Edge e : nedges) // fill with new hash-values
			{
				n.edgeneighs.put(e.getOther(n), e);
				n.nodeneighs.add(e.getOther(n));
			}
		}
	}
	
	public abstract Graph createSubgraph(BitSet nodebits);
	
	public Map<String,Integer> getEdgeSpectrum()
	{
		return edgespectrum;
	}
	
	public Map<String,Integer> getNodeSpectrum()
	{
		return nodespectrum;
	}

	public int hashCode()
	{
		if (hash != 0)
			return hash;
		
		// size
		hash = (nodes.length+1) * (edges.length+1);
		
		int tmp;
		
		// formula
		if (formula != null)
		{		
			tmp = 1;
			for (String s : formula.keySet())
			{
				int primevalue = 1;
				for (char c : s.toCharArray())
					primevalue *= PRIMES[c];
				
				tmp *= formula.get(s) * primevalue;
			}
			hash += tmp;
		}
		
		
		// atomspectrum
		if (nodespectrum != null)
		{
			tmp = 1;
			for (String s : nodespectrum.keySet())
			{
				int primevalue = 1;
				for (char c : s.toCharArray())
					primevalue *= PRIMES[c];
				
				tmp *= nodespectrum.get(s) * primevalue;
			}
			hash += tmp;
		}
		
		// bondspectrum
		if (edgespectrum != null)
		{
			tmp = 1;
			for (String s : edgespectrum.keySet())
			{
				int primevalue = 1;
				for (char c : s.toCharArray())
					primevalue *= PRIMES[c];
				
				tmp *= edgespectrum.get(s) * primevalue;
			}
			hash += tmp;
		}
		
		return this.hash;
	}
	
	public boolean equals(Object other)
	{
		Isomorphism iso = new Isomorphism(this, (Graph)other);
		return iso.VF2();
	}
	
	public static int[] PRIMES = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657};

	public int getMapNum()
	{
		return 0;
	}

}





