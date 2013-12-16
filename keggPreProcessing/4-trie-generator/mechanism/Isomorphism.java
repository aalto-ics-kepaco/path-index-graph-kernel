package mechanism;

import java.util.*;

import mechanism.graphs.*;

public class Isomorphism
{
	/*
	 * Isomorphism algorithm implementation (VF2) for molecular graphs
	 * 
	 */

	class NodePair
	{
		public Node lhs;
		public Node rhs;
		public NodePair(Node lhs, Node rhs)
		{
			this.lhs = lhs;
			this.rhs = rhs;
		}
	}
	
	class StateNode
	{
		public Map<Node,Node> map;
		public BitSet lhsborder;
		public BitSet rhsborder;
		
		StateNode()
		{
			// do nothing
		}
		
		public int size()
		{
			return map.size();
		}
		
		public String toString()
		{
			String s = "";
			for (Node a : map.keySet())
				s += a + " <-> " + map.get(a) + "\n";
			return s;
		}
	}
	
	


	private Graph g1;
	private Graph g2;
	private int itercount = 0;
	private boolean iso = false;
	@SuppressWarnings("unused")
	private Collection<Map<Node,Node>> solutions;
	private int size;
//	private boolean checkBondtype = false;
	
	public Isomorphism(Graph g1, Graph g2)
	{
		this.g1 = g1;
		this.g2 = g2;
		solutions = new ArrayList<Map<Node,Node>>();
		iso = false;
	}
	
	
	public boolean VF2()
	{
		// pre-test for compatibility:
		//  size, formula, atom spectrum and bond spectrum
		//  have to match
		
		// sizes have to match
		if (this.g1.getSize() != this.g2.getSize())
			return false;
		
		if (g1.getEdgeCount() != g2.getEdgeCount())
			return false;
		
		// also spectra's have to match
		if (g1.getFormula() != null && !g1.getFormula().equals(g2.getFormula()))
			return false;
		
		if (g1.getNodeSpectrum() != null && !g1.getNodeSpectrum().equals(g2.getNodeSpectrum()))
			return false;

		if (g1.getEdgeSpectrum() != null && !g1.getEdgeSpectrum().equals(g2.getEdgeSpectrum()))
			return false;
		
		
		size = g1.getSize();
		
		// otherwise continue with algorithm
		StateNode start = new StateNode();
		start.map = new HashMap<Node,Node>();
		start.lhsborder = new BitSet(size);
		start.rhsborder = new BitSet(size);
		match(start);
		
		return iso;
	}
	
	private void match(StateNode s)
	{
		if (iso)
			return;
		
		itercount++;
		
		// full isomorphic mapping found
		if (s.size() == size)
		{
			iso = true;
//			solutions.add( new HashMap<RGAtom,RGAtom>(s.map) ); // shallow copy
			return;
		}
		
		for (NodePair p : candidatepairs(s))
		{
			Node lhs = p.lhs;
			Node rhs = p.rhs;
			
			if (feasible(s, lhs, rhs))
			{
				StateNode ns = new StateNode();
				ns.map = new HashMap<Node,Node>(s.map); // shallow copy, ok because RGAtom's don't ever change
				ns.map.put(lhs, rhs);
				
				ns.lhsborder = (BitSet)s.lhsborder.clone();
				ns.lhsborder.set(lhs.getId(), false);
				ns.rhsborder = (BitSet)s.rhsborder.clone();
				ns.rhsborder.set(rhs.getId(), false);
				
				for (Node ne : lhs.getNodeNeighbors())
					if (ns.lhsborder.get(ne.getId()) == false && !ns.map.containsKey(ne))
						ns.lhsborder.set(ne.getId());
				for (Node ne : rhs.getNodeNeighbors())
					if (ns.rhsborder.get(ne.getId()) == false && !ns.map.containsValue(ne))
						ns.rhsborder.set(ne.getId());

				match(ns);
			}
		}
	}
	
	
	private List<NodePair> candidatepairs(StateNode s)
	{
		List<NodePair> pairs = new ArrayList<NodePair>();
		
		// both sides have border to go through
		if (!s.lhsborder.isEmpty() && !s.rhsborder.isEmpty())
		{
			Node rhs = null;
			for (int i = 0; i < s.rhsborder.size(); i++)
			{
				if (s.rhsborder.get(i) == true)
				{
					rhs = g2.getNodes()[i];
					break;
				}
			}
			
			for (int i = 0; i < s.lhsborder.size(); i++)
			{
				if (s.lhsborder.get(i) == true)
					pairs.add( new NodePair(g1.getNodes()[i],rhs) );
			}
		}
		// either side has emptied its border
		else
		{
			// set of all atoms not mapped on 'rg2'
			Set<Node> rhsleft = new HashSet<Node>( new ArrayList<Node>( Arrays.asList(g2.getNodes())));
			rhsleft.removeAll(s.map.values());
			// pick smallest id atom from them
			int min_id = 100000;
			Node min_atom = null;
			for (Node a : rhsleft)
			{
				if (a.getId() < min_id)
				{
					min_id = a.getId();
					min_atom = a;
				}
			}
		
			for (Node lhs : g1.getNodes())
				if (!s.map.containsKey(lhs))
					pairs.add(new NodePair(lhs,min_atom));
		}
		
		return pairs;
	}
	
	private boolean feasible(StateNode s, Node lhs, Node rhs)
	{
		// the regions of G1 and G2 are divided into three areas:
		// - (MR) mapped region
		// - (BR) adjacent borders to mapped regions
		// - (RR) remote regions
		//
		// the regions don't overlap
		//
		
		// atom symbol has to match
		if (!lhs.getSymbol().equals(rhs.getSymbol()))
			return false;
		
		// check that lhs's mapped neighbors match rhs's mapped neighbors
		Set<Node> lhs_mapped_neighs = new HashSet<Node>();
		for (Node ne : lhs.getNodeNeighbors())
			if (s.map.containsKey(ne))
				lhs_mapped_neighs.add( ne );
		Set<Node> rhs_mapped_neighs = new HashSet<Node>();
		for (Node ne : rhs.getNodeNeighbors())
			if (s.map.containsValue(ne))
				rhs_mapped_neighs.add( ne );
		
		// mapped neighborhoods have to match
		if (lhs_mapped_neighs.size() != rhs_mapped_neighs.size())
			return false;
		
		// check that each atom matches some atom from other side through mapping
		for (Node a : lhs_mapped_neighs)
		{
			if (!rhs_mapped_neighs.contains(s.map.get(a)))
				return false;
			
//			// molecular graph version, check that we have identical bond types
//			if (a.getBond(lhs).getType() != s.map.get(a).getBond(rhs).getType() )
//				return false;
			
			Edge e = a.getEdge(lhs);
			Edge me = s.map.get(a).getEdge(rhs);
			
			// reaction graph version
			if (e.getChangetype() != me.getChangetype() )
				return false;
			
//			if (checkBondtype)
//				if (e.getType() != me.getType() || e.getOldType() != me.getOldType() || e.getNewType() != me.getNewType())
//					return false;
		}
		
		
		// border area sizes have to match
		int lhsneighs = 0;
		for (Node ne : lhs.getNodeNeighbors())
			if (s.lhsborder.get(ne.getId()) == true)
				lhsneighs++;
		int rhsneighs = 0;
		for (Node ne : rhs.getNodeNeighbors())
			if (s.rhsborder.get(ne.getId()) == true)
				rhsneighs++;
		
		if (lhsneighs != rhsneighs)
			return false;
		
		// remote area sizes have to match
		int lhsremotes = 0;
		for (Node ne : lhs.getNodeNeighbors())
			if (s.lhsborder.get(ne.getId()) == false && !s.map.containsKey(ne))
				lhsremotes++;
		int rhsremotes = 0;
		for (Node ne : rhs.getNodeNeighbors())
			if (s.rhsborder.get(ne.getId()) == false && !s.map.containsValue(ne))
				rhsremotes++;
		
		if (lhsremotes != rhsremotes)
			return false;
		
		return true;
	}

}
