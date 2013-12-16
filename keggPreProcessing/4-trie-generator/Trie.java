import java.util.*;

import mechanism.graphs.*;

public class Trie
{
	private TrieNode root;
	private List<TrieNode> nodes;
	
	public Trie()
	{
		nodes = new ArrayList<TrieNode>();
		
		// create root automatically
		root = new TrieNode(null, root); // points to itself
	}
	
	public void addNode(TrieNode v)
	{
		nodes.add(v);
	}
	
	public TrieNode getRoot()
	{
		return root;
	}
	
	public int getSize()
	{
		return nodes.size();
	}
	
	public List<TrieNode> getNodes()
	{
		return nodes;
	}
}

class TrieNode
{
	public String symbol;
	public TrieNode parent;
	public List<TrieNode> children;
	
	public TrieNode(String s, TrieNode p)
	{
		symbol = s;
		parent = p;
		children = new LinkedList<TrieNode>();
	}
	
	public TrieNode AddChild(Node a)
	{
		TrieNode temp = new TrieNode(a.getSymbol(), this);
		children.add(temp);
		return temp;
	}
	public void AddChild(String s)
	{
		children.add(new TrieNode(s, this));
	}
	public void AddChild(TrieNode child)
	{
		children.add(child);
	}
}