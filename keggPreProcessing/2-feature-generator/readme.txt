Generates atom features of mol files. An atom feature is a property of atom and its k-context. Hydrogens are ignored.

Supports generation of

 - atom distribution [e.g. CCO]
 - bond distribution [e.g. 112]
 - Wiener index      [sum of the lengths of the shortest paths between all pairs of atoms in the context]
 - Morgan index      
 - ring memberships  [true/false, whether atom is member of a k-ring]

of context size 'k' around each atom. Uses a message passing algorithm, for citation, use

Markus Heinonen, Sampsa Lappalainen, Taneli Mielikäinen and Juho Rousu
Computing atom mappings for biochemical reactions without subgraphs isomorphism
Journal of Computational Biology 18(1):43-58, 2011
