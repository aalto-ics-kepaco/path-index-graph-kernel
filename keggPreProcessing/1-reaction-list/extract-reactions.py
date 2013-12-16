#!/usr/bin/python
#
# extracts from KEGG all reactions which are valid
# 
# valid reactions contain exact same atoms on both sides of the reaction
# there can be missing mol-files, empty mol-files, etc. as long as atom spectra match
#

import os, optparse, re, sys, time, kegg


REACTIONFILE = "/home/group/icomic/data/kegg/ligand/LATEST/reaction"
MOLFOLDER = "/home/group/icomic/data/kegg/ligand/LATEST/mol/"
KEGGFOLDER = "/home/group/icomic/data/kegg/ligand/LATEST/"
#REACTIONFILE = "/home/mqheinon/temp/reaction"
#MOLFOLDER = "/home/mqheinon/mqheinon/duuni/ec-clustering/tsuda-data/loo-dataset/KEGGCOMPOUND-July2007/"




parser = optparse.OptionParser(usage="%prog [-opts] KEGGFOLDER")
parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Verbose mode")
parser.add_option("-q", "--quiet", dest="quiet", action="store_true", default=False, help="Quiet mode")

(options, args) = parser.parse_args()

if args:
	KEGGFOLDER = args[0]
	REACTIONFILE = KEGGFOLDER + "reaction"
	MOLFOLDER = KEGGFOLDER + "mol/"
	
	kegg.Kegg.MOL_FOLDER = MOLFOLDER
	kegg.Kegg.REACTION_FILE = REACTIONFILE
	



R = kegg.Parsers.parse_reaction(REACTIONFILE)
reaclist = []

# expand equations, remove "+", and replace "<=>" with "p"
# remove (n)'s from equations
for c,r in R.items():
	if "EQUATION" not in r:
		continue
	
	e = r["EQUATION"]
	if "G" in e:
		continue
	
	words = e.split()
	words = [w for w in words if w != "+"]
	cane = []  # canonical e
	
#	if "n" in e:
#		print words
	
	i = 0
	while i < len(words):
		w = words[i]
		
		if w == "2n":
			w = "2*n"
		if w == "2m":
			w = "2*m"
		
		if w.isdigit() or (len(w) < 6 and w != "<=>"):
			cane.append(words[i+1] + "(" + w + ")" )
			i += 2
		else:
			cane.append(w)
			i += 1
	
#	if "n" in e:
#		print cane
	
#	for i in range(len(cane)):
#		w = cane[i]
#		try:
#			x = int(w)
#			cane[i+1] = (cane[i+1] + " ") * x
#			cane[i+1] = cane[i+1].strip()
#		except ValueError:
#			pass
	
#	e = [w for w in cane if not w.isdigit() and w != "+"]
	cat = " ".join(cane)
	
	# formula is recurring
	if "(" in cat:
		n = 2
		m = 2
		
		if "(n-1)" in cat:
			n = 3
		if "(m-1)" in cat:
			m = 3
		
		newe = []
		for w in cane:
			if "(" in w:
#				print c,e,cane,w
				
				multi = eval(w[6:])
				for i in range(multi):
					newe.append( w[0:6] )
			else:
				newe.append(w)
		
		cane = newe
	
	e = " ".join(cane)
	e = e.replace("<=>", "p")
	
	R[c]["EQUATION"] = e
	reaclist.append(c)


finalreacs = []

for c in reaclist:
	if "EQUATION" in R[c]:
		reac = kegg.Reaction(c)
#		print reac, reac.has_hydrogens(), reac.equation, len(reac.subs), len(reac.prods), len(reac.atoms()), "".join([a.symbol for a in reac.atoms()])
		if reac.is_balanced():
			finalreacs.append(c)
			R[c]["REACTION"] = reac
#			print R[c]["EQUATION"]
#			print reac.equation
#			print reac.trim_formula.replace("<=>", "p").replace("+ ","")

for c in finalreacs:
#	print c, len(R[c]["REACTION"].atoms()), R[c]["REACTION"].trim_formula.replace("<=>", "p").replace("+ ","")
	print c, len(R[c]["REACTION"].atoms()), R[c]["EQUATION"]










