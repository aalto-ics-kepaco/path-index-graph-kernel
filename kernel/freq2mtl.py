#!/usr/bin/python
#
#
#
#

import sys, os, glob, string

print " # Converts suffix tree generated .freqs files into matlab sparse matrix files"
print


def compare(a,b):
	alen = len(filter(str.isupper, a.split(" ", 1)[0]))
	blen = len(filter(str.isupper, b.split(" ", 1)[0]))
	
	if alen < blen:
		return -1
	elif alen > blen:
		return 1
	return 0

SUFFIXDIR = "../../results/computed/5-tbwt/"


# open file
filename = SUFFIXDIR + "result-kegg.freqs"

if len(sys.argv) > 1:
	filename = sys.argv[1]

prefix = filename[0: filename.find(".")]

print "reading freqs file", filename

# read file
f = open(filename)
lines = map(str.strip, f.readlines())
f.close()

# sort the file
print "sorting..."
lines.sort(compare)

# read reactions into a dictionary
filelist = SUFFIXDIR + "rgraphlist.txt"

if len(sys.argv) > 2:
	filelist = sys.argv[2]


print "reading rgraphlist file", filelist

# read file
f = open(filelist)
reacs = map(str.strip, f.readlines())
f.close()


for i in range(len(reacs)):
	r = reacs[i][:-4]
	words = r.split("_")
	reacs[i] = words[0] + "_" + words[2] + "_" + words[1] # changes e.g. R00004_0_b to R00004_b_0

# reaction dictionary
reacdict = dict((r, i+1) for i,r in enumerate(reacs))

print "#-debug-clemens-# reacdict: " + str(reacdict)

# start writing sparse matlab matrices
g = open(prefix + ".mtl", "w")
gc = open(prefix + "-core.mtl", "w")
gi = None
gci = None
	
oldlen = -1
newlen = -1

# read line by line
for i,line in enumerate(lines):
	words = line.strip().split()

	newlen = len(filter(str.isupper, words[0]))
	words = words[1:] # clemens: words are graph:frequency items for the i'th reaction

	# changed length
	if oldlen != newlen:
		try:
			gi.close()
		except:
			pass
		try:
			gci.close()
		except:
			pass
		
		print "opening new files", newlen
		gi = open(prefix + "_%d.mtl" % (newlen), "w")
		gi.write("%d %d 0\n" % (20665838,17430)) # header line to indicate size of matrix
		
		gci = open(prefix + "-core_%d.mtl" % (newlen), "w")
		gci.write("%d %d 0\n" % (20665838,17430)) # header line to indicate size of matrix
	
	for w in words:
		r,c = w.split(":")


		print "#-debug-clemens-# reacdict: " + str(reacdict)

		print "#-debug-clemens-# graph r: " + str(r)
		print "#-debug-clemens-# freq c: " + str(c)

		print "#-debug-clemens-# reacdict[r]: " + str(reacdict[r])

		g.write("%d %d %s\n" % (i+1,reacdict[r],c))
		gi.write("%d %d %s\n" % (i+1,reacdict[r],c))

		if "{" in line or "[" in line:
			gc.write("%d %d %s\n" % (i+1,reacdict[r],c))
			gci.write("%d %d %s\n" % (i+1,reacdict[r],c))
	
	oldlen = newlen


gi.close()
g.close()
gci.close()
gc.close()

