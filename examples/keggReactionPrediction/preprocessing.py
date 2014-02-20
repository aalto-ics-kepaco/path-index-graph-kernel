# Runs all preprocessing steps for preprocessing mol files from KEGG for 
# reaction prediction task
#
# @author Clemens Westrup

import argparse
import os, errno, shutil, sys
import subprocess

# constants
outputdir = "results"
molfeaturesdir = outputdir + "/mol-features"
reactiongraphsdir = outputdir + "/reaction-graphs"

# parse command line arguments
def parseargs():
	parser = argparse.ArgumentParser(description='Preprocess kegg files.')
	parser.add_argument('-k', '--keggpath', type=str, dest='keggpath', 
		required=True, help='path to kegg ligand database')
	parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
		help='show verbose information')
	args = parser.parse_args()
	return [args.keggpath, args.verbose]

# mkdir -p
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

# main function
if __name__ == '__main__':
	[keggpath, verbose] = parseargs()

	keggpath = os.path.abspath(keggpath)
	print keggpath

	print ""
	print "--- preprocessing kegg data for path index graph kernel ---"
	print ""

	# remove ending slash on path if existing
	if keggpath.endswith('/'):
		keggpath = keggpath[:-1]

	# show or hide subprocess output
	if (verbose):
		OUTPUT = None
	else:
		OUTPUT = open(os.devnull, 'w')

	# check if output dir exists
	if (os.path.isdir(outputdir)):
		print "Directory \"" + outputdir + "\" for storing results exists already."
		overwrite = raw_input("Press enter to overwrite or crtl+c to abort.")
		shutil.rmtree(outputdir)

	# create temp dir
	print "Writing results in dir: " + outputdir
	mkdir_p(outputdir)

	### reaction list

	# Extract reaction list
	print "Extracting reaction list from kegg, this might take a while..."
	p = subprocess.Popen(
		'python reaction-list/extract-reactions.py -k ' + keggpath 
		+ ' > ' + outputdir + '/kegg-reactions.txt',
		shell=True, stdout=OUTPUT)
	p.wait()
	if p.returncode != 0:
		sys.exit(1)

	### mol-features

	# create output dir for mol-features
	print "Creating dir " + molfeaturesdir
	mkdir_p(molfeaturesdir)
	# create atom features from KEGG mol files
	print "Create atom features from KEGG mol files ..."
	p = subprocess.Popen(
		'python feature-generator/generator2011.py ' + keggpath + '/mol/*'
		+ ' -k all --output-dir ' + molfeaturesdir + '/',
		shell=True, stdout=OUTPUT)
	p.wait()
	if p.returncode != 0:
		sys.exit(1)

	### reaction graphs

	# create output dir reaction graphs
	print "Creating temp dir " + reactiongraphsdir
	mkdir_p(reactiongraphsdir)
	# Create reactiongraph mol files with atommapper
	print "Create reactiongraph mol files with atommapper ..."
	os.chdir("./atommapper/bin");
	p = subprocess.Popen(
		'java Mapper2000 -rgraphs -moldir ' + keggpath + '/mol/ '
		+ '-featdir ../../' + molfeaturesdir + '/ '
		+ '-reacfile ../../' + outputdir + '/kegg-reactions.txt ' 
		+ '-output ../../' + reactiongraphsdir + '/',
		shell=True, stdout=OUTPUT)
	p.wait()
	if p.returncode != 0:
		sys.exit(1)

	print("All done. Reaction graph .mol files can be "
		+ "found in " + reactiongraphsdir)
