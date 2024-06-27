'''
Python wrapper over Jens' AmpliconClassifier
'''

import argparse
import glob
import os
import sys
import subprocess

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Classify AA amplicon type")
	parser.add_argument("-d", "--dir", help="directory containing cycles, graph files", required=True)
	parser.add_argument("--ref", help="Reference genome name used for alignment, one of hg19, GRCh37, or GRCh38", choices=["hg19", "GRCh37", "hg38", "GRCh38"], required=True)
	parser.add_argument("-n", "--name", help="Name of sample", required=True)
	args = parser.parse_args()

	# Identify cycles, graph files
	cycles_files = glob.glob(os.path.join(args.dir,"*_cycles.txt"))
	cycles_files.sort()
	graph_files = glob.glob(os.path.join(args.dir,"*_graph.txt"))
	graph_files.sort()
	n_amplicons = len(cycles_files)
	print (n_amplicons)

	# Run AmpliconClassifier
	ac_path = "/expanse/lustre/projects/csd677/collab/bin/AmpliconClassifier/amplicon_classifier.py"
	for i in range(n_amplicons):
		cycles_path = cycles_files[i]
		graph_path = graph_files[i]
		subprocess.run(["python",ac_path, "--cycles", cycles_path, "--graph", graph_path, "--ref", args.ref])
		print (cycles_path)

