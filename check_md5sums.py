#!python

import argparse
import pathlib
import warnings

def read_md5(file):
	d={}
	with open(file,'r') as f:
		for line in f:
			md5, path = line.strip().split()
			base=pathlib.Path(path).name
			d[base]=md5
	return d

def main():
	parser = argparse.ArgumentParser(description="Compare md5sums written to 2 text files.")
	parser.add_argument('file1')
	parser.add_argument('file2')
	args = parser.parse_args()

	d1 = read_md5(args.file1)
	d2 = read_md5(args.file2)
	all = set(d1.keys()) | set(d2.keys())

	for filename in all:
		if filename not in d1:
			warnings.warn(f"{filename} not present in {args.file1}")
			continue
		if filename not in d2:
			warnings.warn(f"{filename} not present in {args.file2}")
			continue
		try:
			assert(d1[filename] == d2[filename])
		except AssertionError:
			warnings.warn("{filename} md5sums do not match.")
		except:
			raise

if __name__ == "__main__":
	main()
