#!/usr/bin/env python
# coding: UTF-8

from __future__ import print_function
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('dirname', help='directory name')
args = parser.parse_args()

files = glob.glob(args.dirname+'/*.jpg')

outf = open('image_list.txt', 'w')

for f in files:
	outf.write(f + '\n')

outf.close()
