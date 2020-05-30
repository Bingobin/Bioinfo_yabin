#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-09-07 10:13:09
# @Author  : Abel Yabin Liu (abel.yabin.liu@gmail.com)
# @Link    : ${link}
# @Version : $Id$

import os
import sys
import getopt
from collections import defaultdict


def merge_by_fist_line(f1, f2):
	file_merge = defaultdict(dict)
	title = ''
	f1_col_num = 1
	f2_col_num = 1
	for cur_line_number, line in enumerate(open(f1, 'rU')):
		tmp = line.strip().split("\t")
		f1_col_num = len(tmp) - 1
		if cur_line_number == 0:
			title = line.strip()
		else:
			file_merge[tmp[0]]['file1'] = "\t".join(tmp[1:])
	for cur_line_number, line in enumerate(open(f2, 'rU')):
		tmp = line.strip().split("\t")
		f2_col_num = len(tmp) - 1
		if cur_line_number == 0:
			title += "\t" + "\t".join(tmp[1:])
		else:
			file_merge[tmp[0]]['file2'] = "\t".join(tmp[1:])
	print (title)
	for i in sorted(file_merge.keys()):
		content = i
		if 'file1' in file_merge[i].keys():
			content += "\t" + file_merge[i]['file1']
		else:
			content += "\tNull" *  f1_col_num
		if  'file2' in file_merge[i].keys():
			content += "\t" + file_merge[i]['file2']
		else:
			content += "\tNull" * f2_col_num
		print(content)

def usage():
	print ("Usage: python ", sys.argv[0], " --file1(-a)=<file1> --file2(-b)=<file2>")

def main(argv):
	file1 = ''
	file2 = ''
	try:
		opts, args = getopt.getopt(argv, "ha:b:", ["help","file1=","file2="])
		if (len(opts) == 0):
			usage(); sys.exit()
	except getopt.GetoptError:
		usage(); sys.exit()
	for name, value in opts:
		if name in ('-h', '--help'):
			usage(); sys.exit()
		if name in ('-a', '--file1'):
			file1 = value
		if name in ('-b', '--file2'):
			file2 = value
	merge_by_fist_line(file1, file2)

if __name__ == '__main__':
	main(sys.argv[1:])