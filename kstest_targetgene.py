#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-12-29 15:54:59
# @Author  : Abel Yabin Liu
# @Email   : abel.yabin.liu@gmail.com
# @Version : $Id$

import os
import sys
import getopt
import re
import numpy
import math
from scipy.stats import ks_2samp


def cal_mean_var(array):
	array = numpy.array(array)
	sum1 =array.sum()
	array2 = array*array
	sum2 =  array2.sum()
	mean = sum1/len(array)
	var = sum2/len(array)-mean**2
	return(mean, var)


def ks_test(gene, treat, control):
	(t_mean, t_var) = cal_mean_var(treat)
	(c_mean, c_var) = cal_mean_var(control)
	ks_result=str(ks_2samp(treat,control))
#	if(float(re.split(r"=|,|\)", ks_result)[3]) < 0.05):
	if (t_mean!=0):
		fc=math.log(t_mean/c_mean, 2)
		print(gene, t_mean, t_var, c_mean, c_var, sep="\t", end="\t")
		print(fc, re.split(r"=|,|\)", ks_result)[1], re.split(r"=|,|\)", ks_result)[3], treat, control,sep="\t")


def calculate_ks(target_file, fpkm_file):
	sample_index = {}
	sample_fpkm = {}
	for cur_line_number, line in enumerate(open(fpkm_file, 'rU')):
		tmp = line.strip().split("\t")
		if cur_line_number == 0:
			for i in range(2, len(tmp)):
				sample_index[tmp[i]] = i
		else:
			sample_fpkm[tmp[1]] = tmp
#	print(sample_fpkm['WT1'][sample_index['TCGA-AB-3007']])
	for cur_line_number, line in enumerate(open(target_file, 'rU')):
		tmp = line.strip().split("\t")
		if tmp[3] in sample_fpkm.keys():
			variant = []
			control = []
			samp = tmp[-1].split(",")
			for i in range(0,len(samp)):
				if samp[i] in sample_index.keys():
					variant.append(float(sample_fpkm[tmp[3]][sample_index[samp[i]]]))
				else:
					continue
			for s in sample_index.keys():
				if s in samp:
					continue
				else:
					control.append(float(sample_fpkm[tmp[3]][sample_index[s]]))
			if (tmp[3] == "WT1"):
				variant.append(3.1)
				variant.append(2.7)
			if len(variant) >= 2:
				ks_test(tmp[3], variant, control)
		else:
			continue




def usage():
	print ("Usage: python ", sys.argv[0], " --target(-t)=<recurrent_mutation.targetgene> --fpkm(-f)=<fpkm.matrix>")

def main(argv):

	try:
		opts, args = getopt.getopt(argv, "ht:f:", ["help","target","fpkm"])
		if (len(opts) == 0):
			usage()
			sys.exit()
	except getopt.GetoptError:
		usage()
		sys.exit()
	for name, value in opts:
		if name in ('-h', '--help'):
			usage()
			sys.exit()
		if name in ('-t', '--target'):
			target = value
		if name in ('-f', '--fpkm'):
			fpkm = value
	calculate_ks(target, fpkm)


if __name__ == '__main__':
	main(sys.argv[1:])