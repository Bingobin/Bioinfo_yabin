#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import getopt
import re
import numpy
import math
from scipy.stats import ks_2samp
from scipy import stats


def cal_mean_var(array):
	array = numpy.array(array)
	sum1 =array.sum()
	array2 = array*array
	sum2 =  array2.sum()
	mean = sum1/len(array)
	var = sum2/len(array)-mean**2
	return(mean, var)


def ks_test(treat, control):
	(t_mean, t_var) = cal_mean_var(treat)
	(c_mean, c_var) = cal_mean_var(control)
	ks_result=str(ks_2samp(treat,control))
	t_test = str(stats.ttest_ind(treat,control))
	print(t_mean, t_var, c_mean, c_var, sep="\t", end="\t")
	print(ks_result,t_test)


def calculate_ks(file1, file2):
	treatment = []
	control = []
	for cur_line_number, line in enumerate(open(file1, 'rU')):
		treatment.append(float(line.strip()))
	for cur_line_number, line in enumerate(open(file2, 'rU')):
		control.append(float(line.strip()))
	ks_test(treatment, control)



def usage():
	print ("Usage: python ", sys.argv[0], " --treatment(-t)=<data1> --control(-c)=<data2>")

def main(argv):

	try:
		opts, args = getopt.getopt(argv, "ht:c:", ["help","treatment","control"])
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
		if name in ('-t', '--treatment'):
			treatment = value
		if name in ('-c', '--control'):
			control = value
	calculate_ks(treatment, control)


if __name__ == '__main__':
	main(sys.argv[1:])
