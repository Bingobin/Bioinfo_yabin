#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-11-02 20:05:59
# @Author  : Abel Yabin Liu
# @Email   : abel.yabin.liu@gmail.com
# @Version : $Id$

import os
import sys
import getopt
import math


def factorial(n):
	if n == 0 or n == 1:
		return 1
	else:
		return (n*factorial(n-1))

def sumpvalue(number, lamdar):
	p = float(lamdar) ** int(number) / factorial(int(number)) * math.e ** -float(lamdar)
	if int(number) == 0:
		return 0
	else:
		return (p +  sumpvalue(int(number)-1, lamdar))


def usage():
	print ("Usage: python ", sys.argv[0], " --lamdar(-l)=<float> --number(-n)=<integer>")

def main(argv):
	lamdar = 0
	number = 0
	try:
		opts, args = getopt.getopt(argv, "hl:n:", ["help","lamdar=","number="])
		if (len(opts) == 0):
			usage()
			sys.eixt()
	except getopt.GetoptError:
		usage()
		sys.exit()
	for name, value in opts:
		if name in ('-h', '--help'):
			usage()
			sys.exit()
		if name in ('-l', '--lamdar'):
			lamdar = value
		if name in ('-n', '--number'):
			number = value
	print(sumpvalue(number, lamdar))


if __name__ == '__main__':
	main(sys.argv[1:])