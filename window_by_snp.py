#!/usr/bin/env python3
# -*- coding: utf-8 -*-
## python3 ./window_by_snp.py -i Target.Broiler.chr1-33.3col.site.fst -o Target.Broiler.chr1-33.3col.site.fst.20SNP.bed -w 20 -t 0.0001 -f 0:1:2 --header False
import sys
import time
import numpy as np
import argparse
def calculate_snp_wise_top_x_value(in_file, top_x):
	top_x = float(top_x)
	with open(in_file, "r") as fp:
			# print(header)
			if header == 'True':
				next(fp)
			value = []
			while True:
				try:
					line = next(fp)
					while missing in line.strip().split()[column_2]:
							# print(line.strip().split()[column_2])
							line = next(fp)
					value.append(float(line.strip().split()[column_2]))
				except StopIteration:
					break
	value_array = np.array(value)
	sort_index = np.argsort(value_array)
	count = np.size(value_array)
	if decrease:
		threshold_index = int(count * (1 - top_x))
	else:
		threshold_index = int(count * top_x)
	threshold_value = value_array[sort_index[threshold_index]]
	return threshold_value

def calculate_mean(value_in):
	value_array = np.array(value_in)
	count = np.size(value_array)
	mean = np.mean(value_array)
	if decrease:
		count_top = np.size(value_array[value_array > threshold_value])
	else:
		count_top = np.size(value_array[value_array < threshold_value])
	return [str(count), str(mean), str(count_top)]

def output(chrom, pos, value, output_file):
	out = [chrom, pos[0], pos[-1]]
	out.extend(calculate_mean(value))
	# out.extend([str(i) for i in value])
	output_file.write('\t'.join(out) + '\n')

def window_output(in_file, out_file, window):
	window = int(window)
	with open(out_file, "w") as wp:
		out = ["chrom", "start", "end", "count", "mean", "count_top"]
		# value_list = [ "snp"+str(i+1) for i in range(window)]
		# out.extend(value_list)
		wp.write('\t'.join(out)+"\n")
		with open(in_file, "r") as fp:
			# print(header)
			if header:
				next(fp)
			
			line = next(fp)
			chrom = line.strip().split()[column_0]
			
			while True:
				try:
					pos = []
					value = []
					i = 1
					while i <= window:
						while missing in line.strip().split()[column_2]:
							line = next(fp)
						if (not line.strip().split()[column_0] == chrom) and i !=1:
							output(chrom, pos, value, wp)
							chrom = line.strip().split()[column_0]
							break
						else:
							chrom = line.strip().split()[column_0]
							pos.append(line.strip().split()[column_1])
							value.append(float(line.strip().split()[column_2]))
							line = next(fp)
							i += 1
					else:
						output(chrom, pos, value, wp)
						i = 1
				except StopIteration:
					if i !=1:
						output(chrom, pos, value, wp)
					break
				
def final_run():
	parser = argparse.ArgumentParser(
	formatter_class=argparse.RawDescriptionHelpFormatter,
	description="Calculate fst, LSBL, etc using n snps window",
	epilog='''
WindowbySNP (1.0)
@author:     Jingfang SI
@copyright:  2018 China Agricultural University. All rights reserved.
@contact:    sijingfang@foxmail.com
	'''
	)
	dict_True_False_conversion = {'True':True, 'False': False}
	
	parser.add_argument("-f", "--format", required = True, default = "0:1:2",
						help = "specify the format of the input file in which there are three necessary columns(chromosome names, positions, and values). \
						e.g. The default, '0:1:2', indicates that the first column(0) are chromosome names, the second column(1) are positions and the third column(2) are values.")
	parser.add_argument("-i", "--input", required = True, help = "specify the input file")
	parser.add_argument("-o", "--output", required = True, help = "specify the output file")
	parser.add_argument("-w","--window", required = True, help = "specify the window by n snps")
	parser.add_argument("-t","--top", required = True, help = "a number which is specified the top x of the SNP-wise values.")
	parser.add_argument("--header", choices = ["True", "False"], default = "True",
						help = "a logical value indicating whether the file contains the names of the variables as its first line. Default True")
	parser.add_argument("-m","--missing", default = "nan", help = "a character vector of strings which are to be interpreted as missing values. Default nan")
	parser.add_argument("--decrease", choices = ["True", "False"], default = "True", 
						help = "a logical value indicating whether the values are sorted increasing or decreasing. Default True")

	args = parser.parse_args()
	
	# get the column number from '--format'.
	format = args.format
	global column_0
	global column_1
	global column_2
	
	column_0 = int(format.split(':')[0])
	column_1 = int(format.split(':')[1])
	column_2 = int(format.split(':')[2])
	
	in_file = args.input
	out_file = args.output
	window = args.window
	top_x = args.top
	global header
	header = dict_True_False_conversion[args.header]
	global missing
	missing = args.missing
	global decrease
	decrease = dict_True_False_conversion[args.decrease]
	
	# calculate the threshold value of top x in snp-wise
	if decrease:
		print("Calculate the top {0} threshold of SNP-wise. Using the decreasing order.".format(top_x)) 
	else:
		print("Calculate the top {0} threshold of SNP-wise. Using the increasing order.".format(top_x))
	global threshold_value
	threshold_value = calculate_snp_wise_top_x_value(in_file, top_x)
	print("The top {0} threshold value is {1}.".format(top_x, threshold_value)) 
	
	
	# calculate the mean value in window
	print("Calculate the window wise mean values.")
	window_output(in_file, out_file, window)
	print("Finished.")
if __name__ == '__main__':
	final_run()
