#!/usr/bin/python
# not compatible with python 3
__author__		= "Sander Granneman"
__copyright__	= "Copyright 2017"
__version__		= "2.0"
__credits__		= ["Sander Granneman"]
__maintainer__	= "Sander Granneman"
__email__		= "sgrannem@staffmail.ed.ac.uk"
__status__		= "beta"

import sys
import re
from optparse import *
from collections import defaultdict
from pyCRAC.Methods import getfilename
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = ["Arial"]

colordict = {"gray":"#CCCCCC","black":"#000000","yellow":"#CCCC00","red":"#FF0000","cyan":"#00FFFF","blue":"#0000FF"}

def getNucleotideColorInfo(data,column=0,colormap="jet",minvalue=None,maxvalue=None,filter=[]):
	""" processes the file containing the nucleotide positions and the relative reactivity"""
	column = data.columns[column]
	selecteddata = data[column].values
	colorvalues = list()
	colordata = defaultdict(lambda: [0,0,0])
	if colormap == "SHAPE":
		values = [-999,0.0,0.40,0.85]
		for i in selecteddata:
			if i == -999:
				colorvalues.append(colordict["gray"])
			elif (i >= 0.0) and (i < 0.40):
				colorvalues.append(colordict["black"])
			elif (i >= 0.40) and (i < 0.85):
				colorvalues.append(colordict["yellow"])
			elif (i > 0.85):
				colorvalues.append(colordict["red"])
			else:
				sys.stderr.write("Could not find the right color for value %s" % i)
				colorvalues.append(colordict["black"])
	elif colormap == "binary":
		values = [-1,0,1]
		for i in selecteddata:
			if i == 0:
				colorvalues.append(colordict["black"])
			elif i == -1:
				colorvalues.append(colordict["blue"])
			elif i == 1:
				colorvalues.append(colordict["red"])
			else:
				sys.stderr.write("Could not find the right color for value %s" % i)
				colorvalues.append(colordict["black"])
	elif colormap == "DELTASHAPE":
		values = [-999,-0.30,-0.15,0,0.15,0.30]										# -999 values are positions were no reactivity could be calculated due to low coverage or drop-off counts.
		for i in selecteddata:
			if i == -999:
				colorvalues.append(colordict["gray"])
			elif (i <=-0.30):
				colorvalues.append(colordict["blue"])
			elif (i > -0.30) and (i <= -0.15):
				colorvalues.append(colordict["cyan"])
			elif (i > -0.15) and (i < 0.0):
				colorvalues.append(colordict["black"])
			elif (i >= 0.0) and (i < 0.15):
				colorvalues.append(colordict["black"])
			elif (i >= 0.15) and (i < 0.30):
				colorvalues.append(colordict["yellow"])
			elif (i >= 0.30):
				colorvalues.append(colordict["red"])
			else:
				sys.stderr.write("Could not find the right color for value %s" % i)
				colorvalues.append(colordict["black"])
	elif colormap == "BUM_HMM":
		values = [-999,0.0,0.75,0.95]												# -999 values are positions were no reactivity could be calculated due to low coverage or drop-off counts.
		for i in selecteddata:
			if i == -999:
				colorvalues.append(colordict["gray"])
			elif (i >= 0.0) and (i < 0.75):
				colorvalues.append(colordict["black"])
			elif (i >= 0.75) and (i < 0.95):
				colorvalues.append(colordict["yellow"])
			elif (i >= 0.95):
				colorvalues.append(colordict["red"])
			else:
				sys.stderr.write("Could not find the right color for value %s" % i)
				colorvalues.append(colordict["black"])
	else:
		cdict = cm.datad[colormap]
		my_cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)
		nonzerodata = np.array(selecteddata)
		nonzerodata[selecteddata < 0.0] = 0.0
		colorvalues = my_cmap(nonzerodata)
		for nr,i in enumerate(selecteddata,start=1):
			try:
				if i == 0.0:  # if the value is zero, color the nucleotides black.
					colorvalues[nr] = [0.0, 0.0, 0.0, 1]
				if i == -999.0:														# -999 values are positions were no reactivity could be calculated due to low coverage or drop-off counts.
					colorvalues[nr] = [0.5, 0.5, 0.5, 1]
				if filter:
					first,second = filter											# color nucleotides with values in a specified range black
					if i >= first and i <= second:
						colorvalues[nr] = [0.0, 0.0, 0.0, 1]
				if minvalue and i < minvalue: colorvalues[nr] = [0.0, 0.0, 0.0, 1]
				if maxvalue and i > maxvalue: colorvalues[nr] = [0.0, 0.0, 0.0, 1]
			except IndexError:
				break

	colordata.update(dict(zip(data.index,colorvalues)))								# each color list is associated with a specific position.
	return colordata

def colorSVGFile(inputfile,nucinfo,outfile=None):
	""" function for coloring nucleotides in the Ribosome Gallery svg files """
	out = sys.stdout
	if outfile: out = open(outfile,"w")
	count = 0
	with open(inputfile,"r") as svgfile:
		for line in svgfile:
			if re.search("<text id=",line) and re.search(">[AGCTU]</text>",line):
				text_id = re.search("<text id=\"(\d+)\"",line)
				position = int(text_id.group(1))
				try:
					colorstring = nucinfo[position]
					line = re.sub("fill=\"[#0-9a-zA-Z]+\"","fill=\"%s\"" % colorstring,line)
				except IndexError:
					sys.stderr.write("IndexError at position %s\n" % position)
					pass
				out.write(line)
			elif re.search("<text onmouseover=\"mOver\(evt,",line):			# PseudoViewer .svg file.
				basepaired = re.search("<text.*\(([0-9]+)\)",line)
				singlestranded = re.search("<text.*\'([0-9]+)\'",line)
				if basepaired:
					position = int(basepaired.group(1))
				elif singlestranded:
					position = int(singlestranded.group(1))
				result = re.search("(>([AUCG])</text>)",line)
				selection = result.group(1)
				nucleotide = result.group(2)
				try:
					colorstring = nucinfo[position]
					line = line.replace(selection," fill=\"%s\"%s" % (colorstring,selection))
				except IndexError:
					sys.stderr.write("IndexError at position %s\n" % position)
					pass
				out.write(line)
			elif re.search("<text x=",line) and re.search(">[AGCTU]</text>",line):		# VARNA svg files
				count += 1
				position = count
				try:
					colorstring = nucinfo[position]
					line = re.sub("fill=\"[a-zA-Z0-9]+\"","fill=\"%s\"" % colorstring,line)
				except IndexError:
					sys.stderr.write("IndexError at position %s\n" % position)
					pass
				out.write(line)
			else:
				out.write(line)

def printColorValues(inputfile,nucinfo,outfile=None):
	""" prints out the color value for each data point """
	out = sys.stdout
	if outfile: out = open(outfile,"w")
	with open(inputfile,"r") as datafile:
		for line in datafile:
			if line.startswith("#"): continue
			Fld = line.strip().split()
			position = int(Fld[0])
			value = Fld[1]
			colorvalue = nucinfo[(position-1)]
			out.write("%s\t%s\t%s\n" % (position,value,",".join([str(i) for i in colorvalue])))

def main():
	parser = OptionParser(usage="usage: %prog [options] -f svgfile.svg -d datafile.txt -o svgfile.svg -c SHAPE --colorbar --filetype=svg", version="%s" % __version__)
	parser.add_option("-f", dest="file",help="Type the path of the svg, eps or postscript (ps) input file. Note that not all different svg or eps formats are supported", metavar="mysecondarystructure.svg")
	parser.add_option("-d", "--datafile", dest="datafile",help="Type the path of the data file. This should contain two columns, one with the nucleotide positions, the other with the color coding index",metavar="indexfile.txt",default=None)
	parser.add_option("-o", "--outfile", dest="outfile",help="Type the name of your output file. Don't forget to add the .ps extension", metavar="mycoloredpostscriptfile.ps",default=None)
	parser.add_option("-c", "--colormap", dest="colormap",help="to select a specific matplotlib color gradient for coloring. To use the Weeks lab coloring scheme for SHAPE data set colormap to 'SHAPE'. 0 to 0.40 = black, 0.40 to 0.85 = yellow, >0.85 is red. Default is a discrete coloring black, blue, green, yellow, orange, red", metavar="SHAPE",default="jet")
	parser.add_option("--column",dest="column",type="int",help="specify which data column you would like to use. Default is column 2",default=0)
	parser.add_option("--header",dest="header",action="store_true",help="use this flag if your reactivities file have column headers. Default is OFF",default=False)
	parser.add_option("--filter",dest="filter",type="float",nargs=2,help="with this flag you can set a range of values that you do not want to have colored. Expects two values, smallest one should be first",default=None)
	parser.add_option("--max",dest="max",type="float",help="to set a maximum value for your coloring. This allows you to remove outliers",default=None)
	parser.add_option("--min",dest="min",type="float",help="to set a maximum value for your coloring. This allows you to change the lower end of the colorscale",default=None)
	parser.add_option("--printcolorvalues",action="store_true",dest="printcolorvalues",help=SUPPRESS_HELP,default=False)
	(options, args) = parser.parse_args()
	if not options.datafile:
		parser.error("No data file was provided!\n")
	if not options.file:
		parser.error("You forgot to add the path to the original structure (eps,ps or svg) file\n")
	if not options.outfile:
		options.outfile = "%s.svg" % (getfilename(options.datafile))
	columnheaders = None
	if options.header:
		columnheaders = 0
	data = pd.read_table(options.datafile,index_col=0,header=columnheaders)
	nucinfo = getNucleotideColorInfo(data,column=options.column,colormap=options.colormap,minvalue=options.min,maxvalue=options.max,filter=options.filter)
	if options.printcolorvalues:
		colorvaluefile = "%s_colorvalues.txt" % getfilename(options.data)
		printColorValues(options.data,nucinfo,colorvaluefile)
	colorSVGFile(options.file,nucinfo,outfile=options.outfile)

if __name__ == "__main__":
	main()
