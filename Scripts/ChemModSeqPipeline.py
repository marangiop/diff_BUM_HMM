#!/usr/bin/python
__author__		= "Sander Granneman"
__copyright__	= "Copyright 2018"
__version__		= "0.3.0"
__credits__		= ["Sander Granneman"]
__email__		= "sgrannem@staffmail.ed.ac.uk"
__status__		= "beta"

from ruffus import *
from ruffus.cmdline import MESSAGE
from collections import defaultdict
import ruffus.cmdline as cmdline
import subprocess
import glob
import re
import os
import argparse
import sys
import time

DEFAULTADAPT = "AGATCGGAAGAGCACACG"
parser = cmdline.get_argparse(description="ChemModSeqPipeline, a pipeline for processing ChemModSeq data")
parser.add_argument("-f",dest="forwardreads",help="the path to your fastq files with forward reads.",metavar="data_1.fastq",nargs="*",default=None)
parser.add_argument("-r",dest="reversereads",help="the path to your fastq files with reverse reads.",metavar="data_2.fastq",nargs="*",default=None)
parser.add_argument("-g","--gtf",dest="gtf",help="the path to your gtf annotation file",metavar="rRNA.gtf",default=None)
parser.add_argument("-c","--chromosome",dest="chromosome",help="the path to your chromosome length file",metavar="chromosome_lengths.txt",default=None)
#parser.add_argument("-t","--trim",dest="trimnumber",type=int,metavar="1",help="to set the number of nucleotides you want to trim after barcode removal. Default = 1.",default=1)
parser.add_argument("-i","--ignoremutations",dest="ignoremuts",action="store_true",help="to tell the novoalign parser to ignore mutations. Useful when analyzing ChemModSeq data. Default is OFF",default=False)
parser.add_argument("--novoindex",dest="novoindex",help="the path to your novoindex file",metavar="yeast.novoindex",default=None)
parser.add_argument("--name",dest="name",help="provide a single word describing the run. Default is 'analysis' with a time stamp",default="analysis_%s" % time.strftime("%d%m%Y"))
parser.add_argument("-b","--barcodes",dest="barcodes",help="the path to the file containing the list of barcodes. Default is None",metavar="barcodes.txt",default=None)
parser.add_argument("-a","--adapter",dest="adapter",help="provide the name of the file containing the 3' adapter sequence. This is for trimming purposes. Default is None",default=None)
parser.add_argument("-m","--mismatches",dest="mismatches",type=int,help="indicate how many mismatches you allow for demultiplexing. Default is 1",default=1)
parser.add_argument("-p","--processors",dest="processors",type=int,help="indicate how many processors you want to use for analyses. Default is 8",default=8)
args = parser.parse_args()

def getFileBaseName(filename):
	""" removes path and extension from file name """
	return os.path.splitext(os.path.basename(filename))[0]

def getFullFilePath(filename):
	""" grabs the full path name for the file """
	return os.path.abspath(filename)

def getBarcodeInfo(barcodefile):
	return ["%s.fastq" % "_".join(line.strip().split()) for line in open(barcodefile,"r").readlines()]

def decompressFiles(inputfile,outputfile):
	""" decompresses the inputfile """
	cmd = "gzip -d %s" % inputfile
	logger.info(cmd)
	os.system(cmd)

def collateFiles(input):
	""" collates the two fastq input files. Forward read first """
	filedict = defaultdict(lambda: defaultdict(str))
	for filename in input:
		filesearch = re.search(r".+/\w+_(\d)_(\w+).fastq",filename)
		readorientation = filesearch.group(1)
		barcode = filesearch.group(2)
		if barcode == "others": continue
		filedict[barcode][readorientation] = filename
	output = list()
	for barcode in filedict:
		output.append([filedict[barcode]["1"],filedict[barcode]["2"]])
	return output
	
def runFlexBar(inputfiles,outputfiles):
	""" runs Flexbar on the data to remove the adapter sequence from the forward reads """
	outputfilename = "%s/%s_trimmed" % (os.path.join(root_dir,"flexbar_trimmed"),re.search("^.+/([^/]+)_1.(san)?fastq$",inputfiles[0]).group(1))
	cmd = "flexbar -r '%s' -p '%s' -f i1.8 -t '%s' -n 10 -ao 7 -a '%s' --pre-trim-phred 30" % (inputfiles[0],inputfiles[1],outputfilename,args.adapter)
	logger.info(cmd)
	os.system(cmd)
	#for i in outputfiles: os.system("touch '%s'" % i)

def demultiplexSamples(inputfiles,outputfiles):
	""" demultiplexes all the samples """
	os.chdir(os.path.join(root_dir,"demultiplexed"))
	cmd = "pyBarcodeFilter.py -f '%s' -r '%s' -b '%s' -m '%s'" % (inputfiles[0],inputfiles[1],os.path.join(home_dir,args.barcodes),args.mismatches)
	logger.info(cmd)
	os.system(cmd)
	#for i,j in outputfiles: 
	#	 os.system("touch '%s'" % i)
	#	 os.system("touch '%s'" % j)
	os.chdir(root_dir)

def trimAndCollapse(inputfiles,outputfiles):
	""" Removes the last random nucleotide from the forward read, collapses the data and then splits it again into fasta files """
	outputfilename = "%s/%s" % (os.path.join(root_dir,"collapsed"),re.search("(.*)_1.fasta",os.path.basename(outputfiles[0])).group(1))
	cmd = "pyFastqJoiner.py -f '%s' '%s' -c '|' | pyFastqDuplicateRemover.py | pyFastqSplitter.py -c '|' --file_type=fasta -o '%s'" % (inputfiles[0],inputfiles[1],args.trimnumber,outputfilename)
	logger.info(cmd)
	os.system(cmd)
	#for i in outputfiles: 
	#	 os.system("touch '%s'" % i)
	
def alignReads(inputfiles,outputfile):
	""" Runs novoalign on all the collapsed files"""
	cmd = "novoalign -d '%s' -f '%s' '%s' > '%s'" % (args.novoindex,inputfiles[0],inputfiles[1],outputfile)
	logger.info(cmd)
	os.system(cmd)
	#os.system("touch '%s'" % outputfile)

def runPyReadCounters(inputfile,outputfile):
	""" runs pyReadCounters on all the files but only selects properly paired reads"""
	outputfile = re.search("(.*)_count_output_reads.gtf",outputfile).group(1)
	if args.ignoremuts:
		string = "--mutations=nomuts"
	else:
		string = ""
	cmd = "pyReadCounters.py -f '%s' --gtf '%s' -v --rpkm --unique --properpairs -o '%s' '%s'" % (inputfile,args.gtf,outputfile,string)
	logger.info(cmd)
	os.system(cmd)
	#os.system("touch '%s'_count_output_reads.gtf" % outputfile)
	
def makeCoverageSgrFiles(inputfile,outputfiles):
	""" takes the pyReadCounters gtf output files and generates coverage sgr files """
	outputfile = re.search("(.*)_plus_strand_reads.sgr",outputfiles[0]).group(1)
	cmd = "pyGTF2sgr.py --gtf '%s' --zeros --count -v -o '%s' -c '%s'" % (inputfile,outputfile,args.chromosome)
	logger.info(cmd)
	os.system(cmd)
	#for i in outputfiles:
	#	os.system("touch '%s'" % i)
	
def makeDropoffCountsSgrFiles(inputfile,outputfiles):
	""" takes the pyReadCounters gtf output files and generates coverage sgr files """
	outputfile = re.search("(.*)_plus_strand_startpositions.sgr",outputfiles[0]).group(1)
	cmd = "pyGTF2sgr.py --gtf '%s' --zeros --count --min 5 -v -o '%s' -c '%s' --type startpositions" % (inputfile,outputfile,args.chromosome)
	logger.info(cmd)
	os.system(cmd)
	#for i in outputfiles:
	#	os.system("touch '%s'" % i)
	
def makeCoverageBedgraphFiles(inputfile,outputfiles):
	""" takes the pyReadCounters gtf output files and generates bedgraph files for viewing of the data in genome browsers """
	outputfile = re.search("(.*)_plus_strand_reads.bedgraph",outputfiles[0]).group(1)
	cmd = "pyGTF2bedGraph.py --gtf '%s' --count -v --permillion -o '%s' -c '%s'" % (inputfile,outputfile,args.chromosome)
	logger.info(cmd)
	os.system(cmd)
	#for i in outputfiles:
	#	os.system("touch '%s'" % i)

### checking if the correct pyBarcodeFilter version is installed. Need version 3.0 or higher!

result = subprocess.Popen(['pyBarcodeFilter.py', '--version'], stdout=subprocess.PIPE)
out,err = result.communicate()
print float(out)
if float(out) < 3.0:
	sys.stderr.write("To run this script you need to have pyBarcodeFilter version 3.0 or later installed. Please install the latest version of pyCRAC\n")
	exit()	

### start of pipeline

# setting the main working directory
home_dir = os.getcwd()
root_dir = "%s/%s" % (home_dir,args.name)
try:
	os.mkdir(root_dir)
except OSError:
	pass

### setting up logging
if not args.log_file:
	logfile = "%s/%s" % (root_dir,"log_file.txt")
else:
	logfile = "%s/%s" % (root_dir,args.log_file)
logger, logger_mutex = cmdline.setup_logging ("ChemModSeqPipeline",logfile,10)

### starting the pipeline commands

pipeline = Pipeline(name = "ChemModSeqPipeline")

args.forwardreads = [os.path.abspath(i) for i in args.forwardreads]
args.reversereads = [os.path.abspath(i) for i in args.reversereads]
startingfiles = zip(args.forwardreads,args.reversereads)
#print startingfiles

if args.adapter:
	pipeline.transform(task_func = runFlexBar,
						input  = startingfiles,
						filter = formatter("^.+/([^/]+)_1.(san)?fastq$","^.+/([^/]+)_2.(san)?fastq$"),						
						output = ["%s/flexbar_trimmed/{1[0]}_trimmed_1.fastq" % root_dir,"%s/flexbar_trimmed/{1[1]}_trimmed_2.fastq" % root_dir],
					).follows(pipeline.mkdir(os.path.join(root_dir,"flexbar_trimmed")))
	startingfiles = runFlexBar

if args.barcodes:
	pipeline.subdivide(task_func = demultiplexSamples,
						input  = startingfiles,
						filter = formatter("^.+/([^/]+)_1.fastq$","^.+/([^/]+)_2.fastq$"),
						output = [["%s/demultiplexed/{1[0]}_1_%s" % (root_dir,barcodestring),"%s/demultiplexed/{1[1]}_2_%s" % (root_dir,barcodestring)] for barcodestring in getBarcodeInfo(args.barcodes)]
					).follows(pipeline.mkdir(os.path.join(root_dir,"demultiplexed")))

	pipeline.transform(task_func = trimAndCollapse,
						input  = demultiplexSamples,
						filter = formatter("^.+/([^/]+)_trimmed_1_(?P<BARCODE>.*).fastq","^.+/([^/]+)_trimmed_2_(?P<BARCODE>.*).fastq"),						
						output = ["%s/collapsed/{1[0]}_trimmed_{BARCODE[0]}_1.fasta" % root_dir,"%s/collapsed/{1[1]}_trimmed_{BARCODE[1]}_2.fasta" % root_dir],
					).follows(pipeline.mkdir(os.path.join(root_dir,"collapsed"))).follows(demultiplexSamples)

else:
	pipeline.transform(task_func = trimAndCollapse,
						input  = startingfiles,
						filter = formatter("^.+/([^/]+).fastq","^.+/([^/]+).fastq"),						
						output = ["%s/collapsed/{1[0]}.fasta" % root_dir,"%s/collapsed/{1[1]}.fasta" % root_dir],
					).follows(pipeline.mkdir(os.path.join(root_dir,"collapsed")))
					
pipeline.transform( task_func = alignReads,
					input  = trimAndCollapse,
					filter = formatter("^.+/([^/]+)_1.fasta","^.+/([^/]+)_2.fasta"),
					output = "%s/novo_files/{1[0]}.novo" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"novo_files")))

pipeline.transform( task_func = runPyReadCounters,
					input  = alignReads,
					filter = formatter("^.+/([^/]+).novo"),
					output = "%s/pyReadCounters_analyses/{1[0]}_count_output_reads.gtf" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"pyReadCounters_analyses")))

pipeline.subdivide( task_func = makeCoverageSgrFiles,
					input  = runPyReadCounters,
					filter = formatter("^.+/([^/]+)_count_output_reads.gtf"),
					output = ["%s/sgr_files/{1[0]}_plus_strand_reads.sgr" % root_dir,"%s/sgr_files/{1[0]}_minus_strand_reads.sgr" % root_dir]
				).follows(pipeline.mkdir(os.path.join(root_dir,"sgr_files"))).follows(runPyReadCounters)

pipeline.subdivide( task_func = makeDropoffCountsSgrFiles,
					input  = runPyReadCounters,
					filter = formatter("^.+/([^/]+)_count_output_reads.gtf"),
					output = ["%s/sgr_files/{1[0]}_plus_strand_startpositions.sgr" % root_dir,"%s/sgr_files/{1[0]}_minus_strand_startpositions.sgr" % root_dir]
				).follows(runPyReadCounters)
				
pipeline.subdivide( task_func = makeCoverageBedgraphFiles,
					input  = runPyReadCounters,
					filter = formatter("^.+/([^/]+)_count_output_reads.gtf"),
					output = ["%s/bedgraph_files/{1[0]}_plus_strand_reads.bedgraph" % root_dir,"%s/bedgraph_files/{1[0]}_minus_strand_reads.bedgraph" % root_dir]
				).follows(pipeline.mkdir(os.path.join(root_dir,"bedgraph_files"))).follows(runPyReadCounters)

pipeline_run(multiprocess=args.processors,verbose=5)