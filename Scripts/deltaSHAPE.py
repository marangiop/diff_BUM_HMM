#!/usr/bin/env python
#  deltaSHAPE software for detecting meaningful changes in SHAPE reactivity
#
#  - Requires two .map files as input (see README for details)
#  - See the README for required modules, installation, and execution help.
#  - Version 1.0
#  - Copyright Matthew J. Smola 2015

###########################################################################
# GPL statement:                                                          #
#                                                                         #
# This program is free software: you can redistribute it and/or modify    #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# This program is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.   #
###########################################################################


import sys
import os
import argparse
import numpy as np
import scipy.stats as stats
from operator import itemgetter
from itertools import groupby
import warnings
from numpy import nanmean
from numpy import nanstd

def open_map(filename, front, back):
    data = open(filename, "rU")
    datalines = data.readlines()
    data.close()
    if len(datalines[0].split()) < 4:
        sys.exit('ERROR: Input file '+filename+' contains fewer than 4 data columns. Please provide a file in .map format. See the README for details.')
    data = np.array([float(i.split()[1]) for i in datalines])
    errs = np.array([float(i.split()[2]) for i in datalines])
    seq = [str(i.split()[3]) for i in datalines]
    sequence = ''.join(seq)
    # mask 3' primer nucleotides
    for i in range(len(data)-back, len(data)):
        data[i] = -999
        errs[i] = 0
    # mask 5' nucleotides
    for i in range(front):
        data[i] = -999
        errs[i] = 0
    # -999 data points as 'nan', including 5' and 3' masked regions
    for i in range(len(data)):
        if data[i] == -999:
            data[i] = np.nan
            errs[i] = np.nan
    return data, errs, sequence

def smooth(data,err,pad):
    new_data, new_err = [], []
    # eventually we want to exclude no-data nucleotides.
    # create a list ("mask") to store which positions to ignore later.
    mask = []
    for i in range(len(data)):
        if data[i] == -999 or np.isnan(data[i]) == True:
            mask.append(i)
    # you can't center a window at the first nucleotide so mask until a full centered window can be placed
    for i in range(pad):
        new_data.append(np.nan)
        new_err.append(np.nan)
    # proceed by windows to smooth the data
    for i in range(pad, len(data)-pad):
        # use numpy masked array to calculate average without including no-data (nan) nucleotides.
        new_data.append(np.mean(np.ma.MaskedArray([j for j in data[i-pad:i+pad+1]], np.isnan([j for j in data[i-pad:i+pad+1]]))))
        
        # use stats.nanmean to calculate average without including no-data (nan) nucleotides. This causes long_scalars runtime warnings.
        #new_data.append(stats.nanmean([j for j in data[i-pad:i+pad+1] if np.isnan(j) != True]))
        errs = np.array(err[i-pad:i+pad+1])
        squerrs = np.power([j for j in errs if np.isnan(j) != True], 2)
        total = np.sum(squerrs)
        sqrt = np.sqrt(total)
        new_err.append(sqrt/len(data[i-pad:i+pad+1]))
    for i in range(pad):
        new_data.append(np.nan)
        new_err.append(np.nan)
    for i in mask:
        new_data[i] = np.nan
        new_err[i] = np.nan
    return np.array(new_data), np.array(new_err)

def z_factor(data1, data2, err1, err2, factor=1.96):
    z_factors = []
    for i in range(len(data1)):
        if data1[i] == 'nan' or data2[i] == 'nan':
            z_factors.append(float('nan'))
        else:
            top = factor * (err2[i] + err1[i]) #1.645 = 90% confidence interval
            bot = abs(data2[i] - data1[i])
            if bot == 0:
                z_factors.append(float('nan'))
            else:
                z = (1 - (top / bot))
                z_factors.append(z)
    return z_factors

def calc_zScores(diffs):
    mean = nanmean(diffs)
    sigma = nanstd(diffs)
    # calc Z-score
    z_scores = (diffs - mean) / sigma
    return np.array(z_scores)


if __name__ == '__main__':
    
    #############################################
    ## Set up arguments #########################
    #############################################
    
    # pre-parse the arguments so that argparse doesn't interpret negative numbers (for y-min) as option flags.
    for i, arg in enumerate(sys.argv):
      if (arg[0] == '-') and arg[1].isdigit():
          sys.argv[i] = ' ' + arg
    
    # parse the arguments
    parse = argparse.ArgumentParser(
        description="deltaSHAPE computes statistically significant changes in SHAPE-MaP reactivity between two conditions. See README file for further details and file descriptions.", 
        epilog="deltaSHAPE v0.91 by Matt Smola ( matt.smola@gmail.com )",
        add_help=False)
    #parse.negative_number_matcher = _re.compile(r'^-(\d+\.?|\d*\.\d+)([eE][+\-]?\d+)?$')
    
    required = parse.add_argument_group('Required files', 'These files are required in order to run deltaSHAPE analysis.')
    required.add_argument('mapFile1', type=str, help='SHAPE-MaP .map file')
    required.add_argument('mapFile2', type=str, help='SHAPE-MaP .map file, values in this file will be subtracted from those in mapFile1')
    
    data_opt = parse.add_argument_group('Data manipulation', 'Options to specify how SHAPE-MaP data are manipulated and analyzed.')
    data_opt.add_argument('--mask5', type=int, default=0, help="Specify the number of nucleotides at the 5' end to ignore. Default: 0")
    data_opt.add_argument('--mask3', type=int, default=0, help="Specify the number of nucleotides at the 3' end to ignore. Default: 0")
    data_opt.add_argument('-p', '--pad', type=int, default=1, help='Indicate the smoothing window size. Window = 2*pad+1. To turn off smoothing, set PAD = 0. Default: 1')
    data_opt.add_argument('-z', '--Zcoeff', type=float, default=1.96, help='Ajust the Z-factor stringency by changing the equation coefficient. See the README for details. Default: 1.96')
    data_opt.add_argument('-t', '--Zthresh', type=float, default=0, help='Adjust the Z-factor stringency by changing the cutoff threshold. See the README for details. Default: 0')
    data_opt.add_argument('-s', '--SSthresh', type=float, default=1, help='Set the cutoff threshold of standard score filtering. Default: 1.0')
    data_opt.add_argument('-f', '--FindSite', type=str, default='2,3', help='Comma-separated pair of numbers indicating the window pad size and number of required hits when finding binding sites. Default settings look for 3+ nucleotides within a 5-nucleotide window. See the README for details. Default: 2,3')
    
    out_opt = parse.add_argument_group('Output', 'Options specifying plotting and output file details.')
    out_opt.add_argument('-o', '--out', type=str, default="differences.txt", help='Name and location of output file to be written. Default: ./differences.txt')
    out_opt.add_argument('--magrank', action='store_true', help='Sort output file by decreasing deltaSHAPE magnitude. Default: OFF')
    out_opt.add_argument('--all', action='store_true', help='Output data for all nucleotides. Insignificant changes are listed as zero. Default: OFF')
    out_opt.add_argument('--pdf', action='store_true', help='Save plot as PDF. If output file is given, PDF will have same prefix. Default: OFF')
    out_opt.add_argument('--noshow', action='store_true', help='Generate the plot but do not show it. Typically used with --pdf. Default: display plot')
    out_opt.add_argument('--noplot', action='store_true', help='Skip plotting completely. Default: OFF')
    out_opt.add_argument('--dots', action='store_true', help='Plot markers indicating nucleotides that pass Z-factor and standard score filtering. This can get unweildy for large RNAs (>1000). Standard score (open) dots are plotted above Z-factor (filled) dots. Default: OFF')
    out_opt.add_argument('--Zdots', action='store_true', help='Plot markers indicating only nucleotides that pass Z-factor filtering. Default: OFF')
    out_opt.add_argument('--SSdots', action='store_true', help='Plot markers indicating only nucleotides that pass standard score filtering. Default: OFF')
    out_opt.add_argument('--colorfill', action='store_true', help='Highlight deltaSHAPE sites with coloration beneath the plot line for "prettier" figures. Default: OFF')
    out_opt.add_argument('--ymin', type=float, default=-999, help='Set plot y-axis minimum. Default: Determined automatically')
    out_opt.add_argument('--ymax', type=float, default=-999, help='Set plot y-axis maximum. Default: Determined automatically')
    out_opt.add_argument('--xmin', type=float, default=-999, help='Set plot x-axis minimum. Default: Determined automatically')
    out_opt.add_argument('--xmax', type=float, default=-999, help='Set plot x-axis maximum. Default: Determined automatically')
    
    help_opt = parse.add_argument_group('Help')
    help_opt.add_argument('-h', '--help', action="help", help="show this help message and exit")
    args = parse.parse_args()
    #############################################
    
    
    #######################################################
    ## Set up and check analysis parameters ###############
    #######################################################
    
    # front and back refer to how many nucleotides should be masked at the 5' and 3' end, respectively.
    front = args.mask5
    back = args.mask3
    
    # pad specifies how wide the smoothing window will be. Window size = 2*pad + 1
    pad = args.pad
    
    # z_coeff = A in Z-factor = A*(err1 + err2)/abs(data2-data1)
    # z_thresh is the minimum Z-factor required to pass the test.
    z_coeff = args.Zcoeff
    z_thresh = args.Zthresh
    if z_thresh > 1:
        sys.exit('ERROR: Z-factor can never exceed 1. Change -t/--Zthresh value accordingly.')
    
    # ss_thresh is the minimum standard score (Z-score) required to pass the test.
    ss_thresh = args.SSthresh
    
    # site_pad specifies window size used when searching for significant hits. Window size = 2*site_pad + 1.
    # site_min sets the required number of hits within the window.
    # e.g., requiring 3 nts within a 5-nt window has site_pad=2 and site_min=3.
    site_pad = int(args.FindSite.split(',')[0]) # leads to window size used in searching for significant hits
    site_min = int(args.FindSite.split(',')[1]) # set minimum number of hits within window
    if site_pad * 2 + 1 < site_min:
        sys.exit('ERROR: Binding site window size and hit minimum are incompatible.\nDouble-check the -f --FindSite flag or consult the README file.')
    
    # determine how to color the plot
    color=''
    if not args.colorfill:
        color='bar'
    elif args.colorfill:
        color='fill'
    
    # output filenames and prefixes
    outfile = os.path.normpath(args.out)
    if args.pdf:
        #pdf_file = str(os.path.basename(os.path.normpath(args.out)).split('.')[:-1][0])+".pdf"
        pdf_file = str(os.path.normpath(args.out)).split('.')[0]+".pdf"
    
    # check other variables
    if args.ymin > args.ymax:
        sys.exit('ERROR: --ymin must be less than --ymax.')
    if args.xmin > args.xmax:
        sys.exit('ERROR: --xmin must be less than --xmax.')
    #######################################################
    
    
    #######################################################
    ## Run the analysis ###################################
    #######################################################
    
    '''STEP ONE'''
    # open .map files
    data1, err1, seq1 = open_map(args.mapFile1, front, back)
    data2, err2, seq2 = open_map(args.mapFile2, front, back)
    
    #plt.figure("Smoothed Reactivities")
    #plt.plot(range(len(data1)), s_data1, drawstyle='steps-mid', color='red')
    #plt.plot(range(len(data2)), s_data2, drawstyle='steps-mid', color='blue')

    '''STEP TWO'''
    # smooth data and errors
    s_data1, s_err1 = smooth(data1, err1, pad)
    s_data2, s_err2 = smooth(data2, err2, pad)
    
    #plt.figure("Smoothed Reactivities")
    #plt.plot(range(len(data1)), s_data1, drawstyle='steps-mid', color='red')
    #plt.plot(range(len(data2)), s_data2, drawstyle='steps-mid', color='blue')

    '''STEP THREE'''
    # subtract raw and smoothed data
    diff = data1 - data2
    s_diff = smooth(diff, err1, pad)[0]

    '''STEP FOUR'''
    # calculate Z-factors from smoothed data and smoothed errs
    z_factors = z_factor(s_data1, s_data2, s_err1, s_err2, z_coeff)

    '''STEP 5'''
    # calculate Z-scores from difference of smoothed data1 and smoothed data2
    z_scores = calc_zScores(s_diff)

    '''STEP 6'''
    # identify 5-nt windows where 3+ nts are sig. diff.
    # (the window size and number of required hits can be changed with the --FindSite option flags)
    
    sigdiff = []
    for i in range(site_pad, len(diff)-site_pad):
        win = range(i-site_pad ,i+site_pad+1)
        count = 0
        maybes = []
        for j in win:
            if z_factors[j] > z_thresh and np.abs(z_scores[j]) >= ss_thresh:
                count += 1
                maybes.append(j)
        if count >= site_min:
            for k in maybes:
                if k not in sigdiff:
                    sigdiff.append(k)


    '''STEP 7 --- PREPARE FOR PLOTTING & OUTPUTTING'''
    # this is mostly for figuring which regions to highlight in the plot.
    
    pos_consec, neg_consec = [], []
    for k, g in groupby(enumerate([i for i in sigdiff if s_diff[i] >= 0]), lambda(i,x):i-x):
            pos_consec.append(map(itemgetter(1), g))
    for k, g in groupby(enumerate([i for i in sigdiff if s_diff[i] < 0]), lambda(i,x):i-x):
            neg_consec.append(map(itemgetter(1), g))

    pos_shade_bits, pos_x_bits = [], []
    pos_span = []
    data_out = []
    for region in pos_consec:
        pos_shade, pos_x = [], []
        for i in region:
            pos_shade.extend((s_diff[i], s_diff[i]))
            pos_x.extend((i+0.5, i+1.5))
            data_out.append([i+1, seq1[i], s_diff[i], z_factors[i], z_scores[i], s_data1[i], s_data2[i], diff[i], data1[i], data2[i]])
        pos_shade_bits.append(pos_shade)
        pos_x_bits.append(pos_x)
        pos_span.append([region[0]+.5, region[-1]+1.5])
    
    neg_shade_bits, neg_x_bits = [], []
    neg_span = []
    for region in neg_consec:
        neg_shade, neg_x = [], []
        for i in region:
            neg_shade.extend((s_diff[i], s_diff[i]))
            neg_x.extend((i+0.5, i+1.5))
            data_out.append([i+1, seq1[i], s_diff[i], z_factors[i], z_scores[i], s_data1[i], s_data2[i], diff[i], data1[i], data2[i]])
            #outfile.write(outline)
        neg_shade_bits.append(neg_shade)
        neg_x_bits.append(neg_x)
        neg_span.append([region[0]+.5, region[-1]+1.5])
    
    
    '''STEP 8 --- PLOTTING'''
    if not args.noplot:
        
        if args.noshow:
            import matplotlib
            matplotlib.use('Agg')
        
        import matplotlib.pyplot as plt
    
        plt.figure(figsize=(11,4))
        x = range(1,len(s_diff)+1)
        
        plt.plot(x, s_diff, drawstyle='steps-mid', color='black')
        plt.axhline(0, color='black')
    
        # mask primer-binding regions
        plt.axvspan(0,front+0.5, color="grey", alpha=0.25)
        plt.axvspan(len(diff)-back+0.5,len(diff)+0.5, color="grey", alpha=0.25)
    
        # color deltaSHAPE sites
        if color=='fill':
            for i in range(len(pos_shade_bits)):
                plt.fill_between(pos_x_bits[i], pos_shade_bits[i], color='#3EB452')
            for i in range(len(neg_shade_bits)):
                plt.fill_between(neg_x_bits[i], neg_shade_bits[i], color='#7F3B95')
        elif color == 'bar':
            for i in pos_span:
                plt.axvspan(i[0], i[-1], color='#3EB452')
            for i in neg_span:
                plt.axvspan(i[0], i[-1], color='#7F3B95')


        # Are z-factor and standard score dots to be plotted?
        # (this will affect both the plotting and the y-axis limits)
        if args.dots or args.SSdots or args.Zdots:
            dots=True
        else:
            dots=False
        

        # set axes limits
        # default is min=1, max=length of RNA
        if args.xmax == -999:
            args.xmax = len(diff)
        if args.xmin == -999:
            args.xmin = 1
                
        plt.xlim(args.xmin, args.xmax)
        
        # set ymin and ymax automatically from data, or from option flags.
        if args.ymin == -999:
            y_min = min(filter(lambda x: np.isnan(x)==False, s_diff))-0.25
        else:
            y_min = args.ymin

        if args.ymax == -999:
            y_max = max(filter(lambda x: np.isnan(x)==False, s_diff))+0.25
            # adjust y_max if dots are involved.    
            if dots:
                y_max = max(filter(lambda x: np.isnan(x)==False, s_diff))+0.6
        else:
            y_max = args.ymax
        
        plt.ylim(y_min, y_max)
        
        # plot the Z-factor/standard score dots
        if dots:
            for i in range(len(z_scores)):
                if (args.dots or args.SSdots) and abs(z_scores[i]) >= ss_thresh:
                    plt.scatter(i+1, y_max-0.2, marker="s", s=5, color='none', edgecolor='black', zorder=3)
                if (args.dots or args.Zdots) and z_factors[i] >= z_thresh:
                    plt.scatter(i+1, y_max-0.4, marker="o", s=5, color='black', zorder=3)

        # label axes and format ticks
        plt.xlabel("Nucleotide")
        plt.ylabel(r'$\Delta$SHAPE')
        plt.tick_params(which='both', direction='out', top='off', right='off')
        
        # turn off UserWarnings temporarily so that plt.tight_layout() doesn't print a warning to the screen.
        warnings.simplefilter("ignore", UserWarning)
        # set the plot layout
        plt.tight_layout()
        # turn warnings back on in case something terrible happens.
        warnings.resetwarnings()
        
        if args.pdf:
            plt.savefig(pdf_file, format='pdf')
        if not args.noshow:
            plt.show()

    
    '''STEP 9 --- OUTPUT FILE GENERATION'''    
    
    if args.all:
        # get data for all nucleotides not already in data_out, but replace s_diff[i] with zero.
        for i in filter(lambda x: x+1 not in [j[0] for j in data_out], range(len(seq1))):
            data_out.append([i+1, seq1[i], 0, z_factors[i], z_scores[i], s_data1[i], s_data2[i], diff[i], data1[i], data2[i]])
    
    # write the file
    o = open(outfile, 'w')
    # write header
    o.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('Nuc', 'Seq', 'DeltaSHAPE', 'Z-factor', 'Std_Score', 'Smoothed_Data1', 'Smoothed_Data2', 'Unsmoothed_Diff', 'Data1', 'Data2'))
    # sort by decreasing absolute smoothed difference if 
    if args.magrank:
        data_out.sort(reverse=True, key=lambda x: abs(x[2]))
    else:
        data_out.sort(key=lambda x: x[0])
    
    for i in data_out:
        # check for NaN values in columns 3 onward and replace them with -999.
        # I'm sure there's a simpler way to do this but I gave up.
        for j in range(2,len(i)):
            if np.isnan(i[j]):
                i[j]=-999
        
        # write the output
        o.write(('\t').join(map(str, i))+"\n")
    o.close()
