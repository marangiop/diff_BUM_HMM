                             deltaSHAPE

  What is it?
  -----------

  deltaSHAPE [1] presents a statistically rigorous approach to analyzing
  changes in chemical probe reactivity between two conditions. It
  was originally designed for finding protein binding sites using
  SHAPE-MaP data [2], which include estimated measurement errors for
  every nucleotide. deltaSHAPE compares reactivity differences relative
  to the associated errors, and also compares the magnitude of each
  difference relative to every other change in reactivity. Nucleotides
  that exhibit strong, meaningful changes in reactivity are thus
  identified.
  
  Software Requirements
  ---------------------
  
  In order to run deltaSHAPE.py, you must have Python 2.7 installed, as
  well as two additional modules:
  
  matplotlib, freely available from http://www.matplotlib.org
  NumPy, freely available from http://www.numpy.org
  SciPy, freely available from http://www.scipy.org
  
  Installation instructions for each of these modules are provided on
  their respective websites.
  
  Documentation
  -------------
  
  Simple usage: python deltaSHAPE.py <data1.map> <data2.map>
  
  Required input
  
    data1.map, data2.map
    (The .map file format is described below.)
  
  Options
  
  Data handling:
  
    --mask5 INT
      
      Tell deltaSHAPE to ignore INT number of nucleotides at the 5' end of
      the RNA. Often used in combination with --mask3 to exclude primer-
      binding sites from datasets generated with targeted primers [2,3].
      Default: 0
    
    --mask3 INT
    
      Tell deltaSHAPE to ignore INT number of nucleotides at the 3' end of
      the RNA. Often used in combination with --mask5 to exclude primer-
      binding sites from datasets generated with targeted primers [2,3].
      Default: 0

    -p INT, --pad INT

      Indicate the smoothing window size. By convention, SHAPE reactivity
      values are smoothed over a centered 3-nt window prior to calculating
      differences. Window size is calculated as (2 * pad + 1). The default
      smoothing window size is 3 (--pad 1). To eliminate smoothing, set
      --pad to 0.
      Default: 1

    -z FLOAT --Zcoeff FLOAT

      Adjust the Z-factor test stringency. Default setting of 1.96 tests
      whether the 95% confidence intervals of two SHAPE-MaP measurements
      overlap. Can be interpreted as the area under the normal distribution
      curve within Zcoeff standard deviations of the mean. To increase test
      stringency, raise this value. Alternatively, adjust the value of
      -t/--Zthresh.
      Default: 1.96

    -t --Zthresh FLOAT

      Adjust the Z-factor test stringency. Default setting of 0 tests
      whether the ratio of Zcoeff*(err1 + err2):(reactivity1 - reactivity2)
      is greater or equal to 1. Raising Zthresh decreases the maximum
      allowed ratio. Alternatively, adjust the value of -z/--Zcoeff.
      Default: 0

    -s --SSthresh FLOAT

      Set the cutoff threshold used during standard score filtering.
      Default value of 1 indicates that only nucleotides which undergo a
      change in SHAPE reactivity with a magnitude of 1 or more standard
      deviations from the mean will be considered as passing the standard
      score filter. Increase SSthresh to increase stringency, or lower
      SSthresh to decrease stringency.
      Default: 1

    -f --FindSite INT,INT

      Indicate the window size and number of required hits to use when
      determining binding sites. These are provided as a comma-separated
      pair of integers (e.g., 2,3). The first number indicates the window
      pad size where window size = 2 * pad + 1. This "pad" value is
      independent of the pad value used to set the smoothing window size.
      The second number sets the number of nucleotides within a given
      window that must pass both the Z-factor and Standard Score tests.
      The default of 2,3 means that 3 nucleotides within a 5 nucleotide
      window must pass both tests, otherwise they will not be considered
      part of a binding site and will not be reported. To eliminate this
      windowed searching, set -f --FindSite equal to 0,1.
      Default: 2,3

  Output options:

    -o --out STRING

      Provide a filename and location for output file(s) to be written.
      Output file is ordered by nucleotide position, 5' to 3', unless
      otherwise indicated by use of the --magrank option.
      Default: ./differences.txt

    --magrank

      Sort output text file by absolute magnitude of the deltaSHAPE
      differences.
      Default: OFF

    --all

      Output data for all nucleotides, regardless of whether they are
      significantly different between experimental conditions. When this
      option is used, insignificant deltaSHAPE values are reported as 0.
      Default: OFF

    --pdf

      Save a PDF file of the plotted reactivity differences. The PDF file
      will have the same prefix as the output file, by default,
      "differences.pdf". Use of the -o --out flag will change the prefix
      of the PDF file as well. This option is commonly used with the
      --noshow option. Use with the --noplot flag has no effect.
      Default: OFF

    --noshow

      Prevent the plot from showing. Useful when processing many datasets
      or when running on a remote machine without window forwarding.
      Commonly used with the --pdf option. Use with the --noplot flag has
      no effect.
      Default: OFF (the plot is displayed)

    --noplot

      Skip all plotting steps, and only write the output text file. Use of
      this option nullifies any options set with the --pdf, --noshow,
      --dots, --Zdots, --SSdots, --colorbar, --colorfill, --ymin, --ymax,
      --xmin, and --xmax options.
      Default: OFF (plotting is performed)

    --dots

      Plot markers indicating which nucleotides passed both the Z-factor
      and standard score filters. Standard score dots are open and plotted
      the Z-score dots, which are filled. For large RNAs, this can get
      unweildy, especially when using an interactive plotting window.
      However, it is useful for troubleshooting when very many or very few
      nucleotides are flagged as significantly different.
      Default: OFF

    --Zdots

      Exactly like --dots, except only the dots indicating nucleotides
      that pass the Z-factor filter are plotted.
      Default: OFF

    --SSdots

      Exactly like --dots, except only the dots indicating nucleotides
      that pass the Standard Score filter are plotted.
      Default: OFF

    --colorfill

      Significant deltaSHAPE differences are highlighted on the plot with
      colors that fill between the difference line and zero. This produces
      somewhat "prettier" images but can be harder to see for large RNAs
      when scaled down. Default behavior is to highlight the plot with
      colors that span the entire height of the plot.
      Default: OFF

    --ymin FLOAT

      Set the minimum y-axis value.
      Default: Determined automatically

    --ymax FLOAT

      Set the maximum y-axis value
      Default: Determined automatically.

    --xmin FLOAT

      Set the minimum x-axis value.
      Default: Determined automatically

    --xmax FLOAT

      Set the maximum x-axis value
      Default: Determined automatically.

  
  File formats
  ------------
  
  .map files
  
    The .map file is a 4-column tab-delimited format. The columns list:
    
      (1) Nucleotide number (starting at 1)
      (2) SHAPE reactivity
      (3) Estimated standard error of the reactivity measurement
      (4) Nucleotide identity (A, G, C, U)
    
    Positions for which no SHAPE data are available should be included,
    with SHAPE reactivity and standard error values set to -999.

  output text file
  
    The output file created by deltaSHAPE has eight tab-deliminted
    columns:
    
      (1) Nucleotide number (starting at 1)
      (2) Nucleotide identity (A, G, C, U)
      (3) deltaSHAPE value (smoothed difference between conditions 1 & 2)
      (4) Z-factor
      (5) Standard score
      (6) Smoothed reactivity, condition 1
      (7) Smoothed reactivity, condition 2
      (8) Unsmoothed difference between conditions 1 & 2
      (9) Unsmoothed reactivity, condition 1
     (10) Unsmoothed reactivity, condition 2
      
    deltaSHAPE values are the smoothed difference between conditions 1 & 2.
    They are only reported for nucleotides that are identified as having
    strong changes between the two conditions. Values of 0 in the
    deltaSHAPE column are reported when the --all flag is used and all
    nucleotides, regardless of significance, are reported.
  
  Author
  ------
  
  Matt Smola ( matt.smola@gmail.com )
  Weeks Lab, Department of Chemistry
  University of North Carolina at Chapel Hill
  
  
  Copyright
  ---------
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
  References
  ----------
  1. M.J. Smola, J.M. Calabrese, and K.M. Weeks. Detection of RNA-protein
     interactions in living cells with SHAPE. doi:10.1021/acs.biochem.5b00977
  
  2. N.A. Siegfried, et al. RNA motif discovery by SHAPE and mutational
     profiling (SHAPE-MaP). Nat. Methods 11, 959-965 (2014).

  3. M.J. Smola, et al. Selective 2′-hydroxyl acylation analyzed by
     primer extension and mutational profiling (SHAPE-MaP) for direct,
     versatile, and accurate RNA structure analysis. Nat. Protoc. 10,
     1643-1669 (2015).