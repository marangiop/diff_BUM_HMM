# diffBUM_HMM
Bayesian modelling approach for detecting RNA flexibility changes in high-throughput structure probing data

## Background 
DiffBUM-HMM (differential BUM-HMM) is a natural extension of the beta-uniform mixture hidden Markov model developed previously [Selega2017](https://pubmed.ncbi.nlm.nih.gov/27819660/). 

ADD MORE INFO about LDR and LMR

For each experimental condition (e.g. Condition 1 and 2), the null distribution of log (with base e) drop-off rates or log mutation rates (SHAPE-MaP) at each nucleotide position is computed for control samples in order to quantify variability in drop-off or mutation rates observed by chance. Subsequently, coverage-dependent biases are removed by applying a variance stabilization transformation. The different null distributions are computed for different nucleotide patterns to address sequence-dependent bias. Subsequently, per-nucleotide empirical P values are computed for all pairs of treatment and control replicate samples in each condition, by comparing the corresponding log of drop-off rates (LDRs) or log of mutation rates (LMRs, for SHAPE-MaP) to the null distribution. DiffBUM-HMM is run on P values associated with the two independent conditions as observations, leaving out any nucleotides with missing data. The resulting output is a posterior probability of modification for each nucleotide, ranging from 0 to 1. DiffBUM-HMM reports whether nucleotides were unmodified in both conditions, modified in either of the conditions or modified in both conditions.


![img/Figure_1.jpg](img/Figure_1.jpg)

## Reproducing figures from the paper
| Figure | Instructions for raw data analysis | Jupyter Notebook for figure generation |
|   ------------- |-------------        | -------------|
| 2-3  | [Instructions](./Jupyter_notebooks/Figure_2_3/instructions_data_analysis_fig2_3.txt)  | [Notebook](./Jupyter_notebooks/Figure_2_3/Plotting_5'ETS_and_35S_data.ipynb)  |
| 4   |  [Instructions](./Jupyter_notebooks/Figure_4/instructions_data_analysis_fig4.txt)  | TBA     |
| 5A-B   | [Instructions](./Jupyter_notebooks/Figure_5/instructions_data_analysis_fig5.txt)    | [Notebook](./Jupyter_notebooks/Figure_5/Binning_and_smoothing_diffBUM_HMM_signal/notebook_binned_results.ipynb)   |
| 5C   | [Instructions](./Jupyter_notebooks/Figure_5/instructions_data_analysis_fig5.txt)     | [Notebook](./Jupyter_notebooks/Figure_5/Heatmap_diffBUM-HMM_&_deltaSHAPE_with_protein_binding_sites/heatmap.ipynb)    |
| 5D   | [Instructions](./Jupyter_notebooks/Figure_5/instructions_data_analysis_fig5.txt)     | [Notebook](./Jupyter_notebooks/Figure_5/Hypergeometric_test_Xist_bindingsites/notebook_hypergeometric_test.ipynb)  |
| 6   | [Instructions](./Jupyter_notebooks/Figure_6/instructions_data_analysis_fig6.txt)   | [Notebook](./Jupyter_notebooks/Figure_6/Nucleotide_analyses.ipynb)     |
| S1    | [Instructions](./Jupyter_notebooks/Supplementary_Figure_1/instructions_data_analysis_figS1.txt)   |[Notebook](./Jupyter_notebooks/Supplementary_Figure_1/Plotting_pertubation_tests.ipynb)    |
| S3    | [Instructions](./Jupyter_notebooks/Supplementary_Figure_3/instructions_SF3.txt)   | N/A    |
| S4    | [Instructions](./Jupyter_notebooks/Supplementary_Figure_4/instructions_SF4.txt)    | N/A    |
| S5    | [Instructions](./Jupyter_notebooks/Supplementary_Figure_5/instructions_SF5.txt)  | N/A   |

The table above only includes instructions for figures from the paper that have been generated programmatically. 

## Dependencies
The pipeline is built in R. Python and Jupyter (notebook) are needed for performing some of the raw data analysis and figure generation. 

- R 4.0.0 (2020-04-24) (version 3.6.3 2020-02-29 tested working)
- RStudio 1.2.5001 (version 1.1.442 tested working)
- Python 3.7.6 (unless specified otherwise)

## Requirements 
We provide a requirements.txt file listing R and Python packages with correct versions used in the development and benchmarking of the diffBUM-HMM pipeline describe in this study, as reference and also to enable quick installation of the packages. TBA





