# simDominanceModifier

Our simulation model builds off of Porter et al. (2017) and is based on the premise that transcription factors (TF) and binding sites behave according to the thermodynamic and kinetic properties of molecular interactions. 

Our diploid model consists of two unlinked interacting genes, A and B. The alternative alleles of Gene A codes for TFs that have opposite effects on the expression of Gene B, and Gene B’s expression is under additive sexually antagonistic selection. Gene A’s experession is determined by cis-regulatory binding site that can mutate/evolve toward being a better binder of either the male- or female-limited regulatory stimuli, approximating sex-limited hormones. Gene B's expression is determined by cis-regulatory binding site that can mutate/evolve toward being a better binder to either of the TFs expressed by Gene A. Finally, fitness is determined by the expression of Gene B, where females and males benefit from increased and decreased expression, respectively. There is no dominance coefficient; sex-specific dominance reversal for fitness between Gene A's alleles evolves at the genotype-expression level as an emergent property of their fractional occupancy of Gene B’s binding site. Populations are simulated following the standard Wright-Fisher model with discrete, non-overlapping generations undergoing mutation, reproduction and selection. 

Please refer to our publication (Grieshop et al. 2023) and supplementary documents for a detailed description of the model and its parameters.


# Dependencies
* Python 3.10.9 (https://www.python.org/)
* Pandas 1.5.2 (https://pandas.pydata.org/)
* Numpy 1.23.5 (https://numpy.org/)


# Usage

```
usage: simDominanceModifier.py [-h] -n NPROC -p PARFILE

options:
  -h, --help            show this help message and exit
  -n NPROC, --NPROC NPROC
                        Number of processes to run simultaneously
  -p PARFILE, --PARFILE PARFILE
                        CSV file containing parameters for simulations
```

## Input parameter file

The script requires a CSV file containing parameters for each run (`-p PARFILE`). Each row of the parameter file represents an independent simulation. Each row must contain all the parameters listed below. Please refer to `example_parameter_file.csv` for guidance.

| Parameter | Description                                                                                                                                                                                |
|-----------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| DIR       | Output directory.                                                                                                                                                                          |
| NAME      | Output file name.                                                                                                                                                                          |
| REP       | Number of replicate simulations to run.                                                                                                                                                    |
| CONC_D    | Concentration of sex-limited regulatory stimulus, D.                                                                                                                                       |
| GA        | Maximum concentration achievable by the expression of gene A.                                                                                                                              |
| KA        | Stepwise change in the dissociation constant for binding between sex-limited regulatory stimulus and the binding site of gene A..                                                          |
| KB        | Stepwise change in the dissociation constant for binding between the expression product of gene A to the binding site of gene B.                                                           |
| ASAT      | Concentration for the expression of gene A that would saturate the binding site of gene B.                                                                                                 |
| SR        | Selection regime. 0 for sexual antagonism, 1 for sexually concordant selection using the "female" fitness function, 2 for sexually concordant selection using the "male" fitness function. |
| S_F       | Selection coefficient for females.                                                                                                                                                         |
| S_M       | Selection coefficient for males.                                                                                                                                                           |
| GAMMA_F   | Fitness curvature for females.                                                                                                                                                             |
| GAMMA_M   | Fitness curvature for males.                                                                                                                                                               |
| UA        | Mutation rate for binding site of gene A.                                                                                                                                                  |
| UB        | Mutation rate for binding site of gene B.                                                                                                                                                  |
| N         | Population size to simulate.                                                                                                                                                               |
| PROPF     | Proportion of females in the population.                                                                                                                                                   |
| NUMGEN    | Number of generations to run simulation.                                                                                                                                                   |

## Output files
Three CSV files are output for each run of the similation in the directory specificied by `DIR` of the parameter file. The output of each run with be prefixed with `{NAME}.{replicate number}` as specified by `NAME` and `REP` in the parameter file.

### **Summary results of simulation run**: `{NAME}.{replicate number}.general.csv`

The file outputs a variety of summary statistics for the population every 10<sup>th</sup> generation of the simulation. Refer to the table below for descriptions of each column in the file.

| Column        | Description                                                                                               |
|---------------|-----------------------------------------------------------------------------------------------------------|
| GEN           | Generation of simulation.                                                                                 |
| maxFit_F      | Maximum value of absolute fitness in females.                                                             |
| maxFit_M      | Maximum value of absolute fitness in males.                                                               |
| maxFit        | Maximum value of absolute fitness in population.                                                          |
| avgFit_F      | Average value of absolute fitness in females.                                                             |
| avgFit_M      | Average value of absolute fitness in males.                                                               |
| avgFit        | Average value of absolute fitness in population.                                                          |
| avgFit_F00    | Average value of absolute fitness in females homozygous for allelic value 1 at coding site of gene A.     |
| avgFit_F01    | Average value of absolute fitness in females heterozygous at coding site of gene A.                       |
| avgFit_F11    | Average value of absolute fitness in females homozygous for allelic value 2 at coding site of gene A.     |
| avgFit_M00    | Average value of absolute fitness in males homozygous for allelic value 1 at coding site of gene A.       |
| avgFit_M01    | Average value of absolute fitness in males heterozygous at coding site of gene A                          |
| avgFit_M11    | Average value of absolute fitness in males homozygous for allelic value 2 at coding site of gene A.       |
| avgFit_00     | Average value of absolute fitness in individuals homozygous for allelic value 2 at coding site of gene A. |
| avgFit_01     | Average value of absolute fitness in individuals heterozygous at coding site of gene A.                   |
| avgFit_11     | Average value of absolute fitness in individuals homozygous for allelic value 2 at coding site of gene A. |
| avgExpA_F     | Average expression of gene A in females.                                                                  |
| avgExpA_M     | Average expression of gene A in males.                                                                    |
| avgExpA       | Average expression of gene A in population.                                                               |
| varExpA_F     | Variance in expression of gene A in females.                                                              |
| varExpA_M     | Variance in expression of gene A in males.                                                                |
| varExpA       | Variance in expression of gene A in population.                                                           |
| avgStExpB_F   | Average standardized expression of gene B in females.                                                     |
| avgStExpB_M   | Average standardized expression of gene B in males.                                                       |
| avgStExpB     | Average standardized expression of gene B in population.                                                  |
| varStExpB_F   | Variance in standardized expression of gene B in females.                                                 |
| varStExpB_M   | Variance in standardized expression of gene B in males.                                                   |
| varStExpB     | Variance in standardized expression of gene B in population.                                              |
| avgHetExpA0_F | For females heterozygous at the coding site of gene A - average expression for allele with value 1.       |
| avgHetExpA1_F | For females heterozygous at the coding site of gene A - average expression for allele with value 2.       |
| avgHetExpA0_M | For males heterozygous at the coding site of gene A - average expression for allele with value 1.         |
| avgHetExpA1_M | For males heterozygous at the coding site of gene A - average expression for allele with value 2.         |
| freqA_F       | Frequency of allele with value 1 at coding site of gene A across females.                                 |
| freqA_M       | Frequency of allele with value 1 at coding site of gene A across males.                                   |
| freqA         | Frequency of allele with value 1 at coding site of gene A across population.                              |
| gtA00_F       | Frequency of homozygotes with allelic value 1 at coding site of gene A across females.                    |
| gtA01_F       | Frequency of heterozygotes at coding site of gene A across females.                                       |
| gtA11_F       | Frequency of homozygotes with allelic value 2 at coding site of gene A across females.                    |
| gtA00_M       | Frequency of homozygotes with allelic value 1 at coding site of gene A across males.                      |
| gtA01_M       | Frequency of heterozygotes at coding site of gene A across males.                                         |
| gtA11_M       | Frequency of homozygotes with allelic value 2 at coding site of gene A across males.                      |
| gtA00         | Frequency of homozygotes with allelic value 1 at coding site of gene A across individuals.                |
| gtA01         | Frequency of heterozygotes at coding site of gene A across individuals.                                   |
| gtA11         | Frequency of homozygotes with allelic value 2 at coding site of gene A across individuals.                |
| avgAlpha_F    | Average allelic value of binding site of gene A in females.                                               |
| avgAlpha_M    | Average allelic value of binding site of gene A in males.                                                 |
| avgAlpha      | Average allelic value of binding site of gene A in population.                                            |
| avgBeta_F     | Average allelic value of binding site of gene B in females.                                               |
| avgBeta_M     | Average allelic value of binding site of gene B in males.                                                 |
| avgBeta       | Average allelic value of binding site of gene B in population.                                            |
| cAlphaA_F     | Correlation between binding site value of gene A and coding site value of gene A in females.              |
| cAlphaBeta_F  | Correlation between binding site value of gene A and binding site value of gene B in females.             |
| cABeta_F      | Correlation between coding site value of gene A and binding site value of gene B in females.              |
| cAlphaA_M     | Correlation between binding site value of gene A and coding site value of gene A in males.                |
| cAlphaBeta_M  | Correlation between binding site value of gene A and binding site value of gene B in males.               |
| cABeta_M      | Correlation between coding site value of gene A and binding site value of gene B in males.                |

### **Distribution of values for the binding and protein coding sites**: `{NAME}.{replicate number}.dist.csv`

This file outputs the occurence for each allelic value of the binding site at gene A and the binding site of gene B for every 10<sup>th</sup> generation of the simulation. It also outputs the number of females and males binned by their standardized expression levels of Gene B. Refer to the table below for descriptions of each column in the file.

| Column | Description                                                                               |
|--------|-------------------------------------------------------------------------------------------|
| GEN    | Generation of simulation.                                                                 |
| alF{x} | Binding site of gene A: number of alleles in females with value {x}.                      |
| alM{x} | Binding site of gene A: number of alleles in males with value {x}.                        |
| beF{x} | Binding site of gene B: number of alleles in females with value {x}.                      |
| beM{x} | Binding site of gene B: number of alleles in males with value {x}.                        |
| al{x}  | Binding site of gene A: number of alleles in population with value {x}.                   |
| be{x}  | Binding site of gene B: number of alleles in population with value {x}.                   |
| BF{x}  | Standardized expression of gene B: number of females with values between {x} and {x+0.1}. |
| BM{x}  | Standardized expression of gene B: number of males with values between {x} and {x+0.1}.   |

### **Genotype and phenotype of each individual**: `{NAME}.{replicate number}.indiv.csv`

This file outputs properties for each individual at the final generation of the simumlation. Each row of this file represents a separate individual. Refer to the table below for descriptions of each column in the file.

| Column     | Description                                                                                                                                       |
|------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| ID         | Unique identifier of individual.                                                                                                                  |
| SEX        | Sex of individual: `f` for female, `m` for male.                                                                                                  |
| GTA        | Genotye of protein coding region for gene A. 0 and 2 represents homozygosity for allelic value 1 and 2, respectively. 1 represents heterozygotes. |
| A1         | Allelic value at homolog 1 for protein coding region of gene A.                                                                                   |
| A2         | Allelic value at homolog 2 for protein coding region of gene A.                                                                                   |
| ALPHA1     | Allelic value at homolog 1 for binding site of gene A.                                                                                            |
| ALPHA2     | Allelic value at homolog 2 for binding site of gene A.                                                                                            |
| BETA1      | Allelic value at homolog 1 for binding site of gene B.                                                                                            |
| BETA2      | Allelic value at homolog 2 for binding site of gene B.                                                                                            |
| concA1     | Concentration of protein expressed by homolog 1 of gene A.                                                                                        |
| concA2     | Concentration of protein expressed by homolog 2 of gene A.                                                                                        |
| stconB     | Standardized concentration of protein expressed by gene B.                                                                                        |
| absFitness | Absolute fitness of individual.                                                                                                                   |

# References
Porter, AH., Johnson, NA. & Tulchinsky, AY. A New Mechanism for Mendelian Dominance in Regulatory Genetic Pathways: Competitive Binding by Transcription Factors. Genetics 205, 101–112 (2017).

Grieshop, K., Ho, EKH. & Kasimatis, KR. Dominance reversals and the maintenance of genetic variation. *In review.* (2023).
