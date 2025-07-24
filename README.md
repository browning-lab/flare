# flare

The **flare** program uses a set of reference haplotypes
to infer the ancestry of each allele in a set of admixed study samples.
The **flare** program is fast, accurate, and memory-efficient.

Last updated: July 24, 2025  
Current version: 0.5.3

## Contents

* [Installation](#installation)
* [Running flare](#running-flare)
  * [Required parameters](#required-parameters)
  * [Optional parameters](#optional-parameters)
* [Ancestry indices](#ancestry-indices)
* [Output files](#output-files)
  * [Global ancestries file format](#global-ancestries-file-format)
  * [Model file format](#model-file-format)
* [Running flare with multiple chromosomes](#running-flare-with-multiple-chromosomes)
* [Running flare with small or large sample sizes](#running-flare-with-small-or-large-sample-sizes)
* [Change log](#change-log)
* [License](#license)
* [Acknowledgments](#acknowledgements)
* [Citation](#citation)

## Installation

You can download the latest executable file,
[flare.jar](https://faculty.washington.edu/browning/flare.jar),
with the command:

    wget https://faculty.washington.edu/browning/flare.jar

or you can download the source files and create the executable file
with the commands:

    git clone https://github.com/browning-lab/flare.git
    javac -cp flare/src/ flare/src/admix/AdmixMain.java
    jar cfe flare.jar admix/AdmixMain -C flare/src/ ./

[Contents](#contents)

## Running flare

The **flare** program requires Java version 1.11 (or a later version). Use of an
earlier Java version will produce an "Unsupported Class Version" error.

The command:

    java -jar flare.jar

prints a summary of the command line arguments.

To run **flare**, enter the following command:

    java -Xmx[GB]g -jar flare.jar [arguments]

where **[GB]** is the maximum number of gigabytes of memory to use, and
**[arguments]** is a space-separated list of parameter values, each expressed as
**parameter=value**.

The shell script
[run.flare.test](https://raw.githubusercontent.com/browning-lab/flare/master/test/run.flare.test)
will run a test **flare** analysis.

[Contents](#contents)

### Required parameters

The **flare** program has five required parameters. Two of the required
parameters specify
[Variant Call Format](https://faculty.washington.edu/browning/intro-to-vcf.html)
(VCF) files.  A VCF record may have multiple ALT alleles and must
include a genotype (GT) FORMAT subfield. **All genotypes must be phased and have
no missing alleles**.
If a VCF file has unphased or missing genotypes, you can phase the genotypes and
fill in the missing genotypes using the
[Beagle](https://faculty.washington.edu/browning/beagle/beagle.html) program.
Any input file with a name ending in ".gz" is assumed to be gzip-compressed.
Any input VCF file with a name ending in ".bref3" is assumed to be
bref3-compressed.  Software for bref3 compression and decompression can be
downloaded from the
[Beagle web site](https://faculty.washington.edu/browning/beagle/beagle.html).

* **ref=[file]** where **[file]** is the reference VCF file that contains
genotype data for each reference sample. Flare will ignore samples in the
reference VCF file that are not present in the reference panel file
(see the **ref-panel** parameter).

* **ref-panel=[file]** where **[file]** is a reference panel file with two
whitespace-delimited fields per line. The first field is a sample identifier
in the reference VCF file (see the **ref** parameter), and the second field
is the name of the reference panel containing the reference sample.
Flare will ignore samples in the reference VCF file that are not present
in the reference panel file. A reference panel should contain individuals
from the same source population. A reference panel should not normally
contain admixed samples.

* **gt=[file]** where **[file]** is the study VCF file containing genotype
data for admixed study samples whose ancestry is to be inferred.
The **gt-samples** parameter can be used to restrict the analysis to
a subset of samples in the study VCF file. All admixed study samples in an
analysis should be from the same population.

* **map=[file]** where **[file]** is a
[PLINK format genetic map](https://zzz.bwh.harvard.edu/plink/data.shtml#map)
with cM units. Positions of markers that are between genetic map positions are
estimated using linear interpolation. The chromosome identifiers
in the genetic map and the input VCF files must match. HapMap genetic maps
in cM units are available for
[GRCh36](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/),
[GRCh37](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/), and
[GRCh38](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/).

* **out=[string]** where **[string]** is the output filename prefix.

[Contents](#contents)

### Optional parameters

* **array=[true/false]** specifies whether the input data are from a SNP array.
The **min-mac** parameter is ignored if **array=true**.  By default,
flare assumes the input data are sequence data (**default: array=false**).

* **min-maf=[number < 0.5]** specifies the minimum minor allele frequency in
the reference VCF file in order for a marker to be included in the
analysis (**default: min-maf=0.005**). For multi-allelic markers,
the minor allele is the allele with the second-largest frequency.

* **min-mac=[number ≥ 0]** specifies the minimum minor allele count in
the reference VCF file in order for a marker to be included in the
analysis (**default: min-mac=50**).  The **min-mac** parameter is ignored
if **array=true**. If **array=false**, the **min-mac** parameter must be less
than one-half the number of reference haplotypes.  For multi-allelic markers,
the minor allele is the allele with the second-largest frequency.

* **probs=[true/false]** specifies whether posterior ancestry probabilities
are reported (**default: probs=false**). At each marker, the ancestry with
the highest posterior probability for each haplotype is _always_ reported
in the output VCF file.
If **probs=true**, posterior probabilities for each ancestry, haplotype, and
marker will also be reported in the output VCF file.
Reporting posterior probabilities will modestly increase memory use and
computation time and significantly increase the size of the output VCF file.

* **gen=[integer ≥ 1]** specifies the number of generations since
admixture (**default: gen=10**). If **em=true**, the specified **gen** parameter
is an initial value for the **gen** parameter that will be used in the
parameter estimation algorithm. The **gen** parameter is ignored if the
**model** parameter is used.

* **model=[file]** where **[file]** is a whitespace-delimited file containing
model parameters (see [**Model file format**](#model-file-format)). If the
**model** parameter is not used, **flare** will supply a reasonable set of
model parameters (see the [**flare** paper](#citation) for details).
If **em=true** (the default), **flare** will estimate the number of
generations since admixture and the proportion of genotypes with each
ancestry and will replace the values for these two parameters
in the **model** file with their estimated values.
The model parameters used in the analysis are reported in the output
[**model**](#output-files) file.

* **em=[true/false]** specifies whether the number of generations since
admixture and the proportion of genotypes with each ancestry will
be estimated using an iterative expectation maximization (EM) algorithm
prior to inferring local ancestry (**default: em=true**).

* **gt-samples=[file]** (or **gt-samples=^[file]**) where **[file]**
is a text file containing one sample identifier per line.
Only admixed study samples that are present in **[file]** (or absent from
**[file]** if **[file]** is preceeded by **^**) will be analyzed. If the
**gt-samples** parameter is omitted, all admixed study samples will be included
in the analysis. The **gt-samples** parameter filters the study samples,
and the **ref-panel** parameter filters the reference samples.

* **gt-ancestries=[file]** where **[file]** is a whitespace-delimited text file
containing ancestry proportions for a set of samples.  The **gt-ancestries**
file has the same format as the output **global ancestries** file
(see [**Global ancestries file format**](#global-ancestries-file-format)).
If a study sample is present in the first column of the **gt-ancestries** file,
the ancestry proportions in the **gt-ancestries** file will be used
when estimating local ancestry in that sample; otherwise, study-wide
ancestry proportions will be used.

* **excludemarkers=[file]** where [file] is a text file containing markers
(one marker identifier per line) that are to be excluded from the analysis.
A marker identifier can be an identifier from the VCF record ID field, or it
can be a VCF record's CHROM and POS fields separated by a colon
(i.e. "CHROM:POS").

* **nthreads=[integer ≥ 1]** specifies the number of computational threads to
use for the analysis. The default **nthreads** parameter is the number of
CPU cores.  The **nthreads** parameter value is printed in the output **log**
file.

* **seed=[integer]** specifies the seed for random number generation
(**default: seed=-99999**). Repeating an analysis with the same **seed** and
**nthreads** parameters will produce the same local ancestry estimates.

[Contents](#contents)

## Ancestry indices

The **flare** program obtains the list of ancestries from the input files.
If the **model** parameter is specified, **flare** will use the list of
of ancestries at the start of the **model** file. If the **model**
parameter is not specified, **flare** will use the reference panel
names in the **ref-panel** file, in the order each name first appears in the
**ref-panel** file, as the ancestry names.

If the **gt-ancestries** parameter is specified, the list of ancestries in the
header line of the **gt-ancestries** file must have the same ancestries in
the same order as the list of ancestries obtained from the **model** or
**ref-panel** file.

The index of an ancestry in the list of ancestries determines a unique,
non-negative integer associated with each ancestry. The first ancestry in
the list has index 0, the second ancestry in the list has index 1, and so on.

[Contents](#contents)

## Output files

The **flare** program produces four output files: a **log** file, a
**model** file, a **VCF** file, and a **global ancestry** file.

The output **log** file (.log) contains a summary of the analysis.

The output **model** file (.model) contains the model parameters used in the
analysis. The output model file has the same format as the optional input
model file (see [Model file format](#model-file-format)).

The output **VCF** file (.anc.vcf.gz) contains the phased input genotypes and
the estimated local ancestry for each allele. The most probable
ancestry at each marker for a admixed sample's first and second haplotype
are reported in the **AN1** and **AN2** FORMAT subfields.
If [**probs=true**](#optional-parameters), the posterior
ancestry probabilities at each marker for the admixed sample's first and
second haplotypes are reported in the **ANP1** and **ANP2** FORMAT subfields.
The integer that identifies each ancestry is given in the
"##ANCESTRY=<...>" meta-information line.

The output **global ancestries** file (.global.anc.gz) contains the
estimated global ancestry proportions for each admixed sample.
The global ancestry file has the same format as the optional input
**gt-ancestries** file
(see [Global ancestry file format](#global-ancestry-file-format)).

[Contents](#contents)

### Global ancestries file format

A **global ancestries** file is a whitespace-delimited text file containing
the global ancestry proportions for a set of samples. If there are
$A$ ancestries, each line of a **global ancestry** file has $A + 1$ fields.
The first line of a **global ancestry** file is a header line.  The
first field of the header line is `SAMPLE`, and the remaining fields
are the $A$ ancestry identifiers. The first field in each subsequent line
is a sample identifier, and the remaining $A$ fields are the
proportions of the genome having the ancestry given in the corresponding
fields of the header line.  A sample identifier is not permitted to appear
more than once in the first column.

[Contents](#contents)

### Model file format

A **model** file contains model parameters.
The model file can contain comment lines, blank lines, and data lines.
A comment line is a line whose first non-whitespace character is the
'#' character. A blank line contains only whitespace characters.
All other lines are data lines. Data lines contain whitespace-delimited
fields that specify the model parameters.

If there are $A$ ancestries and $P$ reference panels, the model file will
contains $(2A + 5)$ data lines.

* The first data line is the list of $A$ ancestry names.  The first ancestry
in the list has index 1.

* The second data line is the list of the $P$ reference panel names. The
first reference panel in the list has index 1.

* The third data line is the number of generations since admixture.

* The fourth data line is a vector of length $A$ whose $i$-th element
is the proportion of admixed sample genotypes with ancestry $i$.

* The next $A$ data lines contain the first $A \times P$ matrix. The $(i,j)$-th
element of the matrix is the probability that a model state haplotype is in
reference panel $j$ when the model state ancestry is $i$.

* The next $A$ data lines contain the second $A \times P$ matrix.  The $(i,j)$-th
element of the matrix is the probability that a model state haplotype and the
admixed sample haplotype carry different alleles when the model state haplotype
is in reference panel $j$ and the model state ancestry is $i$.

* The final data line is a vector of length $A$ whose
$i$-th element is the the rate of the exponential identity-by-descent segment
cM-length distribution when the most recent common ancestor is pre-admixture
and has ancestry $i$.

It is not normally necessary to use a model file
because **flare** will automatically estimate model parameters by default
(see the [**em**](#optional-parameters) parameter). If you want
to specify the model parameters, the easiest way to ensure that the **model**
file is in the correct format is to run **flare** without the
[**model**](#optional-parameters) parameter, and then modify the values in the
output model file.

[Contents](#contents)

## Running flare with multiple chromosomes

If an input VCF file contains multiple chromosomes, **flare** will estimate model
parameters using data from the first chromosome and use these model parameters
for all subsequent chromosomes in the VCF file.

If you analyze each chromosome separately, you can use the same
model parameters for all chromosomes by analyzing one chromosome, and then
analyze each remaining chromosome with
[**em=false**](#optional-parameters) and the
[**model**](#optional-parameters) parameter set equal to the output
[**model**](#output-files) file from the first chromosome's analysis.

[Contents](#contents)

## Running flare with small or large sample sizes

If the sample size is too small for accurate estimation of model parameters,
you can use the [**model**](#optional-parameters) and
[**em=false**](#optional-parameters) parameters to specify the model
parameters used in the analysis.

If you are analyzing an extremely large number of admixed samples
and need to reduce memory use, you can partition the admixed samples
into subsets and analyze each subset of admixed samples separately (see the
[**gt-samples**](#optional-parameters) parameter).  The inferred
ancestry for a partitioned and a non-partitioned analysis will be the same if
you specify [**em=false**](#optional-parameters) and use the same
[**model**](#optional-parameters), [**seed**](#optional-parameters), and
[**nthreads**](#optional-parameters) parameters for each analysis.

[Contents](#contents)

## Change log

* Add a header line to the output
[global ancestries](#global-ancestries-file-format) file (version 0.5)
* Add an optional **gt-ancestries** command line parameter
(see [Optional parameters](#optional-parameters)) (version 0.5)
* Added [global ancestries](#global-ancestries-file-format) file with 
global ancestry probabilities (version 0.4)
* Add support for bref3-compressed VCF files (version 0.4)

[Contents](#contents)

## License

The **flare** program is licensed under the Apache License, Version 2.0 (the License).
You may obtain a copy of the License from https://www.apache.org/licenses/LICENSE-2.0

[Contents](#contents)

## Acknowledgements

Thank you to Kai Zeng and Thomas Hickman for suggesting several of the
new features in **flare** versions 0.4 and 0.5, and to Thomas Hickman
for creating and sharing implementations of two features: support 
for bref3-compressed VCF files and for user-specified, per-sample 
admixture proportions.

[Contents](#contents)

## Citation

If you use **flare** in a published analysis, please report the program
version printed in the **log** file and cite the article describing
the **flare** method:

> S R Browning, R K Waples, B L Browning (2023). Fast, accurate local ancestry
estimation with FLARE. The American Journal of Human Genetics 110(2):326-335.
doi: http://dx.doi.org/10.1016/j.ajhg.2022.12.010

[Sharon Browning](https://sites.uw.edu/sguy/) developed the **flare** method.  
[Brian Browning](https://faculty.washington.edu/browning) developed the **flare** software.

[Contents](#contents)

