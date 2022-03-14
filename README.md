# Project 2
final project 2 for bimm 143 final

new branch with necessary files for code

new branch with .nb.html and .Rmd files of project 2


## Introduction
Scientific Question: Why do mutations in the BRCA1 gene lead to a 20% or higher risk of developing ovarian cancer compared to mutations in the BRCA2 gene?

BRCA1 & BRCA2:

The breast cancer 1 gene, also known as BRCA1, is a tumor suppressor gene found on chromosome 17 in humans. The breast cancer 2 gene, also known as BRCA2, is another tumor suppressor gene located on chromosome 13 (Mehrgou 2016). Both of these genes produce proteins that are involved in different steps of DNA damage response (DDR) and DNA repair (Roy 2011). There are 1863 amino acids in the BRCA1 protein and 3418 amino acids in the BRCA2 protein (Mehrgou 2016). Mutations in either of the genes can cause genetic instability and lead to a higher risk of developing breast and ovarian cancer, especially those that occur in the germline because they can continue to be passed down throughout generations (Welcsh 2001). It was found that BRCA1 mutation carriers have a 50% lifetime risk of developing ovarian cancer while BRCA2 mutation carriers have a 30% risk (Roy 2011) and that BRCA1 mutations become more common than BRCA2 mutations as the number of ovarian cancer cases within a family increase (Ramus 2009). In the 2000s, PALB2 was established as the physical and functional connection between BRCA1 and BRCA2 (Zhang 2009). 


Scientific Hypothesis: If BRCA1 mutations are at least 20% more likely to lead to ovarian cancer than BRCA2 mutations, then there must be significant genetic and protein structural differences that account for increased risk. 


I will be performing a pairwise sequence alignment to compare the protein sequences of BRCA1 and BRCA2 and then creating a dotplot to visualize the results of the alignment. The results from this will help me to identify any significant genetic differences between the two genes. Sequences will be downloaded from the NCBI database. Significant genetic differences could be used to explain why mutations in BRCA1 are more likely to cause ovarian cancer than mutations in BRCA2. 

I will also be using structural bioinformatics to analyze the protein structures/complexes of BRCA1 and BRCA2. Additionally, I will attempt to determine how and where BRCA1 and BRCA2 interact by searching for binding site residues. The protein structures I analyze will be from the PDB and and I will import them into PyMOL to visualize the protein structural differences and interactions. Identifying the protein structural differences could also reveal why BRCA1 mutations have higher penetrance than BRCA2 mutations. Observing how close their interactions may be can also account for their similar functions and explain why their mutations have comparable effects. 


## Loading In Packages

1. BiocManager: This package allows us to use packages from the Bioconductor project. Packages from this project have functions that can be used for statistical analysis of genomic data in R.

2. Biostrings: This is one of the packages from Bioconductor. It allows for manipulation of biological strings, such as DNA, RNA, and amino acids. I will be using this package to perform a global pairwise alignment of protein sequences. 

3. seqinr: This package is used for retrieving and analyzing biological sequences. I will use seqinr to read in the fasta files of the protein sequences that I download from NCBI. I will also use the dotPlot function from this package to visualize the pairwise sequence alignment results.

4. Bio3D: The Bio3D package can be used for reading and writing biomolecular data types, analyzing sequences and structures, converting and manipulating data, plotting, and more. I will be using this package to read in coordinate files from PDB, look at internal strucure, and identify interacting residues. 

5. knitr: This package provides general-purpose tools for generating reports in R with Literate Programming techniques. It combines combines computing and reporting by incorporating code into texts documents. I will be using the include_graphics function to embed pictures of BRCA1 and BRCA2 complexes from PyMOL.  

6. vembedr: This package is meant for embedding videos into R. Videos from YouTube, Vimeo, Box, and Microsoft stream can be embedded with their respective functions in this package. I will be embedding my PyMOL movie which I uploaded to YouTube. 

```{r}
# if (!require("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")

# BiocManager::install("Biostrings")

library(BiocManager)
library(Biostrings)
library(seqinr)

```
## Performing Bioinformatics Analysis

Pairwise Sequence Alignment: This bioinformatics method is used to compare two sequences by finding the optimal alignment between them. It then measures sequence similarity with an alignment score. High sequence similarity means that the sequences may share a functional, structural, or evolutionary relationship. We will be using the pairwiseAlignment function from the Biostrings package. (https://bioconductor.org/packages/devel/bioc/vignettes/Biostrings/inst/doc/PairwiseAlignments.pdf) 

```{r}

# Change global options so that entire AA sequence prints out - https://statisticsglobe.com/reached-getoption-max-print-omitted-warning-in-r
options(max.print = 10000)

# Make sure bio3d package is detached if later part of code was already run since it also has a read.fasta function
# detach(package:bio3d, unload = TRUE)

# Read in fasta files and designate variables
# brca1_AA: https://www.ncbi.nlm.nih.gov/protein/AAC37594.1?report=fasta
brca1_AA <- read.fasta("BRCA1_AA.fasta", seqtype = "AA")
brca1_AA

# brca2_AA: https://www.ncbi.nlm.nih.gov/protein/AAB07223.1?report=fasta
brca2_AA <- read.fasta("BRCA2_AA.fasta", seqtype = "AA")
brca2_AA

```{r}
# Check class of data
#class(brca1_AA) 
#class(brca2_AA) 
# Convert from list to character strings to prepare for pairwise sequence alignment 
brca1 <- toString(brca1_AA)
brca2 <- toString(brca2_AA)
# Check for successful conversion
#class(brca1) 
#class(brca2) 


# Find optimal global alignment between protein sequences
brca_alignment <- pairwiseAlignment(pattern = brca1, subject = brca2, type = "global")

# Print optimal global alignment and score
brca_alignment

```

Structural Bioinformatics: This method is used to determine, analyze, and predict 3D structures of biomolecules. These biomolecules can be DNA, RNA, or proteins. The information obtained from this method can help us to strengthen our understanding of the relationship between structure and function. We will be using the read.pdb function from the Bio3D package to read in the PDB coordinate files for the crystal structure of the BARD1 BRCT domains (2NTE) and the crystal structure of a PALB2 / BRCA2 complex (3EUZ). Additionally, we will determine if and where the two structures have interacting residues with the binding.site function which is also from Bio3D.
(http://thegrantlab.org/bio3d_v2/html/read.pdb.html, http://thegrantlab.org/bio3d_v2/html/binding.site.html) 

```{r}
# Call bio3d package after pairwise alignment to avoid mixing up read.fasta functions
library(bio3d)

# Read PDB files and print composition summaries
BARD1_BRCT <- read.pdb("2nte")  
BARD1_BRCT 

PALB2_BRCA2 <- read.pdb("3euz")
PALB2_BRCA2

```

```{r}
# Display internal structures 
str(BARD1_BRCT)
str(PALB2_BRCA2)

# This wasn't a necessary function to do for this project. I was just trying it out. 
```

```{r}
# Determine and print interaction(s) between BARD1_BRCT & PALB2_BRCA with a cutoff of 50
# No interacting residues appeared until a cutoff of 42
BSR <- binding.site(BARD1_BRCT, PALB2_BRCA2, cutoff = 50, hydrogens = TRUE, byres = TRUE, verbose = TRUE)
BSR
```

## Plotting The Results
Dotplot: Creating a dotplot is one way to analyze and visualize our results from the pairwise sequence alignment that we performed earlier between the protein sequences of BRCA1 and BRCA2. The BRCA1 sequence is on the x-axis and the BRCA2 sequence is on the y-axis. Similarity between sequences is indicated by a black dot. The dots appears where the sequences match up. Sequences with high similarity would show a distinct line that spans the majority of main diagonal of the dotplot, whereas sequences with low similarity not present a visible diagonal line. 
(https://omicstutorials.com/interpreting-dot-plot-bioinformatics-with-an-example/)

```{r}
# Convert lists into vectors of single characters (necessary for dotplot function)
# https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/do.call 
# I'm assuming there's an easier way to do this part, I just couldn't figure it out. 

do.call("paste", c(brca1_AA, sep = ', ' )) -> BRCA1_vector
do.call("paste", c(brca2_AA, sep = ', ' )) -> BRCA2_vector

# Create and load dotplot 
BRCA1vs2 <- dotPlot(BRCA1_vector, BRCA2_vector)
BRCA1vs2
```

PyMOL: PyMOL is a molecular visualization system. It can import 3D protein structures straight from the PDB and has functions that allow for analyzing molecules. Pictures and movies can be downloaded from here. It provides many options for different ways to look at molecular structure and even how they align with other molecules. Residues, chains, and atoms can even be labeled and the structure sequences are available to look at as well. 
(https://en.wikipedia.org/wiki/PyMOL)


```{r}
# Call knitr package to include images
# https://www.jumpingrivers.com/blog/knitr-include-graphics-external/
library(knitr)
knitr::include_graphics(c("BARD1_BRCT.png","PALB2_BRCA2.png", "PyMOL align.png"))

# Install and call vembedr package to embed video
# install.packages("vembedr")
library("vembedr")
suggest_embed("https://youtu.be/S4B7enDivjM")
embed_youtube("S4B7enDivjM")

```

## Analyzing The Results
Pairwise sequence alignment & dotplot:

The alignment score generated from our sequence alignment was -2500.069, indicating that the BRCA1 and BRCA2 protein sequences had very low or no sequence similarity. There is a positive correlation with a higher score and a better alignment. This makes sense though since the BRCA1 protein sequence had nearly double the amount of amino acids than the BRCA2 protein sequence (3418 AA compared to 1863 AA). The dotplot is also very scattered. There are no apparent diagonal lines. There is actually a lot of empty space on the dotplot which makes sense again going back to how the protein sequences differed greatly in amino acid length. The low sequence similarities between BRCA1 and BRCA2 mean that they do have a significant amount of genetic differences. Although their mutations have similar outcomes, the genetic differences between them may be a factor as to why BRCA1 mutations lead to a significantly higher risk of ovarian cancer. 

Structural bioinformatics & PyMOL: 

From using the structural bioinformatics functions in the Bio3D package, I was able to access information from the PDB files of crystal structures related to BRCA1 and BRCA2, such as the number of models, atoms, chains, and the specific sequences for both structures The BARD1 BRCT domain is part of BRCA1. BRCA1 interacts with PALB2 which in turn interacts with BRCA2. With a cutoff distance of 50, 17 residues were identified as binding residues between the structures. No residues were identified until the cutoff distance surpassed 42. This means that there are interactions between the BARD1/BRCT complex and the PALB2/BRCA2 complex. Since the cutoff distance is a bit high, this could mean there are other molecules involved in the interactions of these complexes. I analyzed the same structures in PyMOL by attempting to align them. PyMOL generated an alignment score of 27 between the two structures. These protein complexes match up better than the protein sequences of BRCA1 and BRCA2 and supports the connection between BRCA1 and BRCA2 proteins. In PyMOL, I was able to show the individual structures as well as their alignment. I made sure to include a scene with their protein interfaces as well. The red lines that appeared between the structures are where residues were aligned. There was no evidence of high sequence or structural similarity between BRCA1 and BRCA2, but they do interact together. The protein structural differences in combination with the genetic differences are significant and can be used to explain at least part of the reason why BRCA1 mutations are more likely to cause ovarian cancer.  
