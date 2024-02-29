Assignment 5 Due February 23, 2024 at 11:59pm

Please provide the exact command line arguments you used to generate your results.
{How to run gene_expression.py}
python3 gene_expression.py raw_counts.txt
-
Question 1:
{How many genes are left after removing genes with zero expression in all samples?}
23259
-
Question 2:
{How many genes are left after removing genes where 20 or more samples have cpm < 1?}
14464
-
Question 3:
{What is the range of library sizes (min, max)?}
(11464513, 20839468)
-
Question 4:
{What is the range of library sizes (min, max) after normalization?}
13671809, 22214610
-
Question 5:
{Compare the two library size bar charts you made. How did the distribution of library sizes change after normalization?}
The distribution between the 
-
{Briefly discuss why it is important to normalize your RNA-seq data.}
It is important to normalize data (as a whole, not just RNA-Seq data) to place the data on a common scale
to account for any abnormalities or outliers in the datasets that could potentially skew potential conclusions made.
For RNA Seq data in particular, normalization accounts for features that could skew the number of reads a particular gene 
has such as overall gene length, GC-content, and sequencing depth [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6171491/].
-
Question 6:
{What are the top ten differentially expressed genes according to your FLD analysis? (Copy and paste your function's output.)}
AK000953: 1.0184876950008257
RAB30: 0.9882430732465937
DBNDD1: 0.9576122203523499
CTDSPL: 0.9151622166434654
GRB14: 0.8851667900187445
YTHDC2: 0.8516774330582303
MYLK4: 0.842607579063623
ABCG1: 0.7465802521032551
FNDC5: 0.739653965436164
BC010186: 0.7372827164007741
-
{Do these genes make sense given the tissue and groups in the experiment?}
Some of these genes make sense given the nature and design of this experiment. In particular, 
GRB14 (weight-loss responsive gene), FNDC5 (fibronectin - released during exercise), MYLK4 (myosin kinase)
appear to be most relevant in their differential expression. 
-
Question 7: 
{Does your result point toward one gene with large effect or many genes with small effects?}
Given that many genes were differentially expressed across both time points (before & after),
this would lead me to reasonably conclude that these genes may be additive with their effects (there is no one gene with a large effect)
-
{Does RNA-seq expression data always give researchers a clear answer?}
No! Sometimes it can lead to more questions, which can be both a positive and negative thing depending on your mentality and stage of the PhD.
In all seriousness, RNA Seq can create more directions from which a researcher can consider the effects of genes that would not have been considered otherwise,
however care must be taken to keep things in biological perspective/relevance
-
Question 8:
{How does the study design of this experiment relate to the assumptions made when studying gene expression data?}
-
Question 9:
{If you were going to spend time and money following up on one of these top ten genes, what would be your candidate and why? (There could be many correct answers.)}
I would follow up on the FNDC5 (fibronectin) or GRB14 (weight-loss gene), given their already verified and demonstrated
function. FNDC5 could be interesting given that its product is produced as as result of exercise and is involved in the production of brown fat [https://www.ncbi.nlm.nih.gov/gene/252995].
GRB14 has been demonstrated to be involved in the receptor kinase pathways associated with insulin, which was the entire aim of the proposed
RNA Seq experiment. 
-