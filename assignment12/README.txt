Bio5488 2024 Spring Assignment 12

Due: Friday April 19th 2024, by 11:59 pm

What to turn in:
1. assignment_12.R
2. README.txt


#=======================================#
Q1: What is the approximate range for the number of sequences observed per ASV?

The approximate range is between 0 and 500 sequence reads per amplicon sequence variant (ASV)

#=======================================#
Q2: Which sample has the highest number of reads?

Naivestool6 has the highest number of reads (~22,400 reads)

#=======================================#
Q3: Is there a statistical difference in the number of reads between the facilities?

Given that the p-value was 0.13 (threshold is <= 0.05), this indicates that the difference in the number of reads between facilities is not statistically significant

#=======================================#
Q4: Do you see any obvious groupings of samples by this analysis? What does it suggest to you? Are there any obvious outliers?

Per the MDS plot, there are 3 clusters:
1. Naivestool 1-3
2. Naivestool 4-6
3. Stool 1-6

From these groups, there are no distinct outliers per cluster

#=======================================#
Q5: Do you observe any phyla that are low in prevalence and abundance in this dataset? 

Actinobacteria demonstrates both low prevalence and abundance in samples per both our generated plot & dataframe as compared to other phyla

#=======================================#
Q6: Are there any obvious differences in the composition of the microbiota between Facility 1 and 2? 

In Facility 2, there appears to be a higher fraction of Bacteriodetes, along with a slight increase in Actinobacteria (though minimal). Other phyla such as Proteobacteria & Verrucomicrobia appear in similar fractions between the two facilities

#=======================================#
Q7: Are there any significant differences in alpha diversity between samples from Facility 1 and 2? 

According to the Shannon Diversity, there is a statistically significant difference in alpha diversity between Facility 1 & 2 (p-value = 0.039); in the other two metrics measured, richness & Pielou's evenness, there was not a statistically significant difference between the two facilities (p-value of 0.12 & 0.62 respectively)

#=======================================#
Q8: Are the samples statistically significantly different by Site? 

PCoA does not illustrate statistical significance, but rather dimension reduction. With this in mind, one can observe that the two facilities group together independently (the two sites are different) with no apparent overlap in the clusters present

#=======================================#
Q9: Taxa enriched in Facility 2 are on the left of the volcano plot, while taxa enriched in Facility 1 are on the right. Do you notice any broad trends? How does this compare with the community composition plots we looked at earlier? 

Facility 2 illustrates an enrichment of Bacteriodetes, while Facility 1 illustrates an enrichment in Firmicutes. These trends are also demonstrated in the earlier generated community composition plot. For all other taxa, they appear to be evenly distributed between the two facilities. 

#=======================================#