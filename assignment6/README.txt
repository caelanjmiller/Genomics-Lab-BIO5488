Assignment 6 Due Mar 01, 2024 at 12am midnight (by end of Friday)

Please provide the exact command line arguments you used to generate your results.
{How to run find_hic_freqs.py}
python3 find_hic_freq.py genomic_bins.tsv matrix_gm.tsv matrix_k562.tsv TADs.csv
-
Question 1:
{Which cell line have merged TADs in this genomic region}
From the heatmaps, one can tell that the GM128718 cell line has merged TADs
-
Question 2:
{Quantitative evidence that the two TADs are merged in one cell line but retained in another}
From the TSV I created it's evident that the inter-TAD contact frequencies for the GM128718 cell line
are significantly higher than in the K562 cell line (~340,000 vs ~50,000), potentially suggesting that
the two TADs are merged in that cell line 
-