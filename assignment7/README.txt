Assignment 7 Due March 22, 2024 at 11:59pm midnight

Please provide your script with comments where necessary to explain your code. 
python3 demultiplex.py call.matrix.csv
-
Question 1:
{Control Count}
19
{Treatment Count}
24
{Non-Determined Count}
7
-
Question 2:
{Control Count}
22
{Treatment Count}
24
{Non-Determined Count}
4
-
Question 3:
{3 strategies to maximize group assignment}
1. For assignment of treatment and control groups, more than 2 CellTags can be utilized to prevent any potential overlap
2. Although it can be precarious, increasing the acceptable Hamming Distance (from 1 to 2 in this exercise for example),
would help to further delineate the different populations (e.g. GTGTAGC vs ATGTTGC - HD of 2 whereas ATGTAGC vs ATGTTGC - HD of 1)
3. Utilize an allowlist for acceptable CellTags between both populations
-