## Week 6 -- 10 points possible

1 + 1 + 1 + 1 + 1 + 1 + 1 + 1 + 1 + 1 + 1 + 1 + 1 = 13 of 10 points possible

1. Given data question: What percentage of reads are valid interactions?

2. Given data question: What constitutes the majority of invalid 3C pairs?/What does it mean?

Dangling pairs were the majority of invalid  3C pairs. I found this in the plot `plotHiCFragment_dCTCF.pdf` which is in the `analysis/pic/dCTCF` subfolder.

3. Script set up to (0.5 pts each)

  * read in data  
  * Filter data based on location  

4. Script set up to log transform the scores

* you don't need to use a `for` loop for this, like the following. `data1['score'] = numpy.log(data1['score'])`. It applies the log to all the scores in the array/along the axis.

5. Script set up to shift the data by subtracting minimum value

* You can find the minimum without a for loop:
  `minimum1 = np.amin(data1['score'])`. It finds the minimum value in the whole array.
* then you can subtract the minimum without a `for` loop: `data1['score'] -= minimum1`. It subtracts the minimum value from each of the scores in the array

6. Script set up to convert sparse data into square matrix

  * again you don't need a for loop for this. You can use the `filtrd` array, or whatever it's name would be after the other steps, something like `mat[filtrd["F1"], filtrd["F2"]] = filtrd["score"]` and this will fill in the whole `mat` array for each combo/row of F1, F2, and score values. You do need to make sure you've subtracted that min_bin or start_bin from F1 and F2 before you do that though

7. Script set up to (0.33 pts each)

  * remove distance dependent signal
  * smooth
  * subtract

8. Turned in the plot of the 3 heatmaps (ddCTCF, dCTCF, and difference) for subset dataset (0.33 pts each ddCTCF/dCTCF/difference)

9. Turned in the plot of the 3 heatmaps (ddCTCF, dCTCF, and difference) for full dataset (0.33 pts each ddCTCF/dCTCF/difference)

10. Heatmap questions (0.33 pts each)

  * See the highlighted difference from the original figure?
  * impact of sequencing depth? --> what you've said seems true as well; for me the biggest noticeable difference is that increased sequencing depth decreases noise/smooths the signal in the data
  * highlighted signal indicates?

Possible Bonus points:

1. Insulation script set up to (0.25 pts each)

  * read in data
  * filter data based on location
  * log transform the data
  * shift the data by subtracting minimum value

2. Insulation script set up to (0.5 pts each)

  * convert sparse data into square matrix
  * find the insulation score by taking mean of 5x5 squares of interactions around target

3. Turned in the plot of the heatmap + insulation scores below (0.5 pts each panel)

Interesting use of the for loop and dictionary storage to process the two matrices! You can use a function to do something similar, and just call the function twice.
