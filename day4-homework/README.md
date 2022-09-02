# QBB2022 - Day 4 - Homework Exercises Submission

**Part A**

```
probs = numpy.around(numpy.arange(0.55, 1.05, 0.05), decimals=2)[::-1]
```

* The numpy.arange function creates an array in which 0.55 is the start value and 1.05 is the stop value (not inclusive of the stop value), and all the intervening entries are steps up from the start value by 0.05, which is the 3rd argument. 
* This array is then fed into the numpy.around function, which rounds it to 2 decimal places (even though the array was already rounded to 2 decimals or lower)
* Finally, the actual probabilities list is set to be the array produced by the around function, but stepping backwards from the final entry, producing the following list:

```
[1.   0.95 0.9  0.85 0.8  0.75 0.7  0.65 0.6  0.55]
```

**Part C**

* In both the corrected and uncorrected charts, power increases with increasing probability (further away from the expected probability of 0.5) and with increasing number of tosses.
* In the corrected chart, an even greater number of tosses and greater difference between the expected and observed probabilities are needed for the test to have a higher power score.

**Part D**

* This study is focused on the biological phenomenon of transmission distortion, which is when a diploid, heterozygous organism passes on its alleles to offspring in proportions other than the expected 50% as predicted by Mendel's Law of Segregation (or, described another way, the organism produces unequal numbers of gametes possessing each of the two alleles).
* Comparison between my simulation and the paper: 
	- "Probability" in my paper is analagous to "Transmission rate" in the paper, as both refer to the phenomenon which is tested deviating from the expected 0.5 probability. In both cases, the power of the test increses as probability/transmission rate deviates further from the expected value of 0.5.
	- Number of tosses and number of sperm are also analogs, referring to the number of "trials" or "test subjects"; again, in both cases, increasing the number increases the power of the statistical test.
	- n_iters in my simulation corresponds to the "1000 independent simulations" run in the study to compute the power for each study design
	- One notable difference is that in my simulation, multiple testing correction decreased power of the tests more dramatically for low-toss number tests regardless of probability, whereas in the study, low-transmission rate studies were more sensitive to losing power when subjected to multiple testing correction.
	- Also, in my simulation, the number of tosses is approximately on a log scale, while in the study, the number of sperm axis is on a linear scale.
	- Both simulations use a binomial test because the aim is to test the statistical significance of a deviation from the expected binomial distribution of observations that can take two values: heads vs. tails, or allele 1 vs. allele 2.