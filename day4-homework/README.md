#QBB2022 - Day 4 - Homework Exercises Submission

**Part A**

```
probs = numpy.around(numpy.arange(0.55, 1.05, 0.05), decimals=2)[::-1]
```

* The numpy.arange function creates an array in which 0.55 is the start value and 1.05 is the stop value (not inclusive of the stop value), and all the intervening entries are steps up from the start value by 0.05, which is the 3rd argument. 
* This array is then fed into the numpy.around function, which rounds it to 2 decimal places
* Finally, the actual probabilities list is set to be the array produced by the around function, but stepping backwards from the final entry, producing the following list:

```
[1.   0.95 0.9  0.85 0.8  0.75 0.7  0.65 0.6  0.55]
```

