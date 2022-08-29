# QBB2022 - Day 1 - Homework Exercises Submission

*Exercise 1*

The error message we see is 'awk: illegal field $(), name "nuc"' for each of A, G, C, and T

To fix it, we need to set the identity of the variable with the -v flag within the awk command

Output from working script:

```
Considering  A
 354 C
1315 G
 358 T
Considering  C
 484 A
 384 G
2113 T
Considering  G
2041 A
 405 C
 485 T
Considering  T
 358 A
1317 C
 386 G
 ```
 
 It makes sense that A would be most commonly replaced with G (and vice versa) and C with T (and vice versa) since A and G are purines and thus similar sizes/shapes, while C and T are pyrimidines
