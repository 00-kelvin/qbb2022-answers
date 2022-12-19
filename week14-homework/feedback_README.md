Excellent work. One small note:

1. For question 3B: when you're estimating the total length of the bins, unfortunately you can't just pipe the output of `grep -v ">" bin.*.fa` to `wc` because when you grep on multiple files at once, it also outputs which file the line came from for each line. Which means you're going to end up with a bunch of extra characters. On top of this, you're counting new line characters as well with this method, although those can be removed with something like `tr -d '\n'`. The for loop you used for 3C would work just fine though! (no points deducted)

Otherwise great job!

(10/10)
