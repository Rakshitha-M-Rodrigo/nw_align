nw_anchored matlab function which aligns sequences using the needleman-wunsch algorithm with the
option of passing in known matches between the two sequences.

You can call the function in matlab by:
nw_anchored('seq1.fa','seq2.fa', [matches.txt])

This will return an alignment score where matches = 1, mismatches = -3, and gaps = -2, 
as well as the alignment of the two sequences.

Note that the 3rd parameter is optional (and currently not fully functional*), and in order for it to work properly
the first sequence listed in the matches text file must be input as seq1 because of the
order.

In this same directory are the HOX and PAX sequences that have been randomly permutated 
10,000 times and then aligned. As well as the alignment.txt which contains the alignments
for both the sequences (and the TITIN sequences) and their scores.

*alignment using known matches is commented out as it's not currently implemented correctly.
