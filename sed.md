```bash
sed 's/functional_gene/#f7f9ad/;s/partial_sequences/#E7C951/;s/lost_gene/#A6325A/;s/pseudogene_pre_stop_codon/#0620b0/;s/frame_shift_as_pseudogene/#4FB3A4/' filename
```
I use sed to replace, and different replacements can do at the same time.


How to replace multiple patterns at once with sed?
Maybe something like this: 
```
sed 's/ab/~~/g; s/bc/ab/g; s/~~/bc/g'
```
Replace ~ with a character that you know won't be in the string.
