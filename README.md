K
=

`K` implements Kondrashov's (1985) infinite-population, infinite-sites model
for genomic mutation, tracking the proportion of the total population having _i_
sites homozygous for a mutation and _j_ sites heterozygous for a mutation.

The major extension provided by `K` is (so far) in allowing for each of the
(_i_,_j_) bins for a given class of mutations to have its own submodel for a
class of mutations.  In the simplest case (and I can't think that I'd want to
be more complex here...) this means that there are two classes of mutations
possible within individual of the population, one containing (_i_,_j_) bins and
one containing (_k_,_l_) bins with the interpretation of _k_ and _l_ being
identical to that for _i_ and _j_, respectively.

This project has lain dormant for nearly 8 years, and is now under active
development.

