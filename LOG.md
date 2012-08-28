2012/08/22
==========

It liiiives....  I'm resurrecting this project because there is still a need
for solid modeling of deleterious mutation, especially at the genomic scale.
Maria Orive and her students have been examining interactions with population
size and patterns are complex... Martin Morgan's work from 2001 still needs
more followup, particularly on the possible dynamics of the mitotic mutation
'sieve' he proposed that favors the fixation of highly deleterious, highly
recessive mutation.  There needs to be more breadth to mutational effects being
modeled, and that is where K comes in.  It can duplicate Kondrashov's original
results from 1985, albeit more efficiently after changes from Charlesworth et
al. 1990, and it can duplicate selective interference from Lande et al. 1994.
My real interest here is in the nested mutation classes, and I'll get that
working, finally, and bang out a little paper.  Then if all goes well I'll move
on to incorporating mitotic mutation and look for the sieve.

**Practical issues related to development that I can recall just now include**

-  --- Make each mutation representation 'K' be a proper C++ class.  During
initial development I dummied that up while sticking with strict C, and the
'translation' to C++ from the original C is basically nonexistent.  Is it worth
doing?  I don't know, but I'll evaluate it.

-  --- Why don't I get proper fitness calculations at mutation rates of 0 and more
importantly 1?  Things were trending nicely at e.g. 0.001 and 0.999, but then
cratered.  This is a persistent issue that has bugged me for years whenever I
remembered back to this project and I never tracked it down to my satisfaction.

-  --- Get nested classes working.  The simplest approach was originally suggested
by Stewart Schultz, and that is to have the most deleterious class be
completely deleterious so that homozygotes need not be represented (they're all
dead).

**Research issues**

-  --- Establish research questions for the initial work.  A natural place to go
is to fix the total mutation rate (say at 0.1/0.2, 0.5, 1.0, 2.0, and 5.0
[going to high mutation rates is something I think is necessary to do]) and
examine the dynamics of inbreeding depression as the relative proportion
contributed by highly and slightly deleterious mutations changes.  How is
purging of highly deleterious mutations affected by slightly deleterious
mutations?  How does this interact with selfing rate, and with mutational
characteristics?



1/21/04
=======
-  --- Figure out the problem with `sum_KArray_n()` values from 1/16/04.


1/16/04
=======
-  -X- Results of test of [6][0][2][0]

The test passed, so at least the behavior is consistent following the changes.

-  -X- Why is `sum_KArray_n()` returning 2?  It's happening immediately after
  `compute_zygotes_n()`.  My guess is, that somehow we're not decreasing the
proportions from the gametes like we're supposed to be doing.  Every other
routine is well-behaved.  SOLVED!!  In `apply_zygotes()` we're dividing the
zygote distribution by 2 because the Male and Female gamete distributions both
sum to 1, and we're combining the M+F and the F+M contributions.  This was
repeated in `apply_zygotes_n()`, but was wrong, because there we're combining the
M[0]+F[0], M[0]+F[1], M[1]+F[0], and M[1]+F[1] contributions; all distributions
sum to 1 here as well, so we must divide the zygote distribution by 4.

-  -X- We should profile the [6][0][2][0] run on finn or adrienne.  Before I
   start poking around in routines I'm not certain of.  I've added Kprof script
to create `Kpg`, an executable that's compiled with `-pg`.  [WAS: Any way to speed
up `compute_gametes_n()`?  It takes a while no matter what the parameter set.  So
let's look into this...  `apply_gametes_n()` is the real culprit... Any way to
speed up `compute_mutation_n()`?] Here's the first few lines of the gprof results
on finn:

~~~~~~
      %   cumulative   self              self     total           
     time   seconds   seconds    calls  ms/call  ms/call  name    
      97.22   519.93   519.93       24 21663.75 21672.50  apply_gametes_full_n
      0.99    525.25     5.32       24   221.67   221.67  apply_zygotes_n
      0.42    527.51     2.26       69    32.75    32.75  copy_KArray_n
      0.27    528.98     1.47 14791887     0.00     0.00  fitness_n
      0.22    530.13     1.15       47    24.47    24.47  sum_KArray_n
      0.21    531.26     1.13 14791887     0.00     0.00  fitness_from_mfitness
      0.18    532.21     0.95       23    41.30   207.66  apply_selection_n
      0.15    532.99     0.78       23    33.91    58.38  normalize_KArray_n
      0.14    533.74     0.75       24    31.25    86.57  cumulative_fitness_n
      ...
~~~~~~

Thus, `apply_gametes_full_n()` is definitely a serious culprit, with all that
self time.  And it's definitely within that routine, because if we look at its
call graph:

    -----------------------------------------------
                  519.93    0.21      24/24          apply_gametes_n [2]
    [3]     97.3  519.93    0.21      24         apply_gametes_full_n [3]
                    0.16    0.03 3516552/3516552     lnbinomial [27]
                    0.02    0.00 1758276/1758276     lnpow_half [34]
    -----------------------------------------------

it's clear that its own code, and not called routines, is responsible for all
that time.  Now that I look at that code, it's got to be the looping, because
the code itself is not time-consuming.  The `is_lethal` code should help
short-circuit some of this, but it doesn't apply here.  The runonce speedup
must have worked very well, or did it... the calls to `lnbinomial` and `lnpow_half`
have the right relationship to each other, but I would have expected `lnpow_half`
to be called 250000 times (50 * 10 * 50 * 10); instead, that indicates that
it's been called 7.03 times that amount.  Could the runonce static flag not be
working?  No, my head is not working.  The way the runonce loops are set up is:

        for (v1=0; v1 <= KN->MI1; v1++) {
            for (v0=0; v0 <= KN->MI0; v0++) {
                for (n1=v1; n1 <= KN->MI1; n1++) {
                    for (n0=v0; n0 <= KN->MI0; n0++) {
                        t1 = lnbinomial(n0, v0);
                        t2 = lnbinomial(n1, v1);
                        t3 = lnpow_half(n0 + n1);
                        t4 = exp(t1 + t2 + t3);
                        p[n0][v0][n1][v1] = t4;
                    }
                }
            }
        }

Note that n1 and n0 do not start at 0, and are bounded by MI0 and MI1, not
MJ-anything.  So, I screwed up... plus debugging proved me wrong.

Anyway, it appears that `apply_gametes_full_n()` doesn't offer much hope for
speedup, because we have to do lots of nested loops, and the interior code is
pretty streamlined.  Maybe we could remove some `if()` statements that produce
debug code.

-  -X- Put the inclusion of `IF_DEBUG` in the above routine on an #if defined().

I added a `TIMECRITICAL_INCLUDEDEBUG[_n]` flag that is not defined by default.
This can be a generally-used check.

I'm going to do some speedups for `apply_mutation[_n]()`.

-  -X- Keep from calling `mut_term()` and `mut_term_n()`.
-  -X- Apply the `is_lethal` speedups here.

Also, 

-  -X- Speedups to `apply_self_progeny_n().`

Added `TIMECRITICAL_INCLUDEDEBUG_n` statement, and added `is_lethal` code.  I
still need to do the TODO announced therein, with better loop control.

-  -X- Apply warnings near `check_normalization_n()` calls.

I made this to be a single general warning at the beginning.

New profiling run doesn't reveal a speedup, at least for the no-mutation,
no-selection case. :-(

-  -X- Results of test of [6][0][2][0]

Success!! So with the mods, the test passed.

I've also modified the Kmake scripts to optimize with -O.

-  -X- Better loop control as described in `apply_self_progeny_n()`.

I've modified both `apply_self_progeny[_n]()` to remove all conditional
statements checking i, j, n, v, etc.  There are no longer any continue or break
conditions, nor os there the if() wrapper around the `s_self_n()` call.  It has
all been accomplished by adjusting the values of the loop endpoints.

-  -X- Results of test of [6][0][2][0]

Success.  This is still not the bread-n-butter test that requires the
production of selfed progeny... so...

-  -X- Do a LSS test with the one-class version, on adrienne.

FAILED!!!  Failed at h=0.  Now to back out the changes and see what it was.

-  -X- First, back out the use of the -O flag, see what that does.

No effect, good.

-  -X- Back out the loop controls on `apply_self_progeny()`.

Success!!  When backing out the v loop controls, at both h=0 and h=0.02.  This
is the problem, then, and backing out the `is_lethal` code in
`apply_self_progeny()` and `is_lethal` code in `apply_mutation()` is not necessary.
Now to narrow down the issue.

-  -X- Remove the (apparently redundant) extra loop check.

FAILED!!!  There may be an interaction with the `is_lethal` code, so

-  -X- Set the -nolethal flag and rerunning K.

FAILED!!!  Plus, it takes a *ton* of time now.  In fact, it was about 100 to
500 times slower.

-  -X- Restore the original loop controls on the v loop.

Test passed after this.

-  -X- Make the appropriate changes to the two-class version.
-  -X- Test the two-class version with [6][0][2][0].

Success!!  At home now, this is my task for the evening.  And it's completed,
it looks good.

Now, I should consider the tasks I need to complete in order to do a good
testing of the two-class model.  It's fast enough, now I need to know that it
makes sense.  I think there is an obvious first attempt -- have one neutral,
no-mutation class, and have the second be a lethal with h=0 and h=0.02, to see
if indeed the marginal distributions are identical to the single-class
versions.

-  -X- Run test with U=1/0, s=1/0, h=0/0.

FAILED!!  See below.

-  --- What is the deal with those wacked `sum_KArray_n()` values?

They seem to be settling down as the number of generations increases, but
honestly I don't get it.  I stopped the test, because a) the values seemed to
be converging to ~40, and b) they're whacked to begin with! 



1/15/04
=======

It's taking a long time to converge on finn (1105 min so far), so I'm going to
try some things to speed it up.  The first order of business, is going to be to
speed up the handling of lethals.  If a mutation class is lethal, then there is
no sense in processing mutant homozygous genotypes.  So, I'll add a "lethals"
flag of some kind to the single-class model and test that against the LSS data.
Critical here is the assumption that `apply_selection()` would already be
removing all the homozygous progeny anyway; maybe there's a way to speed up
`apply_selection()` too, if needed?  Later that day, it seems that
`apply_selection()` doesn't need much speedup, as the only function call within
the loop is to `fitness()` and that just accesses an array value.  It could buy
us some time, but it doesn't seem worth it at the moment.

In `initiate_model_state()`, we check to see if `K->fit_s == 1.0`; if so, set
`K->is_lethal` and print a message if `DEBUG_LETHALS` is set.  The heart of the
changes are in the reproduction routines.  Lethals can only be created by
`apply_self_progeny()` or copied by `apply_apomixis_progeny()`.  In both routines,
we set all classes in which j>0 equal to 0.0, and we don't look at classes in
which j>0.  In `apply_gametes_full()`, we don't look at homozygous mutant adult
classes.  There's no change in `apply_zygotes()`, as outcrossing never created
homozygous mutant adults to begin with.  Finally, in `apply_summed_progeny()`, if
`K->is_lethal` and `K->S[0] >0.0`, then print a warning that normalization is not
expected just prior to the call to `check_normalization()`.  There have been no
changes to `s_self()`, `a_apomixis()` or `o_outcross()`; all implementation details
are in the functions that ultimately call those routines.  The Kondrashov
versions have remained unchanged.

The first cut didn't work, both with h=0 and h=0.02, and I think I know why...
the stats calculations are screwed up because selfing is not producing any
progeny with j>0.  Perhaps we should be producing selfed progeny with j>0, just
not looking at classes with j>0...  That's being tested on finn & adrienne
right now, but maybe i only have to do that when I'm computing the final stats,
not each and every time.  So `apply_self_progeny()` could not create j>0 classes
until told to do so via an argument...

Now that I look at the selection code in `mean_fitness()`, it depends on
`sum_KArray()` for the summed progeny array being sensible... maybe with
`K->is_lethal`, I could always assume that to be 1.0... I'll go back to ignoring
selfed progeny with j>0 and see what values returned from `sum_KArray()`...
`sum_KArray()` is definitely returning values that are <1.0, so that change has
to happen.  The final stats values are non-sensible just as before, maybe what
I'll do is, make the `sum_KArray()` change in `mean_fitness()`, make the change to
`apply_self_progeny()` to only create j>0 classes if the argument says so, and
compute progeny at the very end only.

-  -X- `sum_KArray` change

This was done by creating `mean_fitness_allprogeny()`, making the changes
described above to that function, and calling that from `apply_selection()`.
`mean_fitness()` remains the same, as it is needed for other statistics
computation.

-  -X- `apply_self_progeny()`, `compute_self_progeny()`, reproduction routines.

These changes were done by adding a field `K->createlethal`, the value of this
field is checked everywhere `K->is_lethal` is checked, and the creation is
bypassed only if `K->createlethal==0`.  

-  -X- determine where to make modification so that stats will be sensible.

This change was done by adding an extra section of code to main().  Once
equilibrium is reached, we set `K->createlethal=1` and then do every `compute_*()`
routine except `compute_mutation()`.  The progeny arrays are created as normal,
and then the stats can operate on the progeny arrays at the end of the model
run.

Success!!  Tested on finn and adrienne against the Lande-Schemske-Schultz model
with U=1, s=1, h=0/h=0.02, everything is identical to the previous results
without the speedup.  Now to work this into the multi-class model.

By the way, to turn off this mechanism, just disable the setting of
`K->is_lethal` or `KN->is_lethal[]` in `initial_model_state[_n]()`.

-  -X- Command-line option to disable the `is_lethal` mechanism.

This option is -nolethal.

Two-class model changes:  Note that because of the two-array approach used by
the two-class model, the reproduction routines add proportions to an existing
array, and do *not* set proportions directly in elements of a dedicated array
for each mode of reproduction, as does the one-class model.  So, if we don't
have anything to add, do nothing!  Don't zero the array entries.

-  -X- KN fields
-  -X- `initiate_model_state_n()`
-  -X- `apply_self_progeny_n()`
-  -X- `apply_apomixis_progeny_n()`
-  -X- `apply_gametes_n()`

At the same time as I applied this, I added the `apply_gametes_full_n()` and
`adjust_gametes_n()` mechanism to the two-class model.

-  -X- `mean_fitness_allprogeny_n()`
-  -X- `apply_selection_n()`
-  -X- main() additional calls for stats

I have to review how stats are done here...

-  -X- check on stats generation, period

Also, two-class stats procedure may be wrong because they create progeny
without subjecting them to selection...  That doesn't matter, because the
one-class model looks at progeny fitness in the progeny array; selection
applies to the progeny after they leave the progeny array.  Two-class stats
reconstitute all the progeny from the arrays at equilibrium, and then do
`mean_fitness_n()` on the reconstituted progeny arrays.  As such, their normal
way of doing stats is close to the new lethals way of doing stats.  Mods to
allow for stats generation with lethals merely amount to turning on the
createlethal[] flag for the appropriate mutation class(es).

Test the two-class model on our old favorite, the [6][0][2][0] no-selection
no-mutation model, to see if everything remains consistent...

-  -X- Testing with [6][0][2][0]...  Passed, see 1/16/04.
-  -X- Wow, very strange that `sum_KArray_n()` returns 2.  Why???  Solved, see 1/16/04.
-  -X- Is `mean_fitness` computed correctly with two-class?  It appears so, but this has not been tested. 

Also, we may be able to speed up `apply_selection_n()` in a way that we declined
to for the single-class case.

Regarding mutation class sizes: the results with h=0 have mean_hetloci=176,
with var_hetloci=134.  These should be equal, but they can't be because
K->MI=200 and we're getting underrepresentation of the higher classes.  With
h=0.02, mean_hetloci=49 as does var_hetloci.  For both, mean_homloci=0 and
var_homloci=0, so lethal classes don't need a very high K->MJ; here K->MJ=40
and that seems to be plenty.  I want to do some test runs with different
mutation parameters to see what kinds of values their mean_hetloci etc. have.
The results will help determine what sizes to make the classes for the
two-class runs coming up.

-  --- Mutation class sizes with mildly deleterious mutations (s=0.1, h=0.3).



1/14/04
=======

I should confirm the observation that marginals for the two-class are identical
to the corresponding single-class; that would be expected, because of
independence, and is an important point for the manuscript.  It may be, that
with selection that is *not* true... though I can't think of why that would be,
at the moment.

I need to see if the command line works for two-class.  If it doesn't, I'll fix
it.  It appears that it does, with -U0/-U1, -s0/-s1, -h0/-h1.

I'll start a test run with selection, giving class 0 the same mutation
parameters as I used to verify the single-class example yesterday against the
LSS paper; class 1 will start with all members having 2 heterozygous mutations,
but will otherwise have no mutation or selection.  The arguments are thus -U0 1
-s0 1 -h0 0/0.02 -U1 0 -s1 0 -h1 0.

I had started the run with S=0.01, but it occurred to me that it's a better use
of time to start the run with S=0.99, as the approach to equilibrium with one
class was much more quick with a high selfing rate (30 gens at S=0.99 vs. 255
at S=0.01).  It's also going to dump a 50x10x50x10 mutation genotype array at
the end.  It's already got 56 minutes of CPU time, even at S=0.99...  Now it's
got 119 minutes on finn; so I added the debug flag `DEBUG_GENERATIONS`, that can
be used to follow the number of the current generation, and whether execution
has stopped because the generation cutoff was exceeded.  I will go ahead and
test this on adrienne.  OK, that's running, but it wasn't fine enough
resolution for what I want to see at this point, so I turned on `DEBUG_TRACE1`,
`DEBUG_EQUILIBRIUM` and `DEBUG_TRUNCATE` as well.  I killed the finn and adrienne
runs and restarted them via telnet from anna, as I'm planning on taking the
laptop home for the evening.

**Issues for the future**

-  --- One issue to keep an eye on, is that the number of mutation classes used
   is 200 hets, 40 homs for single-class, and 50 hets, 10 homs for two-class.
This speeds things up but will probably be too small for higher mutation rates
with lower selection; let's see what happens with the above test runs.

-  -X- Now that S=0.99, we seem to be taking a long time in `apply_mutation_n()`,
   is there some way to speed those up?

-  -X- Now that S=0.99, we seem to be taking a long time in
   `apply_self_progeny_n()`, is there some way to speed those up?  Did a bunch of
stuff on 1/16/04.



1/13/04
=======

Single-class confirmed against Lande-Schemske-Schultz paper,
results/test_LSSdata.xls, for U=1, h=0 & h=0.02, s=1, S=0.01 .. 0.99.

I'm computing inbreeding depression using mean fitnesses, *not* the dummied-up
approach that Schultz used.  The Schultz approach is valid at all selfing
rates, but my version is much faster.  The drawback is that it doesn't work if
S=0 or S=1.  For the time being, I have to make my selfing rate endpoints be
just inside the actual endpoints (e.g. 0.001 and 0.999), so that we have
production of selfed and outcrossed progeny for the mean fitness calculations.

It occurred to me that we might be able to speed up the approach to equilibrium
if we do a `-save_savefile` from the previous settings, and use it as a
`-load_savefile` for the next settings.  That didn't work on my initial attempt,
because the model would exit because it detected equilibrium each time.  Maybe
if we allow the model to run a little before checking for equilibrium if we're
loading a file, then that would be the way to do it... anyway, another time.

Tested the nested model starting with 6 heterozygous mutations in class 0 and 2
in class 1.  I had U=0, h=0, s=0, I was just interested in seeing how the
distribution approached equilibrium.  The marginal distributions of each class
at equilibrium appear to be equivalent to an independent non-nested equilibrium
distribution for each 



6/25/03
=======

On plane from St Louis to Miami, 1.5 hrs late

Trying to get back into the spirit of this model Sorting out the problems with
the inbreeding depression computations, starting with the single-locus case.
Working out of the "savenew.txt" saved class proportions file.  Someday, figure
out why there is the initial iteration when loading at-equilibrium class
proportions from that file.

within `stats_inbreeding_depression()`
Debug watch commands
`dump_KArray(K,K->x,10,0,1)`
`dump_KVector1(K,allmgam,10,1)`
`dump_KVector1(K,allfgam,10,1)`

Strange, the allfgam vector appears to be all zeros That's because when it's
set in `apply_gametes()`, it's still f-ing scaled by the outcrossing rate
K->O[g].  I've got to abstract that out of that function.... or should I?  I
want to keep the interface a step-by-step series of logically linked function
calls with minimal intervening statements, and having something like
`apply_gametes(...)`; `mating_system(...)`; `apply_zygotes(...)`.  Anyway, I've
figured it out.  I should have been calling `apply_gametes_stats(...)` for that
first call for the inbreeding depression calcs.  I'll change that, and see what
happens.  Note also the difference in the interface between these, I have to
remind myself what the fi, fj, fg args in the `_stats(...)` version mean.

    void        apply_gametes       (KConfig K, 
                                     KVector1 mgam, KVector1 fgam,
                                     KArray from)
    void        apply_gametes_stats (KConfig K, 
                                     KVector1 mgam, KVector1 fgam,
                                     KScalar fromval,
                                     KInt fi, KInt fj, KInt fg)

Now I remember, the fi, fj, fg args refer to the class from which we generate
the resultant distributions.  So, I'll augment `_stats(...)` to use the whole
dimension (taken from K) if the corresponding arg is <0.  Actually, I didn't do
that.  I split `apply_gametes()` into two functions, `apply_gametes_full()` and
`adjust_gametes()`.  Both have identical interfaces, but `apply_gametes_full()`
creates gamete classes from the full set of source classes, and does no
"post-processing" to the gamete proportions.  The new function `adjust_gametes()`
is where that post-processing occurs, and the only thing it does at the moment
is scale fgam by K->O[g], as used to occur in `apply_gametes()`.  To generate
allfgam and allmgam the `stats_inbreeding_depression()` function calls
`apply_gametes_full()`.  Standard model reproduction calls `apply_gametes()` just
as before.  Now that I've run that, it doesn't work in the short example I've
been loading. Why would that be?  It's a little too late to be worrying about
it now, we're getting ready to land.  It's 12:32 am EST, 9:32 am Chico time,
and I'm running on about 3 hrs sleep.


