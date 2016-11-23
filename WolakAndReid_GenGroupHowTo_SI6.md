*Supporting Information to*
**Accounting for genetic differences among unknown parents in microevolutionary studies: How to include genetic groups in quantitative genetic animal models**
================
Matthew E. Wolak & Jane M. Reid

-   Appendix S6: Construction and implementation of Q and A\*
    -   6.1 A small example pedigree
    -   6.2 Constructing Q
    -   6.3 Constructing A\*
        -   6.3.1 Genetic groups on top or bottom of A\*?
        -   6.3.2 Removing singularities and other problems within A\*
    -   6.4 How to fit genetic groups in MCMCglmm, ASReml-R, ASReml-standalone, & WOMBAT
        -   6.4.1 Tutorial dataset `ggTutorial`
        -   6.4.2 `MCMCglmm`
            -   6.4.2.1 Preparing a pedigree with genetic groups
            -   6.4.2.2 Fixed explicit genetic group effects with Q (from `nadiv`)
            -   6.4.2.3 Random implicit genetic group effects with A\* (from `nadiv`)
        -   6.4.3 `ASReml` in `R`
            -   6.4.3.1 Preparing a pedigree with genetic groups
            -   6.4.3.2 Fixed explicit genetic group effects with Q (from `nadiv`)
            -   6.4.3.3 Fixed implicit genetic group effects with A\* (from `nadiv` or `asreml`)
            -   6.4.3.4 Random implicit genetic group effects with A\* (from `nadiv`)
        -   6.4.4 ASReml Standalone
            -   6.4.4.1 Preparing a pedigree with genetic groups
            -   6.4.4.2 Fixed explicit genetic group effects with Q (from `nadiv`)
            -   6.4.4.3 Fixed implicit genetic group effects with A\* (from `nadiv`)
            -   6.4.4.4 Fixed implicit genetic group effects with A\* (from ASReml)
            -   6.4.4.5 Random implicit genetic group effects with A\* (from ASReml)
        -   6.4.5 WOMBAT
            -   6.4.5.1 Preparing a pedigree with genetic groups
            -   6.4.5.2 Fixed explicit genetic group effects with Q (from `nadiv`)
            -   6.4.5.3 Random implicit genetic group effects with A\* (from `nadiv`)
            -   6.4.5.4 WOMBAT's genetic groups (Random group effects with Q from `nadiv`)
    -   References: Appendix S6
-   Software version information

Appendix S6: Construction and implementation of Q and A\*
=========================================================

Fitting animal models with genetic groups first necessitates formulating the appropriate **Q** and/or <b>A\*</b> matrices. First, we generally describe how **Q** and <b>A\*</b> are formed, with an emphasis on the code used to do this in the `R` package `nadiv` (\>v2.14.0; Wolak 2012). This is followed by explanations of minor changes that can be made to <b>A\*</b> to improve animal model convergence and how to deal with warnings from software regarding singularities and non-positive definite matrices. Next are sections with code that illustrate fitting genetic group animal models in the `MCMCglmm` and `asreml` packages in `R`, the **ASReml** standalone program, and finally the **WOMBAT** standalone program. In each of these four sections, code and examples are provided to fit genetic groups using the `nadiv` package in `R` to create the **Q** and <b>A\*</b> matrices that can then be supplied to these animal model fitting programs (in some cases outside of the `R` environment), followed by examples using existing functionality within each software program (see Table S6.1 for a summary of these programs and their capabilities).

The tutorials for these software programs all use the simulated dataset (`ggTutorial` available in `nadiv` or to download in **WOMBAT** and **ASReml** standalone formats) depicted in Figs 2 and 4 of the main text and described below. In *Appendix S4*, the `R` code in the `nadiv` function `simGG()` that is used to simulate the `ggTutorial` dataset is explained. Details of package and software program versions are available at the end of the supporting information.

**Table S6.1.** Overview of animal model software highlighted in the tutorials. The genetic group capabilities row contains a description of the functionality that is available solely within each software program. The **Form Q** and **Form** <b>A\*</b> rows indicate if each program can create the **Q** or <b>A\*</b> matrices, respectively, from a supplied pedigree. We indicate whether fitting genetic groups explicitly (as fixed effects) or implicitly within the random effects structure (as either fixed or random effects) is recommended (+) or not (-) when using each software program.

|                                | **`MCMCglmm`**                                                           | **`ASReml`** in **`R`**                             | **ASReml**                                                                                                           | **WOMBAT**                                                                                                                                                                                    |
|:-------------------------------|:-------------------------------------------------------------------------|:----------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Genetic group capabilities** | None in `MCMCglmm`, but can use <b>A\*</b> or **Q** created with `nadiv` | The `asreml.Ainverse()` function creates <b>A\*</b> | `!GROUPS` qualifier to the pedigree creates <b>A\*</b> and see the `!LAST` and `!GOFFSET` qualifiers (details below) | Random genetic group effects (with separate variance component from the additive genetic variance) by including `GENGROUPS` in a `SPECIAL` block with a supplied **Q** (created with `nadiv`) |
| **Form Q**                     | No                                                                       | No                                                  | No                                                                                                                   | No                                                                                                                                                                                            |
| **Form** <b>A\*</b>            | No                                                                       | Yes                                                 | Yes                                                                                                                  | No                                                                                                                                                                                            |
| **Explicit approach**          | +                                                                        | +                                                   | +                                                                                                                    | +                                                                                                                                                                                             |
| **Implicit approach**          |                                                                          |                                                     |                                                                                                                      |                                                                                                                                                                                               |
|     *Fixed*                    | -                                                                        | +                                                   | +                                                                                                                    | -                                                                                                                                                                                             |
|     *Random*                   | +                                                                        | +                                                   | +                                                                                                                    | +                                                                                                                                                                                             |
| **Free to use**                | Yes                                                                      | No                                                  | No                                                                                                                   | Yes                                                                                                                                                                                           |
| **Reference**                  | Hadfield 2010                                                            | Butler et al. 2009                                  | Gilmour et al. 2014                                                                                                  | Meyer 2007                                                                                                                                                                                    |

<!---------     ------------    -----------   ------------->
6.1 A small example pedigree
----------------------------

Below we demonstrate how the matrices **Q** and <b>A\*</b> are constructed using an example pedigree from Quaas (1988). A detailed version of the Quaas pedigree is available in the `nadiv` package in `R` as the `Q1988` dataset:

<br>

| id  | dam | sire | damGG | sireGG | phantomDam | phantomSire | group |
|:----|:----|:-----|:------|:-------|:-----------|:------------|:------|
| g1  | NA  | NA   | NA    | NA     | NA         | NA          | NA    |
| g2  | NA  | NA   | NA    | NA     | NA         | NA          | NA    |
| a   | NA  | NA   | g1    | g1     | NA         | NA          | g1    |
| b   | NA  | NA   | g2    | g2     | NA         | NA          | g2    |
| c   | NA  | NA   | g2    | g2     | NA         | NA          | g2    |
| d   | NA  | NA   | g1    | g1     | NA         | NA          | g1    |
| e   | NA  | NA   | g2    | g2     | NA         | NA          | g2    |
| A   | NA  | NA   | g1    | g2     | a          | b           | NA    |
| B   | A   | NA   | A     | g2     | A          | c           | NA    |
| C   | NA  | NA   | g1    | g2     | d          | e           | NA    |
| D   | B   | C    | B     | C      | B          | C           | NA    |

The table above is essentially just a pedigree (first three columns indicate an individual's unique identifying code, dam, and sire). However, two rows have been added at the top of the pedigree for the unique genetic groups and extra columns are included for ease of using all possible pedigree formats accepted by the `nadiv` functions that implement genetic group methods.

In the dataset, "g1" and "g2" are two genetic groups, lower case letters are phantom parent identities, and upper case letters are observed individual identities. The genetic groups and phantom parents are assigned as in Quaas (1988). Additional details regarding the `Q1988` dataset in `nadiv` are in the `R` help file that can be viewed by running the command `?Q1988` in an `R` session.

Below, we will only deal with observed individuals in this dataset, such that the subset of `Q1988` we will use is:

<br>

``` r
Q1988sub <- Q1988[-c(3:7), c("id", "damGG", "sireGG")]
```

| id  | damGG | sireGG |
|:----|:------|:-------|
| g1  | NA    | NA     |
| g2  | NA    | NA     |
| A   | g1    | g2     |
| B   | A     | g2     |
| C   | g1    | g2     |
| D   | B     | C      |

Note, this is not the only format for including genetic group information in functions of the `nadiv` package. To see a full description of alternative formats see the function help files, particularly for the `ggcontrib` (to make **Q**) and `makeAinv` (to make <b>A\*</b>) functions.

<!---------     ------------    -----------   ------------->
6.2 Constructing Q
------------------

As described in the main text, the matrix **Q** contains the fractional contributions from each of *r* genetic groups to each individual's genome, calcluated from the pedigree (see also Robinson 1986; Mrode 2005, ch. 2). Therefore, if there are *n* individuals in the pedigree (e.g., *n*=4 in `Q1988sub`), **Q** will be an *n* row by *r* column matrix. In Henderson's (1976) decomposition of the numerator relationship matrix (**A**=**TDT'**), the lower triangle matrix **T** reflects the contribution of each individual in the pedigree to every other individual. In other words, **T** traces the flow of genes from one generation to another. Therefore, if we add the *r* genetic groups to the top of the pedigree (originally of length *n*) as identities with unknown dam and sire (NA) and fill in genetic group labels as dam and sire identities of every individual with a missing parent (such as how `Q1988sub` is arranged) we can obtain the **T** matrix for a pedigree of length *n*+*r*. The first *r* columns of the **T** matrix contain the fractional contributions from each genetic group. Therefore, **Q** for the *n* individuals in a pedigree is the *r*+1 to *r*+*n* rows and the first *r* columns of **T**.

The `nadiv` function `ggcontrib()` constructs **Q** for a given pedigree with *r* potential genetic group contributions to each of the *n* individuals in the pedigree. For example, **Q** for `Q1988sub` is obtained:

``` r
(Q <- ggcontrib(Q1988sub))
```

    ##      g1    g2
    ## A 0.500 0.500
    ## B 0.250 0.750
    ## C 0.500 0.500
    ## D 0.375 0.625

Each row in **Q** (each entry is a proportional contribution) sums to one and offspring inherit half of each genetic group contribution from each of their parents. Individuals `A` and `C` in `Q1988sub` each have a phantom dam from genetic group "g1" and phantom sire from genetic group "g2". These individuals correspond to the first and third rows of **Q** above. As expected, the two genetic groups (each column in **Q**) contribute equally to the genomes of `A` and `C`.

Individual `B` has observed dam `A` and a phantom sire from genetic group "g2". Row two of **Q** above is a sum of each parents' row (genetic group contributions) divided by two: dam `A` `[0.500 0.500]/2` plus phantom sire "g2" `[0.000 1.000]/2`. In other words, individual `B` gets half of dam `A`'s "g1" group contribution (`0.500/2`) and its phantom sire from "g2" does not contribute to the proportion of `B`'s genome derived from "g1". Similarly, individual `B` gets half of dam `A`'s "g2" group contribution (`0.500/2`) plus half of the phantom sire "g2" contribution (`1.000/2`).

Finally, individual `D` has observed dam `B` and observed sire `C`. Therefore, row four of **Q** above is the sum of each parents' row (genetic group contributions) divided by two: dam `B` `[0.250 0.750]/2` plus sire `C` `[0.500 0.500]/2`.

Implementing an animal model with genetic group effects fitted explicitly as separate fixed regressions requires the **Q** matrix. Because every row of **Q** sums to one, it is not possible to solve an animal model to obtain estimates for every genetic group effect and an overall model intercept. Thus, genetic group effects are not estimable themselves, but can only be estimated when expressed as deviations (or functions) from other group effects (Lynch & Walsh 1998, p.839-841). Therefore, a column/genetic group is chosen as the "reference group" and is assigned a genetic group effect of 0 (see practical implementation of this below). All other group effects will represent deviations from this reference group. For example, we might choose genetic group "g1" in the `Q1988sub` pedigree as our reference group. The estimate of the fixed regression slope from the second column of **Q** equals the "g2" genetic group mean expressed as a deviation from the reference group ("g1") mean of zero. If we had used a pedigree with more than two genetic groups, still setting the first group as the reference, then the genetic group estimates for each group are all deviations from the "g1" reference group mean of zero.

**Q** is also needed to calculate breeding values (**a**) from an animal model that has fitted fixed genetic group effects implicitly within the random effects syntax of the model software. This is because such a model will yield predictions of total additive genetic effects (**u**). From equation (6) in the main text, an individual's breeding value (a<sub>i</sub>) equals the total additive genetic effect (u<sub>i</sub>) minus a weighted sum of genetic group effects contributing to that individual. The weighted sum of genetic group contributions for all individuals can be calculated as the matrix product **Qg**, where **g** is a vector of all genetic group effects and **Q** is calculated above.

<!---------     ------------    -----------   ------------->
6.3 Constructing A\*
--------------------

A major advance leading to the widespread use of animal models comes from Henderson (1976) and subsequent researchers' (e.g., Quaas 1976, 1995; Meuwissen and Luo 1992) work to directly construct the matrix inverse of the additive genetic relatedness matrix (<b>A<sup>-1</sup></b>). Direct methods for obtaining <b>A<sup>-1</sup></b> are necessary, because <b>A<sup>-1</sup></b> and not **A** is used when solving the equations in the animal model. <b>A<sup>-1</sup></b> has a unique structure that can be formed by tracing the flow of alleles from each individual to its offspring. Thus methods to directly construct <b>A<sup>-1</sup></b> rely on the basic premise that the flow of alleles through a pedigree can be traced because each individual inherits, on average, half of its alleles from each parent with some sampling variation that can be modeled based on the Mendelian sampling variance in the population.

To estimate genetic group effects by implicitly fitting the genetic groups within the random effect syntax of the model software, requires an inverse matrix to model the covariance among individual total additive genetic values and group effects based on alleles shared identical by descent. Conveniently, the matrix inverse accounting for individual-individual, group-group, group-individual, and individual-group sharing of alleles can be constructed directly from the pedigree in a very similar manner to <b>A<sup>-1</sup></b> (see also *Appendix S5*; Westell, Quaas & Van Vleck 1988; Quaas 1988). This matrix, <b>A\*</b> (following notation of Quaas 1988), is a compilation of four sub-matrices that are each a mathematical function of **Q**, <b>A<sup>-1</sup></b>, or both (Fig. 3e in main text).

The <b>A\*</b> for any pedigree with genetic groups can be obtained using the `makeAinv()` function in the `nadiv` package. For example, <b>A\*</b> for the `Q1988sub` pedigree is:

``` r
AstarOut <- makeAinv(Q1988sub, ggroups = 2)
(Astar <- AstarOut$Ainv)
```

    ## 6 x 6 sparse Matrix of class "dgCMatrix"
    ##                                                 
    ## A   1.3333333 -0.6666667  .    . -0.5 -0.1666667
    ## B  -0.6666667  1.8333333  0.5 -1  .   -0.6666667
    ## C   .          0.5000000  1.5 -1 -0.5 -0.5000000
    ## D   .         -1.0000000 -1.0  2  .    .        
    ## g1 -0.5000000  .         -0.5  .  0.5  0.5000000
    ## g2 -0.1666667 -0.6666667 -0.5  .  0.5  0.8333333

There are several points to highlight regarding <b>A\*</b>. First, this matrix can be split into four blocks: (1) the upper-left 4-by-4 (*n*x*n*) block (individual-individual covariance structure), (2) the lower-left 2-by-4 (*r*x*n*) block (group-individual covariance structure), (3) the upper-right 4-by-2 (*n*x*r*) block (individual-group covariance structure), and (4) the lower-right 2-by-2 (*r*x*r*) block (group-group covariance structure). Note, the upper-left block of <b>A\*</b> is the same as the <b>A<sup>-1</sup></b> for the `Q1988` pedigree without implicitly including genetic group effects (see also Fig. 3e in main text). For example, <b>A<sup>-1</sup></b> for the Q1988 pedigree without genetic groups can be obtained and compared to <b>A\*</b> above using the same `nadiv` function without invoking the genetic group algorithm [`makeAinv(..., ggroups = NULL, ...)`, the default value of the `ggroups` argument]:

``` r
makeAinv(Q1988[-c(1:7), c("id", "dam", "sire")])$Ainv
```

    ## 4 x 4 sparse Matrix of class "dgCMatrix"
    ##                                
    ## A  1.3333333 -0.6666667  .    .
    ## B -0.6666667  1.8333333  0.5 -1
    ## C  .          0.5000000  1.5 -1
    ## D  .         -1.0000000 -1.0  2

<b>A\*</b> constructed in `nadiv` can be used in various animal model software programs (Table S6.1) to implicitly include genetic group effects within the random effects structure of the model. Although the details vary among software programs (see tutorials in the subsections to *Appendix S6.4* for detailed instructions in each of the highlighted software programs), the general approach is to associate the inverse covariance matrix <b>A\*</b> with the term in the model representing additive genetic effects; almost exactly the same way as the typical <b>A<sup>-1</sup></b> is associated with additive genetic effects in a basic animal model. However, some potential alterations to <b>A\*</b> should be considered before and during the model fitting process.

<!---------     ------------    -----------   ------------->
### 6.3.1 Genetic groups on top or bottom of A\*?

The entries in <b>A\*</b> corresponding to genetic groups can be placed in the last *r* rows and columns (as above) or the first *r* rows and columns. When setting up <b>A\*</b>, switching between these placements is trivial. The consequences for using <b>A\*</b>, however, can be rather large. Solving an animal model requires solving a large set of equations (the Mixed Model Equations or MME). Often, some of these equations may not be unique or in other words two equations may be exactly the same, very similar, or have a linear dependency (such that knowing the solution to one equation means the solution to the other equation is also known). This is a cause of singularities in the MME coefficient matrix of an animal model.

The MME are solved by either moving from the first set of equations to the last set ('top-down') or from the last set to the first set ('bottom-up'). In either case, if two equations are nearly the same then the model will obtain a solution for the first of the pair, but then encounter a problem with the second. Some software packages might automatically deal with this second, non-unique equation. For example, **ASReml** solves the equations from the bottom up and will delete the second equation (equation on top) to continue to solve the MME (Butler et al. 2009, pp. 134-135; Gilmour et al. 2014, pp. 106-107).

Unfortunately, the algorithms to solve the MME are not this simple in practice. Even though **ASReml** states which way the equations are solved, some options allow the program to re-order the equations to simplify solving the model. This re-ordering is based on the contents of each equation without respect to its categorization in the model. In other words, one genetic group equation might be re-ordered so that it is placed next to an equation for a year random effect level while another genetic group equation in the same model may be placed elsewhere with a random effect of individual additive genetic value. Although one cannot be completely certain where the genetic group equations will end up when the software program optimizes the order of equations to most efficiently solve the model, the ability to easily switch where in <b>A\*</b> the genetic groups are positioned has been a strategy that has sometimes worked well in our own personal experience. Switching the placement of the genetic groups from the top to the bottom (or vice versa) of <b>A\*</b> has sometimes solved convergence issues.

However, the standalone version of **ASReml** has the `!LAST` job qualifier (see Table S6.1 and the **ASReml** tutorial in *Appendix S6.4.4*) to consistently re-order the MME such that the genetic groups equations are solved last. This option assumes that the genetic groups are the first *r* rows and columns of <b>A\*</b>, where *r* is the number of levels following the `!LAST` qualifier, and therefore 'on top'. We are not aware of a similar argument that can be implemented for a model using `asreml` in `R`.

**WOMBAT**, also, has several MME ordering strategies that can be implemented. However, these re-order the entire set of equations and the order of certain levels of random effects cannot be specified (Meyer 2007, section 5.3.1). `MCMCglmm` re-orders the MME according to algorithms associated with the underlying `CSparse` functions (Davis 2006; Hadfield 2010) and also cannot designate certain levels of random effects to a certain location in the re-ordered equations. Therefore, it is not clear how specifying genetic groups at the top or bottom of <b>A\*</b> will affect the convergence of any particular model.

The `nadiv` function `makeAinv()`, enables the user to place genetic groups at the top or bottom of <b>A\*</b>, in an attempt to solve genetic group equations first or last. Therefore, if any singularities in the MME coefficient matrix occur then the user can attempt to re-order the genetic group equations in the hope that this will lead to model convergence. The `gOnTop` argument to the `makeAinv()` function switches between placing the genetic groups on top or on bottom of <b>A\*</b>, with the default to place genetic groups at the bottom of <b>A\*</b> [`makeAinv(..., gOnTop = FALSE, ...)`].

<!---------     ------------    -----------   ------------->
### 6.3.2 Removing singularities and other problems within A\*

Fitting genetic group effects implicitly often causes singularities in the coefficient matrix of the MME, that can sometimes be overcome by slight changes to the strategy for grouping (Schaeffer 1991). Software programs may report this issue as a non-positive definite generalized inverse matrix (often abbreviated GIV, GIN, or ginverse), or alternatively stating that a GIV/GIN/ginverse is negative definite, non-positive semi-definite, or ill-conditioned. For example, **ASReml** often deals with these types of singularities according to its own strategies (Butler et al. 2009; Gilmour et al. 2014). If genetic groups are specified in an **ASReml** analysis and a singularity occurs the program will introduce a row to the matrix to overcome this (see *Lagrangian multipliers* in Gilmour et al. 2014, pp. 164-165). Regardless of which software program is used, however, adding values to some of the matrix elements may also remove singularities occurring within the MME in addition to the aforementioned strategy of ordering <b>A\*</b> with genetic groups either on top or at the bottom.

Adding a small number to some diagonals can remove singularities (Schaeffer 1994, 1999; Oikawa and Yasuda 2009; Gilmour 2010). Schaeffer (1994, 1999) advocates adding an identity matrix to the *r*x*r* group-group block of <b>A\*</b> such that all of the diagonals (d<sub>ii</sub>) become 1+d<sub>ii</sub>. This addition to <b>A\*</b> means that the genetic group effects are no longer, conceptually speaking, fixed effects but are random effects with expectations of zero and a variance describing their distribution (see *Appendices S1 & S5*). Further, predicted genetic group effects will be biased toward the expected value of zero when genetic groups are random effects, just like any other random effects, with the amount of bias depending on how much information the data provide to predict genetic group effects (Hadfield et al. 2010). This bias will extend to any functions involving the predicted genetic group effects, specifically the predicted total additive genetic effects (**u**). Although, in many cases the mixed model equations and subsequent variance component estimates change only slightly, thus basic interpretations of results will remain the same (but see cautionary note below; further discussion in Schaeffer 1994). The strategy of adding an identity matrix is desirable if there are many individuals assigned to each group, because it does not re-rank the predicted breeding values and can improve convergence (Schaeffer 1994). The main change to the model when the leading diagonal elements of the group-group portion of <b>A\*</b> are altered is a (co)variance matrix for traits within groups is added to the model equations dealing with the genetic group effects. In the simple case of a single trait, this variance is assumed to be the same across group effects with no covariances between groups (i.e., groups effects are normally distributed with a variance of **I** sigma<sup>2</sup><sub>g</sub>). In practice this variance is often assumed to equal the additive genetic variance for the individual total additive genetic effects (Sullivan 1999, see also **WOMBAT**'s genetic groups in *Appendix S6.4.5.4*).

A further interpretation of this assumption means that a model's estimate of additive genetic variance will also include the magnitudes of the genetic group effects. This may cause problems if the variance in the genetic group effects themselves is much greater than the true additive genetic variance. In such a case, adding a smaller value than 1 to the diagonals (d<sub>ii</sub>) of <b>A\*</b> structures the model to reflect this difference. Using values other than 1 for the addition to the diagonal elements effectively means that the ratio of the variance among genetic group effects to the additive genetic variance is no longer constrained to be 1:1. Thus, when modelling genetic groups as random effects we strongly caution users to carefully consider the sensitivity of the additive genetic variance estimate to the choice of value added to the group-group diagonal elements of <b>A\*</b>. Multiple models fitted with different values added to the diagonal elements are recommended to ensure an unbiased estimate of the additive genetic variance. Schaeffer (1994) provides a detailed discussion and interpretation of different strategies for altering <b>A\*</b>.

Note that **WOMBAT** fits a separate variance component to the genetic group effects (see *Appendix S6.4.5.4*), thus offering a direct way to check that the variance in genetic group effects can be assumed to be the same as the additive genetic variance. **ASReml** provides the `!GOFFSET` pedigree file qualifier to add a value to the group-group diagonal elements when using **ASReml**'s `!GROUPS` pedigree qualifier (see *Appendix S6.4.4.4*). Otherwise, this alteration of <b>A\*</b> is easily accomplished manually after creating <b>A\*</b> with the `makeAinv()` function in the `nadiv` package. Altering the <b>A\*</b> for the `Q1988sub` pedigree above (with genetic groups on the bottom) would proceed as follows and result in the matrix displayed below:

``` r
n <- 4
r <- 2
(ggIdentity <- bdiag(Diagonal(n = n, x = 0), Diagonal(n = r, x = 1)))
```

    ## 6 x 6 sparse Matrix of class "dtCMatrix"
    ##                 
    ## [1,] 0 . . . . .
    ## [2,] . 0 . . . .
    ## [3,] . . 0 . . .
    ## [4,] . . . 0 . .
    ## [5,] . . . . 1 .
    ## [6,] . . . . . 1

``` r
(AstarI <- Astar + ggIdentity)
```

    ## 6 x 6 sparse Matrix of class "dgCMatrix"
    ##                                                 
    ## A   1.3333333 -0.6666667  .    . -0.5 -0.1666667
    ## B  -0.6666667  1.8333333  0.5 -1  .   -0.6666667
    ## C   .          0.5000000  1.5 -1 -0.5 -0.5000000
    ## D   .         -1.0000000 -1.0  2  .    .        
    ## g1 -0.5000000  .         -0.5  .  1.5  0.5000000
    ## g2 -0.1666667 -0.6666667 -0.5  .  0.5  1.8333333

Caution is necessary, however, when dealing with singularities in an animal model. Specifically, model comparison may not be possible if singularities occur. **ASReml** alters the order of equations to remove singularities. Consequently, different terms can cause singularities between runs of the model, where the model may have been slightly modified for the purposes of testing the significance of variance components or computing profile likelihoods of variance components. Singularities caused by different terms will affect the log-likelihood of the model and therefore procedures using likelihood ratio test statistics may no longer be valid (Gilmour 2005; Butler et al. 2009; Gilmour et al. 2014)

<!------------------------------------------------------------->
<!---------     ------------    -----------   ------------->
6.4 How to fit genetic groups in MCMCglmm, ASReml-R, ASReml-standalone, & WOMBAT
--------------------------------------------------------------------------------

<!---------     ------------    -----------   ------------->
### 6.4.1 Tutorial dataset `ggTutorial`

We provide a step-by-step demonstration of how to fit genetic groups in animal models, by modelling genetic group effects either explicitly as separate fixed regressions or implicitly within the random effects, using the different animal model software programs. All demonstrations use the dataset `ggTutorial`, which is available in the `nadiv` package for working in `R`, as supplementary file `ggTutorial.dat` for working in the standalone **ASReml** program, and as supplementary file `ggTutorial.d` for working in the **WOMBAT** program. This dataset was simulated using the `nadiv` package (details in *Appendix S4*) and are the same data as those depicted in Figs 2 and 4 from the main text. The dataset contains 6000 individuals across 15 generations. The simulation specifies a carrying capacity of 400 individuals per generation, 200 mating pairs per generation, and 40 immigrants per generation. The expected difference between founder and immigrant mean breeding values equals `3` and both the environmental and additive genetic variances for both the founder and immigrant groups equal `1`. The data are loaded into `R` by loading the `nadiv` package.

``` r
library(nadiv)
```

The first three columns of `ggTutorial` contain the complete pedigree and all 6000 individuals have a phenotypic record in `p`.

``` r
head(ggTutorial)
```

    ##   id dam sire parAvgU mendel          u           r        p is gen
    ## 1  1  NA   NA      NA     NA -0.3928101 -0.08355464 19.52364  0   1
    ## 2  2  NA   NA      NA     NA -0.6750007  0.14201698 19.46702  0   1
    ## 3  3  NA   NA      NA     NA -0.4416978  1.89195823 21.45026  0   1
    ## 4  4  NA   NA      NA     NA -1.0299707 -1.29820650 17.67182  0   1
    ## 5  5  NA   NA      NA     NA -1.9126222  0.46681834 18.55420  0   1
    ## 6  6  NA   NA      NA     NA -0.9089290 -0.82577884 18.26529  0   1

``` r
str(ggTutorial)
```

    ## 'data.frame':    6000 obs. of  10 variables:
    ##  $ id     : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ dam    : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ sire   : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ parAvgU: num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ mendel : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ u      : num  -0.393 -0.675 -0.442 -1.03 -1.913 ...
    ##  $ r      : num  -0.0836 0.142 1.892 -1.2982 0.4668 ...
    ##  $ p      : num  19.5 19.5 21.5 17.7 18.6 ...
    ##  $ is     : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ gen    : int  1 1 1 1 1 1 1 1 1 1 ...

The rest of the columns contain:

-   `parAvgU` is the average *u* (total additive genetic effect) of each individual's parents

-   `mendel` is the Mendelian sampling deviation from the mid-parent value that is unique to each individual

-   `u` is the total additive genetic effect of each individual's genotype

-   `r` is the residual deviation of each individual

-   `p` is the phenotype for each individual (population expected mean is 20)

-   `is` indicates the immigration status, where focal population residents have `is == 0` and immigrants have `is == 1`

-   `gen` is the generation in which each individual was born

The supplementary files for use in **ASReml** and **WOMBAT** contain a subset of these variables plus three derived variables. The **ASReml** and **WOMBAT** data files contain continuous variables for the coefficient of inbreeding (*f*) and two genetic group genomic contributions (from **Q**). For demonstration purposes, these are added to the basic `ggTutorial` below within each section. Note the **WOMBAT** pedigree and data files have had 10000 added to the individual identity codes to comply with **WOMBAT** data structure requirements.

It is necessary to fit a fixed regression on the coefficient of inbreeding (*f*) to account for inbreeding depression and minimize bias in estimates of quantitative genetic parameters (Kennedy, Schaeffer & Sorensen 1988; Reid & Keller 2010; Wolak & Keller 2014). This is not directly part of the process to fit genetic groups in the animal model, but it is considered necessary when fitting animal models to population data in which inbreeding occurs.

------------------------------------------------------------------------

<!------------------------------------------------------------->
<!---------     ------------    -----------   ------------->
### 6.4.2 `MCMCglmm`

Modelling genetic groups either explicitly as separate fixed regressions or implicitly within the random effects syntax of `MCMCglmm` requires either **Q** or <b>A\*</b> to be created using the `nadiv` package. Below we demonstrate how to create these matrices and fit each in a `MCMCglmm` model.

``` r
library(MCMCglmm)
library(nadiv)
```

To complete the basic preparations of the data for an animal model, all that is necessary is to calculate the coefficients of inbreeding (*f*; to minimize bias from inbreeding depression, see the last paragraph of *Appendix 6.4.1* above). This calculation can be done using the `inverseA()` function of `MCMCglmm`. Consequently, the <b>A<sup>-1</sup></b> matrix is also created, which we will assign as its own object in `R`, because it will be needed in the model fitting genetic group effects explicitly as fixed covariate regressions with **Q** (see *Appendix 6.4.2.2*).

``` r
ainvOut <- inverseA(ggTutorial[, 1:3])
Ainv <- ainvOut$Ainv
ggTutorial$f <- ainvOut$inbreeding
```

`MCMCglmm` facilitates modelling of non-Gaussian response variables, which has greatly extended the application of animal models in evolutionary ecology to a range of non-Gaussian phenotypes. Fitting genetic groups is no different when modelling non-Gaussian response variables. In such models, the genetic group effects (fitted either explicitly as separate fixed regressions or implicitly within the random effects) are estimated on the underlying latent scale.

#### 6.4.2.1 Preparing a pedigree with genetic groups

Here, we provide general instructions to prepare pedigrees. Calculating the matrix of genetic group contributions (**Q**) and the augmented inverse relatedness matrix (<b>A\*</b>) requires a pedigree that indicates groups instead of missing values for parents. In the `ggTutorial` data, we will assign all individuals with unknown parents in the first generation phantom parents from one genetic group (`foc0`). All individuals in subsequent generations that have unknown parents will be assigned phantom parents from a second genetic group (`g1`). Alternatively, this second group could be further divided into genetic groups based on generation. However, the data were simulated such that every immigrant has the same expected mean breeding value, regardless of the generation in which the immigrant was born. Therefore, we will stick to a total of two genetic groups (but, feel free to model more groups!).

First, create an object that will be the pedigree with genetic groups (`ggPed`).

``` r
ggPed <- ggTutorial[, c("id", "dam", "sire", "is", "gen")]
naPar <- which(is.na(ggPed[, 2]))
ggPed$GG <- rep("NA", nrow(ggPed))
  # 'focal' genetic group = "foc0" and 'immigrant' = "g1"
  # obtained by pasting "foc" & "g" with immigrant status "0" or "1", respectively
  ggPed$GG[naPar] <- as.character(ggPed$is[naPar])
  ggPed$GG[ggPed$GG == "0"] <- paste0("foc", ggPed$GG[ggPed$GG == "0"])
  ggPed$GG[ggPed$GG == "1"] <- paste0("g", ggPed$GG[ggPed$GG == "1"])
ggPed[naPar, 2:3] <- ggPed[naPar, "GG"]
```

Note, the approach and format of the pedigree above is different from the `Q1988sub` pedigree earlier - mostly because the format here makes it easier to specify genetic groups in a pedigree for this particular case, but partly to illustrate the flexibility of `nadiv` functions. To be more specific, the `Q1988sub` pedigree contained two extra rows for the two genetic groups. The `ggroups` argument to `makeAinv()` supplied an integer indicating how many rows at the begining of the pedigree contained genetic groups and not individuals. However, in the `ggPed` above no extra rows are added to the pedigree, but in the next few sections the character vector given to the `ggroups` argument in both `ggcontrib()` and `makeAinv()` specifies the name of the unique genetic groups that have been filled in for individuals instead of missing dam and sire identities.

<!---------     ------------    -----------   ------------->
#### 6.4.2.2 Fixed explicit genetic group effects with Q (from `nadiv`)

Fitting genetic group effects explicitly requires the columns of **Q** to be included as separate fixed covariate regressions. The code below creates **Q** for the `ggPed` pedigree and adds the columns (`foc0` and `g1`) as variables in `ggTutorial` so that they can be included in a model.

``` r
Q <- ggcontrib(ggPed[, 1:3], ggroups = c("foc0", "g1"))
ggTutorial <- cbind(ggTutorial, Q)
str(ggTutorial)
```

    ## 'data.frame':    6000 obs. of  13 variables:
    ##  $ id     : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ dam    : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ sire   : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ parAvgU: num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ mendel : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ u      : num  -0.393 -0.675 -0.442 -1.03 -1.913 ...
    ##  $ r      : num  -0.0836 0.142 1.892 -1.2982 0.4668 ...
    ##  $ p      : num  19.5 19.5 21.5 17.7 18.6 ...
    ##  $ is     : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ gen    : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ f      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ foc0   : num  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ g1     : num  0 0 0 0 0 0 0 0 0 0 ...

An important point to make is that, because the pedigree and data for `ggTutorial` are in the same dataframe, the genetic group coefficients from **Q** can be added directly to the data with `ggTutorial <- cbind(ggTutorial, Q)` above. However, whenever a pedigree contains more individuals than the data (or in a different order), an intermediate step is required. Specifically, since the **Q** matrix contains a row for every individual in the pedigree, only the subset of rows in **Q** that correspond to phenotyped individuals in a dataset should be extracted and ordered according to the identities in the data.

The fixed regressions on `foc0` and `g1` can be combined with the standard <b>A<sup>-1</sup></b> to specify an animal model in `MCMCglmm`. Below, we write out the specification of the `MCMCglmm` default priors for the fixed effects in the model. We recommend a normal prior distribution for each genetic group fixed regression coefficient with a mean of 0 and a very large variance. Our fixed effect prior specification below has an element in `mu` and `V` of `prPE$B` to correspond with the model's: intercept, regression on the coefficient of inbreeding (*f*), and regressions on the genetic group contributions (`foc0` and `g1`). Because unique solutions to the model do not exist (see discussion of estimability in *Appendix S6.2*) for the intercept, coefficient of inbreeding (*f*), and genetic group effects (`foc0` and `g1`), we fit the `foc0` genetic group as the last fixed effect. This ensures that the model sets `foc0` as a reference group and estimates the `g0` genetic group effect as a deviation from the mean breeding value in the `foc0` group. The model is written out below, and a saved version is available on Dryad (Wolak & Reid 2016) as supporting information file "./MCMCglmm/ggRegMC.RData" so that it is not necessary to re-run this model to examine the results.

``` r
# Use diffuse normal priors (MCMCglmm defaults) for all fixed effects
# Parameter expanded prior on additive genetic variance
## Scaled F-dsitribution with numerator and denominator df=1
prPEexp <- list(B = list(mu = rep(0, 3), V = diag(3)*1e10), 
         R = list(V = 1, nu = 0.002),
         G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

# NOTE: the saved model could be loaded instead of running it below
## load(file = "./MCMCglmm/ggRegMC.RData") 
ggRegMC <- MCMCglmm(fixed = p ~ f + g1 + foc0,
    random = ~ id,
    ginverse = list(id = Ainv),
    data = ggTutorial,
    prior = prPEexp,
    singular.ok = FALSE,
    pr = TRUE,
    nitt = 28000, burn = 3000, thin = 25)
```

The posterior distributions of genetic group effects are included in the `Sol` object along with the posterior distributions of the other fixed effects. Note the software has dropped the `foc0` group estimate in this model, implying this is the reference group and the `g1` genetic group effects is the estimated deviation from the reference. To report the posterior modes and 95% limits of the highest posterior density:

``` r
with(ggRegMC, cbind(postMode = posterior.mode(Sol[, 1:3]),
    HPDinterval(Sol[, 1:3])))
```

    ##               postMode     lower     upper
    ## (Intercept) 19.7952878 19.664293 19.931933
    ## f            0.3756279 -2.652013  4.529125
    ## g1           3.1375194  2.941707  3.285130

The above posterior mode estimate of 3.138 agrees with the expected difference of `3` between the mean breeding values of the immigrant group (`3`) and the founder population (`0`), specified in the simulation.

The posterior distribution of the breeding values (**a**) can be accessed from `ggRegMC$Sol[, -c(1:3)]`. Because the order of rows in **Q** matches the order of individual identities in `ggTutorial`, the posterior distribution for the total additive genetic effects (**u**) can be calculated as (eqn 5 in main text):

``` r
# Add 0 genetic group effect for the reference group 'foc0' to the 'g1' samples
regPost_gHat <- rbind(0, ggRegMC$Sol[, "g1"])
# For the ith mcmc iteration that was saved, calculate 'u'
regPost_u <- as.mcmc(t(sapply(seq(nrow(ggRegMC$Sol)),
    FUN = function(i){Q %*% regPost_gHat[, i] + ggRegMC$Sol[i, -c(1:3)]})))
```

<!---------     ------------    -----------   ------------->
#### 6.4.2.3 Random implicit genetic group effects with A\* (from `nadiv`)

Fitting genetic group effects implicitly within the random effects portion of the `MCMCglmm` model specification requires the augmented inverse **A** matrix (<b>A\*</b>) to be supplied as a sparse matrix in the `ginverse` argument. The code below creates <b>A\*</b> for the `ggPed` pedigree using the `makeAinv` function in the `nadiv` package.

``` r
Astar <- makeAinv(ggPed[, 1:3], ggroups = c("foc0", "g1"), gOnTop = TRUE)$Ainv
```

In the above code, we directed `makeAinv()` to include the genetic groups on top (top-left block) of <b>A\*</b> using the `gOnTop = TRUE` argument. At this point we could try fitting the model with this version of <b>A\*</b> (i.e., `Astar`). However, the model returns an error that `G-structure 1 is ill-conditioned`. In other words, <b>A\*</b> has at least one eigenvalue that is either negative or zero. Re-defining the groupings can eliminate such a singularity. However, removing the singularity in the matrix while maintaining the current groupings is accomplished by adding a small value to the diagonal elements of the group-group portion of the matrix. Note that such an addition changes the interpretation of the effects in the model, as discussed above (*Appendix 6.3.2* as well as *Appendix S5*), now implying that the genetic group effects are random effects in the model. Below, we add 0.1 to the diagonal element of the <b>A\*</b> which will constrain the variance in genetic group effects to be one tenth the additive genetic variance estimated by the model. A value of 0.1 may bias estimates of additive genetic variance in models of other datasets. Therefore, we advise testing values other than 0.1 (e.g., 1 or 10) and comparing estimates from these alternative models. Note that a value of 1 constrains the variance in genetic group effects to be equal to the additive genetic variance estimated by the model.

``` r
AstarAdd <- Astar
(ggrows <- match(c("foc0", "g1"), dimnames(AstarAdd)[[1]]))
```

    ## [1] 1 2

``` r
diag(AstarAdd)[ggrows] <- diag(AstarAdd)[ggrows] + 0.1
Astar[1:3, 1:3]
```

    ## 3 x 3 sparse Matrix of class "dgCMatrix"
    ##                
    ## foc0 400   . -1
    ## g1     . 520  .
    ## 1     -1   .  1

``` r
AstarAdd[1:3, 1:3]
```

    ## 3 x 3 sparse Matrix of class "dgCMatrix"
    ##                    
    ## foc0 400.1   .   -1
    ## g1     .   520.1  .
    ## 1     -1.0   .    1

The `AstarAdd` sparse inverse matrix can now be associated with the individual identities (`id`) by supplying it as an argument to `ginverse` in the `MCMCglmm` model. We specify the same default prior distribution for the fixed effects below (i.e., the same ones as we wrote out in *Appendix 6.4.2.2*). The prior distribution for the `id` term in the model is the prior used for the additive genetic variance and, in this case, also the same prior for the variance of the genetic group effects. It is unclear how the prior specified for the additive genetic variance will affect the estimate of the genetic group effects or even what prior is to be used in this context. As usual with Bayesian analyses, we recommend priors that are framed within the context of the data at hand, the model being fitted, and any *a priori* belief in the model parameters. Further we suggest testing the sensitivity of any particular prior used for a given model and especially when genetic group effects are treated as random effects by the model.

The model is written out below, but a saved version is available on Dryad (Wolak & Reid 2016) as supporting information file "./MCMCglmm/ggAstarRandMC.RData" so that it is not necessary to re-run this model to examine the results.

``` r
# Parameter expanded prior on additive genetic variance
prPEimp <- list(B = list(mu = rep(0, 2), V = diag(2)*1e10), 
         R = list(V = 1, nu = 0.002),
         G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

# NOTE: the saved model could be loaded instead of running it below
## load(file = "./MCMCglmm/ggAstarRandMC.RData") 
ggAstarRandMC <- MCMCglmm(fixed = p ~ f,
    random = ~ id,
    ginverse = list(id = AstarAdd),
    data = ggTutorial,
    prior = prPEimp,
    pr = TRUE,
    nitt = 53000, burn = 3000, thin = 50)
```

Alternatively, the non-altered matrix `Astar` could be supplied to the `ginverse` argument of `MCMCglmm` to see the types of warnings that led us to fit the model with the matrix `AstarAdd`.

The posterior distributions of genetic group effects are included in the `Sol` object (as long as `pr = TRUE` in the model statement) along with the posterior distributions of the fixed effects. To report the posterior modes and 95% limits of the highest posterior density:

``` r
with(ggAstarRandMC, cbind(postMode = posterior.mode(Sol[, 1:4]),
    HPDinterval(Sol[, 1:4])))
```

    ##              postMode     lower     upper
    ## (Intercept) 22.023119 16.604035 25.597429
    ## f            1.348870 -2.842843  4.164732
    ## id.foc0     -2.164536 -5.801390  3.099577
    ## id.g1        0.786011 -2.742551  6.267101

Note that all fixed effects in this model are estimable and the model returns these estimates along with predicting both of the genetic group effects. The genetic group predictions, themselves, are random effects and thus are deviations from their expected value of 0. The posterior mode of differences between genetic group predictions 3.107 (95% HPD interval:2.945, 3.277) agrees with the expected difference between the mean breeding values of the immigrant group (`3`) and the founder population (`0`) that was specified in the simulation.

The posterior distribution of the total additive genetic effects (**u**) can be accessed from `ggAstarRandMC$Sol[, -c(1:4)]`. Because the order of rows in **Q** matches the order of individual identities in `ggTutorial`, the posterior distribution for the breeding values (**a**) can be calculated as (eqn 6 in main text):

``` r
AstarPost_gHat <- ggAstarRandMC$Sol[, 3:4]
# For the ith mcmc iteration that was saved, calculate 'a'
# Assuming Q has already been calculated as above
AstarPost_a <- as.mcmc(t(sapply(seq(nrow(ggAstarRandMC$Sol)),
    FUN = function(i){ggAstarRandMC$Sol[i, -c(1:4)] -
            Q %*% matrix(AstarPost_gHat[i, ], ncol = 1)})))
```

------------------------------------------------------------------------

<!------------------------------------------------------------->
<!---------     ------------    -----------   ------------->
### 6.4.3 `ASReml` in `R`

To model genetic group effects explicitly as separate fixed regressions in an animal model using the `R` package `asreml` requires **Q** to be created using the `nadiv` package and then included in the `asreml()` model function. Including genetic group effects implicitly within the random effects of an animal model in `asreml()` can be done either using <b>A\*</b> constructed with `nadiv` or with `asreml::asreml.Ainverse()`. All three of these approaches will be demonstrated below using the simulated dataset `ggTutorial`:

``` r
library(asreml)
library(nadiv)
```

To complete the basic preparations of the data for an animal model, all that is necessary is to calculate the coefficients of inbreeding (*f*; to minimize bias from inbreeding depression, see the last paragraph of *Appendix 6.4.1*). This calculation can be done using the `asreml.Ainverse()` function of `asreml`. Consequently, the <b>A<sup>-1</sup></b> matrix is also created, which we will assign as its own `R` object, because it will be needed in the model fitting genetic group effects explicitly as separate fixed covariate regressions with **Q** below (*Appendix 6.4.3.2*).

``` r
ainvOut <- asreml.Ainverse(ggTutorial[, 1:3])
listAinv <- ainvOut$ginv
ggTutorial$f <- ainvOut$inbreeding
## Have to also convert the 'id' variable to a factor for models below
ggTutorial$id <- as.factor(ggTutorial$id)
```

#### 6.4.3.1 Preparing a pedigree with genetic groups

Calculating the matrix of genetic group contributions (**Q**) and the augmented inverse relatedness matrix (<b>A\*</b>) requires a pedigree that has genetic groups indicated instead of missing values for parents. In the `ggTutorial` data, we will assign all individuals with unknown parents in the first generation phantom parents from one genetic group (`foc0`). All individuals in subsequent generations that have unknown parents will be assigned phantom parents from a second genetic group (`g1`). Alternatively, this second group could be further divided into genetic groups based on generation. However, the data were simulated such that every immigrant has the same expected mean breeding value, regardless of the generation in which it was born. Therefore, we will stick to a total of two genetic groups (but, feel free to model more groups!).

First, create an object that will be the pedigree with genetic groups (`ggPed`).

``` r
ggPed <- ggTutorial[, c("id", "dam", "sire", "is", "gen")]
naPar <- which(is.na(ggPed[, 2]))
ggPed$GG <- rep("NA", nrow(ggPed))
  # 'focal' genetic group = "foc0" and 'immigrant' = "g1"
  # obtained by pasting "foc" & "g" with immigrant status "0" or "1", respectively
  ggPed$GG[naPar] <- as.character(ggPed$is[naPar])
  ggPed$GG[ggPed$GG == "0"] <- paste0("foc", ggPed$GG[ggPed$GG == "0"])
  ggPed$GG[ggPed$GG == "1"] <- paste0("g", ggPed$GG[ggPed$GG == "1"])
ggPed[naPar, 2:3] <- ggPed[naPar, "GG"]
# Two rows need to be added for the genetic groups
## Genetic groups will be given the missing value NA instead of parent identities
ggPed <- data.frame(id = c("foc0", "g1", as.character(ggPed$id)),
        dam = c(NA, NA, as.character(ggPed$dam)),
        sire = c(NA, NA, as.character(ggPed$sire)))
head(ggPed)
```

    ##     id  dam sire
    ## 1 foc0 <NA> <NA>
    ## 2   g1 <NA> <NA>
    ## 3    1 foc0 foc0
    ## 4    2 foc0 foc0
    ## 5    3 foc0 foc0
    ## 6    4 foc0 foc0

Note that the `nadiv` function `makeAinv()` used to construct <b>A\*</b> is more flexible than `asreml`'s `asreml.Ainverse()` in the way pedigrees containing genetic groups can be formatted (for more detail on `makeAinv()` formats, see the help file by running the command `?makeAinv` in an `R` session - particularly the examples at the end of the file - or the genetic group pedigree used above in the `MCMCglmm` tutorial). To be more specific, the format required by `asreml.Ainverse()` is similar to the `Q1988sub` pedigree used earlier in the tutorial (*Appendix 6.1*) which contains two extra rows for the two genetic groups.

#### 6.4.3.2 Fixed explicit genetic group effects with Q (from `nadiv`)

Fitting genetic group effects explicitly in `asreml` requires the columns of **Q** to be included as separate fixed covariate regressions. The code below creates **Q** for the `ggPed` pedigree and adds the columns (`foc0` and `g1`) as variables in `ggTutorial` so that they can be included in a model.

``` r
Q <- ggcontrib(ggPed)
ggTutorial <- cbind(ggTutorial, Q)
```

``` r
str(ggTutorial)
```

    ## 'data.frame':    6000 obs. of  13 variables:
    ##  $ id     : Factor w/ 6000 levels "1","2","3","4",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ dam    : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ sire   : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ parAvgU: num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ mendel : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ u      : num  -0.393 -0.675 -0.442 -1.03 -1.913 ...
    ##  $ r      : num  -0.0836 0.142 1.892 -1.2982 0.4668 ...
    ##  $ p      : num  19.5 19.5 21.5 17.7 18.6 ...
    ##  $ is     : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ gen    : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ f      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ foc0   : num  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ g1     : num  0 0 0 0 0 0 0 0 0 0 ...

An important point to make is that, because the pedigree and data for `ggTutorial` have essentially the same sizes and structures (only `ggPed` now has two extra rows for the genetic groups - but these will be dropped automatically from **Q** by `ggcontrib()`), the genetic group coefficients from **Q** can be added directly to the data with `ggTutorial <- cbind(ggTutorial, Q)` above. However, whenever a pedigree contains more individuals than the data (or in a different order), an intermediate step is required. Specifically, since the **Q** matrix contains a row for every individual in the pedigree, only the subset of rows in **Q** that correspond to phenotyped individuals in a dataset should be extracted and ordered according to the identities in the data.

The fixed regressions on `foc0` and `g1` can be combined with the standard <b>A<sup>-1</sup></b> to specify an animal model in `asreml()`. Because unique solutions to the model do not exist (see discussion of estimability in *Appendix S6.2*) for the intercept, coefficient of inbreeding (*f*), and genetic group effects (`foc0` and `g1`), we fit the `foc0` genetic group as the last fixed effect. This ensures that the model sets `foc0` as a reference group and estimates the `g0` genetic group effect as a deviation from the mean breeding value in the `foc0` group.

``` r
ggRegASR <- asreml(fixed = p ~ f + g1 + foc0,
    random = ~ ped(id),
    ginverse = list(id = listAinv),
    data = ggTutorial)
```

The estimates of the genetic group effects are included with the fixed effect estimates. Note the software has set the `foc0` group estimate to zero in this model, implying this is the reference group and the `g1` genetic group effect is the estimated deviation. The variance of the estimated genetic group effect (`g1`) is included with the variances of the fixed effect estimates, which are reported on `asreml`'s underlying transformed scale. To re-scale the variances of the fixed effects back to the phenotypic scale and report them as standard errors:

``` r
with(ggRegASR, cbind(Effect = coefficients$fixed,
    seEffect = sqrt(vcoeff$fixed * sigma2)))
```

    ##                 Effect   seEffect
    ## foc0         0.0000000 0.00000000
    ## g1           3.1167796 0.08929429
    ## f            0.6660162 1.83248941
    ## (Intercept) 19.7976160 0.06579030

The above estimate of 3.117 agrees with the expected difference between mean breeding value in the immigrant group (`3`) and the founder population (`0`) that was specified in the simulation.

The predicted breeding values (**a**) can be accessed from `ggRegASR$coefficients$random`. Because the order of rows in **Q** matches the order of individual identities in `ggTutorial`, the total additive genetic effects (**u**) can be calculated without any further manipulation as (eqn 5 in main text):

``` r
reg_gHat <- matrix(ggRegASR$coefficients$fixed[1:2], ncol = 1)
reg_u <- Q %*% reg_gHat + ggRegASR$coefficients$random
```

Note, if more individuals are in the pedigree than the dataset and/or the identities in the pedigree and data are ordered differently, more coding is necessary to ensure that the predicted genetic effect for a particular individual is added to the correct weighted genetic group contribution. Additionally, if more random effects are included in the model, care has to be taken in order to ensure the correct random effect predictions are extracted for the breeding values (**a**). The `R` function `match()` is particularly helpful in these cases.

#### 6.4.3.3 Fixed implicit genetic group effects with A\* (from `nadiv` or `asreml`)

Fitting fixed effects of genetic groups implicity within the random effects of an animal model requires the augmented inverse relatendess matrix (<b>A\*</b>) to be supplied as a list in the `ginverse` argument. First, we create the sparse matrix of <b>A\*</b> for the `ggPed` pedigree using the `makeAinv` function in the `nadiv` package and then demonstrate how to obtain this same matrix with the `asreml` package.

``` r
listAstar <- makeAinv(ggPed, ggroups = 2, gOnTop = FALSE)$listAinv
```

In the above code, we directed `makeAinv()` to include the genetic groups on bottom (bottom-right block) of <b>A\*</b> using `gOnTop = FALSE`. This was done to improve model convergence (see *Appendix 6.3* and its subsections). In this case, the `asreml` model converges. In other pedigrees, models, and/or genetic group classifications, if singularities occur then sometimes this can be overcome by adding one (or a small value) to the diagonal elements of the group-group portion of the matrix (although see discussion of this in *Appendix S5* and cautions/considerations in *Appendix 6.3.2*). The alteration is demonstrated below along with cautions and interpretations (*Appendix 6.4.3.4*). Note, the easiest approach is to adjust the genetic group diagonals of <b>A\*</b> using the matrix returned by `makeAinv()` (i.e., object `Ainv` in the list), then convert the matrix to the list format required by `asreml`.

Alternatively, we can use `asreml`'s function `asreml.Ainverse()` to obtain a list format of <b>A\*</b>.

``` r
listAstar <- asreml.Ainverse(ggPed, groups = 2)$ginv
```

CAUTION: under the current and previous versions of `asreml` (current version at the end of Supporting Information), quite often `asreml.Ainverse()` simply does not work. It will return an object with a `ginv` entry that is a `data.frame`, but with zero rows. Sometimes, re-running the above command will work - otherwise use `nadiv`'s `makeAinv()` (*Appendix 6.4.3*)! Note that in the above code, `asreml.Ainverse()` includes the genetic groups on top (top-left block) of <b>A\*</b>. Often models fitting this structure of <b>A\*</b> have trouble converging due to this matrix not being positive-definite (see *Appendix 6.3* and its subsections). One strategy to remedy this is to place the genetic groups on the bottom of <b>A\*</b>. Since there is no easy option to specify this in `asreml.Ainverse()` we recommend using `nadiv`'s function `makeAinv()` as demonstrated above.

Regardless of which package was used to create `listAstar`, this sparse inverse matrix in list format has two extra properties ('attributes') that are critical for being able to use it in an `asreml` model. First (and not specific to just genetic group models), all `asreml` ginverse objects must have the "rowNames" attribute to map row/column numbers back to the factor levels in the model variable. Second, specific to models fitting fixed genetic groups implicitly within the random effects the `listAstar` object must have the "geneticGroups" attribute where the first number indicates the number of unique genetic groups. This is essential for `asreml` to correctly calculate the residual degrees of freedom as part of the residual variance and log-likelihood calculations. These attributes are the last two rows of the internal object structure:

``` r
str(listAstar)
```

    ## 'data.frame':    19642 obs. of  3 variables:
    ##  $ row   : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ column: num  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ Ainv  : num  1 1.5 1 1 1.5 1 1 1 1 2 ...
    ##  - attr(*, "rowNames")= chr  "1" "2" "3" "4" ...
    ##  - attr(*, "geneticGroups")= num  2 0

Both of these attributes are discussed in the help file to `makeAinv()` and can be accessed by typing `?makeAinv` in `R`. Now the object `listAstar` can be associated with the individual identities (`id`) in the model by supplying the list to the `ginverse` argument in the `asreml()` function.

``` r
ggAstarASR <- asreml(fixed = p ~ f,
    random = ~ ped(id),
    ginverse = list(id = listAstar),
    data = ggTutorial)
```

The estimate of the genetic group effects are included with the random effect predictions. Since we directed `makeAinv()` to position the genetic groups on the bottom of `Astar`, the genetic group effects are the last two predicted values. The variance of the estimated genetic group effects are included with the variances of the random effect predictions, which are reported on `asreml`'s underlying transformed scale. To re-scale the variances of the fixed effects and genetic group effects back to the phenotypic scale and report these as standard errors along with the estimates themselves:

``` r
# find the location of the genetic groups
(ggrows <- match(c("foc0", "g1"), attr(listAstar, "rowNames")))
```

    ## [1] 6001 6002

``` r
# genetic group and fixed effects (with re-scaled standard errors)
with(ggAstarASR, cbind(Effect = c(coefficients$fixed, coefficients$random[ggrows]),
    seEffect = sqrt(c(vcoeff$fixed, vcoeff$random[ggrows]) * sigma2)))
```

    ##                  Effect   seEffect
    ## f             0.6660162 1.83248941
    ## (Intercept)   0.0000000 0.00000000
    ## ped(id)_foc0 19.7976160 0.06579030
    ## ped(id)_g1   22.9143956 0.06008392

Note the software has set the intercept estimate to zero in this model (see discussion of estimability in *Appendix S6.2*), implying this is the reference to which the `foc0` and `g1` genetic group effects are the estimated deviations. Therefore, the predicted total additive genetic effects (**u**) include the overall phenotypic mean. However, the difference between the genetic group effects 3.117 agrees with the expected difference between the mean breeding values of the immigrant group (`3`) and the founder population (`0`) that was specified in the simulation.

The predicted total additive genetic effects (**u**) can be accessed from `ggAstarASR$coefficients$random[-ggrows]`. Because the order of rows in **Q** matches the order of individual identities in `ggTutorial`, the breeding values (**a**) can be calculated without any further manipulation as (eqn 6 in main text):

``` r
Astar_gHat <- matrix(ggAstarASR$coefficients$random[ggrows], ncol = 1)
Astar_a <- ggAstarASR$coefficients$random[-ggrows] - Q %*% Astar_gHat
```

Note, if more individuals are in the pedigree than the dataset and/or the identities in the pedigree and data are ordered differently, more coding is necessary to ensure that the predicted genetic effect for a particular individual is added to the correct weighted genetic group contribution. Additionally, if more random effects are included in the model, care has to be taken in order to ensure the correct random effect predictions are extracted for the total additive genetic effects (**u**) and used to calculate breeding values (**a**). The `R` function `match()` is particularly helpful in these cases.

#### 6.4.3.4 Random implicit genetic group effects with A\* (from `nadiv`)

Fitting random effects of genetic groups implicity within the random effects of an animal model requires the augmented inverse relatendess matrix (<b>A\*</b>) to be supplied as a list in the `ginverse` argument. First, we create the sparse matrix of <b>A\*</b> for the `ggPed` pedigree using the `makeAinv` function in the `nadiv` package, then make the necessary alterations (see discussion of this in *Appendix S5* and cautions/considerations in *Appendix 6.3.2*) so the model interprets the genetic groups as random effects, convert the sparse matrix to a list, and fit the model.

``` r
Astar <- makeAinv(ggPed, ggroups = 2, gOnTop = FALSE)$Ainv
```

In the above code, we directed `makeAinv()` to include the genetic groups on bottom (bottom-right block) of <b>A\*</b> using `gOnTop = FALSE`. This was done to improve model convergence (see *Appendix 6.3* and its subsections).

Next we add a small value to the diagonal elements of the group-group portion of the matrix taking the same approach as in the `MCMCglmm` tutorial above (*Appendix 6.4.2.3*). As a result of this alteration, the genetic group effects are conceptually random effects in the model. Below, we add 0.1 to the diagonal element of the <b>A\*</b> which will constrain the variance in genetic group effects to be one tenth the additive genetic variance estimated by the model. A value of 0.1 may bias estimates of additive genetic variance in models of other datasets. Therefore, we advise testing values other than 0.1 (e.g., 1 or 10) and comparing estimates from these alternative models. Note that a value of 1 constrains the variance in genetic group effects to be equal to the additive genetic variance estimated by the model.

The easiest approach is to adjust the genetic group diagonals of <b>A\*</b>, and this is why we stored the matrix returned by `makeAinv()` above (i.e., object `Ainv` in the list). Once the additions are made, we can then convert the matrix to the list format required by `asreml`.

``` r
AstarAdd <- Astar
(ggrows <- match(c("foc0", "g1"), dimnames(AstarAdd)[[1]]))
```

    ## [1] 6001 6002

``` r
diag(AstarAdd)[ggrows] <- diag(AstarAdd)[ggrows] + 0.1
Astar[ggrows, ggrows]
```

    ## 2 x 2 sparse Matrix of class "dgCMatrix"
    ##             
    ## foc0 400   .
    ## g1     . 520

``` r
AstarAdd[ggrows, ggrows]
```

    ## 2 x 2 sparse Matrix of class "dgCMatrix"
    ##                 
    ## foc0 400.1   .  
    ## g1     .   520.1

``` r
# convert to the list format using `nadiv::sm2list()`
listAstarAdd <- sm2list(AstarAdd, rownames = rownames(AstarAdd))
# Should be NO "geneticGroups" attribute (or element 1 of this attribute must==0)
str(listAstarAdd)
```

    ## 'data.frame':    19642 obs. of  3 variables:
    ##  $ row   : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ column: num  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ A     : num  1 1.5 1 1 1.5 1 1 1 1 2 ...
    ##  - attr(*, "rowNames")= chr  "1" "2" "3" "4" ...

The object `listAstarAdd` can be associated with the individual identities (`id`) in the model by supplying the list to the `ginverse` argument in the `asreml()` function.

``` r
ggAstarRanASR <- asreml(fixed = p ~ f,
    random = ~ ped(id),
    ginverse = list(id = listAstarAdd),
    data = ggTutorial)
```

The predictions of the genetic group effects are included with the random effect predictions. Since we directed `makeAinv()` to position the genetic groups on the bottom of `Astar`, the genetic group effects are the last two predicted values. The variance of the random effect predictions are reported on `asreml`'s underlying transformed scale. To re-scale the variances of the genetic group predictions back to the phenotypic scale and report these as standard errors along with the predictions themselves and the fixed effect estimates (for comparison to other models):

``` r
# find the location of the genetic groups
(ggrows <- match(c("foc0", "g1"), attr(listAstarAdd, "rowNames")))
```

    ## [1] 6001 6002

``` r
# genetic group and fixed effects (with re-scaled standard errors)
with(ggAstarRanASR, cbind(Effect = c(coefficients$fixed, coefficients$random[ggrows]),
    seEffect = sqrt(c(vcoeff$fixed, vcoeff$random[ggrows]) * sigma2)))
```

    ##                  Effect seEffect
    ## f             0.6656207 1.832498
    ## (Intercept)  21.3560577 2.320737
    ## ped(id)_foc0 -1.5578210 2.320741
    ## ped(id)_g1    1.5578209 2.320741

Note the intercept in this model is estimable and so the model returns this estimate along with predicting both of the genetic group effects. The genetic group predictions, themselves, are random effects and thus are deviations from their expected value of 0. The difference between the genetic groups 3.116 agrees with the expected difference between the mean breeding values of the immigrant group (`3`) and the founder population (`0`) that was specified in the simulation.

The predicted total additive genetic effects (**u**) can be accessed from `ggAstarRanASR$coefficients$random[-ggrows]`. Here, the prediced total additive genetic effects include the genetic group effects expressed as deviations from 0. Because the order of rows in **Q** matches the order of individual identities in `ggTutorial`, the breeding values (**a**) can be calculated without any further manipulation as (eqn 6 in main text):

``` r
AstarRan_gHat <- matrix(ggAstarRanASR$coefficients$random[ggrows], ncol = 1)
AstarRan_a <- ggAstarRanASR$coefficients$random[-ggrows] - Q %*% AstarRan_gHat
```

Note, if more individuals are in the pedigree than the dataset and/or the identities in the pedigree and data are ordered differently, more coding is necessary to ensure that the predicted genetic effect for a particular individual is added to the correct weighted genetic group contribution. Additionally, if more random effects are included in the model, care has to be taken in order to ensure the correct random effect predictions are extracted for the total additive genetic effects (**u**) and used to calculate breeding values (**a**). The `R` function `match()` is particularly helpful in these cases.

------------------------------------------------------------------------

<!------------------------------------------------------------->
<!---------     ------------    -----------   ------------->
### 6.4.4 ASReml Standalone

Here we cover four approaches to fit genetic groups in an **ASReml** analysis. The first two use `nadiv` functions in `R` to create **Q** and <b>A\*</b>, which are each separately included in **ASReml** models. The third approach uses **ASReml** arguments and qualifiers to rely entirely upon the **ASReml** capabilities for fitting genetic groups. Whereas the first three approaches are used to fit genetic groups as fixed effects, the final approach demonstrates how to fit random effects of genetic groups.

The preparation of the `ggTutorial` data and requisite outputs from `nadiv` in `R` will first be demonstrated before detailing the code used in the **ASReml** model fitting. However, the basic files necessary to fit the **ASReml** models without any preparation in `R` are available in the accompanying supporting files on Dryad (Wolak & Reid 2016). The content of these files is summarized as:

-   `ggTutorial.dat` contains the `ggTutorial` data. This file contains space separated columns from a subset of columns in the original `data.frame` in the following order: `id`, `p`, `is`, `gen`, `f`, `foc0`, and `g1`. This is needed for all models fitted.

-   `reg.ped` contains the pedigree and is analogous to the first three columns of `ggTutorial`. This is needed to fit a model with genetic group effects explicitly as separate fixed regressions in the model.

-   `nadivAstar.giv` and `nadivAstar.txt` are needed to fit genetic group effects implicitly within the random effects of the model using `nadiv`'s `makeAinv()` in `R` to create <b>A\*</b> (contained in `nadivAstar.giv`). `nadivAstar.txt` provides the row names for `nadivAstar.giv` and is needed to associate levels in `nadivAstar.giv` to identities in the data `ggTutorial.dat`. Further explanation is provided below where these files are created.

-   `asremlGG.ped` contains the genetic group pedigree. It is the same as `reg.ped` except two rows are added for each of the two genetic groups and individuals do not have missing values, but genetic groups, specified for missing dams and sires. This is needed to fit genetic group effects implicitly within the random effects using **ASReml**'s own capabilities.

-   `ggReg.as`, `ggNadivAstarFxd.as`, `ggAsremlAstarFxd.as`, and `ggAsremlRan.as` are included within the folders with similar names. Each specifies a different model to run.

To prepare the data for an animal model, all that is necessary is to calculate the coefficients of inbreeding (*f*; to minimize bias from inbreeding depression, see the last paragraph of *Appendix 6.4.1*). This calculation can be done using the `makeAinv()` function in the `nadiv` package.

``` r
library(nadiv)
```

``` r
ggTutorial$f <- makeAinv(ggTutorial[, 1:3])$f
```

#### 6.4.4.1 Preparing a pedigree with genetic groups

Calculating the matrix of genetic group contributions (**Q**) and the augmented inverse relatedness matrix (<b>A\*</b>) requires a pedigree that has genetic groups indicated instead of missing values for parents. In the `ggTutorial` data, we will assign all individuals with unknown parents in the first generation phantom parents from one genetic group (`foc0`). All individuals in subsequent generations that have unknown parents will be assigned phantom parents from a second genetic group (`g1`). Alternatively, this second group could be further divided into genetic groups based on generation. However, the data were simulated such that every immigrant has the same expected mean breeding value, regardless of the generation in which it was born. Therefore, we will stick to a total of two genetic groups (but, feel free to model more groups!).

First, create an object that will be the pedigree with genetic groups (`ggPed`).

``` r
ggPed <- ggTutorial[, c("id", "dam", "sire", "is", "gen")]
naPar <- which(is.na(ggPed[, 2]))
ggPed$GG <- rep("NA", nrow(ggPed))
  # 'focal' genetic group = "foc0" and 'immigrant' = "g1"
  # obtained by pasting "foc" & "g" with immigrant status "0" or "1", respectively
  ggPed$GG[naPar] <- as.character(ggPed$is[naPar])
  ggPed$GG[ggPed$GG == "0"] <- paste0("foc", ggPed$GG[ggPed$GG == "0"])
  ggPed$GG[ggPed$GG == "1"] <- paste0("g", ggPed$GG[ggPed$GG == "1"])
ggPed[naPar, 2:3] <- ggPed[naPar, "GG"]
# Two rows need to be added for the genetic groups
## Genetic groups will be given the missing value 0 instead of parent identities
ggPed <- data.frame(id = c("foc0", "g1", as.character(ggPed$id)),
        dam = c(0, 0, as.character(ggPed$dam)),
        sire = c(0, 0, as.character(ggPed$sire)))
head(ggPed)
```

    ##     id  dam sire
    ## 1 foc0    0    0
    ## 2   g1    0    0
    ## 3    1 foc0 foc0
    ## 4    2 foc0 foc0
    ## 5    3 foc0 foc0
    ## 6    4 foc0 foc0

#### 6.4.4.2 Fixed explicit genetic group effects with Q (from `nadiv`)

Fitting genetic group effects explicitly requires the columns of **Q** to be included as separate fixed covariate regressions. The code below creates **Q** for the `ggPed` pedigree.

<!--- #FIXME  ---   ---->
``` r
Q <- ggcontrib(ggPed[, 1:3])
```

    ## Warning in numPed(pedigree): Zero in the dam column interpreted as a
    ## missing parent

    ## Warning in numPed(pedigree): Zero in the sire column interpreted as a
    ## missing parent

<!---  #END FIXME  ---- ----    ---->
Note that the two warnings above (regarding zeroes being interpreted as missing dams/sires) are to be expected and should not be cause for any concern in this particular pedigree/example.

Now we add the columns of **Q** (`foc0` and `g1`) as variables in `ggTutorial` so that they can be included in a model.

``` r
ggTutorial <- cbind(ggTutorial, Q)
```

``` r
str(ggTutorial)
```

    ## 'data.frame':    6000 obs. of  13 variables:
    ##  $ id     : Factor w/ 6000 levels "1","2","3","4",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ dam    : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ sire   : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ parAvgU: num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ mendel : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ u      : num  -0.393 -0.675 -0.442 -1.03 -1.913 ...
    ##  $ r      : num  -0.0836 0.142 1.892 -1.2982 0.4668 ...
    ##  $ p      : num  19.5 19.5 21.5 17.7 18.6 ...
    ##  $ is     : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ gen    : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ f      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ foc0   : num  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ g1     : num  0 0 0 0 0 0 0 0 0 0 ...

An important point to make is that, because the pedigree and data for `ggTutorial` have essentially the same sizes and structures (only `ggPed` now has two extra rows for the genetic groups - but these will be dropped automatically from **Q** by `ggcontrib()`), the genetic group coefficients from **Q** can be added directly to the data with `ggTutorial <- cbind(ggTutorial, Q)` above. However, whenever a pedigree contains more individuals than the data (or in a different order), an intermediate step is required. Specifically, since the **Q** matrix contains a row for every individual in the pedigree, only the subset of rows in **Q** that correspond to phenotyped individuals in a dataset should be extracted and ordered according to the identities in the data.

Now we need to write the pedigrees and portion of the data to files in formats that **ASReml** can use. Note, **ASReml** expects pedigree columns ordered ID, Sire, Dam and we will use '0' to denote missing values (instead of NA).

``` r
regPed <- ggTutorial[, 1:3]
regPed[is.na(regPed[, 2]), 2] <- 0
regPed[is.na(regPed[, 3]), 3] <- 0
```

``` r
write.table(regPed[, c(1,3,2)], "./asreml/ggReg/reg.ped",
    row.names = FALSE, quote = FALSE)
asData <- ggTutorial[, -c(2:7)]
write.table(asData, "./asreml/ggReg/ggTutorial.dat",
    row.names = FALSE, quote = FALSE)
```

The fixed regressions on `foc0` and `g1` can be combined with the standard <b>A<sup>-1</sup></b> to specify an animal model in **ASReml** (see also the "./asreml/ggReg" folder in the supporting files). Because unique solutions to the model do not exist (see discussion of estimability in *Appendix S6.2*) for the intercept, coefficient of inbreeding (*f*), and genetic group effects (`foc0` and `g1`), we fit the `foc0` genetic group as the last fixed effect. This ensures that the model sets `foc0` as a reference group and estimates the `g0` genetic group effect as a deviation from the mean breeding value in the `foc0` group.

    Fixed regression for Explicit genetic groups
     id !P
     p
     is !I 0 1
     gen 15
     f
     foc0
     g1
    reg.ped !ALPHA !SKIP 1 !MAKE !GIV !DIAG
    ggTutorial.dat !SKIP 1
    p ~ mu f g1 foc0 !r id

The estimate of the genetic group effects and standard errors are reported as the genetic group fixed effect regression coefficients and standard errors in the "ggReg.sln" file. These can be read back into `R`:

``` r
regPred <- read.table("./asreml/ggReg/ggReg.sln", header = TRUE)
```

so that we can see the fixed effect estimates and their standard errors:

``` r
regPred[1:4, ]
```

    ##   Model_Term Level Effect seEffect
    ## 1       foc0     1  0.000  0.00000
    ## 2         g1     1  3.117  0.08929
    ## 3          f     1  0.666  1.83200
    ## 4         mu     1 19.800  0.06579

Note the software has set the `foc0` group estimate to zero in this model, implying this is the reference group and the `g1` genetic group effect is the estimated deviation from this reference. The above estimate of 3.117 agrees with the expected difference between mean breeding value in the immigrant group (`3`) and the founder population (`0`) that was specified in the simulation.

The predicted breeding values (**a**) are also reported in the "ggReg.sln" file (see now `R` object `regPred`) and follow below the fixed effect estimates. Because the order of rows in **Q** matches the order of individual identities in `ggTutorial`, the total additive genetic effects (**u**) can be calculated without any further manipulation as (eqn 5 in main text):

``` r
reg_gHat <- matrix(regPred[1:2, "Effect"], ncol = 1)
reg_u <- Q %*% reg_gHat + regPred[-c(1:4), "Effect"]
```

Note, if more individuals are in the pedigree than the dataset and/or the identities in the pedigree and data are ordered differently, more coding is necessary to ensure that the predicted genetic effect for a particular individual is added to the correct weighted genetic group contribution. Additionally, if more random effects are included in the model, care has to be taken in order to ensure the correct random effect predictions are extracted for the breeding values (**a**). The `R` function `match()` is particularly helpful in these cases.

#### 6.4.4.3 Fixed implicit genetic group effects with A\* (from `nadiv`)

Genetic group effects can be fit implicitly within the random effects when the augmented inverse relatedness matrix <b>A\*</b> has been supplied as a list to **ASReml**. The code below creates <b>A\*</b> for the `ggPed` pedigree and saves it as a list in the format **ASReml** requires.

First, we will calculate **Q** and include it in the dataset for consistency with other versions of this data.

<!--- #FIXME  ---   ---->
``` r
Q <- ggcontrib(ggPed[, 1:3])
```

    ## Warning in numPed(pedigree): Zero in the dam column interpreted as a
    ## missing parent

    ## Warning in numPed(pedigree): Zero in the sire column interpreted as a
    ## missing parent

<!---  #END FIXME  ---- ----    ---->
Note that the two warnings above (regarding zeroes being interpreted as missing dams/sires) are to be expected and should not be cause for any concern in this particular pedigree/example.

Add the columns of **Q** (`foc0` and `g1`) as variables in `ggTutorial`.

``` r
ggTutorial <- cbind(ggTutorial, Q)
```

``` r
str(ggTutorial)
```

    ## 'data.frame':    6000 obs. of  13 variables:
    ##  $ id     : Factor w/ 6000 levels "1","2","3","4",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ dam    : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ sire   : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ parAvgU: num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ mendel : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ u      : num  -0.393 -0.675 -0.442 -1.03 -1.913 ...
    ##  $ r      : num  -0.0836 0.142 1.892 -1.2982 0.4668 ...
    ##  $ p      : num  19.5 19.5 21.5 17.7 18.6 ...
    ##  $ is     : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ gen    : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ f      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ foc0   : num  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ g1     : num  0 0 0 0 0 0 0 0 0 0 ...

An important point to make is that, because the pedigree and data for `ggTutorial` have essentially the same sizes and structures (only `ggPed` now has two extra rows for the genetic groups - but these will be dropped automatically from **Q** by `ggcontrib()`), the genetic group coefficients from **Q** can be added directly to the data with `ggTutorial <- cbind(ggTutorial, Q)` above. However, whenever a pedigree contains more individuals than the data (or in a different order), an intermediate step is required. Specifically, since the **Q** matrix contains a row for every individual in the pedigree, only the subset of rows in **Q** that correspond to phenotyped individuals in a dataset should be extracted and ordered according to the identities in the data.

Now we create the sparse matrix of <b>A\*</b> for the `ggPed` pedigree using the `makeAinv` function in the `nadiv` package.

``` r
listAstar <- makeAinv(ggPed, ggroups = 2, gOnTop = TRUE)$listAinv
```

In the above code, we directed `makeAinv()` to include the genetic groups on top (top-left block) of <b>A\*</b> using `gOnTop = TRUE`. This is necessary so that **ASReml** solves the equations in an order that improves model convergence (see *Appendix 6.3* and its subsections). In this case, the **ASReml** model converges. In other pedigrees, models, and/or genetic group classifications, if singularities occur then sometimes this can be overcome by adding one (or a small value) to the diagonal elements of the group-group portion of the matrix (although see discussion of this in *Appendix S5* and cautions/considerations in *Appendix 6.3.2*). The alteration is demonstrated above in the `asreml` tutorial along with cautions and interpretations (*Appendix 6.4.3.4*). The resulting list could then be given to the standalone **ASReml** program much the same as we do here (see also *Appendix 6.4.4.5*). Note, the easiest approach is to adjust the genetic group diagonals of <b>A\*</b> using the matrix returned by `makeAinv()` (i.e., object `Ainv` in the list), then convert the matrix to the list format required by **ASReml**.

Now we need to write the data (`ggTutorial`) and list of the <b>A\*</b> (`listAstar`) to a file for **ASReml**. Further, we need to write to a file all of the identity codes in the pedigree that match to the rows and columns in `listAstar`. This enables **ASReml** to match levels in the generalized inverse (row and column numbers of `listAstar`) to observations in `ggTutorial`. Critically, we need to include the "geneticGroups" attribute of `listAstar` as the first line in the file containing this list.

``` r
# Write the nadiv Astar ginv list
write.table(listAstar,
    file = "./asreml/ggNadivAstarFxd/nadivAstar.giv",
    col.names = FALSE, row.names = FALSE)
  fConn <- file("./asreml/ggNadivAstarFxd/nadivAstar.giv", "r+")
  Lines <- readLines(fConn)
  # Add the !GROUPSDF *n* qualifier as the first line of the file
  ngrps <- attr(listAstar, "geneticGroups")[1]
  writeLines(c(paste0("!GROUPSDF ", ngrps), Lines), con = fConn)
  close(fConn) 

# Create mapping between order of identities in listAstar and their unique codes
write.table(ggPed[, 1], file = "./asreml/ggNadivAstarFxd/nadivAstar.txt",
    col.names = FALSE, row.names = FALSE, quote = FALSE)

# Write the data
asData <- ggTutorial[, -c(2:7)]
write.table(asData, "./asreml/ggNadivAstarFxd/ggTutorial.dat",
    row.names = FALSE, quote = FALSE)
```

Now we can run the animal model in **ASReml** using the supplied general inverse matrix `listAstar` in the file "nadivAstar.giv" (see also the "./asreml/ggNadivAstarFxd" folder in the supporting files):

    nadiv's Astar for fixed genetic group effects implicitly within the random effects
     id !A !L nadivAstar.txt
     p
     is !I 0 1
     gen 15
     f
     foc0
     g1
    nadivAstar.giv
    ggTutorial.dat !SKIP 1
    !LAST id 2
    p ~ mu f !r giv1(id)

Here, genetic group effects are fitted in the model using `nadiv`'s generalized inverse by including `giv1(id)` in the random effects syntax statement. The line below the data file, but before the model statement, (`!LAST id 2`) tells **ASReml** to fit the first two equations associated with the `id` variable (the genetic groups) last. This is done to help the model to converge (see *Appendix 6.3.1*).

The genetic group predictions and standard errors are reported in the "ggNadivAstarFxd.sln" file. These can be read back into `R`:

``` r
nadivAstarFxdPred <- read.table("./asreml/ggNadivAstarFxd/ggNadivAstarFxd.sln", header = TRUE)
```

so that we can see the fixed effect estimates and their standard errors along with the genetic group predictions:

``` r
nadivAstarFxdPred[1:4, ]
```

    ##   Model_Term Level Effect seEffect
    ## 1          f     1  0.666  1.83200
    ## 2         mu     1  0.000  0.00000
    ## 3   giv1(id)  foc0 19.800  0.06579
    ## 4   giv1(id)    g1 22.910  0.06008

Note the software has set the intercept (`mu`) estimate to zero in this model (see discussion of estimability in *Appendix S6.2*), implying this is the reference to which the `foc0` and `g1` genetic group effects are the estimated deviations. Therefore, the predicted total additive genetic effects (**u**) include the overall phenotypic mean. However, the difference between the genetic group effects 3.11 agrees with the expected difference between the mean breeding values of the immigrant group (`3`) and the founder population (`0`) that was specified in the simulation.

The predicted total additive genetic effects (**u**) are also reported in the "ggNadivAstarFxd.sln" file (see now `R` object `nadivAstarFxdPred`) and follow below the genetic group predictions. Because the order of rows in **Q** matches the order of individual identities in `ggTutorial`, the predicted breeding values (**a**) can be calculated without any further manipulation as (eqn 6 in main text):

``` r
AstarFxd_gHat <- matrix(nadivAstarFxdPred[3:4, "Effect"], ncol = 1)
AstarFxd_a <- nadivAstarFxdPred[-c(1:4), "Effect"] - Q %*% AstarFxd_gHat
```

Note, if more individuals are in the pedigree than the dataset and/or the identities in the pedigree and data are ordered differently, more coding is necessary to ensure that the predicted genetic effect for a particular individual is added to the correct weighted genetic group contribution. Additionally, if more random effects are included in the model, care has to be taken in order to ensure the correct random effect predictions are extracted for the total additive genetic effects (**u**) and used to calculate breeding values (**a**). The `R` function `match()` is particularly helpful in these cases.

#### 6.4.4.4 Fixed implicit genetic group effects with A\* (from ASReml)

**ASReml** can fit genetic group effects implicitly within the random effects using its own function to construct <b>A\*</b>. First, we construct the data and pedigree files in `R`.

Calculate **Q** and include it in the dataset for consistency with other versions of this data.

<!--- #FIXME  ---   ---->
``` r
Q <- ggcontrib(ggPed[, 1:3])
```

    ## Warning in numPed(pedigree): Zero in the dam column interpreted as a
    ## missing parent

    ## Warning in numPed(pedigree): Zero in the sire column interpreted as a
    ## missing parent

<!---  #END FIXME  ---- ----    ---->
Note that the two warnings above (regarding zeroes being interpreted as missing dams/sires) are to be expected and should not be cause for any concern in this particular pedigree/example.

Add the columns of **Q** (`foc0` and `g1`) as variables in `ggTutorial`.

``` r
ggTutorial <- cbind(ggTutorial, Q)
```

``` r
str(ggTutorial)
```

    ## 'data.frame':    6000 obs. of  13 variables:
    ##  $ id     : Factor w/ 6000 levels "1","2","3","4",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ dam    : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ sire   : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ parAvgU: num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ mendel : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ u      : num  -0.393 -0.675 -0.442 -1.03 -1.913 ...
    ##  $ r      : num  -0.0836 0.142 1.892 -1.2982 0.4668 ...
    ##  $ p      : num  19.5 19.5 21.5 17.7 18.6 ...
    ##  $ is     : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ gen    : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ f      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ foc0   : num  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ g1     : num  0 0 0 0 0 0 0 0 0 0 ...

An important point to make is that, because the pedigree and data for `ggTutorial` have essentially the same sizes and structures (only `ggPed` now has two extra rows for the genetic groups - but these will be dropped automatically from **Q** by `ggcontrib()`), the genetic group coefficients from **Q** can be added directly to the data with `ggTutorial <- cbind(ggTutorial, Q)` above. However, whenever a pedigree contains more individuals than the data (or in a different order), an intermediate step is required. Specifically, since the **Q** matrix contains a row for every individual in the pedigree, only the subset of rows in **Q** that correspond to phenotyped individuals in a dataset should be extracted and ordered according to the identities in the data.

Write the pedigree containing genetic groups and the relevant portion of the data to files in formats that **ASReml** can use. Note, **ASReml** expects pedigree columns ordered ID, Sire, Dam.

``` r
# Write the pedigree
write.table(ggPed[, c(1,3,2)], "./asreml/ggAsremlAstarFxd/asremlGG.ped",
    row.names = FALSE, quote = FALSE)

# Write the data
asData <- ggTutorial[, -c(2:7)]
write.table(asData, "./asreml/ggAsremlAstarFxd/ggTutorial.dat",
    row.names = FALSE, quote = FALSE)
```

The model is specified in **ASReml** (see also the "./asreml/ggAsremlAstarFxd" folder in the supporting files):

    ASReml's Astar for fixed genetic groups implicitly within the random effects
     id !P
     p
     is !I 0 1
     gen 15
     f
     foc0
     g1
    asremlGG.ped !ALPHA !SKIP 1 !MAKE !GIV !DIAG !GROUPS 2
    ggTutorial.dat !SKIP 1
    !LAST id 2
    p ~ mu f !r id

Here, genetic groups are specified as the top two lines of the pedigree file `asremlGG.ped` using `!GROUPS 2` as an additional qualifier to the pedigree file line. **ASReml** then creates <b>A\*</b> from the pedigree and saves this as a list in the file "ggAsremlAstarFxd\_A.giv" (when the `!GIV` qualifier is included in the pedigree line). The line below the data file, but before the model statement, (`!LAST id 2`) tells **ASReml** to fit the first two equations associated with the `id` variable (the genetic groups) last. This is done to help the model to converge (see *Appendix 6.3.1*).

The genetic group predictions and standard errors are reported in the "ggAsremlAstarFxd.sln" file. These can be read back into `R`:

``` r
asremlAstarFxdPred <- read.table("./asreml/ggAsremlAstarFxd/ggAsremlAstarFxd.sln", header = TRUE)
```

so that we can see the fixed effect estimates and their standard errors along with the genetic group predictions:

``` r
asremlAstarFxdPred[1:4, ]
```

    ##   Model_Term Level Effect seEffect
    ## 1          f     1  0.666  1.83200
    ## 2         mu     1  0.000  0.00000
    ## 3         id  foc0 19.800  0.06579
    ## 4         id    g1 22.910  0.06008

Note the software has set the intercept (`mu`) estimate to zero in this model (see discussion of estimability in *Appendix S6.2*), implying this is the reference to which the `foc0` and `g1` genetic group effects are the estimated deviations. Therefore, the predicted total additive genetic effects (**u**) include the overall phenotypic mean. However, the difference between the genetic group effects 3.11 agrees with the expected difference between the mean breeding values of the immigrant group (`3`) and the founder population (`0`) that was specified in the simulation.

The predicted total additive genetic effects (**u**) are also reported in the "ggAsremlAstarFxd.sln" file (see now `R` object `asremlAstarFxdPred`) and follow below the genetic group predictions. Because the order of rows in **Q** matches the order of individual identities in `ggTutorial`, the predicted breeding values (**a**) can be calculated without any further manipulation as (eqn 6 in main text):

``` r
AstarFxd_gHat <- matrix(asremlAstarFxdPred[3:4, "Effect"], ncol = 1)
AstarFxd_a <- asremlAstarFxdPred[-c(1:4), "Effect"] - Q %*% AstarFxd_gHat
```

Note, if more individuals are in the pedigree than the dataset and/or the identities in the pedigree and data are ordered differently, more coding is necessary to ensure that the predicted genetic effect for a particular individual is added to the correct weighted genetic group contribution. Additionally, if more random effects are included in the model, care has to be taken in order to ensure the correct random effect predictions are extracted for the total additive genetic effects (**u**) and used to calculate breeding values (**a**). The `R` function `match()` is particularly helpful in these cases.

#### 6.4.4.5 Random implicit genetic group effects with A\* (from ASReml)

Fitting random effects of genetic groups implicity within the random effects of an animal model requires a modified augmented inverse relatendess matrix (<b>A\*</b>) to be created. This can be accomplished in `R` with `nadiv` such that the object is then passed to **ASReml** or the built-in functions of **ASReml** can do this while fitting the model. Creating the modified matrix in `R` using `nadiv` has been demonstrated in *Appendix 6.4.3.4* and including such a matrix in an **ASReml** analysis is demonstrated in *Appendix 6.4.4.3*. Here, we demonstrate how to use **ASReml** to create a modified <b>A\*</b> for the `ggPed` pedigree that has the necessary alterations (see discussion of this in *Appendix S5* and cautions/considerations in *Appendix 6.3.2*) so the model interprets the genetic groups as random effects and fit the model all in a single execution of **ASReml**.

First, we construct the data and pedigree files in `R`. Calculate **Q** and include it in the dataset for consistency with other versions of this data.

<!--- #FIXME  ---   ---->
``` r
Q <- ggcontrib(ggPed[, 1:3])
```

    ## Warning in numPed(pedigree): Zero in the dam column interpreted as a
    ## missing parent

    ## Warning in numPed(pedigree): Zero in the sire column interpreted as a
    ## missing parent

<!---  #END FIXME  ---- ----    ---->
Note that the two warnings above (regarding zeroes being interpreted as missing dams/sires) are to be expected and should not be cause for any concern in this particular pedigree/example.

Add the columns of **Q** (`foc0` and `g1`) as variables in `ggTutorial`.

``` r
ggTutorial <- cbind(ggTutorial, Q)
```

``` r
str(ggTutorial)
```

    ## 'data.frame':    6000 obs. of  13 variables:
    ##  $ id     : Factor w/ 6000 levels "1","2","3","4",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ dam    : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ sire   : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ parAvgU: num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ mendel : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ u      : num  -0.393 -0.675 -0.442 -1.03 -1.913 ...
    ##  $ r      : num  -0.0836 0.142 1.892 -1.2982 0.4668 ...
    ##  $ p      : num  19.5 19.5 21.5 17.7 18.6 ...
    ##  $ is     : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ gen    : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ f      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ foc0   : num  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ g1     : num  0 0 0 0 0 0 0 0 0 0 ...

An important point to make is that, because the pedigree and data for `ggTutorial` have essentially the same sizes and structures (only `ggPed` now has two extra rows for the genetic groups - but these will be dropped automatically from **Q** by `ggcontrib()`), the genetic group coefficients from **Q** can be added directly to the data with `ggTutorial <- cbind(ggTutorial, Q)` above. However, whenever a pedigree contains more individuals than the data (or in a different order), an intermediate step is required. Specifically, since the **Q** matrix contains a row for every individual in the pedigree, only the subset of rows in **Q** that correspond to phenotyped individuals in a dataset should be extracted and ordered according to the identities in the data.

Write the pedigree containing genetic groups and the relevant portion of the data to files in formats that **ASReml** can use. Note, **ASReml** expects pedigree columns ordered ID, Sire, Dam.

``` r
# Write the pedigree
write.table(ggPed[, c(1,3,2)], "./asreml/ggAsremlRan/asremlGG.ped",
    row.names = FALSE, quote = FALSE)

# Write the data
asData <- ggTutorial[, -c(2:7)]
write.table(asData, "./asreml/ggAsremlRan/ggTutorial.dat",
    row.names = FALSE, quote = FALSE)
```

The model is specified in **ASReml** (see also the "./asreml/ggAsremlRan" folder in the supporting files):

    ASReml's Astar for random genetic groups implicitly within the random effects
     id !P
     p
     is !I 0 1
     gen 15
     f
     foc0
     g1
    asremlGG.ped !ALPHA !SKIP 1 !MAKE !GIV !DIAG !GROUPS 2 !GOFFSET 0.1
    ggTutorial.dat !SKIP 1
    !LAST id 2
    p ~ mu f !r id

Here, genetic groups are specified as the top two lines of the pedigree file `asremlGG.ped` using `!GROUPS 2` as an additional qualifier to the pedigree file line. Further, the `!GOFFSET 0.1` qualifier to the pedigree line instructs **ASReml** to add `0.1` to the diagonal elements of the group-group equations in <b>A\*</b> (see discussion of this in *Appendix S5*, cautions/considerations in *Appendix 6.3.2*, and other examples in *Appendix 6.4.2.3* and *Appendix 6.4.3.4*). As a result of this alteration, the genetic group effects are conceptually random effects in the model. The addition of 0.1 to the diagonal element of the <b>A\*</b> constrains the variance in genetic group effects to be one tenth the additive genetic variance estimated by the model. A value of 0.1 may bias estimates of additive genetic variance in models of other datasets. Therefore, we advise testing values other than 0.1 (e.g., 1 or 10) and comparing estimates from these alternative models. Note that a value of 1 constrains the variance in genetic group effects to be equal to the additive genetic variance estimated by the model.

**ASReml** then creates <b>A\*</b> from the pedigree and saves this as a list in the file "ggAsremlRan\_A.giv" (when the `!GIV` qualifier is included in the pedigree line). The line below the data file, but before the model statement, (`!LAST id 2`) tells **ASReml** to fit the first two equations associated with the `id` variable (the genetic groups) last. This is done to help the model to converge (see *Appendix 6.3.1*).

The genetic group predictions and standard errors are reported in the "ggAsremlRan.sln" file. These can be read back into `R`:

``` r
asremlRanPred <- read.table("./asreml/ggAsremlRan/ggAsremlRan.sln", header = TRUE)
```

so that we can see the fixed effect estimates and their standard errors along with the genetic group predictions:

``` r
asremlRanPred[1:4, ]
```

    ##   Model_Term Level  Effect seEffect
    ## 1          f     1  0.6656    1.832
    ## 2         mu     1 21.3600    2.321
    ## 3         id  foc0 -1.5580    2.321
    ## 4         id    g1  1.5580    2.321

Note the intercept in this model (`mu`) is estimable and so the model returns this estimate along with predicting both of the genetic group effects. The genetic group predictions, themselves, are random effects and thus are deviations from their expected value of 0. The difference between the genetic groups 3.116 agrees with the expected difference between the mean breeding values of the immigrant group (`3`) and the founder population (`0`) that was specified in the simulation.

The predicted total additive genetic effects (**u**) are also reported in the "ggAsremlRan.sln" file (see now `R` object `asremlRanPred`) and follow below the genetic group predictions. Here, the prediced total additive genetic effects include the genetic group effects expressed as deviations from 0. Because the order of rows in **Q** matches the order of individual identities in `ggTutorial`, the predicted breeding values (**a**) can be calculated without any further manipulation as (eqn 6 in main text):

``` r
AstarRan_gHat <- matrix(asremlRanPred[3:4, "Effect"], ncol = 1)
AstarRan_a <- asremlRanPred[-c(1:4), "Effect"] - Q %*% AstarRan_gHat
```

Note, if more individuals are in the pedigree than the dataset and/or the identities in the pedigree and data are ordered differently, more coding is necessary to ensure that the predicted genetic effect for a particular individual is added to the correct weighted genetic group contribution. Additionally, if more random effects are included in the model, care has to be taken in order to ensure the correct random effect predictions are extracted for the total additive genetic effects (**u**) and used to calculate breeding values (**a**). The `R` function `match()` is particularly helpful in these cases.

------------------------------------------------------------------------

<!------------------------------------------------------------->
<!---------     ------------    -----------   ------------->
### 6.4.5 WOMBAT

Here are three approaches to fit genetic groups in a **WOMBAT** analysis. All three require some preparation using the `nadiv` functions in `R`. However, only the explicit regression on the genetic group covariate includes genetic groups as fixed effects. The second and third examples below include genetic groups as random effects in the model. It should be noted that **WOMBAT** sets up the mixed model slightly differently than the previous software programs. Notably, **WOMBAT** fits response variables and covariates as deviations from the mean response and mean of each covariate, respectively. Thus, regression coefficients are equivalent, but have slightly different estimates/interpretations from the estimates in the above tutorials. The preparation of the `ggTutorial` data and requisite outputs from `nadiv` in `R` will first be demonstrated before detailing the code used in the **WOMBAT** model fitting. However, the basic files necessary to fit the **WOMBAT** models without any preparation in `R` are available in the accompanying supporting files on Dryad (Wolak & Reid 2016). The content of these files is:

-   `ggTutorial.d` contains the re-coded `ggTutorial` data. A value of 10000 has been added to the `id` to allow unique integer codes for the genetic groups (**WOMBAT** requires integer identities). This file contains space separated columns from a subset of the columns in the original `data.frame` in the following order: `id`, `p`, `is`, `gen`, `f`, `foc0`, and `g1`. This is needed for all three types of models fitted.

-   `wombatPed.d` contains the pedigree and is analogous to the first three columns of `ggTutorial`

-   `ggReg.par`, `ggNadivAstarRan.par`, and `ggWombatRan.par` are included within each folder and each specifies a separate model to run.

-   `id.gin` contains an inverse covariance matrix, in list format. This is either the `nadiv` created <b>A\*</b> or <b>A<sup>-1</sup></b>, depending on the analysis.

-   `id.codes` contains the mapping of running integers to the particular name of individuals in `ggTutorial.d`. For **WOMBAT**'s implementation of a genetic group analysis, this file must contain an extra column (dummy column in these analyses) that contains an integer code specifying if a genotype is available, as well as extra columns with the **Q** matrix (which must be constructed within `nadiv` as there is no function to do this in **WOMBAT**).

To prepare the data for an animal model, all that is necessary is to calculate the coefficients of inbreeding (*f*; to minimize bias from inbreeding depression, see the last paragraph of *Appendix 6.4.1*). This calculation can be done using the `makeAinv()` function in the `nadiv` package.

``` r
library(nadiv)
```

``` r
ggTutorial$f <- makeAinv(ggTutorial[, 1:3])$f
```

#### 6.4.5.1 Preparing a pedigree with genetic groups

Calculating the matrix of genetic group contributions (**Q**) and the augmented inverse relatedness matrix (<b>A\*</b>) requires a pedigree that has genetic groups indicated instead of missing values for parents. In the `ggTutorial` data, we will assign all individuals with unknown parents in the first generation phantom parents from one genetic group (`foc0`). All individuals in subsequent generations that have unknown parents will be assigned phantom parents from a second genetic group (`g1`). Alternatively, this second group could be further divided into genetic groups based on generation. However, the data were simulated such that every immigrant has the same expected mean breeding value, regardless of the generation in which it was born. Therefore, we will stick to a total of two genetic groups (but, feel free to model more than 2 groups!).

Create an object that will be the pedigree with genetic groups (`ggPed`).

``` r
ggPed <- ggTutorial[, c("id", "dam", "sire", "is", "gen")]
naPar <- which(is.na(ggPed[, 2]))
ggPed$GG <- rep("NA", nrow(ggPed))
  # 'focal' genetic group = "foc0" and 'immigrant' = "g1"
  # obtained by pasting "foc" & "g" with immigrant status "0" or "1", respectively
  ggPed$GG[naPar] <- as.character(ggPed$is[naPar])
  ggPed$GG[ggPed$GG == "0"] <- paste0("foc", ggPed$GG[ggPed$GG == "0"])
  ggPed$GG[ggPed$GG == "1"] <- paste0("g", ggPed$GG[ggPed$GG == "1"])
ggPed[naPar, 2:3] <- ggPed[naPar, "GG"]
```

Note, the approach and format of the pedigree below is different from the `Q1988sub` pedigree (*Appendix 6.1*) - mostly because the format below makes it easier to specify genetic groups in a pedigree for this particular case, but partly to illustrate the flexibility of `nadiv` functions. To be more specific, the `Q1988sub` pedigree contained two extra rows for the two genetic groups. The `ggroups` argument to `makeAinv()` supplied an integer indicating how many rows at the begining of the pedigree contained genetic groups and not individuals. However, below no extra rows are added to the pedigree and the character vector given to the `ggroups` argument in both `ggcontrib()` and `makeAinv()` specifies the unique genetic groups.

#### 6.4.5.2 Fixed explicit genetic group effects with Q (from `nadiv`)

Fitting genetic group effects explicitly requires the columns of **Q** to be included as separate fixed covariate regressions. The code below creates **Q** for the `ggPed` pedigree and adds the columns (`foc0` and `g1`) as variables in `ggTutorial` so that they can be included in a model.

``` r
Q <- ggcontrib(ggPed[, 1:3], ggroups = c("foc0", "g1"))
```

``` r
ggTutorial <- cbind(ggTutorial, Q)
```

``` r
str(ggTutorial)
```

    ## 'data.frame':    6000 obs. of  13 variables:
    ##  $ id     : Factor w/ 6000 levels "1","2","3","4",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ dam    : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ sire   : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ parAvgU: num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ mendel : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ u      : num  -0.393 -0.675 -0.442 -1.03 -1.913 ...
    ##  $ r      : num  -0.0836 0.142 1.892 -1.2982 0.4668 ...
    ##  $ p      : num  19.5 19.5 21.5 17.7 18.6 ...
    ##  $ is     : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ gen    : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ f      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ foc0   : num  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ g1     : num  0 0 0 0 0 0 0 0 0 0 ...

An important point to make is that, because the pedigree and data for `ggTutorial` have essentially the same sizes and structures, the genetic group coefficients from **Q** can be added directly to the data with `ggTutorial <- cbind(ggTutorial, Q)` above. However, whenever a pedigree contains more individuals than the data (or in a different order), an intermediate step is required. Specifically, since the **Q** matrix contains a row for every individual in the pedigree, only the subset of rows in **Q** that correspond to phenotyped individuals in a dataset should be extracted and ordered according to the identities in the data.

Now we need to write the pedigree and portion of the data to files in formats that **WOMBAT** can use. Also, for consistency with analyses below, we will add a constant value of 10000 to the integer `id`s.

``` r
IDaddDF <- ggTutorial[, -c(4:7)]
IDaddDF[, 1] <- as.integer(as.character(IDaddDF[, 1])) + 10000
IDaddDF[, 2] <- as.integer(as.character(IDaddDF[, 2])) + 10000
IDaddDF[, 3] <- as.integer(as.character(IDaddDF[, 3])) + 10000

IDaddPed <- IDaddDF[, 1:3]
IDaddPed[is.na(IDaddPed[, 2]), 2] <- 0
IDaddPed[is.na(IDaddPed[, 3]), 3] <- 0
```

``` r
write.table(IDaddDF[, -c(2:3)], "./wombat/ggTutorial.d", col.names = FALSE,
    row.names = FALSE)
write.table(IDaddPed, "./wombat/wombatPed.d", col.names = FALSE,
    row.names = FALSE)
```

The fixed regressions on `foc0` and `g1` can be combined with the standard <b>A<sup>-1</sup></b> to specify an animal model in **WOMBAT** (see also the "./wombat/ggReg" folder in the supporting files). Because unique solutions to the model do not exist (see discussion of estimability in *Appendix S6.2*) for the coefficient of inbreeding (*f*) and genetic group effects (`foc0` and `g1`), we fit the `foc0` genetic group as the last fixed effect. This ensures that the model sets `foc0` as a reference group and estimates the `g0` genetic group effect as a deviation from the mean breeding value in the `foc0` group.

    COM Explicit genetic groups fitted as Fixed regressions
    PED ../wombatPed.d
    DAT ../ggTutorial.d
      id 6000
      p
      is 2
      gen 15
      f
      foc0
      g1
    end
    ANAL UNI
    MODEL
     COV f(1)
     COV g1(1)
     COV foc0(1)
     RAN id NRM
     TR p
    END
    VAR id 1
     0.5
    VAR error 1
     0.5
    SPECIAL
      COVZER f(1) FIT
      COVZER foc0(1) FIT
      COVZER g1(1) FIT
    END

Note the last `SPECIAL` section, which is necessary because both the coefficient of inbreeding (*f*) and the genetic group contributions (`foc0` and `g1`) have values of zero that are meaningful (i.e., these zeroes do not indicate missing values).

Estimates of the genetic group effects are reported in the "./wombat/ggReg/FixSolutions.out" file. Note the software has set the `foc0` group estimate to zero in this model, implying this is the reference group and the `g1` genetic group effect is the estimated deviation. The estimate of 3.098 agrees with the expected difference between mean breeding value in the immigrant group (`3`) and the founder population (`0`) that was specified in the simulation.

The predicted breeding values (**a**) are also reported in the "RnSoln\_id.dat" file.

#### 6.4.5.3 Random implicit genetic group effects with A\* (from `nadiv`)

Fitting genetic group effects implicitly within the random effects requires the list format of a modified <b>A\*</b> to be supplied to **WOMBAT**. Currently, there is no option to incorporate the number of genetic groups fitted implicitly when the software calculates the residual degrees of freedom. Residual degrees of freedom that incorporate the number of genetic groups are needed to calculate the residual variance and the log-likelihood for a model in which genetic groups are considered fixed effects (e.g., see how the residual degrees of freedom are adjusted in an `asreml` analysis in *Appendix 6.4.3.4*). However, the number of genetic groups are not considered in these same calculations when genetic groups are considered random effects. Therefore, we demonstrate fitting genetic groups implicitly when groups are treated as random effects.

Fitting random effects of genetic groups implicity within the random effects of an animal model requires the augmented inverse relatendess matrix (<b>A\*</b>) to be supplied as a list. First, we create the sparse matrix of <b>A\*</b> for the `ggPed` pedigree using the `makeAinv` function in the `nadiv` package in `R`, then make the necessary alterations (see discussion of this in *Appendix S5* and cautions/considerations in *Appendix 6.3.2*) so the model interprets the genetic groups as random effects, convert the sparse matrix to a list, and fit the model.

The code below creates <b>A\*</b> for the `ggPed` pedigree, the data file, and `id.codes` so that the genetic group and identities in <b>A\*</b> can be associated with identities in the data file. We will also calculate **Q** and include it in the dataset for consistency with other versions of this data.

``` r
Q <- ggcontrib(ggPed[, 1:3], ggroups = c("foc0", "g1"))
ggTutorial <- cbind(ggTutorial, Q)
```

An important point to make is that, because the pedigree and data for `ggTutorial` have essentially the same sizes and structures (only `ggPed` now has two extra rows for the genetic groups - but these will be dropped automatically from **Q** by `ggcontrib()`), the genetic group coefficients from **Q** can be added directly to the data with `ggTutorial <- cbind(ggTutorial, Q)` above. However, whenever a pedigree contains more individuals than the data (or in a different order), an intermediate step is required. Specifically, since the **Q** matrix contains a row for every individual in the pedigree, only the subset of rows in **Q** that correspond to phenotyped individuals in a dataset should be extracted and ordered according to the identities in the data.

Next create <b>A\*</b>:

``` r
Astar <- makeAinv(ggPed[, 1:3], ggroups = c("foc0", "g1"), gOnTop = TRUE)$Ainv
```

Next we add a small value to the diagonal elements of the group-group portion of the matrix. As a result of this alteration, the genetic group effects are conceptually random effects in the model. Below, we add 0.1 to the diagonal element of the <b>A\*</b> which will constrain the variance in genetic group effects to be one tenth the additive genetic variance estimated by the model. A value of 0.1 may bias estimates of additive genetic variance in models of other datasets. Therefore, we advise testing values other than 0.1 (e.g., 1 or 10) and comparing estimates from these alternative models. Note that a value of 1 constrains the variance in genetic group effects to be equal to the additive genetic variance estimated by the model.

The easiest approach is to adjust the genetic group diagonals of <b>A\*</b>, and this is why we stored the matrix returned by `makeAinv()` above (i.e., object `Ainv` in the list). Once the additions are made, we can then convert the matrix to the list format required by `asreml`.

``` r
AstarAdd <- Astar
(ggrows <- match(c("foc0", "g1"), dimnames(AstarAdd)[[1]]))
```

    ## [1] 1 2

``` r
diag(AstarAdd)[ggrows] <- diag(AstarAdd)[ggrows] + 0.1
Astar[ggrows, ggrows]
```

    ## 2 x 2 sparse Matrix of class "dgCMatrix"
    ##             
    ## foc0 400   .
    ## g1     . 520

``` r
AstarAdd[ggrows, ggrows]
```

    ## 2 x 2 sparse Matrix of class "dgCMatrix"
    ##                 
    ## foc0 400.1   .  
    ## g1     .   520.1

``` r
# convert to the list format using `nadiv::sm2list()`
listAstarAdd <- sm2list(AstarAdd, rownames = rownames(AstarAdd))
```

Now we need to write the data, modified <b>A\*</b> list, and identity codes to files in formats that **WOMBAT** can use. Also, for consistency with analyses below, we will add a constant value of 10000 to the integer `id`s.

``` r
IDaddDF <- ggTutorial[, -c(2:7)]
IDaddDF[, 1] <- as.integer(as.character(IDaddDF[, 1])) + 10000
write.table(IDaddDF[, -c(2:3)], "./wombat/ggTutorial.d", col.names = FALSE,
    row.names = FALSE)

# Write the general inverse of A* to file, note the re-ordering of first two columns
## First, calculate the log determinant (see note below)
AstarAddLogdet <- -1 * determinant(AstarAdd, logarithm = TRUE)$modulus[1]
write.table(listAstarAdd[, c(2,1,3)], "./wombat/ggNadivAstarRan/id.gin",
    col.names = FALSE, row.names = FALSE)
  fConn <- file("./wombat/ggNadivAstarRan/id.gin", "r+")
  Lines <- readLines(fConn)
# Add the log-determinant as the first line of the file
  writeLines(c(AstarAddLogdet, Lines), con = fConn)
  close(fConn) 

# Create mapping between order of identities and their unique codes
## including 2 genetic groups
id.codes <- cbind(seq(1, max(listAstarAdd[, 1:2]), 1), c(1, 2, 10000 + seq(nrow(ggTutorial))))
write.table(id.codes, "./wombat/ggNadivAstarRan/id.codes", col.names = FALSE,
    row.names = FALSE)
```

Note above that we added a step to calculate the log-determinant of the modified <b>A\*</b>. **WOMBAT** requires this as the first line of a `.gin` (generalized inverse) file (note this is the log-determinant of the non-inverse matrix, which can be calculated from the inverse).

The animal model in **WOMBAT** using the supplied general inverse matrix <b>A\*</b> (see also the "./wombat/ggNadivAstarRan" folder in the supporting files):

    COM Random genetic groups fitted implicitly using Astar from nadiv
    DAT ../ggTutorial.d
      id 6002
      p
      is 2
      gen 15
      f
      foc0
      g1
    end
    ANAL UNI
    MODEL
     COV  f(1)
     RAN  id    GIN
     TR   p
    END
    VAR id 1
     0.5
    VAR error 1
     0.5
    SPECIAL
      COVZER f(1) FIT
    END

Note the last `SPECIAL` section, which is necessary because the coefficient of inbreeding (*f*) has values of zero that are meaningful (i.e., these zeroes do not indicate missing values).

The genetic group effect predictions are reported as the first two values in the "./wombat/ggNadivAstarRan/RnSoln\_id.dat" file. The genetic group predictions are random effects and thus are deviations from their expected value of zero. The difference between the genetic groups 3.117 agrees with the expected difference between the mean breeding values of the immigrant group (`3`) and the founder population (`0`) that was specified in the simulation.

The predicted total additive genetic effects (**u**) are also reported in the "./wombat/ggNadivAstarRan/RnSoln\_id.dat" file. Here, the prediced total additive genetic effects include the genetic group effects expressed as deviations from 0.

#### 6.4.5.4 WOMBAT's genetic groups (Random group effects with Q from `nadiv`)

Fitting genetic group effects using **WOMBAT**'s own functionality still requires calculating and supplying **Q**. Although **Q** is supplied, genetic groups are included implicitly within the random effects. Further, the genetic group effects themselves are treated as random effects with a variance among group effects that is estimated by the model.

**WOMBAT** does not provide a function to calculate **Q**, so we calculate it using `nadiv`'s `ggcontrib()` function and provide <b>A<sup>-1</sup></b> as a generalized inverse. **WOMBAT** can be used to create <b>A<sup>-1</sup></b>, but here we will use the `makeAinv()` function in `nadiv`. For more details, see example 18 supplied with **WOMBAT**, particularly the pdf file entitled "RNote\_WOMBATS2Step.pdf" that is included within the Example 18 folder and the "wombat.par" file in the folder "~/Example18/C".

The code below creates <b>A<sup>-1</sup></b> for the `ggPed` pedigree and creates the data file and `id.codes` that contains **Q** (we also include the columns of **Q** in the data file for consistency with other versions of this data set).

``` r
Q <- ggcontrib(ggPed[, 1:3], ggroups = c("foc0", "g1"))
```

``` r
ggTutorial <- cbind(ggTutorial, Q)
```

``` r
AinvOut <- makeAinv(ggTutorial[, 1:3], det = TRUE)
listAstar <- sm2list(AinvOut$Ainv)
AinvLogdet <- AinvOut$logDet
```

Note above we told `makeAinv()` to calculate the log-determinant. **WOMBAT** requires this as the first line of a `.gin` (generalized inverse) file (note this is the log-determinant of the non-inverse matrix **A**, which `nadiv` calculates when forming the inverse - here <b>A<sup>-1</sup></b>). Now we need to write the data, <b>A<sup>-1</sup></b> list, and identity codes to files in formats that **WOMBAT** can use. Also, for consistency with other analyses, we will add a constant value of 10000 to the integer `id`s.

``` r
IDaddDF <- ggTutorial[, -c(2:7)]
IDaddDF[, 1] <- as.integer(as.character(IDaddDF[, 1])) + 10000
write.table(IDaddDF[, -c(2:3)], "./wombat/ggTutorial.d", col.names = FALSE,
    row.names = FALSE)

# Write the general inverse of A^-1 to file, note the re-ordering of first two columns
write.table(listAinv[, c(2,1,3)], "./wombat/ggWombatRan/id.gin",
    col.names = FALSE, row.names = FALSE)
  fConn <- file("./wombat/ggWombatRan/id.gin", "r+")
  Lines <- readLines(fConn)
# Add the log-determinant as the first line of the file
  writeLines(c(AinvLogdet, Lines), con = fConn)
  close(fConn) 

# Create mapping between order of identities and their unique codes
## including 2 genetic groups
## also add "dummy" column to indicate no genotypes available (1)
## Add Q matrix as final columns
id.codes <- cbind(seq(1, 6000), 10000 + seq(nrow(ggTutorial)), rep(1, 6000), round(10000*Q, 0))
write.table(id.codes, "./wombat/ggWombatRan/id.codes", col.names = FALSE,
    row.names = FALSE)
```

Note that **WOMBAT** wants the genetic group contributions in the `id.codes` file as an integer. Therefore, we have multiplied **Q** by a large number (10000 used above) and then rounded these to the nearest whole number. This scaling factor of 10000 will need to be indicated in the parameter file ("ggWombatRan.par"). The animal model in **WOMBAT** is then (see also the "./wombat/ggWombatRan" folder in the supporting files):

    COM WOMBATs genetic groups (random genetic group effects with Q from `nadiv`)
    DAT ../ggTutorial.d
      id 6000
      p
      is 2
      gen 15
      f
      foc0
      g1
    end
    ANAL UNI
    MODEL
      COV   f(1)
      SUBJ  id
      RAN   id     GIN
      RAN   ggrps  IDE
      TR    p
    END
    VAR id 1
     0.5
    VAR ggrps 1
     1.0
    VAR error 1
     0.5
    SPECIAL
      COVZER f(1) FIT
    END
    SPECIAL
      GENGROUPS ggrps 2 10000
    END

Here, genetic groups are specified as their own random effect in the model as `ggrps` with a diagonal covariance matrix specified (`IDE`). Thus, genetic group effects are conceptually treated as random effects (*Appendix S1*) in the model. Note the `SPECIAL` sections, which are necessary because the coefficient of inbreeding (*f*) has values of zero that are meaningful (i.e., these zeroes do not indicate missing values) and for fitting genetic groups. The genetic group `SPECIAL` section indicates the scaling factor by which we multiplied elements in **Q** that were supplied in the file "./wombat/ggWombatRan/id.codes".

The variance estimates in "./wombat/ggWombatRan/SumEstimates.out" match the simulated values for the `id` and residual terms. The estimated variance for the genetic group term (`ggrps`) is approximately `2`. All other simulated values are recovered in this analysis (model estimates match simulated values).

The genetic group effect predictions are reported in the "./wombat/ggWombatRan/RnSoln\_ggrps.dat" file. The genetic group predictions are random effects and thus are deviations from their expected value of zero. The difference between the genetic groups 3.111 agrees with the expected difference between the mean breeding values of the immigrant group (`3`) and the founder population (`0`) that was specified in the simulation.

The predicted total additive genetic effects (**u**) are reported in the "./wombat/ggWombatRan/RnSoln\_id.dat" file. Here, the prediced total additive genetic effects include the genetic group effects expressed as deviations from 0.

References: Appendix S6
-----------------------

Butler, D.G., Cullis, B.R., Gilmour, A.R. & Gogel, B.J. (2009) asreml: asreml() fits the linear mixed model. R package version 3.0. www.vsni.co.uk.

Davis, T.A. (2006) Direct Methods for Sparse Linear Systems. SIAM, Philadelphia.

Gilmour, A.R. (2005) Developments in utilizing pedigrees in genetic analysis within ASReml. Proceedings of the Association for the Advancement of Animal Breeding and Genetics, pp. 290–293.

Gilmour, A.R. (2010) Handling non positive definite relationship matrices in mixed models. 9th World Congress on Genetics Applied to Livestock Production, pp. 1–4. Leipzig, Germany.

Gilmour, A.R., Gogel, B.J., Cullis, B.R., Welham, S.J. & Thompson, R. (2014) ASReml 4.1 user guide.

Hadfield, J.D. (2010) MCMC methods for multi-response generalized linear mixed models: the MCMCglmm R package. Journal of Statistical Software, 33, 1–22.

Hadfield, J.D., Wilson, A.J., Garant, D., Sheldon, B.C. & Kruuk, L.E.B. (2010) The misuse of BLUP in ecology and evolution. American Naturalist, 175, 116–125.

Henderson, C.R. (1976) A simple method for computing the inverse of a numerator relationship matrix used in prediction of breeding values. Biometrics, 32, 69–83.

Kennedy, B.W., Schaeffer, L.R. & Sorensen, D.A. (1988) Genetic properties of animal models. Journal of Dairy Science, 71, 17–26.

Lynch, M. & Walsh, B. (1998) Genetics and Analysis of Quantitative Traits. Sinauer Associates, Incorporated, Sunderland, MA.

Meuwissen, T.H.E. & Luo, Z. (1992) Computing inbreeding coefficients in large populations. Genetics, Selection, Evolution, 24, 305–313.

Meyer, K. (2007) WOMBAT - a tool for mixed model analyses in quantitative genetics by restricted maximum likelihood (REML). Journal of Zhejiang University. Science. B, 8, 815–821.

Mrode, R.A. (2005) Linear Models for the Prediction of Animal Breeding Values, 2nd ed. CABI Publishing, Cambridge, MA.

Oikawa, T. & Yasuda, K. (2009) Inclusion of genetically identical animals to a numerator relationship matrix and modification of its inverse. Genetics, Selection, Evolution, 41, 25.

Quaas, R.L. (1976) Computing the diagonal elements and inverse of a large numerator relationship matrix. Biometrics, 32, 949–953.

Quaas, R.L. (1988) Additive genetic model with groups and relationships. Journal of Dairy Science, 71, 1338–1345.

Quaas, R.L. (1995) Fx algorithms. An unpublished note.

Reid, J.M. & Keller, L.F. (2010) Correlated inbreeding among relatives: occurrence, magnitude, and implications. Evolution, 64, 973–985.

Robinson, G.K. (1986) Group effects and computing strategies for models for estimating breeding values. Journal of Dairy Science, 69, 3106–3111.

Schaeffer, L.R. (1991) C. R. Henderson: Contributions to predicting genetic merit. Journal of Dairy Science, 74, 4052–4066.

Schaeffer, L.R. (1994) Multiple-country comparison of dairy sires. Journal of Dairy Science, 77, 2671–2678.

Schaeffer, L.R. (1999) Phantom parents. An unpublished note.

Sullivan, P.G. (1999) Estimating Genetic Variances and Covariances for Models with Genetic Groups. <http://cgil.uoguelph.ca/dcbgc/Agenda9909/agenda9909.htm>.

Westell, R.A., Quaas, R.L. & Van Vleck, L.D. (1988) Genetic groups in an animal model. Journal of Dairy Science, 71, 1310–1318.

Wolak, M.E. (2012) nadiv: an R package to create relatedness matrices for estimating non-additive genetic variances in animal models. Methods in Ecology and Evolution, 3, 792–796.

Wolak, M.E. & Keller, L.F. (2014) Dominance genetic variance and inbreeding in natural populations. Quantitative Genetics in the Wild (eds A. Charmantier, D. Garant & L.E.B. Kruuk), pp. 104–127. Oxford University Press, Oxford.

Wolak, M.E. & Reid, J.M. (2016) Data from: Accounting for genetic differences among unknown parents in microevolutionary studies: How to include genetic groups in quantitative genetic animal models. Journal of Animal Ecology, Dryad Digital Repository, <http://dx.doi.org/10.5061/dryad.jf7cr>.

<!----------------------------------------------------------------------------->
<!----------------------------------------------------------------------------->
<!----------------------------------------------------------------------------->
Software version information
============================

The following software versions were used to run all of the code above:

``` r
R.version.string
```

    ## [1] "R version 3.3.2 (2016-10-31)"

``` r
packageVersion("nadiv")
```

    ## [1] '2.14.3.1'

``` r
packageVersion("MCMCglmm")
```

    ## [1] '2.22.1'

``` r
packageVersion("asreml")
```

    ## [1] '3.0'

**WOMBAT:** Linux 64-bit version 12 June 2015, compiled using the `ifort` compiler with MKL library (threaded)

**ASReml:** ASReml 4.10 [28 Dec 2014] lr [18 Mar 2015], from asreml-4.1.0.1051-lr.64.tgz
