# How to include genetic groups in quantitative genetic animal models

This repository contains a version controlled, editable, and commentable [Appendix S6](WolakAndReid_GenGroupHowTo_SI6/WolakAndReid_GenGroupHowTo_SI6.md) of the Supporting Information to *Accounting for genetic differences among unknown parents in microevolutionary studies: How to include genetic groups in quantitative genetic animal models* ([Wolak & Reid. 2017. Journal of Animal Ecology](http://onlinelibrary.wiley.com/doi/10.1111/1365-2656.12597/full)). Appendix S6 is the heart of the tutorial part of the Supporting Information (SuppInfo), containing the code and discussion of just how to practically go about fitting the models discussed in the manuscript.

The manuscript and SuppInfo were both thoroughly peer reviewed, however, these published versions represent a snapshot in time for the methods discussed. Methods and code are never static, particularly as software changes over time. This repository provides an alternative resource that can be changed while also tracking and keeping a history of every change made along the way.


## Quick inventory of the repository

   - **WolakAndReid_GenGroupHowTo_SI6.md** and **WolakAndReid_GenGroupHowTo_SI6.html** are the files to be viewed

   - **WolakAndReid_GenGroupHowTo_SI6.Rmd** is the actual file that creates the SuppInfo. Any further edits or changes will occur in this file. This is rendered using the `rmarkdown` [CRAN package](https://CRAN.R-project.org/package=rmarkdown) (on a command line: `R -e "rmarkdown::render('WolakAndReid_GenGroupHowTo_SI6.Rmd')"`) to create the `.md` document to be viewed on GitHub. 

     - Any changes should be done to the `.Rmd` file. With each new version of this file, it should be re-rendered.

   - **WolakAndReid_GenGroupHowTo_SuppInfo.pdf** is the original SuppInfo accompanying the journal article.

   - Folders **asreml**, **ASRemlR**, **MCMCglmm**, and **wombat** contain the model/outputs for each of these software programs.

   - **MCMCglmm_modelDiet.R** is how I compressed the original `MCMCglmm` saved model objects (available from [dryad](http://www.datadryad.org/resource/doi:10.5061/dryad.jf7cr))


## Changes

For ease of reference, significant changes to be noted below. Tag with commits or issues, where appropriate.

### Major

### Minor
    - 22 November 2016: Initialize files for GitHub
      - Formatting to remove LaTex etc. within the `.Rmd` file that were used to create a nice looking *pdf*
      - Compress `MCMCglmm` models (see *./MCMCglmm_modelDiet.R*)
