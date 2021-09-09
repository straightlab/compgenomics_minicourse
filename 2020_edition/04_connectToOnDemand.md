# Start a remote Rstudio session on Sherlock

We can run Rstudio on sherlock using an onDemand session that can run on your web browser. 

Before we do start Rstudio, we need to download some R libraries we'll need. It is best to do this in a command-line version of R.

## Running R from the command line

R is not by default accessible on Sherlock but needs to be loaded first. It is part of a "module", just like a lot of software on Sherlock (https://www.sherlock.stanford.edu/docs/software/list/)

On a compute node terminal, load the R module and the start R with
```bash
# load module
ml R/3.6

#start R
R
```

This should open up an R console. In R, install the following packages: rmarkdown, dplyr, tidyr, ggplot2 with

```R
install.packages(c('rmarkdown', 'dplyr', 'tidyr', 'ggplot2'))
```
Once this is done, quit R with
```R
quit()
```
## RStudio

Start an R studio session on sherlock by navigating to this link:
https://login.sherlock.stanford.edu/pun/sys/dashboard/batch_connect/sessions

1. Click on the RStudio session on the left side of the page. Change the following parameters:
    - R version to 3.6.1
    - partition to astraigh
    - cpus to 1
    - runtime to 2 hours
2. Click Launch
3. Click 'Connect to R studio' when it loads (this will take a few minutes to load)
***This will open a new browser window. If it does not load within a minute or so, try these steps again on a private browsing window***


## Open the template notebook

We have put together an R notebook to plot the results from Kracken. Using the bash terminal, copy this file into your notebooks directory that you created in notebook 01. We use the `cp` command which is formatted as `cp <source_path> <target_path>`

```bash
# don't forget to modify teamCKO to your own folder name
cp /scratch/groups/astraigh/genomics_minicourse/shared/notebooks/mock_analysis.Rmd /scratch/groups/astraigh/genomics_minicourse/teamCKO/notebooks
```

Back to Rstudio, Click File>OpenFile. Copy the path to your notebooks dir into the navigation bar on the top: `/scratch/groups/astraigh/genomics_minicourse/teamCKO/notebooks/`. Select `mock_analysis.Rmd`


Here you can run the pre-generated code by selecting a cell that contains the code you want to run and hit command + shift + enter


