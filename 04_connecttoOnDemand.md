# Start a remote Rstudio session on Sherlock

We can run Rstudio on sherlock using an onDeman session that can run on your web browser. 

We have put together an R script to plot the results from Kracken. Copy this file into your notebooks directory that you created in notebook 01. We use the `cp` command which is formatted as `cp <source_path> <target_path>`

```bash
# don't forget to modify teamCKO to your own folder name
cp /scratch/groups/astraigh/genomics_minicourse/shared/notebooks/mock_analysis.Rmd /scratch/groups/astraigh/genomics_minicourse/teamCKO/notebooks
```


Start an R studio session on sherlock by navigating to this link:
https://login.sherlock.stanford.edu/pun/sys/dashboard/batch_connect/sessions

1. Click on the RStudio session on the left side of the page. Change the following parameters:
     -R version to 3.6.1
     -partition to astraigh
     -cpus to 1
     -runtime to 2 hours
     
2.Click Launch

3.Click 'Connect to R studio' when it loads (this will take a few minutes to load)
***This will open a new browser window. If it does not load within a minute or so, try these steps again on a private browsing window***

4.Click File>OpenFile. Copy the path to your notebooks dir into the navigation bar on the top: `/scratch/groups/astraigh/genomics_minicourse/teamCKO/notebooks/`. Select `mock_analysis.Rmd`


Here you can run the pre-generated code by selecting a cell that contains the code you want to run and hit command + shift + enter


