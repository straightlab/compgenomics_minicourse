Now that have analyzed our data using Kracken, we can plot the results using R. We can run R on sherlock using an interative session that can run on your web browser. 

We have put together an R script to plot the results from Kracken. Copy this file into your notebooks directory that you created in notebook 01
cp copies a file or directory from one location to another. It is formatted as:
cp source_path target_path

cp /scratch/groups/astraigh/genomics_minicourse/shared/notebooks/mock_analysis.Rmd /scratch/groups/astraigh/genomics_minicourse/<Charles>/notebooks



Activate an R studio session on sherlock by navigating to this link:
https://login.sherlock.stanford.edu/pun/sys/dashboard/batch_connect/sessions

Click on the RStudio session on the left side of the page
-Change the following parameters:
-R version to 3.6.1
-partition to astraigh
-cpus to 1
-runtime to 2 hours

Click Launch

Click 'Connect to R studio' when it loads (this will take a few minutes to load)
***This will open a new browser window. If it does not load within a minute or so, try these steps again on a private browsing window***

Click File>OpenFile
copy the path to your notebooks dir into the navigation bar on the top:
/scratch/groups/astraigh/genomics_minicourse/kelsey/notebooks/

select mock_analysis.Rmd

click open

Here you can run the pre-generated code by selecting a cell that contains the code you want to run and hit command + shift + enter


