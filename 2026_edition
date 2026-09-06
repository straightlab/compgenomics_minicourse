# Connection to Sherlock and familiarization with the linux shell 

We will do all of our bioinformatics analysis on the Stanford High Performance Computer Cluster Sherlock (www.sherlock.stanford.edu). In this part, we will illustrate the basic workflow for connecting to the cluster and getting everything ready for a computational project. The steps are:
- ssh into the cluster
- set up our working directory
- create a persistent session on the cluster to come back to 
- request computational resources
- load necessary packgages


At the end of this notebook you will know how to :
- navigate the linux filesystem with `ls`, `cd`, `mkdir`, `pwd`
- set up a nicely organized folder for the project
- use a Gnu Screen session
- use SLURM to access to computational resources


## Connection to Sherlock
In your browser go to <login.sherlock.stanford.edu> , then click on `>Shell`
Alternatively, you can ssh into sherlock using a terminal app on your computer. On a Mac, you can use the native terminal app Term. Open up Term (Command-space term), then
```
ssh <username>@login.sherlock.stanford.edu
```


## Setting up our workspace for the project.
Go to the group minicourse directory
```
cd $GROUP_SCRATCH/biochem_minicourse_2025
```
Note that `$GROUP_SCRATCH` is a bash variable, which contains the path to a default "scratch" folder for our lab. You can see the content of this variable with
```
echo $GROUP_SCRATCH
```

Make a new folder for yourself, and move to that folder. For example for team Straight
```
mkdir <your_dir_name>
cd <your_dir_name>
```

This folder is currently empty. You can list the content of the current folder with 
```
ls 
```

You can always see your current location in the filesystem with  
```
pwd
```

Now let's create a few subfolders to organize our work. We want our project directory (the team directory in this case) to look like this
```text
<your_dir_name>
├── data
│   └── samples
|   └── external
└──tmp
```

- The `tmp` folder is a scratch folder which we will use to put temporary files (some programs create temporary files and require a dedicated folder for these files).
- The `data` folder will contain analysis for specific datasets, arranged into subfolders.

Make these directories with 
```bash
mkdir data

mkdir data/samples data/external

mkdir tmp
```
Verify the tree structure with 
```bash
tree .
```

Finally, please grant read/write/execute access to your directory to everyone in the group (this will be important after the minicourse, for when we need to delete the files we generate today). Any directories you make in Sherlock will automatically grant read/execute access to everyone, but only the creator of the directory will have write access.  Read access means someone can "read" the file, execute access means someone can run a program with the file, and write access means someone can edit (or delete) the file. When we type:

```bash
ls -l .
```

This will display the permissions for the directory you are currently in, which should look like this: drwxr-sr-x. "d" stands for directory, and that's followed by the read/write/execute (rwx or rws) permissions for the directory owner, the group members (of astraigh), and anyone who uses Sherlock (in that order). 

To change permissions on a directory or file, you need to use the "chmod" command. This command is then followed by "-R" to indicate we want to change permissions recursively (for every file/directory nested within a certain directory). Then, the numerical value below is an "octal value" that represents what permissions you want to grant to the directory owner, group, and everyone. Check out this article for more details on where the numbers come from: https://www.redhat.com/en/blog/manage-permissions. Type the following: 

```bash
chmod -R 775 $GROUP_SCRATCH/biochem_minicourse_2025/<your_directory>
```

Now, when you type 
```bash
ls -l .
```

The permissions should read: drwxrwsr-x

## Running persistent / recoverable sessions with Gnu Screen

What happens if I am in the middle of some task on Sherlock and I lose internet connection or close my computer? To avoid having to go back to square one, we need to set up a persistent bash session. The standard way to do this is using a window manager such as GNU Screen.

```bash
screen -S genomics
```

This creates a new session called genomics. Let's demonstrate what it does by putting a mark in our terminal, leaving and then coming back

```bash
# leave a mark
echo "I am here"
# keep track of the node number in the prompt (for example sh01-ln03)
#the screen session will only be accessible from this login node

# close your terminal, then relogin into sherlock
ssh sunetid@login.sherlock.stanford.edu

# if you're not assigned the same login node as before, connect to the original one. If it's the same skip this step
ssh sh01-ln03

# get back to you persistent session
screen -x
```

## Getting computational resources on Sherlock

When we ssh into sherlock, we automatically obtain a shell on a "login node". So far, we've ran all of our commands on a login node. You can see that from the look of the prompt `[username@sh01-ln03 login ~]`. Login nodes are not meant to do any heavy computations. You can think of them as entry points into Sherlock, and places to do very simple operations, like parsing/creating directories, viewing files, etc...  

For anything more computationally heavy than that, you need to request dedicated resources (in fact, if you run a command on a login node that takes too long to complete or requests too much memory, it will automatically get aborted). The way to ask for resources on Sherlock is through the resource manager and job scheduler Slurm (https://slurm.schedmd.com/documentation.html). There are several ways to ask for resources with Slurm depending on whether you want to run an interactive job (where you keep entering commands), or want to submit a job that will run without your input when resources become available. This is explained in "Running jobs on Sherlock" <https://www.sherlock.stanford.edu/docs/user-guide/running-jobs/>.

For this bootcamp, we are going to request 2 cpus each for 3h. 
`salloc -p astraigh --time=03:00:00 --cpus-per-task=2`

You should quickly get a prompt that looks like that `salloc: Pending job allocation <ID>` 
Your command line will then read: [username@sh02-09n13 /scratch/groups/astraigh/biochem_minicourse_2025/<your_dir_name>]$
This shell is running on a dedicated computational node (here sh02-09n13). Within this shell, you'll have access to 2 CPUS and 32 GB of RAM.

This computational node is part of the partition `astraigh` (specificed by the `-p astraigh` flag) which is reserved for our lab. There are other partitions you can use, refer to the Sherlock documentation.

## Loading packages for our subsequent analysis

We preinstalled some bioinformatics tools we're going to use during this bootcamp inside a conda environment. Conda is a unix package manager that helps maintain your software packages in order, and an environment is an isolated set of software packages that are required for a task or computational pipeline.

First add the conda path to your bash profile with these two commands:

```
echo -e "PATH=/home/groups/astraigh/miniconda3_py38_4.9.2/bin:$PATH" >> ~/.bashrc
export PATH=/home/groups/astraigh/miniconda3_py38_4.9.2/bin:$PATH
```

Load this environment with
```
conda activate minicourse
#or you can try:
source activate minicourse
```
(alternatively conda activate /home/groups/astraigh/miniconda3_py38_4.9.2/envs/minicourse)

## Managing your SLURM jobs

You can see the jobs currently running on the cluster with 
```
squeue -u $USER
```

You can cancel a running job with 
```
scancel <JOBID>
```
where the JOBID is obtained from `squeue`

## Additional resources

- Running jobs on Sherlock : https://www.sherlock.stanford.edu/docs/user-guide/running-jobs/
- advanced connection options : https://www.sherlock.stanford.edu/docs/advanced-topics/connection/
- GNU screen : https://www.howtoforge.com/linux_screen
