# Connection to Sherlock and familiarization with the linux shell 

We will do all of our bioinformatics analysis on the Stanford High Performance Computer Cluster Sherlock (www.sherlock.stanford.edu). In this part, we will illustrate the basic workflow for connecting to the cluster and getting everything ready for a computational project. The steps are:
- ssh into to the cluster
- setup up our working directory
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
cd $PI_SCRATCH/genomics_minicourse
```
Note that `$PI_SCRATCH` is a bash variable, which contains the path to a default "scratch" folder for our lab. You can see the content of this variable with
```
echo $PI_SCRATCH
```

Make a new folder for yourself, and move to that folder. For example for team CKO
```
mkdir -p teamCKO
cd teamCKO
```

The `-p` flag is usefull in case the folder aready exists.

This folder is currently empty. You can list the content of the current folder `.` with 
```
ls .
```

You can always see your current location in the filesystem with 
```
pwd
```

Now let's create a few subfolders to organize our work. We want our project directory (the team directory in that case) to look like this
```text
teamCKO
├── data
│   └── woyke_mockcommunity
└── notebooks
```

The data folder will contain analysis for specific datasets, arranged into subfolders. The first dataset we will look at is a mock bacterial community from Woyke et al. so we'll prepare a folder for it.

Make these directories with 
```bash
mkdir -p data/woyke_mockcommunity notebooks
```
Verify the tree scructure with 
```bash
tree .
```

## Running persistent / recoverable sessions with Gnu Screen

What happens if I am in the middle of some task on Sherlock and I loose internet connection or close my computer? To avoid having to back to square one, we need to set up a persistent bash session. The standard way to do this is using a window manager such as GNU Screen.

```bash
screen -S genomics
```

This creates a new session called genomics. Let's demomstrate what it does by putting a mark in our terminal, leaving and coming back

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
screen -r
```

## Getting computational resources on Sherlock

When we ssh into sherlock, we automatically obtain a shell on a "login node". So far, we've ran all of our commands on a login node. You can see that from the look of the prompt `[username@sh01-ln03 login ~]` Login nodes are not meant to do any heavy computations. You can think of them as entry points into Sherlock, and places to do very simple operations, like parsing/creating directories, viewing files, etc...  

For anything more computationally heavy than that, you need to request dedicated resources (in fact, if you run a command on a login node that takes too long to complete or requests too much memory, it will automatically get aborted). The way to ask for resources on Sherlock is through the resource manager and job scheduler Slurm (https://slurm.schedmd.com/documentation.html). There are several way to ask for resources with Slurm depending on whether you want to run an interactive job (where you keep entering commands), or want to submit a job that will run without your input when resources become available. This is explained in "Running jobs on Sherlock" found in the Additional Resources section below ``

For this bootcamp, we are going to request 2 cpus each for 3h. 
`srun -p astraigh --time=03:00:00 --cpus-per-task=2 --pty bash`

The last part of this command (--pty bash) is to create an interactive job, more specically we are requesting a bash shell.

You should quickly a prompt that looks like that `srun: job <ID> queued and waiting for resources` 
Your command line will then read: [username@sh02-09n13 /scratch/groups/astraigh/genomics_minicourse/username_dir]$
This shell is running on a dedicated computational node (here sh02-09n13). Within this shell, you'll have access to 2 CPUS and 32gB of RAM.

This computational node is part of the parition `astraigh` (specificed by the `-p astraigh` flag) which is reserved for our lab. There are other partitions you can use, see ().

## Loading packages for our subsequent analysis

We preinistalled some bioinformatics tools we're going to use during this bootcamp inside an anaconda environment. 

First add the anaconda path to your bash profile

```
echo -e "PATH=/share/PI/astraigh/miniconda3/bin:$PATH" >>~/.bash_profile
```

Load this environment with
```
source activate bootcamp
```


## Additional resources

- Running jobs on Sherlock : https://www.sherlock.stanford.edu/docs/user-guide/running-jobs/
- advanced connection options : https://www.sherlock.stanford.edu/docs/advanced-topics/connection/
- GNU screen : https://www.howtoforge.com/linux_screen




