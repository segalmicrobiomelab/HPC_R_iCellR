#!/bin/bash
#SBATCH --job-name=R.work.icellr					# Job name
#SBATCH --mail-type=BEGIN,END,FAIL			       		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=benjamin.wu@nyumc.org				# Where to send mail
#SBATCH --partition=cpu_short						# GPU 
#SBATCH --ntasks=8	             					# CPU (total) each node = 40 CPUs
#SBATCH --nodes=2							# node
#SBATCH --mem=300gb                     				# Job memory request
#SBATCH --output=serial_test_%j.log   					# Standard output and error log
#SBATCH --time=6:00:00							# Time allotted 


### November 4 2021
### Loading R 
### module installing

module load r/4.1.1
module load pandoc/2.2.3.2
module load slurm

##### /gpfs/scratch/wub02/projects/msq.singlecell.chronicmurine.true/2021-04-30/
##### /gpfs/scratch/wub02/projects/msq.singlecell.chronicmurine.true/

cd /gpfs/scratch/wub02/projects/msq.singlecell.chronicmurine.true/

##### Creating Assembly/Classifiers for 2021.04
##### /gpfs/scratch/wub02/projects/milan.elastase.run/gg_13_8_otus

# /gpfs/home/wub02/R/

Rscript --vanilla /gpfs/scratch/wub02/projects/msq.singlecell.chronicmurine.true/code/R/211110_icellr.R
