#!/bin/sh

#SBATCH --job-name=nightlights
#SBATCH --output=slurmlogs/nightlights%j.out
#SBATCH --error=slurmlogs/nightlights%j.err
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20g
#SBATCH --mail-type=END,FAIL

module load geos/3.7.0
module load gdal/2.3.2
module load udunits/2.2
module load proj/4.9
module load R

Rscript calls.R > calls.log
