#!/bin/sh
# 
# A tool for quickly creating and submitting array jobs for submission to the cluster.
# This script will run the same simulation with multiple parameter sets. 
# The simulation must read a single .csv file with one header and one row as an input,
# which would usually be called as follows:  ./mysimulation input.csv
# Note this script runs c++ simulations by default compiled with g++. The arrayjob template contains the line
# module load gcc. If you are using python or something else, you will want to load the corresponding module here
# by replacing this line.
# 
# This script is used to apply the simulation to multiple input files, the rows of which 
# are stored in an input list csv file. 
# Each row of the input list csv contains the input file contents for a single run.
# This code splits each row of the input list csv into input.csv files in separate subdirectories labelled from
# 1 to N, where N is the number of rows in the input list excluding the header.
# Then an array job submission script is created and called to submit these jobs to the cluster on our private queue.
# 
# Syntax for calling this script is as follows:
# ./MakeArrayJob.sh path_to_my_simulation path_to_my_list_of_inputs cluster_walltime_in_hrs cluster_walltime_in_mins
# 
# e.g. ./MakeArrayJob.sh ./simulation ./inputlist.csv 0 30
# The above would run the simulation over the rows in inputlist.csv as an array job with walltime of 30 mins
#
# Best to call this function in the directory containing inputlist.csv. The simulation can be elsewhere.
# Jordan Juritz, Jan 2021

sim=$1
inputlist=$2
walltimehrs=$3
walltimemins=$4
walltimehrspost=$5
walltimeminspost=$6

# Get current directory
maindirectory=$(pwd)

# Get absolute path to the simulation
sim="$(cd "$(dirname "$sim")"; pwd)/$(basename "$sim")"

# Get absolute path to the inputlist
inputlist="$(cd "$(dirname "$inputlist")"; pwd)/$(basename "$inputlist")"

echo "Generating array job script to run $sim on the rows of $inputlist..."

# Get date
dt=$(date '+%d/%m/%Y %H:%M:%S');

# Get number of parameter sets. 
numparams=$(($(wc -l < $inputlist )-1 ))

# Get the header of the list of input files. The header will be copied into each input file.
header=$(awk 'NR == 1' $inputlist )

echo "Creating $numparams subdirectories and input files..." 
# For each distinct input file...
for ((i=1;i<=$numparams;i++));
do 	
# If the directory with label i doesn't exist
	if [ ! -d "./$i" ]
	then
		# make the directory and create an input file in it, with the header
		mkdir "./$i"
		echo $header >> ./$i/input.csv
		mkdir "./$i/figures"
	fi
done

# Run through the list of input files and copy the relevant line into each input.csv in each subdirectory.
awk -F, -vOFS="," 'NR>1{ 
print $0 >> "./"(NR-1)"/input.csv"
}' $inputlist

echo "Building array job script and requesting ${walltimehrs}h ${walltimemins}m wall-time..."  
# Now we build the arrayjob.sh file which is used for submission
# First if the hrs or minutes requested are a single character, append a 0 before them.
if [ ${walltimehrs} -lt 10 ]; then 
walltimehrs=0$walltimehrs 
fi
if [ ${walltimemins} -lt 10 ]; then 
walltimemins=0$walltimemins 
fi

# The following constructs the ArrayJob.sh script which is used for submission
# By default we request 1gb of memory 1 node and 1 cpu
# A timeout command is called with the execution line.
# The timeout stops the code 2 mins before the walltime limit is reached.
# The code is executed within each parameter's sub directory, so input files stay with the corresponding outputs.
{
        echo "#!/bin/sh"
        echo "#PBS -lwalltime="$walltimehrs":"$walltimemins":00"
        echo "#PBS -lselect=1:ncpus=1:mem=1gb"
        echo "#PBS -J 1-"$numparams
        echo ""
        echo "module load gcc"  
        echo "maindirectory="$maindirectory
        echo "subdir=""$""{maindirectory}/""$""PBS_ARRAY_INDEX/"
        echo "cd ""$""subdir"
        echo 'runcommand="'$sim" ""$"'{subdir}/input.csv"'
	echo "timeout "$(( $walltimemins+60*$walltimehrs -2))"m ""$""runcommand"
} >> ArrayJob.sh

#{
#    echo "#!/bin/sh"
#    echo "#PBS -lwalltime="$walltimehrspost":"$walltimeminspost":00"
#    echo "#PBS -lselect=1:ncpus=1:mem=1gb"
#    echo ""
#    echo "module load anaconda3/personal"
#    echo "maindirectory="$maindirectory
#    echo "cd ""$""maindirectory"
#    echo 'runcommand="python main.py"'
#    echo "timeout "$(( $walltimeminspost+60*$walltimehrspost -2))"m ""$""runcommand"
#} >> PostJob.sh

echo "Submitting array job..."
# Finally we submit the array job to our private queue
firstjob=$(qsub ${maindirectory}/ArrayJob.sh)
#qsub -W depend=afterany:$firstjob $maindirectory/PostJob.sh

echo "Done. type qstat -t to check the status of your array job"
