
from __future__ import with_statement
from __future__ import absolute_import
from subprocess import check_output
import re
import sys
from io import open

def submitJob_local(index, commandExecutable):
   
	RUN_FILENAME = 'myrun'
	JOB_NAME = 'Sm1Fe12_-{}'.format(index)

	# Step 1
	# SBATCH -t 06:00:00
	# module load intel/mkl-11.2.3 vasp/vasp-5.4.4 mpi/impi-5.0.3
	# srun --ntasks=1 --nodes=1 {3} & 

	myrun_content = '''#!/bin/sh
#SBATCH -J {0}
#SBATCH -o log.output
#SBATCH -e log.error
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
. /opt/intel/oneapi/setvars.sh --force

env
{1}
'''.format(JOB_NAME, commandExecutable)


	with open(RUN_FILENAME, 'w') as fp:
		fp.write(myrun_content)

	# Step 2
	output = check_output('sbatch {0}'.format(RUN_FILENAME), 
		shell=True, universal_newlines=True)


	# Step 3
	# Here we parse job ID from the output of previous command
	jobNumber = int(output.split(' ')[3])
	return jobNumber


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', dest='index', type=int)
	parser.add_argument('-c', dest='commandExecutable', type=str)
	args = parser.parse_args()

	jobNumber = submitJob_local(index=args.index, commandExecutable=args.commandExecutable)
	print('<CALLRESULT>')
	print(int(jobNumber))

	