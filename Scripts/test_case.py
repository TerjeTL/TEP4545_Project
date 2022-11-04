#import numpy as np
#import subprocess

from subprocess import call # To call bash scripts from this file

CASE_DIR = "./../FoamProject/blockMeshVortexShedding/data-export-test/" 

## SET-UP ##

NUM_CPU = 6

############

# Fully clean up before making a new run
rc = call([CASE_DIR + "cleanAll.sh", "-y"])

# Run the test case
rc += call([CASE_DIR + "runAll.sh", "-y", "-p", str(NUM_CPU)])

# Export data from run
name = "test_data_export.csv"
rc += call([CASE_DIR + "exportResults.sh", "-n", name])


print("done from python!")

print("exit code:", rc)
