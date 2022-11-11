#import numpy as np
#import subprocess

from subprocess import call # To call bash scripts from this file

CASE_DIR = "./../FoamProject/blockMeshVortexShedding/data-export-test/" 

## SET-UP ##

NUM_CPU = 6

############

def Case:
    def __init__(self, case_name, length_to_thickness_ratio, cell_density):
        self.case_name = case_name
        self.length_to_thickness_ratio = length_to_thickness_ratio
        self.cell_density = cell_density

    def run(self):
         # Fully clean up before making a new run
        rc = call([CASE_DIR + "cleanAll.sh", "-y"])

        # Run the test case
        rc += call([CASE_DIR + "runAll.sh", "-y", "-p", str(NUM_CPU)])

        rc += call([CASE_DIR + "exportResults.sh", "-n", self.case_name])

        return rc

if __name__ == '__main__':

    test_case = Case("test_data_case", 5, 100)
    test_case.run()

    sys.exit(main())  # next section explains the use of sys.exit