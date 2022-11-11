#import numpy as np
#import subprocess

from subprocess import call # To call bash scripts from this file

CASE_DIR = "./../FoamProject/blockMeshVortexShedding/data-export-test/" 

## SET-UP ##

NUM_CPU = 6

############

class Case:
    def __init__(self, case_name, length_to_thickness_ratio, cell_density, wall_distance):
        self.case_name = case_name
        self.length_to_thickness_ratio = length_to_thickness_ratio
        self.cell_density = cell_density
        self.wall_distance = wall_distance

    def run(self):
         # Fully clean up before making a new run
        rc = call([CASE_DIR + "cleanAll.sh", "-y"])

        # Run the test case
        rc += call([CASE_DIR + "runAll.sh", "-y", "-p", str(NUM_CPU)])

        rc += call([CASE_DIR + "exportResults.sh", "-n", self.case_name])

        return rc

def main():
    for i in range(5):
        case_name = "test_case" + str(i)
        test_case = Case(case_name, 5, 100, 10)
        test_case.run()

if __name__ == '__main__':
    main()
