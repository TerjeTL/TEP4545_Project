#import numpy as np
#import subprocess

from operator import lt
from subprocess import call # To call bash scripts from this file
import numpy as np

CASE_DIR = "./../FoamProject/blockMeshVortexShedding/templateCase/" 

## SET-UP ##

NUM_CPU = 6

############

class Case:
    def __init__(self, case_name, time_step, length_to_thickness_ratio, cell_density, wall_distance):
        self.case_name = case_name
        self.length_to_thickness_ratio = length_to_thickness_ratio
        self.cell_density = cell_density
        self.wall_distance = wall_distance
        self.time_step = time_step

    def run(self):
         # Fully clean up before making a new run
        rc = call([CASE_DIR + "cleanAll.sh", "-y"])

        # Run the test case
        rc += call([CASE_DIR + "runAll.sh", "-y", "-p", str(NUM_CPU), "-wh", str(self.wall_distance),
            "-dt", str(self.time_step), "-lt", str(self.length_to_thickness_ratio),
            "-dh", str(self.cell_density) ])

        rc += call([CASE_DIR + "exportResults.sh", "-a", "-f", "export_data", "-d", self.case_name])

        return rc

def main():
    lt_list = [5.0, 6.0, 8.0]
    wd_list = [3.0, 1.0]

    for i in range(len(lt_list)):
        for j in range(len(wd_list)):
            case_name = "LT"+str(lt_list[i])+"/"+"WD"+str(wd_list[j])
            test_case = Case(case_name, 1.2e-4, lt_list[i], 2.25e3, wd_list[j])
            test_case.run()

if __name__ == '__main__':
    main()
