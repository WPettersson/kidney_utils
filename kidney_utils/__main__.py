"""If called as a module, print the solution as a string for each file in the command line."""

import sys
from .readers import read_solution, read_instance, ReadException

def print_things(files):
    """Print each thing."""
    for filename in files:
        try:
            solution = read_solution(filename)
            print(solution.to_string())
        except ReadException:
            # Probably not a solution.
            try:
                instance = read_instance(filename)
                print(instance.to_string())
            except ReadException:
                # Uuh .. no idea
                print("Could not do anything with %s" % filename)

if __name__ == "__main__":
    print_things(sys.argv[1:])
