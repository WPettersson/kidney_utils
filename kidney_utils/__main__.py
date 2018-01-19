"""If called as a module, print the solution as a string for each file in the command line."""

import sys
from .readers import read_solution

def print_solutions(files):
    """Print each solution."""
    for filename in files:
        solution = read_solution(filename)
        print(solution.to_string())

if __name__ == "__main__":
    print_solutions(sys.argv[1:])
