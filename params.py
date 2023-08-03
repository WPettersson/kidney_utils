"""Print some parameters of input graphs."""

import sys
import os
from multiprocessing import Pool
from kidney_utils.readers import read_instance

def calc_params(filename):
    graph = read_instance(filename)
    if graph:
        density = graph.density()
        diameter = 0 #  graph.diameter()
        radius = 0 #  graph.radius()
        size = graph.size()
        edges = graph.edge_count()
        cycles = graph.find_cycles()
        three_cycles = len([x for x in cycles if len(x) == 3])
        two_cycles = len([x for x in cycles if len(x) == 2])
        name = os.path.basename(filename)
        tw = graph.approx_treewidth()
        scc = graph.strongly_connected_components()
        scc_sizes = [len(s) for s in scc]
        return (name, size, edges, density, two_cycles, three_cycles, diameter, radius, tw, scc_sizes)
    return None


with Pool() as pool:
    results = [pool.apply_async(calc_params, (arg.rstrip(),)) for arg in sys.argv[1:]]
    print("name,size,edge_count,density,two_cycles,three_cycles,diameter,radius,approx tree-width, strongly connected component sizes")
    for res in results:
        out = res.get()
        if out is not None:
            (name, size, edges, density, two_cycles, three_cycles, diameter, radius, tw, scc_sizes) = out
            print("%s,%s,%s,%.3f,%.3f,%.3f" % (name, size, edges, density, two_cycles, three_cycles, diameter, radius, tw, scc_sizes))
        else:
            print("Error on some file?")
