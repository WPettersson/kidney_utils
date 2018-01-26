"""Print some parameters of input graphs."""

import sys
import os
from kidney_utils.readers import read_instance


print("name,size,edge_count,density,diameter,radius")
for filename in sys.argv[1:]:
    filename = filename.rstrip()
    graph = read_instance(filename)
    if graph:
        density = graph.density()
        diameter = graph.diameter()
        radius = graph.radius()
        size = graph.size()
        edges = graph.edge_count()
        name = os.path.basename(filename)
        print("%s,%s,%s,%.3f,%.3f,%.3f" % (name, size, edges, density, diameter, radius))
