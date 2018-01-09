"""Print some parameters of input graphs."""

import sys
import os
from kidney_utils.readers import read_xml, read_json

for filename in sys.argv[1:]:
    filename = filename.rstrip()
    graph = None
    if filename[-4:] == ".xml":
        graph = read_xml(filename)
    elif filename[-5:] == ".json":
        graph = read_json(filename)
    if graph:
        density = graph.density()
        diameter = graph.diameter()
        radius = graph.radius()
        name = os.path.basename(filename)
        print("%s,%.3f,%.3f,%.3f" % (name, density, diameter, radius))
