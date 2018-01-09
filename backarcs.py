"""Print density of input graphs."""

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
        print("%s: %.2f%%" % (os.path.basename(filename), 100*graph.backarc()))
