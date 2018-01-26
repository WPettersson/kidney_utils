"""An instance of a kidney exchange problem.
"""

from kidney_utils.graph import Graph

class Instance(Graph):
    """An instance of a kidney exchange problem. Currently, just a graph
    really.
    """

    def to_string(self):
        """Create some sort of string representation of this instance.
        """
        string = "KidneyExchange instance with %d vertices\n" % self.size()
        string += "Treewidth <= %d\n" % self.approx_treewidth()
        return string
