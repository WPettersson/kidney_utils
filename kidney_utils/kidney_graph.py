"""A directed graph representing kidney patients and donors."""


class Vertex(object):
    """A vertex in a directed graph."""
    def __init__(self, name):
        self._name = name
        self._leaving = []
        self._entering = []

    def add_out_edge(self, edge):
        """Add an edge leaving this vertex."""
        self._leaving.append(edge)

    def add_in_edge(self, edge):
        """Add an edge entering this vertex."""
        self._entering.append(edge)


class Edge(object):
    """An edge in a directed graph."""
    def __init__(self, v1, v2, weight=1):
        self._v1 = v1
        self._v2 = v2
        self._weight = weight

    def __str__(self):
        name = "%s->%s" % (self._v1, self._v2)
        if self._weight != 1:
            name += " (%.2f)" % (self._weight)
        return name

    def __repr__(self):
        return self.__str__()


class Graph(object):
    """A directed graph representing kidney patients and donors."""
    def __init__(self):
        self._edges = []
        self._vertices = {}

    def add_edge(self, v1, v2, weight=1):
        """Add an edge from v1 to v2 with weight."""
        edge = Edge(v1, v2, weight)
        self._edges.append(edge)
        if v1 not in self._vertices:
            self._vertices[v1] = Vertex(v1)
        if v2 not in self._vertices:
            self._vertices[v2] = Vertex(v2)
        self._vertices[v1].add_out_edge(edge)
        self._vertices[v2].add_in_edge(edge)

    def get_density(self):
        """Get the density of the graph."""
        num_verts = len(self._vertices)
        if num_verts < 2:
            return 0
        return len(self._edges)/(num_verts*(num_verts-1))
