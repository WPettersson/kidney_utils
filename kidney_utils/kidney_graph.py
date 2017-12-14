"""A directed graph representing kidney patients and donors."""

from Queue import Queue
import snap


class Vertex(object):
    """A vertex in a directed graph."""
    def __init__(self, index):
        if not type(index) is int:
            raise TypeError()
        self._name = "V%d" % index
        self._index = index
        self._leaving = []
        self._entering = []

    def add_out_edge(self, edge):
        """Add an edge leaving this vertex."""
        self._leaving.append(edge)

    def add_in_edge(self, edge):
        """Add an edge entering this vertex."""
        self._entering.append(edge)

    def edges_in(self):
        """The edges that leave this vertex."""
        return self._entering

    def edges_out(self):
        """The edges that leave this vertex."""
        return self._leaving

    def neighbours_out(self):
        """Return the list of neighbours when leaving this vertex."""
        return [edge.tail() for edge in self.edges_out()]

    def neighbours_in(self):
        """Return the list of neighbours which entering this vertex."""
        return [edge.head() for edge in self.edges_in()]

    def index(self):
        """The (integer) index of this vertex"""
        return self._index

    def __str__(self):
        return self._name

    def __repr__(self):
        return self.__str__()


class Edge(object):
    """An edge in a directed graph."""
    def __init__(self, v1, v2, weight=1):
        self._v1 = v1
        self._v2 = v2
        self._weight = weight

    def head(self):
        """The start of the edge."""
        return self._v1

    def tail(self):
        """The end of the edge."""
        return self._v2

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
        self._snap_graph = snap.TNGraph.New()
        self._eccentricity = None
        self._shortest_paths = None

    def add_edge(self, vert_1, vert_2, weight=1):
        """Add an edge from vert_1 to vert_2 with weight."""
        # Mark eccentricity and shortest-paths as not known
        self._eccentricity = None
        self._shortest_paths = None
        if vert_1 not in self._vertices:
            self._vertices[vert_1] = Vertex(vert_1)
            self._snap_graph.AddNode(vert_1)
        if vert_2 not in self._vertices:
            self._vertices[vert_2] = Vertex(vert_2)
            self._snap_graph.AddNode(vert_2)
        edge = Edge(self._vertices[vert_1], self._vertices[vert_2], weight)
        self._edges.append(edge)
        self._vertices[vert_1].add_out_edge(edge)
        self._vertices[vert_2].add_in_edge(edge)
        self._snap_graph.AddEdge(vert_1, vert_2)

    def vertex(self, index):
        """Get a Vertex from an index."""
        return self._vertices[index]

    def calculate_shortest_paths(self):
        """For each pair of distinct vertices u and v, calculate the shortest
        path between u and v."""
        size = len(self._vertices)
        self._shortest_paths = [[[] for x in range(size)] for x in range(size)]
        for destination in self._vertices.values():
            self.calculate_shortest_paths_to(destination)

    def calculate_shortest_paths_to(self, dest):
        """Calculate all shortest paths to dest."""
        to_check = Queue()
        d_ind = dest.index()
        # The neighbours have one hop to this vertex
        for near in dest.neighbours_in():
            self._shortest_paths[near.index()][d_ind] = [dest]
            for to_do in near.neighbours_in():
                to_check.put((to_do, near))
        while not to_check.empty():
            (source, next_vert) = to_check.get()
            s_ind = source.index()
            if self._shortest_paths[s_ind][d_ind]:
                continue
            path = [next_vert]
            path.extend(self._shortest_paths[next_vert.index()][d_ind])
            self._shortest_paths[s_ind][d_ind] = path
            for to_do in source.neighbours_in():
                to_check.put((to_do, source))

    def calculate_eccentricity(self):
        """Calculate the eccentricity of each vertex. For a vertex v, this is
        the maximum of shortest_path(v,u) over all other vertices u."""
        if not self._shortest_paths:
            self.calculate_shortest_paths()
        self._eccentricity = [max([len(path) for path in paths]) for paths in
                              self._shortest_paths]

    def density(self):
        """Get the density of the graph."""
        num_verts = len(self._vertices)
        if num_verts < 2:
            return 0
        return float(len(self._edges))/(num_verts*(num_verts-1))

    def out_degree_dist(self):
        """Get the out-degree distribution of the graph"""
        snap_results = snap.TIntPrV()
        snap.GetOutDegCnt(self._snap_graph, snap_results)
        return [[x.GetVal1(), x.GetVal2()] for x in snap_results]

    def in_degree_dist(self):
        """Get the in-degree distribution of the graph"""
        snap_results = snap.TIntPrV()
        snap.GetInDegCnt(self._snap_graph, snap_results)
        return [[x.GetVal1(), x.GetVal2()] for x in snap_results]

    def diameter(self):
        """Get the diameter of the graph"""
        if not self._eccentricity:
            self.calculate_eccentricity()
        return max(self._eccentricity)

    def radius(self):
        """Get the radius of the graph"""
        if not self._eccentricity:
            self.calculate_eccentricity()
        return min(self._eccentricity)

    def __str__(self):
        return "Graph on %d nodes" % len(self._vertices)
