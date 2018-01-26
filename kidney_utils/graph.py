"""A directed graph representing kidney patients and donors."""

try:
    from queue import Queue
except ImportError:
    from Queue import Queue

import community
import networkx
from progressbar import ProgressBar, Bar, Percentage

class Vertex(object):
    """A vertex in a directed graph."""
    def __init__(self, desc, index):
        self._desc = desc
        self._leaving = []
        self._entering = []
        self._index = index

    def desc(self):
        """Return a description of this vertex. This is the key used by the
        graph to refer to this vertex.
        """
        return self._desc

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

    def index(self):
        """The (integer) index of this vertex in the graph."""
        return self._index

    def neighbours_out(self):
        """Return the list of neighbours when leaving this vertex."""
        return [edge.head() for edge in self.edges_out()]

    def neighbours_in(self):
        """Return the list of neighbours which entering this vertex."""
        return [edge.tail() for edge in self.edges_in()]

    def all_neighbours(self):
        """Return all neighbours, for instance if treating this directed graph
        as an undirected graph.
        """
        return self.neighbours_in() + self.neighbours_out()

    def __str__(self):
        return "V%s" % (self._desc)

    def __repr__(self):
        return self.__str__()


class Edge(object):
    """An edge in a directed graph."""
    def __init__(self, v1, v2, weight=1):
        self._v1 = v1
        self._v2 = v2
        self._weight = weight

    def tail(self):
        """The start of the edge."""
        return self._v1

    def head(self):
        """The end of the edge."""
        return self._v2

    def weight(self):
        """The weight of the edge."""
        return self._weight

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
        self._vertices_list = []
        self._eccentricity = None
        self._shortest_paths = None
        self._nxgraph = networkx.Graph()

    def size(self):
        """Size, aka number of vertices."""
        return len(self._vertices)

    def add_vertex(self, vert):
        """Add a vertex to the graph."""
        if vert not in self._vertices.keys():
            self._vertices[vert] = Vertex(vert, len(self._vertices_list))
            self._vertices_list.append(self._vertices[vert])

    def add_edge(self, vert_1, vert_2, weight=1):
        """Add an edge from vert_1 to vert_2 with weight."""
        # Mark eccentricity and shortest-paths as not known
        self._eccentricity = None
        self._shortest_paths = None
        self.add_vertex(vert_1)
        self.add_vertex(vert_2)
        edge = Edge(self._vertices[vert_1], self._vertices[vert_2], weight)
        self._edges.append(edge)
        self._vertices[vert_1].add_out_edge(edge)
        self._vertices[vert_2].add_in_edge(edge)
        if weight != 1:
            self._nxgraph.add_edge(vert_1, vert_2, weight=weight)
        else:
            self._nxgraph.add_edge(vert_1, vert_2)

    def vertex(self, index):
        """Get a Vertex from an index."""
        return self._vertices[index]

    def vertex_list(self):
        """The list of Vertex objects in this graph"""
        return self._vertices_list

    def edge_count(self):
        """Number of edges in the graph."""
        return len(self._edges)

    def edge_list(self):
        """List of Edge objects in this graph"""
        return self._edges

    def calculate_shortest_paths(self, quiet=True):
        """For each pair of distinct vertices u and v, calculate the shortest
        path between u and v."""
        if self._shortest_paths:
            return
        size = len(self._vertices)
        self._shortest_paths = [[[] for _ in range(size)] for __ in range(size)]
        count = len(self._vertices)
        if not quiet:
            pbar = ProgressBar(widgets=[Bar(), Percentage()], maxval=count)
            pbar.start()
            count = 0
        for destination in self._vertices.values():
            self.calculate_shortest_paths_to(destination)
            if not quiet:
                count += 1
                pbar.update(count)
        if not quiet:
            pbar.finish()

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
        #TODO
        pass

    def in_degree_dist(self):
        """Get the in-degree distribution of the graph"""
        #TODO
        pass

    def diameter(self):
        """Get the diameter of the graph"""
        if not self._eccentricity:
            self.calculate_eccentricity()
        return max(self._eccentricity)

    def radius(self, ignore_loops=True):
        """Get the radius of the graph"""
        if not self._eccentricity:
            self.calculate_eccentricity()
        if ignore_loops:
            return min([dist for dist in self._eccentricity if dist != 0])
        return min(self._eccentricity)

    def hop_plot(self):
        """Calculate how much of the graph can be reached in a given number of
        hops. The returned values are the average of the reach for each
        individual vertex.
        :return: A list of tuples (a, b) where a is the number of hops,
        and b is the proportion of the graph that can be reached in that many
        hops.
        """
        per_vertex = self.hop_plot_vertices()
        # Make sure the list corresponding to each vertex is the same length
        longest = max([len(x) for x in per_vertex])
        results = [0] * longest
        for vert in per_vertex:
            for (distance, proportion) in vert:
                results[distance-1] += proportion
            for dist in range(len(vert), longest):
                results[dist] += 1.0
        results = [result/len(self._vertices) for result in results]
        return results

    def hop_plot_vertices(self):
        """Calculate how much of the graph can be reached in a given number of
        hops, for each vertex.

        :return: A list of list of tuples. One list of tuples is returned for
        each vertex. The list contains tuples (a, b) where a is the number of
        hops and b is the proportion of the graph that can be reached from this
        vertex in a hops.
        """
        results = []
        num_verts = len(self._vertices)
        self.calculate_shortest_paths()
        for index, vert in enumerate(self._vertices):
            verts_done = 0
            distance = 1
            result_here = []
            paths = self._shortest_paths[index]
            while verts_done < num_verts:
                verts_done += len([x for x in paths if len(x) == distance])
                result_here.append((distance, float(verts_done)/num_verts))
                distance += 1
            results.append(result_here)
        return results

    def best_partition(self):
        """Get the best possible partition of this graph into clusters. Uses
        the louvain method described in Fast unfolding of communities in large
        networks, Vincent D Blondel, Jean-Loup Guillaume, Renaud Lambiotte,
        Renaud Lefebvre, Journal of Statistical Mechanics: Theory and
        Experiment 2008(10), P10008 (12pp) and the python package from
        https://github.com/taynaud/python-louvain/
        """
        return community.best_partition(self._nxgraph)

    def adjacency(self):
        """Return the adjacency matrix of this graph."""
        matrix = []
        for _, vert in self._vertices.items():
            here = [0.] * len(self._vertices)
            for edge in vert.edges_out():
                here[edge.head().index()] += edge.weight()
            for edge in vert.edges_in():
                here[edge.tail().index()] += edge.weight()
            matrix.append(here)
        return matrix

    def backarc(self):
        """Return the proportion of edges that have a back-arc.
        Note that a value of 1 indicates that every edge has a back-arc.
        """
        count = 0
        edges = []
        for edge in self._edges:
            here = [edge.tail(), edge.head()]
            edges.append(here)
            if [edge.head(), edge.tail()] in edges:
                count += 1
        return (count / 2) / len(self._edges)

    def approx_treewidth(self):
        """Computes an approximation of the treewidth of this graph, using the
        Greedy Fill-In algorithm.
        """
        decomp = Graph()
        # Make decomp a copy of self. We modify this to find the treewidth.
        # Then we throw it away, even though it's almost a complete tree
        # decomposition, because I don't yet need that bit.
        for edge in self._edges:
            decomp.add_edge(edge.tail().desc(), edge.head().desc(), edge.weight())
        for vert in self.vertex_list():
            decomp.add_vertex(vert.desc())
        used = [False] * decomp.size()
        width = 0
        def next_vert():
            """Get the next vertex to add to a bag, according to Greedy Fill-In
            """
            least_added = -1
            to_return = None
            for vert in decomp.vertex_list():
                if used[vert.index()]:
                    continue
                num_added = 0
                bag_size = 1
                neighbours = vert.all_neighbours()
                for vert_a in neighbours:
                    if used[vert_a.index()]:
                        continue
                    bag_size += 1
                    for vert_b in neighbours:
                        if used[vert_b.index()] or vert_b in vert_a.all_neighbours():
                            continue
                        num_added += 1
                if least_added == -1 or num_added < least_added:
                    least_added = num_added
                    to_return = vert
                    this_bag_size = bag_size
            return to_return, this_bag_size
        for _ in range(self.size()):
            best_choice, bag_size = next_vert()
            if bag_size > width + 1:
                width = bag_size - 1
            used[best_choice.index()] = True
            for vert in decomp.vertex_list():
                if not used[vert.index()] and vert in best_choice.all_neighbours():
                    for other in decomp.vertex_list():
                        if vert == other:
                            continue
                        if used[other.index()] or other not in best_choice.all_neighbours():
                            continue
                        if other not in vert.all_neighbours():
                            decomp.add_edge(vert.desc(), other.desc())
        return width

    def __str__(self):
        return "Graph on %d nodes" % len(self._vertices)
