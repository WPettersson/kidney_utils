"""Test cases for the Graph class."""

from kidney_utils.graph import Graph
from nose.tools import eq_


def test_treewidth_complete_graphs():
    """Test treewidth algorithm, which is approximate and not exact.
    """

    def test_kn(size):
        """Test on complete graphs."""
        graph = Graph()
        for one in range(size):
            for two in range(one + 1, size):
                graph.add_edge(one, two)
        eq_(size-1, graph.approx_treewidth())
    for size in range(2, 6):
        test_kn(size)


def test_treewidth():
    """Test treewidth algorithm, which is approximate and not exact.
    """
    graph = Graph()
    for one, two in [(1, 2), (2, 3), (3, 4), (2, 4), (1, 5)]:
        graph.add_edge(one, two)
    eq_(2, graph.approx_treewidth())
    graph = Graph()
    for one, two in [(1, 2), (2, 3), (3, 4), (2, 4), (1, 5), (5, 6), (5, 7),
                     (5, 8), (6, 7), (6, 8), (7, 8)]:
        graph.add_edge(one, two)
    eq_(3, graph.approx_treewidth())

def test_strong_connected_component():
    """Test the strongly connected components implementation.
    """
    graph = Graph()
    for one, two in [(1, 2), (2, 3), (3, 1)]:
        graph.add_edge(one, two)
    scc = graph.strongly_connected_components()
    eq_(1, len(scc))
    eq_(3, len(scc[0]))

    graph = Graph()
    for one, two in [(1, 2), (2, 3), (3, 4), (4, 3)]:
        graph.add_edge(one, two)
    scc = graph.strongly_connected_components()
    eq_(3, len(scc))
    # Tarjan's is deterministic, so we should be able to ignore the ordering
    # here, but if this test fails, check the order of the lengths
    eq_([2, 1, 1], [len(c) for c in scc])

def test_cycles():
    """Test the enumeration of simple cycle finding in directed graphs."""
    graph = Graph()
    for one, two in [(1, 2), (2, 3), (3, 1)]:
        graph.add_edge(one, two)
    cycles = list(graph.find_cycles())
    eq_(len(cycles), 1)
    eq_(cycles[0], [1, 2, 3])

def test_groups():
    """Test the detection of similar groups of vertices."""
    graph = Graph()
    for one, two in [(1, 2), (2, 3), (1, 4), (4, 3), (3, 1)]:
        graph.add_edge(one, two)
    groups = graph.group()
    eq_(len(groups), 3)
