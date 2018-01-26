"""Test cases for the Graph class."""

from kidney_utils.graph import Graph
from nose.tools import eq_


def test_treewidth_complete_graphs():
    """Test treewidth algorithm, which is approximate and not exact.
    """
    def test_kn(n):
        graph = Graph()
        for one in range(n):
            for two in range(one + 1, n):
                graph.add_edge(one, two)
        eq_(n-1, graph.approx_treewidth())
    for n in range(2, 6):
        test_kn(n)


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
