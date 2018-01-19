"""Test cases for reading instance files."""

from kidney_utils.readers import read_instance, read_solution
from nose.tools import eq_, ok_

def almost_equal(float_a, float_b, tolerance=1e-6):
    """Check if two floating point numbers are almost-equal."""
    return abs(float_a - float_b) < tolerance

def test_read_test_one():
    """Run tests on test1.json"""
    instance = read_instance("tests/testfiles/test1.json")
    eq_(instance.size(), 4)
    eq_(instance.radius(), 2)
    eq_(instance.diameter(), 3)
    hop_plot = [0.3125, 0.75, 1.0]
    instance_hop_plot = instance.hop_plot()
    for index, value in enumerate(hop_plot):
        ok_(value, instance_hop_plot[index])

def test_read_test_two():
    """Run tests on test2.json"""
    instance = read_instance("tests/testfiles/test2.json")
    eq_(instance.size(), 6)
    eq_(instance.radius(), 3)
    eq_(instance.diameter(), 5)
    hop_plot = [0.27777777777777773, 0.6666666666666666, 0.8333333333333334, 0.9444444444444445, 1.0]
    instance_hop_plot = instance.hop_plot()
    for index, value in enumerate(hop_plot):
        ok_(value, instance_hop_plot[index])

def test_read_test_three():
    """Run tests on test3.json"""
    instance = read_instance("tests/testfiles/test3.json")
    eq_(instance.size(), 3)
    eq_(instance.radius(), 2)
    eq_(instance.diameter(), 3)
    hop_plot = [0.4444444444444444, 0.8888888888888888, 1.0]
    instance_hop_plot = instance.hop_plot()
    for index, value in enumerate(hop_plot):
        ok_(value, instance_hop_plot[index])

def test_read_solution_one():
    """Run tests on test1.xml"""
    solution = read_solution("tests/testfiles/test1.xml")
    eq_(solution.size(), 2)
    eq_(len(solution.two_ways()), 1)
    eq_(len(solution.three_ways()), 0)
    eq_(len(solution.three_ways_with_backarcs()), 0)
    eq_(len(solution.short_chains()), 0)
    eq_(len(solution.long_chains()), 0)

def test_read_solution_two():
    """Run tests on test2.xml"""
    solution = read_solution("tests/testfiles/test2.xml")
    eq_(solution.size(), 6)
    eq_(len(solution.two_ways()), 3)
    eq_(len(solution.three_ways()), 0)
    eq_(len(solution.three_ways_with_backarcs()), 0)
    eq_(len(solution.short_chains()), 0)
    eq_(len(solution.long_chains()), 0)

def test_read_solution_three():
    """Run tests on test1.xml"""
    solution = read_solution("tests/testfiles/test3.xml")
    eq_(solution.size(), 3)
    eq_(len(solution.two_ways()), 0)
    eq_(len(solution.three_ways()), 1)
    eq_(len(solution.three_ways_with_backarcs()), 1)
    eq_(len(solution.short_chains()), 0)
    eq_(len(solution.long_chains()), 0)
