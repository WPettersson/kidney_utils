"""Handle solutions to the kidney exchange problem, aka maintain the list of
exchanges.
"""

class ExchangeExistsException(Exception):
    """Raised if we are replacing a cycle or chain that already exists in a
    solution.  Solutions are meant to be immutable.
    """
    pass

class Exchange(object):
    """Represents a set of exchanges between donor-patient pairs.
    """

    def __init__(self, pairs, weight=0):
        self._pairs = pairs
        self._weight = weight

    def __len__(self):
        return len(self._pairs)

class Cycle(Exchange):
    """Represents a cycle of exchanges."""

    def __init__(self, pairs, weight=0, backarc=False):
        super().__init__(pairs, weight)
        self._backarc = backarc

    def backarc(self):
        """Does this cycle have a backarc"""
        return self._backarc

class Chain(Exchange):
    """A chain of exchanges initiated by an altruistic donor."""

    def __len__(self):
        # The last pair is just a recipient, not donating
        return super().__len__() - 1

class Solution(object):
    """A set of kidney exchanges, as well as some metadata.
    """
    def __init__(self):
        """Basic constructor."""
        self._desc = None
        self._cycles = {}
        self._chains = {}

    def set_description(self, desc):
        """Set the description. This usually confers which optimality
        objectives were used to find the solution.
        """
        self._desc = desc

    def add_cycle(self, cycle_id, pairs, weight=0, backarc=False):
        """Adds a cycle of exchanges to the solution.

        :param cycle_id: The ID of the cycle
        :param pairs: A list of (donor, patient) pairs.
        :param weight: The total weight of the cycle, if known.
        :param backarc: Does this cycle have a backarc
        """
        if cycle_id in self._cycles:
            raise ExchangeExistsException()
        cycle = Cycle(pairs, weight, backarc)
        self._cycles[cycle_id] = cycle

    def add_chain(self, chain_id, pairs, weight=0):
        """Adds a chain to the solution.

        :param chain_id: The ID of the chain
        :param pairs: A list of (donor, patient) pairs. The first donor must be
        the altruistic donor, the last donor should be None.
        :param weight: The total weight of the chain, if known.
        """
        if chain_id in self._chains:
            raise ExchangeExistsException()
        chain = Chain(pairs, weight)
        self._chains[chain_id] = chain

    def size(self):
        """Count the number of kidney exchanges.

        :return: The number of exchanges.
        """
        return (sum([len(x) for x in self._cycles.values()]) +
                sum([len(x) for x in self._chains.values()]))

    def _n_way(self, length):
        """Return a list of n-way exchanges."""
        return [x for x in self._cycles.values() if len(x) == length]

    def three_ways(self):
        """Return a list of three-way exchanges.

        :return: A list of Cycle objects.
        """
        return self._n_way(3)

    def two_ways(self):
        """Return a list of two-way exchanges.

        :return: A list of Cycle objects.
        """
        return self._n_way(2)

    def _n_chains(self, length):
        """Return a list of chains of length n."""
        return [x for x in self._chains.values() if len(x) == length]

    def short_chains(self):
        """Return a list of short chains.

        :return: A list of short Chain objects.
        """
        return self._n_chains(2)

    def long_chains(self):
        """Return a list of long chains.

        :return: A list of long Chain objects.
        """
        return self._n_chains(3)

    def to_string(self):
        """Produce a human-readable summary of the exchange.

        :return: A summary.
        """
        string = str(self._desc) + "\n"
        string += "%d total transplants\n" % self.size()
        string += "%d two-way exchanges\n" % len(self.two_ways())
        string += "%d three-way exchanges\n" % len(self.three_ways())
        string += "%d short chains\n" % len(self.short_chains())
        string += "%d long chains\n" % len(self.long_chains())
        return string
