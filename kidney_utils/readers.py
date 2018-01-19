"""Utilities to read and inspect files containing kidney exchange data."""

import json
from defusedxml import ElementTree as ET
from .kidney_graph import Graph
from .solutions import Solution

class ReadException(Exception):
    """An error occurred when reading a file."""
    pass

def read_instance(filename):
    """Read an instance of a kidney exchange problem.

    :param filename: The file containing the instance. It must end in either
    .json or .xml, depending on the input format.
    :returns graph: The Graph object which represents this instance.
    :raises ReadException: Raised if the input file is not recognised or
    malformed.
    """
    if filename[-5:] == ".json":
        return read_json_instance(filename)
    if filename[-4:] == ".xml":
        return read_xml_instance(filename)
    raise ReadException("Unknown file format")

def read_solution(filename):
    """Reads a solution from a file.

    :param filename: The file containing the solution. It must end in either
    .json or TODO, depending on the input format.
    :returns solution: The Solution
    :raises ReadException: Raised if the input file is not recognised or
    malformed.
    """
    if filename[-4:] == ".xml":
        return read_xml_solution(filename)
    raise ReadException("Unknown file format")

def read_xml_instance(filename):
    """Read an XML file containing kidney exchange data.
    """
    xml_data = ET.parse(filename)
    donors = xml_data.getroot()#['data']
    graph = Graph()
    for donor in donors:
        sources = donor.find("sources")
        if not sources:
            donor_id = donor.attrib["donor_id"]
            source = "alt_%s" % donor_id
        else:
            sources = list(sources)
            if len(sources) != 1:
                raise ReadException("Only donors with exactly 1 source are supported")
            source = sources[0].text
        matches = donor.find("matches")
        if matches:
            for match in matches:
                target = match.find("recipient").text
                score = float(match.find("score").text)
                graph.add_edge(source, target, score)
    return graph


def read_json_instance(filename):
    """Read a JSON file containing kidney exchange data.
    """
    with open(filename, "r") as infile:
        json_data = json.load(infile)
    graph = Graph()
    for index, donor in json_data["data"].items():
        if "altruistic" in donor:
            source = "alt_%s" % (index)
            for match in donor["matches"]:
                target = match["recipient"]
                score = float(match["score"])
                graph.add_edge(source, target, score)
        else:
            if len(donor["sources"]) != 1:
                raise ReadException("Only donors with exactly 1 source are supported")
            source = donor["sources"][0]
            for match in donor["matches"]:
                target = match["recipient"]
                score = float(match["score"])
                graph.add_edge(source, target, score)
    return graph

def read_xml_solution(filename):
    """Reads a solution from a file.

    :param filename: The file containing the solution.
    :returns solution: The Solution
    :raises ReadException: Raised if the input file is malformed.
    """
    xml_data = ET.parse(filename)
    root = xml_data.getroot()#['data']
    solution = Solution()
    def get_xml_text(element):
        """Get the text within an element"""
        return "".join([x for x in element.itertext()])
    output = root.find("output")
    cycles = {}
    for cycle in output.find("all_cycles"):
        cycle_id = int(cycle.get('id'))
        weight = float(cycle.get('weight'))
        backarc = int(cycle.get('backarcs', 0))
        altruistic = bool(cycle.get('altruistic', None))
        pairs = []
        for pair in cycle.findall('pair'):
            try:
                patient = int(get_xml_text(pair.find("p")))
            except AttributeError:
                # Altruistic, no "first" patient
                patient = -1
            donor = int(get_xml_text(pair.find("d")))
            pairs.append((donor, patient))
        cycles[cycle_id] = (weight, backarc, altruistic, pairs)
    sol = output.find('exchange_data').find('entry')
    solution.set_description(get_xml_text(sol.find('description')))
    for cycle in sol.find('exchanges'):
        cycle_id = int(get_xml_text(cycle))
        (weight, backarc, altruistic, pairs) = cycles[cycle_id]
        if altruistic:
            solution.add_chain(chain_id=cycle_id, pairs=pairs)
        else:
            solution.add_cycle(cycle_id=cycle_id, pairs=pairs, weight=weight,
                               backarc=backarc)
    return solution
