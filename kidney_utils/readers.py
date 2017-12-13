"""Utilities to read and inspect files containing kidney exchange data."""

from defusedxml import ElementTree as ET
import json
from .kidney_graph import Graph


class ReadException(Exception):
    """An error occurred when reading a file."""
    pass


def read_xml(filename):
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
                raise ReadException()
            source = sources[0].text
        matches = donor.find("matches")
        if matches:
            for match in matches:
                target = match.find("recipient").text
                score = float(match.find("score").text)
                graph.add_edge(source, target, score)
    return graph


def read_json(filename):
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
                raise ReadException()
            source = donor["sources"][0]
            for match in donor["matches"]:
                target = match["recipient"]
                score = float(match["score"])
                graph.add_edge(source, target, score)
    return graph
