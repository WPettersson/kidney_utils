"""Basic install script."""

from setuptools import setup, find_packages

setup(
    name="kidney_utils",
    author="William Pettersson",
    author_email="william.pettersson@glasgow.ac.uk",
    version="0.1",
    license="GPL",
    keywords="kidney kidney_exchange",
    packages=find_packages(),
    install_requires=["defusedxml", "progressbar33", "networkx",
                      "python-louvain"]
)
