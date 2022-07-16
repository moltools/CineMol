#!/usr/bin/env python3
from matplotlib.pyplot import install_repl_displayhook
from setuptools import setup, find_packages

import cinnamol.config as config 


setup (
    name="cinnamol",
    version=config.__version__,
    author="David Meijer",
    author_email="david.meijer@wur.nl",
    description=config.__description__,
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=["pyvista"]
)
