#!/usr/bin/env python
"""MultiQC_YC is a plugin for MultiQC containing customized modules and
templates."""

import site
import sys

from setuptools import find_packages, setup

version = "0.0.1"
site.ENABLE_USER_SITE = "--user" in sys.argv[1:]

setup(
    name="multiqc_yc",
    version=version,
    author="Chang Ye",
    description="MultiQC plugin for interal use at yc",
    long_description=__doc__,
    keywords="bioinformatics",
    url="https://github.com/y9c/MultiQC_readsStats",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["multiqc>=1.10"],
    entry_points={
        "multiqc.modules.v1": [
            "readsStats = multiqc.modules.readsStats:MultiqcModule",
        ],
        "multiqc.templates.v1": [
            "yc = multiqc_yc.templates.yc",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Environment :: Web Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
