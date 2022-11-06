#!/usr/bin/env python
"""MultiQC_HeLab is a plugin for MultiQC containing customized modules and
templates."""

import site
import sys

from setuptools import find_packages, setup

version = "0.0.1"
site.ENABLE_USER_SITE = "--user" in sys.argv[1:]

setup(
    name="multiqc_helab",
    version=version,
    author="Matthias De Smet",
    author_email="11850640+matthdsm@users.noreply.github.com",
    description="MultiQC plugin for interal use at helab",
    long_description=__doc__,
    keywords="bioinformatics",
    url="https://github.com/CenterForMedicalGeneticsGhent/MultiQC_helab",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["multiqc>=1.10", "pandas"],
    entry_points={
        "multiqc.hooks.v1": [
            "config_loaded = multiqc_helab.multiqc_helab:update_config",
        ],
        "multiqc.modules.v1": [
            "fastq_screen_fork = multiqc_helab.modules.fastq_screen_fork.fastq_screen_fork:MultiqcModule",
            "sampletracking = multiqc_helab.modules.sampletracking.sampletracking:MultiqcModule",
            "picard_demultiplex = multiqc_helab.modules.picard_demultiplex.demultiplex:MultiqcModule",
        ],
        "multiqc.templates.v1": [
            "helab = multiqc_helab.templates.helab",
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
