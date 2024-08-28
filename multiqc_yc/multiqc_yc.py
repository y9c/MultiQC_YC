#!/usr/bin/env python

import logging

from multiqc.utils import config
from pkg_resources import get_distribution

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_yc_version = get_distribution("multiqc_yc").version
log.info(
    "Running MultiQC (Chang Y forked) Plugin v{}".format(config.multiqc_yc_version)
)


def update_config() -> None:
    """Update MultiQC config object.

    * Update module order
    * Disable unnecessary modules to avoid duplicate data
    * Update search patterns
    """

    log.info("Multiqc (Chang Y forked) - Updating config")
    # Add module to module order
    config.module_order.extend([])

    # Move module to the top
    config.top_modules.extend([])

    # Disable module to avoid duplicate data
    disabled_modules = []
    for module in disabled_modules:
        del config.avail_modules[module]

    # Update search pattern
    update_search_patterns()

    # Update fn cleanup
    update_fn_cleanup()


def update_search_patterns() -> None:
    """
    Update search patterns
    Overwrite default search pattern and set 'shared' to false to avoid running the module from core mqc
    :return: None
    """
    # readsStats
    if "readsStats" not in config.sp:
        # add new key
        config.update_dict(
            config.sp,
            {
                "readsStats/cutadapt": {
                    "contents": "Cutadapt report",
                    "fn": "*.json",
                    "num_lines": 10,
                    "shared": False,
                },
                "readsStats/bowtie2": {
                    "contents": "reads; of these:",
                    "fn": "*_genes.report",
                    "num_lines": 1,
                    "shared": False,
                },
                "readsStats/star": {
                    "contents": "       Mapping speed, Million of reads per hour |",
                    "num_lines": 4,
                    "shared": False,
                },
                "readsStats/dedup": {
                    "contents": "total (QC-passed reads + QC-failed reads)",
                    "num_lines": 1,
                    "shared": False,
                },
            },
        )


def update_fn_cleanup() -> None:
    """
    Update filename cleanup patterns
    :return: None
    """
    config.fn_clean_exts.extend(
        [
            "_duplicate_metrics",
            "_samtools_stats",
            "_samtools_idxstats",
            ".unaligned",
        ]
    )
