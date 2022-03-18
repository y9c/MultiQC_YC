#!/usr/bin/env python

import logging

from multiqc.utils import config
from pkg_resources import get_distribution

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_helab_version = get_distribution("multiqc_helab").version
log.info("Running MultiQC HeLab Plugin v{}".format(config.multiqc_helab_version))


def update_config() -> None:
    """Update MultiQC config object.

    * Update module order
    * Disable unnecessary modules to avoid duplicate data
    * Update search patterns
    """

    log.debug("HeLab - Updating config")
    # Add module to module order
    config.module_order.extend(
        [
            {"sampletracking": {"module_tag": ["DNA", "RNA"]}},
            {"demultiplex": {"module_tag": ["DNA", "RNA", "Demultiplex"]}},
        ]
    )

    # Move module to the top
    config.top_modules.extend(["sampletracking", "demultiplex"])

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
    ## Sampletracking
    if "sampletracking/crosscheckfingerprints" not in config.sp:
        # add new key
        config.update_dict(
            config.sp,
            {"sampletracking/crosscheckfingerprints": {"contents": "CrosscheckFingerprints", "shared": False}},
        )
        # overwrite old key
        config.update_dict(
            config.sp, {"picard/crosscheckfingerprints": {"fn": "nonexistent", "shared": False}},
        )


def update_fn_cleanup() -> None:
    """
    Update filename cleanup patterns
    :return: None
    """
    config.fn_clean_exts.extend(["_duplicate_metrics", "_samtools_stats", "_samtools_idxstats", ".unaligned"])
