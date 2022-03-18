""" MultiQC submodule to parse output from Picard CrosscheckFingerprints """

import logging
import re
from collections import OrderedDict
from distutils.util import strtobool
from itertools import chain

import pandas as pd
from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import heatmap, table

# Initialize the logger
log = logging.getLogger("multiqc")

# This is a subset, the rest of the fields are self descriptive
FIELD_DESCRIPTIONS = {
    "LEFT_SAMPLE": "The name of the left sample.",
    "LEFT_GROUP_VALUE": "The name of the left data-type group.",
    "RIGHT_SAMPLE": "The name of the right sample.",
    "RIGHT_GROUP_VALUE": "The name of the right data-type group.",
    "RESULT": "The categorical result of comparing the calculated LOD score against the threshold.",
    "DATA_TYPE": "The datatype used for the comparison.",
    "LOD_SCORE": "Log10 of the probability that the samples come from the same individual.",
    "LOD_SCORE_TUMOR_NORMAL": "LOD score with the assumption that Left is a Tumor.",
    "LOD_SCORE_NORMAL_TUMOR": "LOD score with the assumption that Right is a Tumor.",
    "LOD_THRESHOLD": "The LOD threshold used for this pairwise comparison.",
    "TUMOR_AWARENESS": "Whether or not this pairwise comparison was flagged for tumor awareness",
}


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        super(MultiqcModule, self).__init__(
            name="Sampletracking",
            anchor="sampletracking",
            info=" - Sampletracking QC based on output from Picard CrosscheckFingerprints",
        )

        report_count = self.parse_reports()
        if report_count > 0:
            log.info(f"Found {report_count} Sampletracking reports")

    def parse_reports(self):
        """
        Find Picard CrosscheckFingerprints reports and parse their data.
        """

        self.picard_CrosscheckFingerprints_df = None
        report_dfs = []
        for f in self.find_log_files("sampletracking/crosscheckfingerprints", filehandles=True):
            report_dfs.append(pd.read_csv(f["f"], sep="\t", comment="#"))
        if len(report_dfs) > 0:
            self.picard_CrosscheckFingerprints_df = pd.concat(report_dfs, ignore_index=True)
            self.fingerprints_table()
            self.lod_heatmap()
        return len(report_dfs)

    def lod_heatmap(self):
        cats, hm_data = self._lod_heatmap_data()

        if len(hm_data) > 0:
            self.add_section(
                name="LOD heatmap",
                anchor="sampletracking-lod-heatmap",
                description="Heatmap of LOD score per comparison",
                plot=heatmap.plot(
                    hm_data,
                    xcats=cats,
                    pconfig={
                        "id": "sampletracking-lod-heatmap-plot",
                        "title": "Sampletracking: LOD Heatmap",
                        "min": -5,
                        "max": 5,
                        "square": False,
                    },
                ),
            )

    def _lod_heatmap_data(self):
        hm_df = self.picard_CrosscheckFingerprints_df.pivot(
            index="LEFT_SAMPLE", columns="RIGHT_SAMPLE", values="LOD_SCORE"
        )
        cats = hm_df.index.tolist()
        hm_data = hm_df.values.tolist()

        return [cats, hm_data]

    def fingerprints_table(self):
        # clean data
        ## remove expexted mismatches
        table_df = self.picard_CrosscheckFingerprints_df[
            self.picard_CrosscheckFingerprints_df["RESULT"] != "EXPECTED_MISMATCH"
        ]

        ## drop rows that contain samples that should be ignored
        for idx, row in table_df.iterrows():
            if self.is_ignore_sample(row["LEFT_SAMPLE"]) or self.is_ignore_sample(row["RIGHT_SAMPLE"]):
                table_df.drop(idx, inplace=True)

        table_data = table_df.to_dict(orient="index")

        if len(table_df) > 0:
            self.add_section(
                name="Crosscheck Fingerprints",
                anchor="sampletracking-crosscheckfingerprints",
                description="Pairwise identity checking between samples.",
                helptext="""
                Checks that all data in the set of input files comes from the same individual, based on the selected group granularity.
                Expected Mismatches have been omitted for clarity.
                """,
                plot=table.plot(
                    table_data,
                    self._get_table_headers(table_data),
                    {
                        "namespace": "sampletracking",
                        "id": "sampletracking-crosscheckfingerprints-table",
                        "table_title": "Sampletracking: Crosscheck Fingerprints",
                        "save_file": True,
                        "col1_header": "ID",
                        "no_beeswarm": True,
                    },
                ),
            )

    def _take_till(self, iterator, fn):
        """Take from an iterator till `fn` returns false.
    
        Returns the iterator with the value that caused false at the front, and all the lines skipped till then as a list.
        """
        headers = []
        try:
            val = next(iterator)
            while fn(val):
                headers.append(val)
                val = next(iterator)
        except StopIteration:
            return ()

        return (chain([val], iterator), headers)

    def _parse_cli(self, line):
        """Parse the Picard CLI invocation that is stored in the header section of the file."""
        tumor_awareness_regex = r"CALCULATE_TUMOR_AWARE_RESULTS(\s|=)(\w+)"
        lod_threshold_regex = r"LOD_THRESHOLD(\s|=)(\S+)"

        tumor_awareness = None
        lod_threshold = None

        tumor_awareness_match = re.search(tumor_awareness_regex, line)
        if tumor_awareness_match is not None:
            tumor_awareness = strtobool(tumor_awareness_match.group(2))

        lod_threshold_match = re.search(lod_threshold_regex, line)
        if lod_threshold_match is not None:
            lod_threshold = float(lod_threshold_match.group(2))

        return (tumor_awareness, lod_threshold)

    def _get_table_headers(self, data):
        """Create the headers config"""

        crosscheckfingerprints_table_cols = [
            "RESULT",
            "LOD_SCORE",
        ]
        crosscheckfingerprints_table_cols_hidden = [
            "LOD_THRESHOLD",
            "LEFT_RUN_BARCODE",
            "LEFT_LANE",
            "LEFT_MOLECULAR_BARCODE_SEQUENCE",
            "LEFT_LIBRARY",
            "LEFT_FILE",
            "RIGHT_RUN_BARCODE",
            "RIGHT_LANE",
            "RIGHT_MOLECULAR_BARCODE_SEQUENCE",
            "RIGHT_LIBRARY",
            "RIGHT_FILE",
            "DATA_TYPE",
        ]

        # Allow customisation from the MultiQC config
        sampletracking_config = getattr(config, "sampletracking_config", {})
        crosscheckfingerprints_table_cols = sampletracking_config.get(
            "CrosscheckFingerprints_table_cols", crosscheckfingerprints_table_cols
        )
        crosscheckfingerprints_table_cols_hidden = sampletracking_config.get(
            "CrosscheckFingerprints_table_cols_hidden", crosscheckfingerprints_table_cols_hidden,
        )

        # Add Left and Right Sample names / groups, keeping it as minimal as possible
        sample_group_are_same = (
            lambda x: x["LEFT_SAMPLE"] == x["LEFT_GROUP_VALUE"] and x["RIGHT_SAMPLE"] == x["RIGHT_GROUP_VALUE"]
        )
        if all(sample_group_are_same(values) for values in data.values()):
            crosscheckfingerprints_table_cols = ["LEFT_SAMPLE", "RIGHT_SAMPLE",] + crosscheckfingerprints_table_cols
            crosscheckfingerprints_table_cols_hidden += [
                "LEFT_GROUP_VALUE",
                "RIGHT_GROUP_VALUE",
            ]
        else:
            crosscheckfingerprints_table_cols = [
                "LEFT_SAMPLE",
                "LEFT_GROUP_VALUE",
                "RIGHT_SAMPLE",
                "RIGHT_GROUP_VALUE",
            ] + crosscheckfingerprints_table_cols

        headers = OrderedDict()
        for h in FIELD_DESCRIPTIONS:
            # Skip anything not set to visible
            if h not in crosscheckfingerprints_table_cols:
                continue

            # Set up the configuration for the column
            h_title = h.replace("_", " ").strip().lower().capitalize().replace("Lod", "LOD")
            headers[h] = {
                "title": h_title,
                "description": FIELD_DESCRIPTIONS.get(h),
                "namespace": "CrosscheckFingerprints",
                "scale": False,
            }

            # Rename Result to be a longer string so the table formats more nicely
            if h == "RESULT":
                headers[h]["title"] = "Result"
                headers[h]["cond_formatting_rules"] = {
                    "pass": [{"s_contains": "EXPECTED_"}],
                    "warn": [{"s_eq": "INCONCLUSIVE"}],
                    "fail": [{"s_contains": "UNEXPECTED_"}],
                }

            # Add appropriate colors for LOD scores
            if h.startswith("LOD"):
                headers[h]["scale"] = "RdYlGn"
                headers[h]["shared_key"] = "LOD"
                headers[h]["bars_zero_centrepoint"] = True

            if h in crosscheckfingerprints_table_cols_hidden:
                headers[h]["hidden"] = True

        return headers
