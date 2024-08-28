#!/usr/bin/env python

"""MultiQC module to parse output from Reads Statistics."""

import json
import logging
from collections import Counter, OrderedDict, defaultdict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Reads Stats",
            anchor="reads_stats",
            href="https://github.com/y9c/pseudoU-BIDseq",
        )

        # Find and load any FastQ Screen reports
        self.reads_data = defaultdict(Counter)
        self.computed_data = dict()
        self.num_orgs = 0
        for f in self.find_log_files("readsStats/cutadapt", filehandles=True):
            # combine runs
            n, _ = f["s_name"].removesuffix("_cutadapt").rsplit("_run", 1)
            self.reads_data[n] += Counter(self.parse_cutadapt(f))
            self.add_data_source(f)

        for f in self.find_log_files("readsStats/bowtie2", filehandles=True):
            n, _ = f["s_name"].removesuffix("_genes").rsplit("_run", 1)
            self.reads_data[n] += Counter(self.parse_bowtie2(f))
            self.add_data_source(f)
        for f in self.find_log_files("readsStats/star", filehandles=True):
            n, _ = f["s_name"].removesuffix("_genome").rsplit("_run", 1)
            self.reads_data[n] += Counter(self.parse_star(f))
            self.add_data_source(f)
        for f in self.find_log_files("readsStats/dedup", filehandles=True):
            parsed_data = self.parse_dedup(f)
            n = f["s_name"].removesuffix("_dedup")
            if n.endswith("_genes"):
                self.reads_data[n.removesuffix("_genes")] += Counter(
                    {k + "_genes": v for k, v in parsed_data.items()}
                )
            if n.endswith("_genome"):
                self.reads_data[n.removesuffix("_genome")] += Counter(
                    {k + "_genome": v for k, v in parsed_data.items()}
                )
            self.add_data_source(f)

        for s_name, c in self.reads_data.items():
            d = dict()
            # Samples with 0 reads and reports with some skipped sections might be missing things here
            d["total_sequences"] = c["before_trimming"]
            d["trimming_discarded_sequences"] = (
                c["before_trimming"] - c["after_trimming"]
            )
            d["trimming_discarded_percentage"] = (
                d["trimming_discarded_sequences"] / d["total_sequences"]
            )
            d["mapping_to_genes_sequences"] = max(
                c["before_genes"] - c["after_genes"], 1
            )
            d["dedup_to_genes_sequences"] = c["after_dedup_genes"]
            d["duplicate_to_genes_sequences"] = (
                d["mapping_to_genes_sequences"] - c["after_dedup_genes"]
            )
            d["duplicate_to_genes_percentage"] = (
                d["duplicate_to_genes_sequences"] / d["mapping_to_genes_sequences"]
            )
            d["dedup_to_genes_percentage"] = (
                d["dedup_to_genes_sequences"] / d["mapping_to_genes_sequences"]
            )
            d["mapping_to_genes_percentage"] = (
                d["mapping_to_genes_sequences"] / c["after_trimming"]
            )
            d["mapping_to_genome_sequences"] = max(
                c["before_genome"] - c["after_genome"], 1
            )
            d["mapping_to_genome_percentage"] = (
                d["mapping_to_genome_sequences"] / c["after_trimming"]
            )
            d["dedup_to_genome_sequences"] = c["after_dedup_genome"]
            d["duplicate_to_genome_sequences"] = (
                d["mapping_to_genome_sequences"] - c["after_dedup_genome"]
            )
            d["duplicate_to_genome_percentage"] = (
                d["duplicate_to_genome_sequences"] / d["mapping_to_genome_sequences"]
            )
            d["dedup_to_genome_percentage"] = (
                d["dedup_to_genome_sequences"] / d["mapping_to_genome_sequences"]
            )
            d["unmapped_sequence"] = (
                c["after_trimming"]
                - d["mapping_to_genome_sequences"]
                - d["mapping_to_genes_sequences"]
            )
            self.computed_data[s_name] = d

        # Add plot
        self.add_section(
            name="Reads Statistics (barplot)",
            anchor="Reads Statistics",
            description="""
            This plot shows the number of reads after each step. For this customized pipeline only, this is not a general purpose module.
            You can chose click the figure legend to show or hide some kinds of reads. For example, if you only select `duplicate to genome` and `dedup to genome`, you can have the duplication ratio stats.
            """,
            plot=self.readstats_simple_plot(),
        )
        # Add to the general statistics table
        self.fastqc_general_stats()

    def parse_cutadapt(self, f):
        """Parse cutadapt data from json.
        {
          "read_counts": {
            "input": 3946922,
            "filtered": {
              "too_short": 177599,
              "too_long": null,
              "too_many_n": null,
              "too_many_expected_errors": null,
              "casava_filtered": null,
              "discard_trimmed": null,
              "discard_untrimmed": null
            },
            "output": 58176,
            ...
        """
        parsed_data = dict()

        json_data = json.load(f["f"])
        parsed_data["before_trimming"] = json_data["read_counts"]["input"]
        parsed_data["after_trimming"] = json_data["read_counts"]["output"]
        return parsed_data

    def parse_bowtie2(self, f):
        """Parse bowtie2 data."""
        parsed_data = dict()
        n_mapped = 0
        for l in f["f"]:
            if "reads; of these:" in l:
                parsed_data["before_genes"] = int(l.split(" ")[0])
            elif " aligned exactly 1 time" in l or " aligned >1 times" in l:
                n_mapped += int(l.split("(")[0])
        parsed_data["after_genes"] = parsed_data["before_genes"] - n_mapped
        return parsed_data

    def parse_star(self, f):
        """Parse STAR data."""
        parsed_data = dict()
        n_mapped = 0
        for l in f["f"]:
            if "Number of input reads |" in l:
                parsed_data["before_genome"] = int(l.split("|")[-1])
            elif (
                "Uniquely mapped reads number |" in l
                or "Number of reads mapped to multiple loci |" in l
                or "Number of chimeric reads |" in l
            ):
                n_mapped += int(l.split("|")[-1])
        parsed_data["after_genome"] = parsed_data["before_genome"] - n_mapped
        return parsed_data

    def parse_dedup(self, f):
        """Parse dedup data."""
        parsed_data = dict()
        for l in f["f"]:
            if l.strip().split("\t")[-1] == "primary":
                parsed_data["after_dedup"] = int(l.split("\t")[0])
        return parsed_data

    def readstats_simple_plot(self):
        """Makes a simple bar plot with summed alignment counts for each
        species, stacked."""

        data = self.computed_data
        # First, sum the different types of alignment counts
        # Sort the categories by the total read counts
        cats = OrderedDict()

        pconfig = {
            "id": "reads_stats_bar",
            "title": "Reads Statistics (barplot)",
            "cpswitch_c_active": False,
            "ylab": "Percentages",
        }

        cats = {
            "trimming_discarded_sequences": {
                "name": "trimming discarded",
                "color": "red",
            },
            "unmapped_sequence": {
                "name": "Unmapped discarded",
                "color": "dimgray",
            },
            "duplicate_to_genes_sequences": {
                "name": "duplicate to genes",
                "color": "darkblue",
            },
            "duplicate_to_genome_sequences": {
                "name": "duplicate to genome",
                "color": "darkgreen",
            },
            "dedup_to_genes_sequences": {
                "name": "dedup to genes",
                "color": "skyblue",
            },
            "dedup_to_genome_sequences": {
                "name": "dedup to genome",
                "color": "lightgreen",
            },
        }

        return bargraph.plot(data, cats, pconfig)

    def fastqc_general_stats(self):
        """Add some single-number stats to the basic statistics table at the
        top of the report."""

        # Prep the data
        data = self.computed_data

        headers = OrderedDict()
        headers["total_sequences"] = {
            "title": "# total Sequences",
            "description": "Total number of raw reads sequenced for this sample",
            "min": 0,
            "scale": "Set1",
            "format": "{:,.0f}",
        }
        headers["trimming_discarded_sequences"] = {
            "title": "# trimming discarded",
            "description": "Number of discarded reads in trimming step (include adapter dimer, low QC, etc.)",
            "min": 0,
            "scale": "Reds",
            "format": "{:,.0f}",
            "hidden": True,
        }
        headers["trimming_discarded_percentage"] = {
            "title": "% trimming discarded",
            "description": "Percentage of discarded reads in trimming step (include adapter dimer, low QC, etc.)",
            "min": 0,
            "max": 1,
            "scale": "Reds",
            "format": "{:.2%}",
        }
        headers["mapping_to_genes_percentage"] = {
            "title": "% mapping to genes",
            "description": "Percentage of reads that aligned into genes reference in mapping step",
            "min": 0,
            "max": 1,
            "scale": "Blues",
            "format": "{:.2%}",
        }
        headers["mapping_to_genes_sequences"] = {
            "title": "# mapping to genes",
            "description": "Number of reads that aligned into genes reference in mapping step",
            "min": 0,
            "scale": "Blues",
            "format": "{:,.0f}",
            "hidden": True,
        }
        headers["dedup_to_genes_percentage"] = {
            "title": "% dedup to genes",
            "description": "Percentage of reads after duplicate removal step vs. total reads that aligned into genes reference",
            "min": 0,
            "max": 1,
            "scale": "Blues",
            "format": "{:.2%}",
        }
        headers["dedup_to_genes_sequences"] = {
            "title": "# usable reads for genes",
            "description": "Number of usable reads for detect modification in genes",
            "min": 0,
            "scale": "Blues",
            "format": "{:,.0f}",
        }
        headers["mapping_to_genome_percentage"] = {
            "title": "% mapping to genome",
            "description": "Percentage of reads that aligned into genome reference in mapping step",
            "min": 0,
            "max": 1,
            "scale": "Greens",
            "format": "{:.2%}",
            # "hidden": True,
        }
        headers["mapping_to_genome_sequences"] = {
            "title": "# mapping to genome",
            "description": "Number of reads that aligned into genome reference in mapping step",
            "min": 0,
            "scale": "Greens",
            "format": "{:,.0f}",
            "hidden": True,
        }
        headers["dedup_to_genome_percentage"] = {
            "title": "% dedup to genome",
            "description": "Percentage of reads after duplicate removal step vs. total reads that aligned into genome reference",
            "min": 0,
            "max": 1,
            "scale": "Greens",
            "format": "{:.2%}",
        }
        headers["dedup_to_genome_sequences"] = {
            "title": "# usable reads for genome",
            "description": "Number of usable reads for detect modification in genome",
            "min": 0,
            "scale": "Greens",
            "format": "{:,.0f}",
        }
        self.general_stats_addcols(data, headers)
