# A plugin for MultiQC containing customized modules and templates

## Installation

This plugin can be installed using the following methods

- using `pip`:

```bash
pip install --upgrade --force-reinstall git+https://github.com/y9c/MultiQC_HeLab.git
```

## Usage

### Modules

| Name           | Description                                                                                                                                                                      |
| -------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Sampletracking | Parse Picard [CrosscheckFingerprints](https://gatk.broadinstitute.org/hc/en-us/articles/360057441151-CrosscheckFingerprints-Picard-) output and format sensible tables and plots |

### Templates

| Name  | Description                                                                                                                        |
| ----- | ---------------------------------------------------------------------------------------------------------------------------------- |
| healb | template with custom logo's and affiliate links. To enable this template, add the `-t/--template helab` option to the command line |
