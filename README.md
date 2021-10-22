# DNAssembly

*Yet another DNA manipulation package, except this one does exactly what I want*

DNAssembly is a framework that lets you do your cloning as you would in the wet lab, but on your computer. In place of
colorless liquids in test tubes, DNAssembly will take your plasmid files in your favorite format (Genbank, Snapgene)
and run them through an assembly protocol to produce fully annotated cloning products. All input sequences are verified
to avoid common cloning pitfalls.

DNAssembly is initially implemented to perform Modular Golden Gate Cloning as described in "A Highly Characterized Yeast
Toolkit for Modular, Multipart Assembly" (Lee et al. 2013).

# Dependencies
 - networkx
 - biopython
 - pytest


# Required Benchling Environmental Variables:
- BENCHLING_URL: Benchling URL (e.g. "https://{tenant_name}.benchling.com/api")
- BENCHLING_API_VERSION: Version of Benchling API (e.g. "v2")
- BENCHLING_API_KEY: Benchling user API tenant key for for querying REST API endpoints
- SEQ_NAMING_STRAT: Benchling DNA plasmid sequence naming strategy 
- SEQ_REGISTRY_ID: Benchling DNA plasmid sequence registry ID code
- SEQ_SCHEMA_ID: Benchling DNA plasmid sequence schema ID code
- CARB: TBD (not sure what this represents)
- KAN: TBD (not sure what this represents)
- PART_NAMING_STRAT: Benchling DNA part naming strategy 
- PART_REGISTRY_ID: Benchling DNA part registry ID code
- PART_SCHEMA_ID: Benchling DNA part schema ID code


# Loading environmental variables:
### By command line
```
export <env_variable>=<env_variable_value>
```
### By file
```
set -o allexport
source <config_file_path>
set +o allexport
```


# Installation:
### Using root repository directory
```
pip install .
```
### In another project using pip CLI
```
pip install git+https://github.com/outpace-bio/rearray-python@[branch-name or commit-id]
```
### In another project using requirements.txt file
```
git+ssh://git@github.com/outpace-bio/dnassembly.git@[branch-name or commit-id]#egg=rearray
```

# Getting Started
1. Install package as outlined above
2. Load up required benchling environmental variables
3. (optional) Run unit tests
4. Package is now ready to use!


# Unit Testing
### Executing the tests
In the project directory execute the following commands:
1. Install repository using root directory (outlined in installation section) 
2. In command line, run: `pytest tests`

### Generating code coverages
1. Run unit tests with coverage package `coverage run -m pytest tests`
2. Generate report and view in standard output `coverage report -m`
3. Generate report in html format for easier viewing `coverage html`

### Testing Road Map
Utils:
- [x] annotation
- [ ] benchlingAPI
- [ ] conversion
- [ ] utils
