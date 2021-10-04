# DNAssembly

*Yet another DNA manipulation package, except this one does exactly what I want*

DNAssembly is a framework that lets you do your cloning as you would in the wet lab, but on your computer. In place of
colorless liquids in test tubes, DNAssembly will take your plasmid files in your favorite format (Genbank, Snapgene)
and run them through an assembly protocol to produce fully annotated cloning products. All input sequences are verified
to avoid common cloning pitfalls.

DNAssembly is initially implemented to perform Modular Golden Gate Cloning as described in "A Highly Characterized Yeast
Toolkit for Modular, Multipart Assembly" (Lee et al. 2013).

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

# Running Unit Testing

In the project directory execute the following commands:
```
pip install .
pytest tests
```

## Coverage

Utils:
- [x] annotation
- [ ] benchlingAPI
- [ ] conversion
- [ ] utils
