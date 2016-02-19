# ensembl-annotation

## Running perltidy
Before commiting to this repository please run the following has been run to ensure the Perl source files are formatted correctly.
```
find . -name '*.p[l|m]' -exec perltidy -pro=/path/to/ensembl-annotation/perltidyrc  -b {} \;
```
