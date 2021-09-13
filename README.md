## The PR2 primer database
---
[![DOI](https://zenodo.org/badge/191439796.svg)](https://zenodo.org/badge/latestdoi/191439796)


An [interactive database](https://app.pr2-primers.org/) of eukaryotic rRNA primers and primer sets for metabarcoding studies compiled from the literature. 

### Database structure

* **Primers**. Primers have been mapped when possible onto the reference SSU sequence for _[Saccharomyces cerevisiae](http://apollo.chemistry.gatech.edu/RibosomeGallery/eukarya/S%20cerevisiae/SSU/index.html)_ ([FU970071](https://www.ncbi.nlm.nih.gov/nuccore/FU970071), 1799 nucleotides, first nucleotide mrearked as 1).  

* **Primer sets**. Primer sets for 18S rRNA have been tested against the [eukaryotic PR2 database](https://pr2-database.org/) version 4.12.0 as well as [Silva Seed release 132](https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v132.tgz) and results can be displayed interactively.

### Panels

* **About**: Basic information 
* **Primers**: table with download
* **Primer sets**: table with download
* **Amplification - overview**: Give for all primer sets tested % of sequences amplified and amplicon size
* **Amplification - detail**: For one primer set tested, detail for different taxonomic levels (kingdom, supergroup, division, class)
* **Test your primer/probe**: Test a single primer or probe against PR2 version 4.12.0 and Silva seed 132
* **Test your primer set**: Test your primer set against PR2 version 4.12.0 and Silva seed 132

### Errors

Please report errors or primer sets not listed here in the [Issues page of the PR2 primer database](https://github.com/pr2database/pr2-primers/issues).


### Citation

Vaulot, D., Mahé, F., Bass, D., & Geisen, S. (2021). [pr2-primer : An 18S rRNA primer database for protists](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13465). Molecular Ecology Resources, in press. DOI: 10.1111/1755-0998.13465

### Resources
* Website: https://app.pr2-primers.org/
* Docker: https://hub.docker.com/repository/docker/vaulot/pr2-primers
* Source code: https://github.com/pr2database/pr2-primers
* Script: https://pr2database.github.io/pr2-primers/PR2_Primers.html

### Maintainer
* Daniel Vaulot: vaulot@gmail.com

### Contributors

* Stefan Geisen:  stefan.geisen@wur.nl
* Fréderic Mahé: frederic.mahe@cirad.fr
* David Bass: david.bass@cefas.co.uk

### Versions

1.1.0 - 2021-05-17
* Panel to test individual primer or probes.
* Add new primers and primer sets.

1.0.4 - 2021-01-26
* Add time out of 10 min: after 10 min of inactivity close browser window and stop app.

1.0.2 - 2020-12-25
* Add non-Metazoa primer and primer sets from Bass and del Campo 2020

1.0.1 - 2020-12-14
* Fix problem with V9 for bacteria

1.0.0 - 2020-11-11
* Initial version
