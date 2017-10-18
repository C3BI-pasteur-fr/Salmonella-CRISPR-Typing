# Development of Salmonella CRISPR tool as a command line interface in Python

CRISPR polymorphism is a powerful tool to subtype Salmonella strains and is now used in routine for epidemiological investigations.
This tool gets a CRISPR profile by identifying the presence of known spacers and direct repeats (DRs) in a given sequence based on a catalogue.

This tool is a reimplemntation of a former tool in Perl/CGI developed by G. Guigon

-------------------------

## Installation

#### Requirements

This tool uses [Biopython](http://biopython.org/) library.

To install, run the following command:

```bash
pip install git+https://gitlab.pasteur.fr/kehillio/salmonella-CRISPR.git 
```

## How does it work ?

Salmonella crispr tool allows the discovery of spacers and DRs from either a PCR product of a full assembly.


```bash
crispr your_file.fa
```

Please, read the help message for advanced features on the tool.

```bash
crispr --help
```

## References

Fabre L1, Zhang J, Guigon G, Le Hello S, Guibert V, Accou-Demartin M, de Romans S, Lim C, Roux C, Passet V, Diancourt L, Guibourdenche M, Issenhuth-Jeanjean S, Achtman M, Brisse S, Sola C, Weill FX. CRISPR typing and subtyping for improved laboratory surveillance of Salmonella infections. PLoS One. 2012;7(5):e36995. DOI:[10.1371/journal.pone.0036995](http://doi.org/10.1371/journal.pone.0036995)