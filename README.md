# Development of Salmonella CRISPR Typing tool as a command line interface in Python

CRISPR polymorphism is a powerful tool to subtype Salmonella strains and is now used in routine for epidemiological investigations.
This tool gets a CRISPR profile by identifying the presence of known spacers and direct repeats (DRs) in a given sequence based on a catalogue.

This tool is a reimplemntation of a former tool in Perl/CGI developed by G. Guigon

-------------------------

## Installation

#### Requirements

This tool is developed in Python3 and uses [Biopython](http://biopython.org/) library.

We recommand the use of a virtual environment to install and run the tool.

To install, run the following commands:

```bash
python3 -m venv crispr_typing
source crispr_typing/bin/activate
pip install git+https://gitlab.pasteur.fr/kehillio/salmonella-CRISPR.git 
```

To check if the installation was successful, try to display the help menu:

```bash
crispr_typing -h
```

## How does it work ?

Salmonella CRISPR Typing tool allows the discovery of spacers and DRs from either a PCR product of a full assembly.


```bash
crispr_typing your_file.fa
```

Please, read the help message for advanced features on the tool.

```bash
crispr_typing --help
```

## References

Fabre L, Zhang J, Guigon G, Le Hello S, Guibert V, Accou-Demartin M, de Romans S, Lim C, Roux C, Passet V, Diancourt L, Guibourdenche M, Issenhuth-Jeanjean S, Achtman M, Brisse S, Sola C, Weill FX. CRISPR typing and subtyping for improved laboratory surveillance of Salmonella infections. PLoS One. 2012;7(5):e36995. DOI:[10.1371/journal.pone.0036995](http://doi.org/10.1371/journal.pone.0036995)
