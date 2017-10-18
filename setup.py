from setuptools import setup
import sys, os

exec(open('salmonella_crispr/version.py').read())

setup(name="salmonella_crispr",
        version=__version__,
        description='This tool gets a CRISPR profile by identifying the presence of known spacers and direct repeats (DRs) in a given sequence based on a catalogue. This tool is a reimplemntation of a former Perl tool in developed by G. Guigon',
        author='Kenzo-Hugo Hillion and Laetitia Fabre',
        author_email='kehillio@pasteur.fr and laetitia.fabre@pasteur.fr',
        license='',
        keywords = ['salmonella','crispr','fasta'],
        install_requires=['biopython'],
        packages=["salmonella_crispr"],
        package_data={
        'salmonella_crispr': ['data/spacers_Salmonella.fa'],
        },
        entry_points={'console_scripts':['crispr=salmonella_crispr.main:run']},
        classifiers=[
            'Development Status :: 4 - Beta',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Operating System :: OS Independent',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'Environment :: Console',
            ],
        )
