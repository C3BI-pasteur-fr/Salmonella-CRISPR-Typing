from setuptools import setup
import sys, os

exec(open('salmonella_crispr/version.py').read())

setup(name="salmonella_crispr",
        version=__version__,
        description='',
        author='Kenzo-Hugo Hillion and Herve Menager',
        author_email='kehillio@pasteur.fr and hmenager@pasteur.fr',
        license='',
        keywords = ['salmonella','crispr','fasta'],
        install_requires=['biopython'],
        packages=["salmonella_crispr"],
        package_data={
        'salmonella_crispr': ['data/spacers_Salmonella.fa'],
        },
        entry_points={'console_scripts':['crispr=salmonella_crispr.main:run']},
        classifiers=[
            'Development Status :: 1 - Planning',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Operating System :: OS Independent',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'Environment :: Console',
            ],
        )
