from setuptools import setup, find_packages

setup(
    name='proteinmpnn',
    version='1.0.1',
    packages=find_packages(),
    install_requires=[
        # List your project dependencies here
    ],
    entry_points={'console_scripts': [
        "proteinmpnn=proteinmpnn.protein_mpnn_run:main",
    ],},
    # Metadata
    author='Justas Dauparas',
    author_email='',
    description='A little convenience fork of the proteinmpnn package by Justas Dauparas et al.',
    url='https://github.com/dauparas/ProteinMPNN',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)