# Install Submodules 'Pysam' and 'BioPython'

all:
	cd biopython; python setup.py install
	cd pysam; python setup.py install
