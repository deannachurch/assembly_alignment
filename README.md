# assembly_alignment

Preliminary assembly-assembly alignment analysis

# Features
* Creates a stats table with 1 row per assembly sequence with numbers per alignment discrepancy type
* Makes graphs, if the assembly has chromosomes
* Creates per assembly bed files for each discrepancy type

# Future work
* Update stats to add the top ten by size for each even type

# Install
Ensure you're using Python 2.7 and have [bedtools](http://bedtools.readthedocs.org/en/latest/) installed, then:

```
$ git clone https://github.com/deannachurch/assembly_alignment
$ cd assembly_alignment

# Work around a pysam issue, https://github.com/pysam-developers/pysam/issues/247
$ export HTSLIB_CONFIGURE_OPTIONS="--disable-libcurl"

# Create a virtualenv and activate it
$ virtualenv env
$ source env/bin/activate

# Install dependencies
$ pip install -r requirements.txt
```

# Run
Create a config file specifying path to files and files to create

