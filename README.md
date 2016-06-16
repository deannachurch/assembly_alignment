# assembly_alignment

Preliminary assembly-assembly alignment analysis

Currently supports:
* creating a stats table with 1 row per assembly sequence with numbers per alignment discrepancy type
* if the assembly has chromosomes, make graphs.
* Creates per assembly bed files for each discrepancy type

To do:
* update stats to add the top ten by size for each even type

To run:
Create a config file specifying path to files and files to create

15 Jun 2015 update
created NCBI_Nonstandard_Output branch.
Not merging with master right now. This branch is for dealing with non-standard NCBI assembly alignment output and
is full of kludges that may or may not generally work
Can't just ignore sequence report as I use this to instantiate sequence objects. Write a script to create this.



Dependencies:

* from __future__ import division
* import csv
* from collections import defaultdict
* import collections
* import pybedtools
* from pybedtools import BedTool
* import subprocess
* import numpy as np
* import matplotlib.pylab as plt
* import os
* import logging
* import logging.handlers
* from logging import config
* import datetime
* import yaml
* import errno
* import argparse
