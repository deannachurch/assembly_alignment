# assembly_alignment

Preliminary assembly-assembly alignment analysis

Currently supports:
* creating a stats table with 1 row per assembly sequence with numbers per alignment discrepancy type
* if the assembly has chromosomes, make graphs.

To do: 
* create bed files for each alignmnet discrepancy type
* update stats to add the top ten by size for each even type

To run: 
Create a config file specifying path to files and files to create

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
* import seaborn as sns
* import os
* import logging
* import logging.handlers
* from logging import config
* import datetime
* import yaml
* import errno
* import argparse
