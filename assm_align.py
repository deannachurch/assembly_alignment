#!/usr/bin/env python
############################################
#
#  Process assm-assm alignment files to get stats
#  For now using a config to point to relevant files (rather than parameters)
#  
#  Deanna M. Church
#  23 Dec, 2015
#
#############################################
from __future__ import division
import csv
from collections import defaultdict
import pybedtools
import subprocess
import numpy as np
import matplotlib.pylab as plt
#import seaborn as sns
import os
import logging
import logging.handlers
from logging import config
import datetime
import yaml
import errno

def parseSeqRep(fi, assm_name, assm_acc):
	assm_list=[]
	try:
		with open(fi, 'r') as infile:
			data=csv.reader(infile, delimiter="\t")
			for line in data:
				if line[0].startswith("# Assembly Name:"):
					if assm_name not in line[0]:
						logging.critical("Wrong report, assembly name mismatch: %s\t%s" % (assm_name, line[0]))
			 	if not line[0].startswith("#"):
			 		rec = Seq()
			 		rec.assm=assm_name
			 		rec.name=line[0]
			 		rec.ref_acc=assm_acc
			 		rec.length=int(line[8])
			 		rec.role=line[1]
			 		rec.assm_unit=line[7]
			 		assm_list.append(rec)

	except IOError:
		print "Can't open %s" % fi
	return assm_list

class Seq(object):
	#set up attributes you want to track about the sequence
	def __init__(self):
		self.assm=""
		self.name=""#sequence name, e.g. 1
		self.ref_acc=""#sequence refseq accession.version
		self.length=0#sequence length
		self.role=""#seq role ('assembled-molecule, 'alt-scaffold', 'unlocalized-scaffold', 'unplaced-scaffold')
		self.assm_unit=""#assembly unit seq is in 

def main():
	#set up logging
	try:
		os.makedirs("log")
	except OSError as e:
		if e.errno != errno.EEXIST:
			print "ERROR: assm_align.py: log directory failure: %s" % e.errno
			sys.exit(1)
	date=datetime.datetime.now().strftime("%Y-%m-%d_%H%M")
	log_file="log/assm_align_%s.log" % (date)
	config_dict=yaml.load(open("resources/assm_align_log_cfg.yml", 'r'))
	config_dict['handlers']['logfile']['filename']=log_file
	logging.config.dictConfig(config_dict)
	logger=logging.getLogger()
	logger.info("================assm_align.py started: log file=%s================" % log_file)
	
	##read config file and get file parameters
	cfg_dict=yaml.load(open("resources/assm_align_cfg.yml", 'r'))
	assm1=cfg_dict['input_files']['assm1']
	assm2=cfg_dict['input_files']['assm2']
	##create sequence objects
	assm1_list=[]
	assm2_list=[]
	assm1_list=parseSeqRep(assm1['seq_rpt'], assm1['name'], assm1['acc'])
	logging.info("Read %s, sequences: %d" % (assm1['name'], len(assm1_list)))
	assm2_list=parseSeqRep(assm2['seq_rpt'], assm2['name'], assm2['acc'])
	logging.info("Read %s, sequences: %d" % (assm2['name'], len(assm2_list)))

	



if __name__=="__main__":
	main()