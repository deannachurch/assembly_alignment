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


def parseAlignReport(fi, assm_name):
	#data structures
	gap_no_hit = defaultdict(int)
	ungap_no_hit = defaultdict(int)
	no_hit_loc = defaultdict(list)
	try:
		with open(fi, 'r') as infile:
			data=csv.reader(infile, delimiter="\t")
			for line in data:
				#skip empty lines
				if (len(line)==0):
					#skip blank lines
					continue
				elif (line[0].startswith("#")):
					#skip header lines
					continue
				elif (len(line) == 1):
					#metadata
					if line[0].startswith('Query Assembly Name:'):
						query_asm_name = line[0].split(':')[1].lstrip()
						if not query_asm_name == assm_name:
							logging.critical("Assembly name mismatch: %s\t%s" % (query_asm_name, assm_name))
							sys.exit(3)
					if line[0].startswith('Sequence Name:'):
						seq_name = line[0].split(':')[1].lstrip()
				else:
					#data lines
					#print line
					data_type=line[0]
					start = int(line[1])
					end = int(line[2])
					gap_len = int(line[3])
					ungap_len = int(line[4])
					loc="%s:%d-%d" % (seq_name, start, end)
					if data_type == "NoHit":
						gap_no_hit[seq_name] += gap_len
						ungap_no_hit[seq_name] += ungap_len
						no_hit_loc[seq_name].append(loc)
				
	except IOError:
		logging.critical("Can't open %s" % fi)
		sys.exit(2)

	print gap_no_hit

def parseSeqRep(fi, assm_name, assm_acc):
	#parse the NCBI seq report to get names, roles, assm-units, etc
	assm_list={}
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
			 		assm_list[line[0]]=rec

	except IOError:
		logging.critical("Can't open %s" % fi)
		sys.exit(1)

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

	def set_nohit(self, length, loc_list):
		self.no_hit_len=length
		self.no_hit_list=loc_list


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
	assm1_dict=parseSeqRep(assm1['seq_rpt'], assm1['name'], assm1['acc'])
	logging.info("Read %s, sequences: %d" % (assm1['name'], len(assm1_dict)))
	assm2_dict=parseSeqRep(assm2['seq_rpt'], assm2['name'], assm2['acc'])
	logging.info("Read %s, sequences: %d" % (assm2['name'], len(assm2_dict)))
	##parse alignment report
	parseAlignReport(assm1['align_rpt'], assm1['name'])

	



if __name__=="__main__":
	main()