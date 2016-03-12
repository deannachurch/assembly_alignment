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
import collections
import pybedtools
from pybedtools import BedTool
import subprocess
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import os
import sys
import logging
import logging.handlers
from logging import config
import datetime
import yaml
import errno
import argparse

def getLength(uniq_loc_list):
	#get length of intervals in a bedtool
	tot_len=defaultdict(int)
	for loc in uniq_loc_list:
		chrom=loc[0]
		tot_len[chrom]+= loc[3]

	return tot_len

def mergeLoc(loc_list):
	#For data types where there can be redundancy, we need to merge the locs and then get the lengths
	uniq_loc_list=[]
	tot_len=0
	loc_str=""
	for loc in loc_list:
		(chrom, pre_loc)=loc.split(":")
		(start, end)=pre_loc.split("-")
		if not loc_str:
			loc_str="%s\t%s\t%s" % (chrom, start, end)
		else:
			loc_str="%s\n%s\t%s\t%s" % (loc_str, chrom, start, end)
	non_uniq=BedTool(loc_str, from_string=True)
	uniq_int=non_uniq.merge()
	for inter in uniq_int:
		loc=(inter.chrom, inter.start, inter.stop, inter.stop-inter.start)
		uniq_loc_list.append(loc)
	return uniq_loc_list

def sort_list(unsort_list):
	sort_seq_list=[]
	try:
		sort_seq_list=sorted(unsort_list, key=lambda item: (int(item.partition(' ')[0]) if item[0].isdigit() else float('inf'), item))
	except:
		sort_seq_list=sorted(unsort_list)
	return sort_seq_list

def parseAlignReport(fi, assm_name, obj_dict):
	#data structures
	gap_no_hit = defaultdict(int)
	ungap_no_hit = defaultdict(int)
	no_hit_loc = defaultdict(list)
	sp_loc=defaultdict(list)
	sp_only_loc=defaultdict(list)
	inv_loc=defaultdict(list)
	mix_loc=defaultdict(list)
	offchrom_loc=defaultdict(list)
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
					start = int(line[1])-1
					end = int(line[2])
					gap_len = int(line[3])
					ungap_len = int(line[5])
					loc="%s:%d-%d" % (seq_name, start, end)
					#can count no hit locs because there should be no overlap
					if data_type == "NoHit":
						gap_no_hit[seq_name] += gap_len
						ungap_no_hit[seq_name] += ungap_len
						no_hit_loc[seq_name].append(loc)
					#all other locs, need to get locs and then uniquify
					elif data_type == "Inv":
						inv_loc[seq_name].append(loc)
					elif data_type == "Mix":
						mix_loc[seq_name].append(loc)
					elif data_type == "SP":
						sp_loc[seq_name].append(loc)
					elif data_type == "SP Only":
						sp_only_loc[seq_name].append(loc)

	except IOError:
		logging.critical("Can't open %s" % fi)
		sys.exit(2)
	#set alignment attribute information for seq_object
	for seq in obj_dict:
		if not seq in gap_no_hit:
			logging.debug("Seq has no NoHit %s: %s" % (assm_name, seq))
			obj_dict[seq].set_nohit(0, 0, [])
		else:
			uniq_no_hit_loc=mergeLoc(no_hit_loc[seq])
			obj_dict[seq].set_nohit(gap_no_hit[seq], ungap_no_hit[seq], uniq_no_hit_loc)
		if not seq in sp_loc:
			logging.debug("Seq has no SP %s: %s" % (assm_name, seq))
			obj_dict[seq].set_sp(0, [])
		else:
			sp_uniq_loc=mergeLoc(sp_loc[seq])
			sp_len=getLength(sp_uniq_loc)
			obj_dict[seq].set_sp(sp_len[seq], sp_uniq_loc)
		if not seq in sp_only_loc:
			logging.debug("Seq has no SP Only %s: %s" % (assm_name, seq))
			obj_dict[seq].set_sp_only(0, [])
		else:
			sp_only_uniq_loc=mergeLoc(sp_only_loc[seq])
			sp_only_len=getLength(sp_only_uniq_loc)
			obj_dict[seq].set_sp_only(sp_only_len[seq], sp_only_uniq_loc)
		if not seq in inv_loc:
			logging.debug("Seq has no INV data: %s: %s" % (assm_name, seq))
			obj_dict[seq].set_inv(0,[])
		else:
			inv_uniq_loc=mergeLoc(inv_loc[seq])
			inv_len=getLength(inv_uniq_loc)
			obj_dict[seq].set_inv(inv_len[seq], inv_uniq_loc)
		if not seq in mix_loc:
			logging.debug("Seq has no Mix data: %s: %s" % (assm_name, seq))
			obj_dict[seq].set_mix(0,[])
		else:
			mix_uniq_loc=mergeLoc(mix_loc[seq])
			mix_len=getLength(mix_uniq_loc)
			obj_dict[seq].set_mix(mix_len[seq], mix_uniq_loc)


def parseSeqRep(fi, assm_name, assm_acc, chrom_list, exclude_mt):
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
			 		if exclude_mt == True:
			 			if line[1] == 'assembled-molecule' and line[7] == "Primary Assembly":
			 				chrom_list.append(line[0])
			 		else:
			 			if line[1] == 'assembled-molecule':
			 				chrom_list.append(line[0])

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

	def set_nohit(self, gap_len, ungap_len, loc_list):
		self.nohit_len=gap_len
		self.ungap_nohit_len=ungap_len
		self.no_hit_list=loc_list

	def set_sp(self, sp_l, sp_list):
		self.sp_len=sp_l
		self.sp_list=sp_list

	def set_sp_only(self, sp_only_l, sp_only_list):
		self.sp_only_len=sp_only_l
		self.sp_only_list=sp_only_list

	def set_inv(self, inv_l, inv_loc_l):
		self.inv_len=inv_l
		self.inv_loc_list=inv_loc_l

	def set_mix(self, mix_l, mix_loc_l):
		self.mix_len=mix_l
		self.mix_loc_list=mix_loc_l

def makeBed(out_file, assm_dict, data_type):
	out_str=out_file.split(".")[0].split("/")[1]
	out=open(out_file, 'w')
	out.write("track name=%s\n" % out_str)
	seq_list=assm_dict.keys()
	#sort list alphanumerically, it looks like the complex sort barfs on non-GRC assemblies- trying work around.
	#I should do this in the class so I don't have to do it twice in the code.
	sort_seq_list=sort_list(seq_list)
	loc_list=[]
	if data_type == "nohit":
		for seq in sort_seq_list:
			loc_list.extend(assm_dict[seq].no_hit_list)
	elif data_type == "collapse":
		for seq in sort_seq_list:
			loc_list.extend(assm_dict[seq].sp_list)
	elif data_type == "expand":
		for seq in sort_seq_list:
			loc_list.extend(assm_dict[seq].sp_only_list)
	elif data_type == "inv":
		for seq in sort_seq_list:
			loc_list.extend(assm_dict[seq].inv_loc_list)
	elif data_type == "mix":
		for seq in sort_seq_list:
			loc_list.extend(assm_dict[seq].mix_loc_list)
	else:
		logging.error("Unknown data type, abandoning bed: %s" % data_type)

	for loc in loc_list:
		out.write("%s\t%d\t%d\n" % (loc[0], loc[1], loc[2]))
	out.close()

def writeStats(fh, assm1, assm2, assm_dict):
	##stats header
	date=datetime.datetime.now().strftime("%Y-%m-%d")
	fh.write("##%s vs %s assembly alignment report\n##%s\n" % (assm1, assm2, date))
	fh.write("##Overall stats\n")
	fh.write("#Sequence\tNoHit\tUnGap_NoHit\tCollapse(SP)\tExpansion(SP Only)\tInversion\tMix\n")
	##by chromosome
	nohit_tot=0
	ungap_nohit_tot=0
	coll_tot=0
	exp_tot=0
	inv_tot=0
	mix_tot=0
	seq_list=assm_dict.keys()
	#sort list alphanumerically, it looks like the complex sort barfs on non-GRC assemblies- trying work around.
	sort_seq_list=sort_list(seq_list)
	for seq in sort_seq_list:
		obj=assm_dict[seq]
		nohit_tot += obj.nohit_len
		ungap_nohit_tot += obj.ungap_nohit_len
		coll_tot += obj.sp_len
		exp_tot += obj.sp_only_len
		inv_tot += obj.inv_len
		mix_tot += obj.mix_len
		fh.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % (seq, obj.nohit_len, obj.ungap_nohit_len, obj.sp_len, obj.sp_only_len, obj.inv_len, obj.mix_len))
	fh.write("total\t%d\t%d\t%d\t%d\t%d\t%d\n" % (nohit_tot, ungap_nohit_tot, coll_tot, exp_tot, inv_tot, mix_tot))
	##top ten for each category (add later)

def makeBarGraph(assm1_chrom_list, assm1_dict, assm1_name, assm2_chrom_list, assm2_dict, assm2_name, out_fi, data_type):
	#check that chrom lists are the same- they should be but better to check
	err=0
	if [item for item in assm1_chrom_list if item in assm2_chrom_list]:
		pass
	else:
		logging.critical("%s not in both lists" % item)
		err += 1
	if err>0:
		logging.error("Chromosome lists not the same, not making graphic")
		return
	#set up lists for graphing
	assm1_list=[]
	assm2_list=[]
	title=""
	for seq in assm1_chrom_list:
		if data_type == "collapse":
			assm1_list.append(assm1_dict[seq].sp_len)
			assm2_list.append(assm2_dict[seq].sp_len)
			title="Collapse sequence: %s and %s" % (assm1_name, assm2_name)
		elif data_type == "expand":
			assm1_list.append(assm1_dict[seq].sp_only_len)
			assm2_list.append(assm2_dict[seq].sp_only_len)
			title="Expanded sequence: %s and %s" % (assm1_name, assm2_name)
		elif data_type == "no_hit":
			assm1_list.append(assm1_dict[seq].nohit_len)
			assm2_list.append(assm2_dict[seq].nohit_len)
			title="NoHit sequence: %s and %s" % (assm1_name, assm2_name)
		elif data_type == "ungap_nohit":
			assm1_list.append(assm1_dict[seq].ungap_nohit_len)
			assm2_list.append(assm2_dict[seq].ungap_nohit_len)
			title="Ungapped NoHit sequence: %s and %s" % (assm1_name, assm2_name)
		else:
			logging.error("Unknown data type, not making image: %s" % data_type)
			return
	#set up plot
	sns.set_style("ticks")
	sns.set_context("talk")
	plt.figure(figsize=(20,10), dpi=100)
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_scientific(False)
	X = np.arange(len(assm1_chrom_list))
	width=0.5
	y1_label=assm1_name
	y2_label=assm2_name
	plt.bar(X, assm1_list, width=0.5, facecolor='seagreen', edgecolor='none', label=y1_label)
	plt.bar(X+0.5, assm2_list, width=0.5, facecolor="blue", edgecolor='none', label=y2_label)
	plt.xticks(np.arange(len(assm1_list))+0.5, assm1_chrom_list, ha='center', size='22')
	plt.yticks(size="22")
	plt.xlabel('sequences', size='36')
	plt.ylabel('number of bases', size='36')
	plt.title(title, size='36')
	plt.legend(loc='upper left', prop={'size':24})
	sns.despine(top=True, right=True)
	plt.savefig(out_fi, dpi=100)

def main():
	#parse parameters
	parser = argparse.ArgumentParser(description="assm_align.py: process NCBI assm-assm alignments (also needs sequence report files)")
	parser.add_argument("--config", dest='cfg_file', help="path to config file (default is resources/assm_align_cfg.yml)")
	args = parser.parse_args()
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
	cfg_file=args.cfg_file
	if not cfg_file:
		cfg_file="resources/assm_align_cfg.yml"
	cfg_dict=yaml.load(open(cfg_file, 'r'))
	assm1=cfg_dict['input_files']['assm1']
	assm2=cfg_dict['input_files']['assm2']
	exclude_mt=cfg_dict['params']['exclude_mt']
	##create sequence objects
	#list to get sequence order of 'chrom' correct
	assm1_chrom_list=[]
	assm2_chrom_list=[]
	assm1_dict=parseSeqRep(assm1['seq_rpt'], assm1['name'], assm1['acc'], assm1_chrom_list, exclude_mt)
	logging.info("Read %s, sequences: %d, chromosomes: %d" % (assm1['name'], len(assm1_dict), len(assm1_chrom_list)))
	assm2_dict=parseSeqRep(assm2['seq_rpt'], assm2['name'], assm2['acc'], assm2_chrom_list, exclude_mt)
	logging.info("Read %s, sequences: %d, chromosomes: %d" % (assm2['name'], len(assm2_dict), len(assm2_chrom_list)))
	##parse alignment report
	logging.info("Processing %s" % assm1['name'])
	parseAlignReport(assm1['align_rpt'], assm1['name'], assm1_dict)
	logging.info("Procssing %s" % assm2['name'])
	parseAlignReport(assm2['align_rpt'], assm2['name'], assm2_dict)
	##produce stats
	assm1_stats_out=cfg_dict['output_files']['assm1']['stats']
	if not assm1_stats_out == False:
		logging.info("Writing stats file: %s" % assm1_stats_out)
		stats1_out=open(assm1_stats_out, 'w')
		writeStats(stats1_out, assm1['name'], assm2['name'], assm1_dict)
		stats1_out.close()
	assm2_stats_out=cfg_dict['output_files']['assm2']['stats']
	if not assm2_stats_out == False:
		logging.info("Writing stats file: %s" % assm2_stats_out)
		stats2_out=open(assm2_stats_out, 'w')
		writeStats(stats2_out, assm2['name'], assm1['name'], assm2_dict)
		stats2_out.close()
	##produce bed files if desired
	##I'm sure there is a better way, but brute forcing it now
	if cfg_dict['params']['make_bed'] == True:
		assm1_nohit_bed=cfg_dict['output_files']['assm1']['no_hit_bed']
		assm1_collapse_bed=cfg_dict['output_files']['assm1']['collapse_bed']
		assm1_expand_bed=cfg_dict['output_files']['assm1']['expand_bed']
		assm1_inv_bed=cfg_dict['output_files']['assm1']['inv_bed']
		assm1_mix_bed=cfg_dict['output_files']['assm1']['mix_bed']
		logging.info("Making assembly 1 beds")
		makeBed(assm1_nohit_bed, assm1_dict, "nohit")
		makeBed(assm1_collapse_bed, assm1_dict, "collapse")
		makeBed(assm1_expand_bed, assm1_dict, "expand")
		makeBed(assm1_inv_bed, assm1_dict, "inv")
		makeBed(assm1_mix_bed, assm1_dict, "mix")
		assm2_nohit_bed=cfg_dict['output_files']['assm2']['no_hit_bed']
		assm2_collapse_bed=cfg_dict['output_files']['assm2']['collapse_bed']
		assm2_expand_bed=cfg_dict['output_files']['assm2']['expand_bed']
		assm2_inv_bed=cfg_dict['output_files']['assm2']['inv_bed']
		assm2_mix_bed=cfg_dict['output_files']['assm2']['mix_bed']
		logging.info("Making assembly 2 beds")
		makeBed(assm2_nohit_bed, assm2_dict, "nohit")
		makeBed(assm2_collapse_bed, assm2_dict, "collapse")
		makeBed(assm2_expand_bed, assm2_dict, "expand")
		makeBed(assm2_inv_bed, assm2_dict, "inv")
		makeBed(assm2_mix_bed, assm2_dict, "mix")

	##produce graphs is desired, and only if chromosomes are available
	if len(assm1_chrom_list)>0 and len(assm2_chrom_list) >0:
		logging.info("Starting image production as both assemblies have chromosomes")
		collapse_img=cfg_dict['output_files']['comp_img']['both_collapse']
		expand_img=cfg_dict['output_files']['comp_img']['both_expand']
		nohit_img=cfg_dict['output_files']['comp_img']['both_nohit']
		ungap_nohit_img=cfg_dict['output_files']['comp_img']['both_ungap_nohit']
		if not collapse_img == False:
			logging.info("Making collapse image")
			makeBarGraph(assm1_chrom_list, assm1_dict, assm1['name'], assm2_chrom_list, assm2_dict, assm2['name'], collapse_img, "collapse")
		if not expand_img == False:
			logging.info("Making expansion image")
			makeBarGraph(assm1_chrom_list, assm1_dict, assm1['name'], assm2_chrom_list, assm2_dict, assm2['name'], expand_img, "expand")
		if not nohit_img == False:
			logging.info("Making no hit image")
			makeBarGraph(assm1_chrom_list, assm1_dict, assm1['name'], assm2_chrom_list, assm2_dict, assm2['name'], nohit_img, "no_hit")
		if not ungap_nohit_img == False:
			logging.info("Making ungap no hit image")
			makeBarGraph(assm1_chrom_list, assm1_dict, assm1['name'], assm2_chrom_list, assm2_dict, assm2['name'], ungap_nohit_img, "ungap_nohit")


if __name__=="__main__":
	main()
