import csv
import argparse
from Bio import SeqIO
from collections import defaultdict

def getLen(fi):
    len_hash={}
    for seq_record in SeqIO.parse(fi, "fasta"):
        len_hash[seq_record.id]= len(seq_record)
    return len_hash

def main():
    parser = argparse.ArgumentParser(description="make_seq_report.py: create seq report for assembly when NCBI omits this")
	parser.add_argument("--fasta", dest="fasta_file", help="path to fasta file")
    parser.add_argument("--asm_name", dest="asm_name", help="assembly name")
    parser.add_argument("--outdir", dest="outdir", help="output directory")
	args = parser.parse_args()
    asm_name=args.asm_name
    fasta=args.fasta
    outdir=args.outdir
    #get sequence lengths
    seq_len={}
    seq_len=getLen(fasta)
    #create output file
    out_file="%s/%s.assembly.txt" % (outdir, asm_name)
    out=open(out_file, 'w')
    out.write("# Assembly Name: %s\n" % asm_name)
    out.write("# Fake assembly report to replace missing NCBI one\n")
    out.write("#Assumptions: Scaffold level assembly, only primary unit\n")
    out.write("#Sequence-name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n")
    for seq in seq_len:
        out.write("%s\unlocalized-scaffold\tNA\tNA\tNA\tNA\tNA\tPrimary Assembly\t%d\tNA\n" % (seq, seq_len[seq]))
    out.close()


if __name__=="__main__":
	main()
