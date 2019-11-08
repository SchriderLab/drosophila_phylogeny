#!/usr/bin/env python3
import sys, argparse, os
from Bio import SeqIO
from Bio.Seq import Seq


def get_cds(file):
   records=SeqIO.to_dict(SeqIO.parse(file, "fasta"))
   for r in records.items():
       # r is tuple of seq name [0] and seq [1]
       cds_l=len(r[1].translate(to_stop=True))*3
       records[r[0]]=r[1][:cds_l]
   SeqIO.write(records.values(),file+".cds", "fasta")

def main():
        parser = argparse.ArgumentParser(description='Extract CDS from BUSCO DNA using Protein information')
        parser.add_argument( '-d',help = "DNA fasta file", dest = 'DNA')
        parser.set_defaults(input = None)
        args = parser.parse_args()
        #if args.input == None:
                #print("Error! No fasta file provided for argument -d")
                #sys.exit( -1 )
        get_cds(args.DNA)


if __name__ == "__main__":
    main()
