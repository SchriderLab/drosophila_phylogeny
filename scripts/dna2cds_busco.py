#!/usr/bin/env python3
import sys, argparse, os
from Bio import SeqIO
from Bio.Seq import Seq


def get_cds(file):
   records=SeqIO.to_dict(SeqIO.parse(file, "fasta"))
   for r in records.items():
       # r is tuple of seq name [0] and seq [1]
       #Check all reading frames
       cds_l0=len(r[1].translate(to_stop=True))*3
       cds_l1=len(r[1][1:].translate(to_stop=True))*3
       cds_l2=len(r[1][2:].translate(to_stop=True))*3
       lengths=[cds_l0,cds_l1,cds_l2] 
       first_index=lengths.index(max(lengths))
       last_index=max(lengths)  
       records[r[0]]=r[1][first_index:last_index]
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
