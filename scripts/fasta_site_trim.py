#!/usr/bin/env python 

from Bio import AlignIO,SeqIO
import sys, argparse
#print("ATTENTION:If codon alignmemnt,  proceed with caution ")
def trim(input,base_allow):
    alignment = AlignIO.read(input, "fasta")
    for column in range(0,alignment.get_alignment_length()):
        if len(alignment[:,int(column)])-alignment[:,int(column)].count("-") >= int(base_allow):
            if 'alignment_trimmed' in vars():
                x=int(column)
                y=column+1
                alignment_trimmed+=alignment[:,x:y]
            else:
                x=int(column)
                y=column+1
                alignment_trimmed=alignment[:,x:y]
    return(alignment_trimmed)


def main():
        parser = argparse.ArgumentParser(description='Remove sites with < N bases allowed from FASTA alignment')
        parser.add_argument( '--Nbase', help = "Remove sites with < N bases",dest='NBASE')
        parser.add_argument( '--input',help = "A fasta file", dest = 'FASTA')
        parser.set_defaults( list = None, input = None)
        args = parser.parse_args()
        if args.FASTA == None:
                print("Error! No fasta file provided for argument --input")
                sys.exit( -1 )
        tr=trim(args.FASTA,args.NBASE)
        SeqIO.write(tr,args.FASTA+".tr", "fasta")


if __name__ == "__main__":
    main()  
