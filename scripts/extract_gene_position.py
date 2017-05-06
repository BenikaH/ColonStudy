from __future__ import print_function
import sys

def get_positions(filename, outputfile):
    with open(filename, 'r') as f:
       with open(str(outputfile), 'w') as fo:
           print('gene\tchr\tstart\tend',file=fo)
           for line in f:
                line = line.rstrip()
                if 'transcript_name' not in line and '#!' not in line:
                    fields = line.split('\t')
                    chromosome = fields[0]
                    start = fields[3]
                    end = fields[4]
                    gene_name = fields[-1].split('"')[1]
                    print(gene_name + '\t' + chromosome + '\t' + start + '\t' + end, file=fo)
def main():
    get_positions(sys.argv[1], sys.argv[2])   

if __name__ == '__main__':
    main()
