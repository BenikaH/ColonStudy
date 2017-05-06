#TODO:Fix last 10
from __future__ import print_function
import sys

def get_mirna_names(filename):
    names = []
    with open(filename, 'r') as f:
        for line in f:
            fields = line.split('\t')
            names.append(fields[0])
            
    return names

def get_positions(names, filename, outputfile):
    with open(filename, 'r') as f:
        with open(str(outputfile), 'w') as fo:
            print('miRNA\tchr\tstart\tend',file=fo)
            for line in f:
                if line.startswith('#') or 'primary-transcript' in line:
                    continue
                line = line.rstrip()
                fields = line.split('\t')
                chromosome = fields[0]
                start = fields[3]
                end = fields[4]
                name = fields[-1].split('=')[-2].split(';')[-2]
                id = fields[-1].split('=')[1].split(';')[0].split('_')[-1]
                #Extracts mosts of miRNA positions, still needs to be tweaked
                if name in names:
                    print(name + '\t' + chromosome + '\t' + start + '\t' + end, file=fo)
                else:
                    if (name + '-3p') in names:
                        print(name + '-3p' + '\t' + chromosome + '\t' + start + '\t' + end, file=fo)
                    if (name + '-5p') in names:
                        print(name + '-5p'+ '\t' + chromosome + '\t' + start + '\t' + end, file=fo)
                    nameArray = name.split('-')
                    if len(nameArray)==4 and (nameArray[0]+'-'+nameArray[1]+'-'+nameArray[2]+'-'+id+'-'+nameArray[3]) in names:
                        print(nameArray[0]+'-'+nameArray[1]+'-'+nameArray[2]+'-'+id+'-'+nameArray[3]+ '\t' + chromosome + '\t' + start + '\t' + end, file=fo)
                    if len(nameArray)==3 and (nameArray[0]+'-'+nameArray[1]+'-'+nameArray[2]+'-'+id+'-3p') in names:
                        print(nameArray[0]+'-'+nameArray[1]+'-'+nameArray[2]+'-'+id+'-3p'+ '\t' + chromosome + '\t' + start + '\t' + end, file=fo)
                    if len(nameArray)==3 and (nameArray[0]+'-'+nameArray[1]+'-'+nameArray[2]+'-'+id+'-5p') in names:
                        print(nameArray[0]+'-'+nameArray[1]+'-'+nameArray[2]+'-'+id+'-5p'+ '\t' + chromosome + '\t' + start + '\t' + end, file=fo)


def main():
    miRNAnames = get_mirna_names(sys.argv[1])
    get_positions(miRNAnames, sys.argv[2], sys.argv[3])

if __name__ == '__main__':
    main()

