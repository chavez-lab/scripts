import gzip
import argparse
import gzip

# Parse input params
parser = argparse.ArgumentParser()
parser.add_argument('file',help='significances.txt.gz file')
parser.add_argument('q', type=float, help='qvalue threshold')
parser.add_argument('resolution',type=int, help='Resolution that FitHiC was run at.')
args = parser.parse_args()

def fithic2bedpe(file,q,resolution):
    '''
    Converts fithic format to bedpe.
    fithic format:
    chr1  fragmentMid1  chr2  fragmentMid2  contactCount   p-value  q-value   bias1     bias2  ExpCC
    bedpe format:
    chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2
    '''
    with gzip.open(file,'rb') as f:
        f.readline()
        for line in f:
            line=line.strip().split(b'\t')
            #(chr1,fragmentMid1,chr2,fragmentMid2,contactCount,pval,qval,bias1,bias2,ExpCC) = line
            line = list(map(lambda x:x.decode("utf-8"),line))
            if line[6] == 'nan' or float(line[6]) > q:
                continue
            start1 = int(line[1]) - resolution//2
            end1 = int(line[1]) + resolution//2
            start2 = int(line[3]) - resolution//2
            end2 = int(line[3]) + resolution//2
            print('\t'.join(map(str,[line[0],start1,end1,line[2],start2,end2,"",line[6]])))

fithic2bedpe(args.file,args.q,args.resolution)