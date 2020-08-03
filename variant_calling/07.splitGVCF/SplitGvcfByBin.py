from optparse import OptionParser
import csv
import sys
import gzip
import copy

def returnRef ( chr = 'chr1A.1', loc = 1) :
    d = {}
    with open('break.ref', mode='r') as f:
        f_csv = csv.reader(f, delimiter=',')
        for row in f_csv:
            d[row[0]] = row[1]
        return d[chr + '_' + str(loc)]

def IntTo2Alpha ( n ) :
    return chr(int(n/26)%26+97) + chr(n%26+97)

def SplitChrPosByBin (input_fn, bin_size, chr_col = 0, pos_col = 1, ref_col = 3, alt_col = 4, endpos_col = 7,
                      output_prefix = "", output_suffix = ".vcf.gz") :
    try:
        if input_fn :
            IN = gzip.open(input_fn, 'r')
        else :
            IN = sys.stdin
        #
    except IOError:
        sys.stderr.write( "[Error] file %s cannot be open." % input_fn)
        exit(-1)
    #
    line = IN.readline().decode()
    header_str = ""

    while line[0] == '#' :
        header_str = header_str + line
        line = IN.readline().decode()
    #
    ID = 0
    BinStart = 0
    BinEnd = bin_size
    #
    cur_chr = line.strip().split("\t")[chr_col]
    #
    OUT = open(output_prefix + cur_chr + "." + IntTo2Alpha(ID)  + output_suffix, "w")
    OUT.write(header_str)

    while line:
        tokens = line.strip().split("\t")
        alt = tokens[alt_col]
        chr = tokens[chr_col]
        pos = int( tokens[pos_col] )
        if alt == '<NON_REF>' :
            pos_end = int( tokens[endpos_col].split("=")[1] )

            # if chr == cur_chr :
            # for now all gvcf files have same chr, so pass.

            if pos_end <= BinEnd :
                OUT.write(line)
            else :
                if pos > BinEnd :
                    OUT.close()
                else :
                    tokens_last = copy.deepcopy(tokens)
                    tokens_first = copy.deepcopy(tokens)
                    tokens_last[endpos_col] = 'END=' + str(BinEnd)
                    OUT.write('\t'.join(tokens_last))
                    OUT.close()
                    #
                    tokens_first[pos_col] = str(BinEnd + 1)
                    tokens_first[ref_col] = returnRef(chr, BinEnd + 1)
                    line = '\t'.join(tokens_first) + '\n'
                #
                ID += 1
                BinStart = BinEnd
                BinEnd += bin_size
                #
                OUT = open(output_prefix + chr + "." + IntTo2Alpha(ID)  + output_suffix, "w")
                OUT.write(header_str)
                OUT.write(line)
        else :
            if pos > BinEnd :
                OUT.close()
                ID += 1
                BinStart = BinEnd
                BinEnd += bin_size
                #
                OUT = open(output_prefix + chr + "." + IntTo2Alpha(ID)  + output_suffix, "w")
                OUT.write(header_str)
            #
            OUT.write(line)
        #
        line = IN.readline().decode()



def main():
    usage = "Usage: SplitChrPosByBin -i <input> -p <prefix[.chr.]> -s <[.chr.]suffix>\n" \
            "Description: Split the files by each chromosomes. \n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2017-09-30"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="INPUT", default=None,
                      help="Input file, CGmap or ATCGmap foramt, "
                           "use STDIN when not specified."
                           "(gzipped if end with \'gz\').", metavar="FILE")
    parser.add_option("-B", dest="bin_size", default=1000000,
                      help="Size of Bin [default: %default]", metavar="INT")
    parser.add_option("-C", dest="chr_col", default=1,
                      help="The column of chromosome in input file [default: %default]",
                      metavar="INT")
    parser.add_option("-P", dest="pos_col", default=2,
                      help="The column of position in input file [default: %default]",
                      metavar="INT")
    parser.add_option("-R", dest="ref_col", default=4,
                      help="The column of ref allel in input file [default: %default]",
                      metavar="INT")
    parser.add_option("-A", dest="alt_col", default=5,
                      help="The column of alt allel in input file [default: %default]",
                      metavar="INT")
    parser.add_option("-E", dest="endpos_col", default=8,
                      help="The column of end position in input file [default: %default]",
                      metavar="INT")
    parser.add_option("-p", dest="output_prefix", default="",
                      help="The prefix for output file", metavar="STRING")
    parser.add_option("-s", dest="output_suffix", default="g.vcf",
                      help="The suffix for output file "
                           "(gzipped if end with \'gz\').", metavar="STRING")
    #
    (options, args) = parser.parse_args()
    #
    if options.output_suffix == "" :
        output_suffix = ""
    else :
        output_suffix = "." + options.output_suffix
    #
    if options.output_prefix == "" :
        output_prefix = ""
    else :
        output_prefix = options.output_prefix + "."
    #

    SplitChrPosByBin(options.INPUT, int(options.bin_size),
                     int(options.chr_col)-1, int(options.pos_col)-1, int(options.ref_col)-1, int(options.alt_col)-1, int(options.endpos_col)-1,
                     output_prefix, output_suffix )

#
# ===========================================
if __name__ == "__main__":
    main()