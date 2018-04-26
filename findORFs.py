#!/usr/bin/env python3
# Name: Dara Rouholiman (drouholi)
# Group Members: None
'''
findORFs finds potential genes in a genome. It is an equivalent tool as the
ORF finder at NCBI database. The program accepts DNA sequences as an input
from any FASTA file and outputs the results in a text file.
User is able to change the options (longestGene : True or False, Start codons :
{ATG, TTG, GTG}, minGene = any integer)

Usage Example 
Input: tass2.fa
Output :tass2ORFdata-ATG-300.txt
tass2 NODE_159_length_75728_cov_97.549133
+1 57166..61908  4743
-1  8192..11422  3231
+2 65963..69004  3042
-3 14589..16862  2274
-2  2968.. 4872  1905
+1 64093..65952  1860
-3    30.. 1694  1665
'''

########################################################################
# CommandLine
########################################################################
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'findORFs finds potential genes just like NCBIs ORF finder software. ', 
                                             epilog = 'Start codon, stop codon, and length are variable', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        self.parser.add_argument('outFile', action = 'store', help='output file name') 
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=True, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= range(0, 2000), action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', nargs='?', help='start Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Main
# Here is the main program
# 
#
########################################################################
   

def main(myCommandLine=None):
    '''
    Implements the Usage exception handler that can be raised from anywhere in process.  

    '''
    if myCommandLine is None:
        myCommandLine = CommandLine([ 'vulgaris.fasta',
                                      'vulgarisORFdata-ATG-100.txt',
                                      '--longestGene',
                                      '--start=ATG',
                                      '--minGene=100'])
    else :
        myCommandLine = CommandLine(myCommandLine)
###### replace the code between comments.
        # myCommandLine.args.inFile has the input file name
        myCommandLine.args.inFile = 'vulgaris.fasta'
        # myCommandLine.args.outFile has the output file name
        myCommandLine.args.outFile = 'vulgarisORFdata-ATG-100.txt'
        # myCommandLine.args.longestGene is True if only the longest Gene is desired
        # myCommandLine.args.start is a list of start codons
        myCommandLine.args.start = ['ATG', 'GTG','CTG', 'TTG']
        # myCommandLine.args.minGene is the minimum Gene length to include
        myCommandLine.args.minGene = 100
######
        
    import sequenceAnalysis as x #importing sequenceAnalysis module with the new name 'x'
    #taking the fastafile (tass2.fa) through FastAreader
    tass = x.FastAreader (myCommandLine.args.inFile)
    output = myCommandLine.args.outFile #define the output file
    # making the output and preparing to write on it:
    with open(output, 'w') as ORFData:
        #getting the header from the file by using readFasta() method from the FastAreader class. 
        for head, seq in tass.readFasta() :
            print (head, file = ORFData) #note for the print statement to write (,file = the name) is essential
            orf = x.ORF() # running the object through ORF class, initializing the lists
            #using the ORFfinder method in the class to get the list of the lists of the ORFs
            geneData = orf.ORFfinder(seq)
            #let's combine all the 6 lists together
            dataList = geneData[0] + geneData[1] + geneData[2] + geneData[3] + geneData[4] + geneData[5]
            #sort the data 
            dataList.sort (key= lambda length: length[3], reverse=True)
            #order the genes by frame,startposition, endposition, and their length
            for genes in dataList:
                frame = genes[0] #frame is the first list
                start = genes[1] #start is the 2nd list
                stop = genes[2]  #stop is the 3rd list
                length = genes[3] #length is the 4th list
                #print statement, note (, file = theName) is essential to write on outputfile
                print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(frame, start, stop, length), file = ORFData)  

    if myCommandLine.args.longestGene:
        print ('longestGene is', str(myCommandLine.args.longestGene) )
    else :
        pass
#######
    
if __name__ == "__main__":
    main()


