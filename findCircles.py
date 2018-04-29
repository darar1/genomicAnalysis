#!/usr/bin/env python3
# Name: Dara Rouholiman (drouholi)
# Group Members: None
"""
 findCircles finds contigs that have been treated as unsable reads by the assembler. It reads off psl files, Not original sequencing
 data.Basically, it finds Pairednames (which has either /1 or /2 at the end) but in different strands. Hence, finds read-pair that
 supports a circular hypothesis.findCircles reads from stdin and writes to stdout. Simply run this program from the Terminal:
 python3 findCircles.py <EEV14-Vf.filtered.psl>  candidates.txt
 
Example:
 Input:<EEV14-Vf.filtered.psl>
 Output:candidates.txt
 
NODE_190026_length_10610_cov_28.822809 46
NODE_451325_length_14821_cov_22.767965 33
NODE_333359_length_18863_cov_27.166304 24
NODE_38071_length_10300_cov_21.263689 14
NODE_340060_length_11744_cov_14.423620 12
NODE_326825_length_36906_cov_33.828186 8
NODE_444083_length_13906_cov_16.088594 7
NODE_290471_length_11815_cov_13.465341 7
NODE_219343_length_14447_cov_14.828892 5
NODE_348878_length_14034_cov_14.812812 4
NODE_198000_length_53530_cov_17.185019 3
NODE_504247_length_13693_cov_12.981158 2
NODE_481727_length_12120_cov_37.470791 2

"""
import sys
class PSLreader :
    '''
    Class to provide reading of a file containing psl alignments
    formatted sequences:
    object instantiation:
    myPSLreader = PSLreader(<file name>):
 
    object attributes:
    fname: the initial file name
 
    methods:
    readPSL() : reads psl file, yielding those alignments that are within the first or last
                1000 nt

    readPSLpairs() : yields psl pairs that support a circular hypothesis 
        
    Author: David Bernick
    Date: May 12, 2013
    '''
    
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        
        self.fname = fname
            
    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
 
    def readPSL (self):
        '''
        using filename given in init, returns each filtered psl records
        that contain alignments that are within the terminal 1000nt of
        the target. Incomplete psl records are discarded.
        If filename was not provided, stdin is used.
 
        This method selects for alignments that could may be part of a
        circle.
 
        Illumina pairs aligned to the top strand would have read1(+) and read2(-).
        For the bottoms trand, read1(-) and read2(+).
 
        For potential circularity,
        these are the conditions that can support circularity:
        read1(+) near the 3' terminus
        read1(-) near the 5' terminus
        read2(-) near the 5' terminus
        read2(+) near the 3' terminus
 
        so...
        any read(+) near the 3', or
        any read(-) near the 5'
    
        '''
    
        nearEnd = 1000   # this constant determines "near the end"
        with self.doOpen() as fileH:
            
            for line in fileH:
                pslList = line.split()
                if len(pslList) < 17:
                    continue
                tSize = int(pslList[14])
                tStart = int(pslList[15])
                strand = str(pslList[8])    #Strand (+/-) is located at the 9th list
                contig = str(pslList[13])   #Contig/Node is located at the 14th list
                ReadPair = str(pslList[9])  #ReadPair/Query is the 10th list
                
                if strand.startswith('+') and (tSize - tStart > nearEnd):
                    continue
                elif strand.startswith('-') and (tStart > nearEnd):
                    continue
                
                #yield line
                
                #Getting the ReadPair, contig, and strand defined above instead of the whole line
                yield ReadPair, contig, strand
                
    def readPSLpairs (self):
        """
        Yields psl pairs that support the circular contig hypothesis.Finding Same Contigs (1 and 2s) but at different
        strands by using the yield of the readPSL() method and the "Pairwise Group Collect" as the template. 

        """
        thisGroup = None                            # making an arbitrary object equal none
        thisGroupList = [[], []]                    # making a list of 2 lists 
        
        for ReadPair, node, strand in self.readPSL(): # getting the 3 yielded objects from the readPSL() method 
            read = ReadPair.split("/")                # getting rid of the "/" in the ReadPair
            Read = read[0]                            # defineing Read as only the first part of the ReadPair
            num = int(read[1]) -1                     # num is just the number in ReadPair
            if Read != thisGroup:                                               #when Read is not thisGroup(none)
                #Just like is "Pairwise Group Collect"
                for leftContig, leftStrand in thisGroupList[0]:                 #left items in the first list
                    
                    for rightContig, rightStrand in thisGroupList[1]:           #right items in the 2nd list
                        
                        if leftContig ==  rightContig:                          #when the left item (contig) is equal to the right ite(contig)
                            
                            if leftStrand != rightStrand:                       #when the left item (2nd item, Strand) is not equal to the right item (strand)
                                
                                yield leftContig                                #yield either left or right equal item (contig)

                thisGroup = Read                                #Now, defining/equalling thisGroup as the Read (ReadPair minus the number)
                thisGroupList = [[], []]                        #empty the list of the two list 
            thisGroupList[num].append([node, strand])           #adding the the contig and the stand to their right list
            #Doing the forloop again to fill up the lists this time
        for leftContig, leftStrand in thisGroupList[0]:                     #left items in the first list
            for rightContig, rightStrand in thisGroupList[1]:               #right items in the 2nd list
                if leftContig == rightContig:                               #when the left item (contig) is equal to the right ite(contig)
                    if leftStrand != rightStrand:                           #when the left item (2nd item, Strand) is not equal to the right item (strand)
                        yield leftContig                                    #yield either left or right equal item (contig)
            

def main():
    """
    Defining the the PSL file, running it through PSLreader class.
    counting the pairreads that appear more than 1
    Outputing/printing the resulted contigs

    """
    
    #Opening the psl file through PSLreader class
    #*stdin is expected to be used: USE the terminal with name of the input psl file*
    myPSLreader = PSLreader('') # can either have the name of the file here or not and just use stdin in the Terminal
    test = myPSLreader.readPSLpairs()                # getting the readPSLpair() method to work on our file 


    contigCounter = {}                               # make a dictionary for contigs and put the pairreads that apppear more than 1 in that as the value.
    for contig in test:                      
        try:
            contigCounter[contig] += 1               # add one for every appearance 
        except KeyError:
            contigCounter[contig] = 2                # KeyError exception for the count of 2s, hence the program stops at two and doesn't count the 1s
    #sort the output and the count by decreasing value.
    for key, value in sorted(contigCounter.items(), key=lambda x : x[1], reverse = True):  
        print('%s %s' % (key, value)) #Print the Node and the count
main()
