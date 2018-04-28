#!/usr/bin/env python3
# Name: Dara Rouholiman (drouholi)
# Group Members: None
"""
Compares GC content, aaComposition and relative codon of 2 separate genomes.
This program is very similar to genomeAnalyzer which calculates composition stats
on only one genome.

Required Module: sequenceAnalyzer.py

Input: 2 different genomic fasta files

output:
sequence lengths: Genome1 = 2.21 Mb  Genome2 = 2.21 Mb

GC contents:      Genome1 = 55.7%    Genome2 = 55.7%

      Genome1                          Genome2

AAA : K  17.4 ( 12810)       AAA : K  17.4 ( 12810)
AAC : N  17.6 ( 12962)       AAC : N  17.6 ( 12962)
AAG : K  35.6 ( 26206)       AAG : K  35.6 ( 26206)
AAU : N   6.0 (  4423)       AAU : N   6.0 (  4423)
ACA : T   9.9 (  7308)       ACA : T   9.9 (  7308)
ACC : T  12.6 (  9253)       ACC : T  12.6 (  9253)
ACG : T  13.1 (  9655)       ACG : T  13.1 (  9655)
ACU : T   8.2 (  6029)       ACU : T   8.2 (  6029)
.
.
.

"""
def main():
    """
    Function imports necessary classes from sequenceAnalysis.py module, uses
    FastAreader class to open and read the fasta files. the main() defines
    two object through the NucParams class, uses the methods in NucParams class
    to get the counts of the characters (nucleotide, codon, amino acid).
    """
    #imports FastAreader and NucParams classes from sequenceAnalyzer module
    from sequenceAnalysis import FastAreader, NucParams
    #read the 1st fasta file by FastAreader
    genome1 = FastAreader('vulgaris.fasta')
    #read the 2nd fasta file by FastAreader
    genome2 = FastAreader('Shewanella.fasta')
    #Makes an object called seq1 form NucParams class
    seq1 = NucParams('')
    #Makes an object called seq2 form NucParams class
    seq2 = NucParams('')
    # for head, sequence in the 1st genomeFasta file
    for headx, seqx in genome1.readFasta():
        #use addSequence method to update the dictionaries
        seq1.addSequence(seqx)
    
        #define GC content for the 1st genome
        GC1 = ((seq1.nucDic.get('G') + seq1.nucDic.get('C'))\
                                  
                                  /(seq1.nucCount())*100)

    print(seqx)
    # for head, sequence in the 1st genomeFasta file
    for heady, seqy in genome2.readFasta():
        #use addSequence method to update the dictionaries
        seq2.addSequence(seqy)
        #define GC content for the 2nd genome
        GC2 = ((seq2.nucDic.get('G') + seq2.nucDic.get('C'))\
                                  
                                  /(seq2.nucCount())*100)

### PRINT STATEMENTS ###
    #prints the sequence length of the two genomes
    print ("Genome1 = %s" %headx)
    print ("Genome2 = %s" %heady)
    print("sequence lengths: Genome1 = %.2f Mb  Genome2 = %.2f Mb"\
          %((seq1.nucCount()/(1000000)),(seq2.nucCount()/(1000000))))
    #blank line
    print("")
    #prints GC contents of the two genomes
    print ("GC contents:      Genome1 = %.1f%%    Genome2 = %.1f%%" \
           %(GC1, GC2))
    #blank line
    print("")
    #prints genome1 and genome2 for organization
    print("      Genome1                          Genome2")
    #blank line
    print ("")
    """
    readCounter = {}
    for nuc in seq1:
        read = seq[nuc:nuc+150]
        if read in seq2:
            readCounter[read] +=1
    print (readCounter)
    """
    #rnaCodonTable dictionary from the NucParams class is used for relative codon usage and amino acid composition 
    #sorts the rnaCodonTable and gets access to keys and values
    #do try:, exceptKeyError: format to address when we reach three stop codon with '-' value
    for keys, values in sorted(seq1.rnaCodonTable.items()) and sorted(seq2.rnaCodonTable.items()): #getting keys(codon)and values(aminoacid) from the rnaCodonTable
        try:

            xxx = keys # defining the key as the 3letter codon
            A = values # defining the values as one letter Amino acid
            D1 = seq1.codonDic[xxx] #for the 1st genome, getting the count of the codons from the codon Dictionary 
            aaCount1 = seq1.aaDic[A] #for the 1st genome, getting the amino acid count from the aminoacid dictionary
            F1 = ((D1/aaCount1)*1000) # for the 1st genome, frequency is codon count over amino acid count times 100
            D2 = seq2.codonDic[xxx] # for the 2nd genome, getting the count of the codons from the codon Dictionary 
            aaCount2 = seq2.aaDic[A]# for the 2nd genome, getting the amino acid count from the aminoacid dictionary
            F2 = ((D2/aaCount2)*1000) #relative frequencyfor the 2nd genome
        except KeyError:
            continue
        
        #prints two genomes codon usage right next to each other.    
        print ("%s : %s %5.1f (%6d)       %s : %s %5.1f (%6d)" % (xxx, A, F1, D1, xxx, A, F2, D2))
    

main()
