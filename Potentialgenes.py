#!/usr/bin/env python3
# Name: Dara Rouholiman (drouholi)
# Group Members: None
"""
Finds ORFs and Translates them to protein sequences. Which will be useful
when running blastp and looking for genes coding region in the genome.
input: genome.fasta
output: ORFs in form of proteins

##Pseudocode:
    1.import SeqIO and use the read() function.
    2. Use fastAreader from sequenceAnalysis or SeqIO.read from Biopython to read through the sequence
    3. Set the Input/output files with 'open method' and create a file handle or use the command line I/O
    4. Use a for loop to run the Orf_finding in everyframe and then make an object to transcribe and another one to translate 
    or in biopython use a for loop to find genes and translate as well
    5. Make one list of translated genes, In Biopython record 

"""
def findgenes():
    """
    Finds ORFs and Translates them to protein using modules from
    biopython 1.61 package.
    """
    from Bio import SeqIO
    record = SeqIO.read('Shewanella.fasta', 'fasta') # Must tell SEqIO what format is being read
    #print(record.seq)
    # The bacterial, Archaeal and Plant Plastid Code, and chloroplast proteins! Can be changed otherwise to table 4, etc
    table = 11
    min_pro_len = 100 #min Orf is set at 100 as a standard approach, can be changed
    #Using the Seq Object's split method to get a list of all the possible ORF translations in the 6 reading frames
    #We are counting the frames from the 5' end of each strand
    #the following code doesn't keep track of the locations in terms of the protein counting
    for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        for frame in range(3):
            for pro in nuc[frame:].translate(table).split("*"):
                if len(pro) >= min_pro_len:
                    print (">%s " % (pro)) # each record begins with >(greater sign) 

findgenes()
