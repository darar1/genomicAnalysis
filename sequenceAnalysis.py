#!/usr/bin/env python3
# Name: Dara Rouholiman (drouholi)
# Group Members: None
"""
The Module provides the composition of the genome to calculate statistics on it.

Classes:
    NucParams - provides essential methods for any sequence analysis
    ProteinParams- calculates physical/chemical properties of a protein sequence
    FastAreader- provides reading of a file containing one of more FASTA(written by David Bernick)


"""
class NucParams:
    """
    Class to provide counts and composition of a sequence (DNA,RNA, or a protein)

    methods:
    addSequence(): accepts additional seqeunces, presumed to start in frame 1
    aaComposition(): returns a dictionary of counts of the aminoacids
    nucComposition(): returns a dictionary of counts of valid nucleotides
    codonComposition(): returns a dictionary of counts of codons
    nucCount(): returns an integer of sum of the nucleotides' count

    """
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C', # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C', # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '*', 'UGA': '*', # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '*', 'UGG': 'W', # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R', # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S', # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G', # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G' # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}
 
    def __init__ (self, seq):
        """
        accepts a string of RNA or DNA or protein sequence
        initializes 3 zeroed dictionaries(amino acid, codon, nucleotide)
        to differ between composition of them .
        """
        #Amino acid dictionary with zero value foe each key
        self.aaDic = {
            'A': 0, 'G': 0, 'M': 0, 'S': 0, 'C': 0,
            'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
            'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
            'W': 0, 'F': 0, 'L': 0, 'R': 0, 'Y': 0
            }
        #Codon Dictionary with zero value for each key
        self.codonDic = {
        # U
        'UUU': 0, 'UCU': 0, 'UAU': 0, 'UGU': 0, # UxU
        'UUC': 0, 'UCC': 0, 'UAC': 0, 'UGC': 0, # UxC
        'UUA': 0, 'UCA': 0, 'UAA': 0, 'UGA': 0, # UxA
        'UUG': 0, 'UCG': 0, 'UAG': 0, 'UGG': 0, # UxG
        # C
        'CUU': 0, 'CCU': 0, 'CAU': 0, 'CGU': 0, # CxU
        'CUC': 0, 'CCC': 0, 'CAC': 0, 'CGC': 0, # CxC
        'CUA': 0, 'CCA': 0, 'CAA': 0, 'CGA': 0, # CxA
        'CUG': 0, 'CCG': 0, 'CAG': 0, 'CGG': 0, # CxG
        # A
        'AUU': 0, 'ACU': 0, 'AAU': 0, 'AGU': 0, # AxU
        'AUC': 0, 'ACC': 0, 'AAC': 0, 'AGC': 0, # AxC
        'AUA': 0, 'ACA': 0, 'AAA': 0, 'AGA': 0, # AxA
        'AUG': 0, 'ACG': 0, 'AAG': 0, 'AGG': 0, # AxG
        # G
        'GUU': 0, 'GCU': 0, 'GAU': 0, 'GGU': 0, # GxU
        'GUC': 0, 'GCC': 0, 'GAC': 0, 'GGC': 0, # GxC
        'GUA': 0, 'GCA': 0, 'GAA': 0, 'GGA': 0, # GxA
        'GUG': 0, 'GCG': 0, 'GAG': 0, 'GGG': 0 # GxG
        }
        #Nucleotide Dictionary with zero value for each key
        self.nucDic = {'A': 0, 'C': 0, 'T': 0, 'U': 0, 'G': 0, 'N':0}
        
    
    def addSequence (self, string):
        """
        accepts additional sequences and presumed to start at frame 1.
        In other words, updates each dictionary defined at __init__ with
        the right count. 
        """
        l = ''.join(string).split() #gets rid of whitespace
        self.string = ''.join(l).upper() #makes sure all characters are uppercase
        #for each nucleotide in the cleaned sequence
        for nuc in self.string: 
            if nuc in self.nucDic:
                self.nucDic[nuc] +=1    #add a count of 1 for each present nucleotide to its dictionary
                if nuc ==0: #for case when there is no U/T just count 1 (or raise a exception KeyError)
                    self.nucDic[nuc] +=1
                
                
        #note the sequence could be either RNA or DNA, choosing RNA to work with
        rnaSeq = self.string.replace('T','U') #translating to RNA for amino acid count
        #for each character(base) in every other 3 position
        for base in range(0,len(rnaSeq),3): 
            codon = rnaSeq[base:base+3] #A codon is 3 character(base) long in a RNA sequence
            if codon in self.codonDic: #if codon is present in the dicionary 
                self.codonDic[codon] +=1 #add a count of one to its value in codon Dictionary
                amino = NucParams.rnaCodonTable[codon] #amino acid is the codons value in rnaCodonTable
                for amino in self.aaDic: #for a present amino acid: 
                    try:
                    
                        self.aaDic[amino] +=1 #add a value of 1 to the key in amino acid Dictionary
                
                    except KeyError: #in case of '-' where the stop codons are not defined right
                        
                        continue #just ignore their present and move on
        
    def aaComposition(self):
        """
        returns a dictionary of counts over the valid 20amino acids.
        It decodes the codon to amino acids first then counts.
        
        """
        return self.aaDic
    
        
    def nucComposition(self):
        
        """
         returns a dictionary of counts of valid nucleotides found in the sequence.
         RNA nucleotides count as RNA nucleotides
         DNA nucleotides count as DNA nucleotides
         Any N bases is also counted, but invalid nucleotides are ignored
        """
        return self.nucDic
         

    def codonComposition(self):
        
         
        """
        returns a dictionary of counts of codons of valid nucleotides in form
        of RNA nucleotied (all U's), even when the given sequence is DNA!!
        """
        
        return self.codonDic
            

         
    def nucCount(self):
        
        """
         Returns an integer which is the sum of every valid nucleotide found in the
         sequence. equals == sum over the nucleotide composition dictionary
        """
        return sum(self.nucDic.values())


class ProteinParam :
    """
    This class calculates the physical and chemical properties of any protein sequence.

    initialized: a dictionary of aminoacids from the given protein string where
                 unexpected characters (B,J,O,U,X,Z) are ignored and the correct characters
                 (A, C, D, E, F, G,H, I, L, K, M, N, P, Q, R, S, T, V, Y, W) are guaranteed
                 to be uppercase. 
    methods: aaCount(), PI(), aaComposition(), charge(), molarExtinction(), mass Extinction()
             molecularweight()
    """
    
    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like

    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

 
    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189}
    
    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}
 
    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

# the __init__ method requires a protein string to be provided, either as a
# string or list of strings that will be concatenated
    def __init__ (self, protein):
        """
       a dictionary of aminoacids from the given protein string where
       unexpected characters (B,J,O,U,X,Z) are ignored and the correct characters
       (A, C, D, E, F, G,H, I, L, K, M, N, P, Q, R, S, T, V, Y, W) are guaranteed
        to be uppercase.
        """
        l = ''.join(protein).split()
        protString = ''.join(l).upper()
        protSeq = protString.replace('B','').replace('J', '').replace('O','').replace('U','').replace('X','').replace('Z','').replace('',',').split(',')
        self.cleanPro = protSeq[1:-1] #the first and last character of the list are empty
        self.aaDic = dict.fromkeys(self.cleanPro)
        
          
        

    def aaCount (self):
        """
        Returns a single integer count of valid aminoacids characters only.
        valid amino acid characters are :
        A, C, D, E, F, G,H, I, L, K, M, N, P, Q, R, S, T, V, Y, W.
        
        """
        aaNumber = 0 # setting the aa total numbers to zero
        #a list of valid amino acid charaters
        validaa = ["A","C","D","E","F","G","H","I","L","K","M","N","P","Q","R","S","T","V","Y","W"]
        for aa in protein: #for every aa characters in the protein string
            aaNumber += 1     #adds 1 to the count of amino acids(aaNumber)
            if aa not in validaa: 
                aaNumber -= 1 #If the character is not one of the valids it removes one from the count
                
        return (aaNumber) 
 
    def pI (self):
        """
        Estimates the theoretical isoelectric point by finding the pH(from the charge method).
        
        """
        pI = 0.0 #setting the PI at 0
        # while the pI is from 0 to less/equal than 14  (the range of pH)
        while pI <= 14:
            charge = ProteinParam.charge(protein, pI) # defining charge by charge method
            if charge <= 0: #If charge is less than 0 
                break #stops
            else:
                pI += .001 #else: just add.001 decimal and retry to see if the value is negative
            
        return pI
        
 
    def aaComposition (self) :
        """
        Returns a dictionary keyed by single letter amino acid code, and valued by integer count of valid amino acid characters in the protein
        sequence. Valid amino acids are: A, C, D, E, F, G,H, I, L, K, M, N, P, Q, R, S, T, V, Y, W. 
        """
        #getting the clean list of aminoacids of protein (same code as the _init_)
        l = ''.join(protein).split() 
        protString = ''.join(l).upper()
        protSeq = protString.replace('B','').replace('J', '').replace('O','').replace('U','').replace('X','').replace('Z','').replace('',',').split(',')
        cleanPro = protSeq[1:-1]
        #Making a dictionary with the aa codes as keys and their counts as values 
        dr = dict() #empty dictionary.
        aa = 'ACDEFGHILKMNPQRSTVYW' #valid aa's
        for c in aa: #for characters in valid amino acids
            dr[c] = 0 #add 0 values to the keys(all valid aa's)
        for c in cleanPro: #for characters in the list of cleaned protein
            dr[c] += 1 #add 1 to the count for every present character in the sequence
            
  
        return (dr) 
        
        
 
    def charge (self, pH):
        """
        Calculates the net charge on the protein at a specific pH. It uses the pKa of each charged
        amino acids and the Nterminus and Cterminus. With a required parameter named pH. 
        """
        #the dictionary created with the count of each aa in the sequence
        aaDic = ProteinParam.aaComposition(protein)

        # Making an empty list to sum all the charges of each character in
        chargeList = []
        #Adding charge of the individual aa  and Nter and Cterminus.
        chargeList.append(-1. / (1+(10**(ProteinParam.aaCterm-pH))))
        chargeList.append(-1.*aaDic['D'] / (1+(10**(ProteinParam.aa2chargeNeg['D']-pH))))
        chargeList.append(-1.*aaDic['E'] / (1+(10**(ProteinParam.aa2chargeNeg['E']-pH))))
        chargeList.append(-1.*aaDic['C'] / (1+(10**(ProteinParam.aa2chargeNeg['C']-pH))))
        chargeList.append(-1.*aaDic['Y'] / (1+(10**(ProteinParam.aa2chargeNeg['Y']-pH))))
        chargeList.append(1.*aaDic['H'] / (1+(10**(pH-ProteinParam.aa2chargePos['H']))))
        chargeList.append( 1.*aaDic['K'] / (1+(10**(pH-ProteinParam.aa2chargePos['K']))))
        chargeList.append( 1.*aaDic['R'] / (1+(10**(pH-ProteinParam.aa2chargePos['R']))))
        chargeList.append( 1. / (1+(10**(pH-ProteinParam.aaNterm))))
        pH_sum = sum(chargeList) #summing the charges
        return pH_sum
        
 
    def molarExtinction (self):
        """
        Returns how much light a protein absorbs at  280nm wavelength.
        """
        #the dictionary created with the count of each aa in the sequence
        aaDic = ProteinParam.aaComposition(protein) 
        #Get the amount of Y,W,and C's we have.
        Y = aaDic.get('Y')
        W = aaDic.get('W')
        C = aaDic.get('C')
        #to get Molar Extinction, multiply Y,W,C with their respective extiction coeff at 280nm
        mE = (Y*ProteinParam.aa2abs280['Y'])+(C*ProteinParam.aa2abs280['C'])+(W*ProteinParam.aa2abs280['W'])       
        return int(mE)
 
    def massExtinction (self):
        """
        calculates the Mass extinction coefficient from the Molar Extinction coefficient by dividing
        by the molecularWeight of the corresponding protein.
        
        Author: David Bernick
        """
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0
 
    def molecularWeight (self):
        """
        Calculates the molecular weight (MW) of the protein sequence. Sums up the weights of each
        valid amino acid in the sequence and excluding the waters that are released with peptide bond.
        """
        #getting the cleaned protein (all valid and uppercase characters)
        l = ''.join(protein).split()
        protString = ''.join(l).upper()
        protSeq = protString.replace('B','').replace('J', '').replace('O','').replace('U','').replace('X','').replace('Z','').replace('',',').split(',')
        cleanPro = protSeq[1:-1]
    
        results = [] #empty list to do the summation 
        for aa in cleanPro: #for amino acids in the clean protein
            aaWeight = ProteinParam.aa2mw.get(aa) #measure the weight of each amino acid from the aa2mw table
            results.append(aaWeight) #append the weights to the list 
            results_sum = sum(results) #sum the numbers in the list
            molecularW = results_sum - (ProteinParam.mwH2O*(ProteinParam.aaCount(protein)-1)) #exclude the water
            
            
                                     
        return (molecularW)

    
 
 
class FastAreader :
    '''
    Class to provide reading of a file containing one or more FASTA
    formatted sequences:
    object instantiation:
    FastAreader(<file name>):
 
    object attributes:
    fname: the initial file name
 
    methods:
    readFasta() : returns header and sequence as strings.
    Author: David Bernick
    Date: April 19, 2013
    '''
    def __init__ (self, fname):
        '''contructor: saves attribute fname '''
        self.fname = fname
 
    def readFasta (self):
        '''
        using filename given in init, returns each included FastA record
        as 2 strings - header and sequence.
        whitespace is removed, no adjustment is made to sequence contents.
        The initial '>' is removed from the header.
        '''
        header = ''
        sequence = ''
        
        with open(self.fname) as fileH:
            # initialize return containers
            header = ''
            sequence = ''
 
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()
 
            # header is saved, get the rest of the sequence
            # up until the next header is found
            # then yield the results and wait for the next call.
            # next call will resume at the yield point
            # which is where we have the next header
            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        # final header and sequence will be seen with an end of file
        # with clause will terminate, so we do the final yield of the data
        yield header,sequence
 
# presumed object instantiation and example usage
# myReader = FastAreader ('testTiny.fa');
# for head, seq in myReader.readFasta() :
#     print (head,seq)
