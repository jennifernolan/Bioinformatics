'''
    Jennifer Nolan
    C16517636
    Assignment porgram for part 5
'''

def main():
    #name of file being read from
    Filename = 'DNA_Seq.fasta'
    
    #open specified file to read its contents
    Fp1 = open(Filename, 'r')
    Fp2 = open(Filename, 'r') # opened a second time to get description line
    
    #read the entire file without the descriptor line
    Data = Fp1.readlines()[1:]
    #read the description line
    Description = Fp2.readlines(1)
    
    #format the description line
    Descript = ''.join(Description)
    
    descript = Descript.split('\n')
    
    description = ''.join(descript)
    
    #remove end of line markers and convert DNA sequence into continuous sequence
    DataString = ''.join(Data)
    
    DnaSeq = DataString.split('\n')
    
    Dna = ''.join(DnaSeq)
    
    #get the reverse of the DNA sequence
    ReverseSeq = reverse(Dna)
    
    #get the compliment of the original DNA sequence using the reversed sequence
    compliment = DnaCompliment(ReverseSeq)
    
    #get reading frames 
    posDNA, posRF, negComp, negRF = ReadingFrames(Dna, compliment)
    
    #convert to Amino Acid sequence
    AminoAcidSeq, AAArray = DnaToAA(posDNA, negComp)
    
    #get orfs and start and stop positions within amino acid sequence
    start, stop, orf = startStopCodon(AminoAcidSeq, AAArray)
    
    #display file information
    print("This program will read the contents of the following file:")
    print(Filename)
    print("\nThe following is the descriptor line of the file:")
    print(description)
    print("\nThe folowing is the DNA split into its reading frames:")
    for i in range(0, len(posDNA)):
        print(posRF[i])
        for j in range(0, len(posDNA[i])):
            print(posDNA[i][j])
    for i in range(0, len(negComp)):
        print(negRF[i])
        for j in range(0, len(negComp[i])):
            print(negComp[i][j])
    
    #display orf information
    ORFInfo(start, stop, orf, AminoAcidSeq)
    
#function that returns the reverse of the original DNA sequence
def reverse(DNA):
    return DNA[::-1]      
  
#function that returns the compliment of the original DNA sequence by using the reversed sequence
def DnaCompliment(seq):
    #dictionary to hold the compliments of each nucleotide
    complimentDict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    DnaComp = list(seq)
    
    #search for the compliment of each element of the reversed sequence using the dictionary declared above
    for index in range(len(DnaComp)):
        DnaComp[index] = complimentDict[DnaComp[index]]
    
    #return the compliment as a string
    return ''.join(DnaComp) 
    
    
def ReadingFrames(DNA, Comp):
    #list containing the 3 positive reading frames for the original DNA sequence
    positivekey = ["+1", "+2", "+3"]
    posRF = []
    
    #divide the original DNA sequence into 3 and add to a reading frame array/list
    for n in range(0, 3):
        posRF.append([DNA[n::]])
    
    #list containing the 3 negative reading frames for the compliment DNA sequence
    negativekey = ["-1", "-2", "-3"]
    negRF = []
    
    #divide the compliment DNA sequence into 3 and add to a reading frame array/list
    for n in range(0, 3):
        negRF.append([Comp[n::]])
    
    return posRF, positivekey, negRF, negativekey

#function to convert DNA reading frames into amino acids   
def DnaToAA(DNA, Comp):
    AminoAcidList = []
    AminoAcidSeq = ''
    
    # DNA <-> AA translation table: CodonTable
    CodonTable = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        }
    
    DNACodons = []
    codons = []
    
    i = 0 
    
    #divide original DNA sequence into codons, groups of three, for each reading frame
    while i < len(DNA):
        for n in range(0, len(DNA[i])):
            Frames = []
            Frames = DNA[i][n]
            for j in range(0, len(Frames), 3):
                codons.append(Frames[j:j+3])
        DNACodons.append(codons)
        codons = []
        i = i + 1
    
    i = 0 
    
    #divide compliment DNA sequence into codons, groups of three, for each reading frame
    while i < len(Comp):
        for n in range(0, len(Comp[i])):
            Frames = []
            Frames = Comp[i][n]
            for j in range(0, len(Frames), 3):
                codons.append(Frames[j:j+3])
        DNACodons.append(codons)
        codons = []
        i = i + 1
    
    AA = []
    AAArray = []
    
    #get amino acid sequence and also the amino acid array for each frame -> so that when getting ORFs later there is no run into next amino acid from another frame
    for i in range(0, len(DNACodons)):
        for j in range(0, len(DNACodons[i])):
            if DNACodons[i][j] in CodonTable:
                AminoAcid = CodonTable[DNACodons[i][j]]
                AminoAcidList.append(AminoAcid)
                AA.append(AminoAcid)
                AminoAcidSeq = ''.join(AminoAcidList)
        AAArray.append(AA)
        AA = []
    
    return  AminoAcidSeq, AAArray
    
#function to get all start and stop positions for Amino Acid sequence and find all ORFs in amino acid sequence 
def startStopCodon(AASeq, AAArray):
    
    AAStart = []
    AAStop = []
    ORFSeq = []
    ORF = []
    
    #go through each amino acid array and find if it has an ORF
    for i in range(0, len(AAArray)):
        j = 0
        while j < len(AAArray[i]):
            #once a start amino acid is found continue adding to ORF array until either a stop amino acid is found or the end of the array/frame is reached
            if AAArray[i][j] == 'M':
                ORFSeq.append(AAArray[i][j])
                #move to the position after the start amino acid
                k = j + 1
                while AAArray[i][k] != '*' and k + 1 != len(AAArray[i]):
                    #keep adding the amino acids to the ORF until either a stop amino acid is reached or the end of the amino acid array/frame is reached
                    ORFSeq.append(AAArray[i][k])
                    k = k + 1 
                #add the found ORF to an array that stores all the ORFs found
                ORF.append(ORFSeq)
                ORFSeq = []
                #ensures that the next loop starts from the value after the stop position
                j = k
            else:
                j = j + 1
                
                
    i = 0   
    #go through the amino acid sequence and find the start and stop positions
    while i < len(AASeq):
        #once a start amino acid is found add position to start array and once stop amino acid found add position to stop array
        if AASeq[i] == 'M':
            AAStart.append(i + 1)
            j = i + 1
            while AASeq[j] != '*' and j < len(AASeq):
                j = j + 1
            m = j
            AAStop.append(m + 1)
            i = m
        else: 
            i = i + 1
    
    return AAStart, AAStop, ORF
  
#function to display all the ORF information gathered from this program  
def ORFInfo(start, stop, orf, AminoAcidSeq):
    print('\nThe following is the amino acid sequence of the file:')
    print(AminoAcidSeq)
    
    count = 1
    for i in range(0, len(start)):
        if len(orf[i]) > 10:
            print("\nORF: " + str(count) + "; The ORF start is: " + str(start[i]) + "; The ORF end is: " + str(stop[i]) + "; The length is, in amino acid: " + str(len(orf[i])) + "; The length is, in nucleotide: " + str(len(orf[i]) * 3 + 3))
            print('\n')
            orfstring = ''.join(orf[i])
            print(orfstring)
            count = count + 1
        
main()

'''
    TEST PLAN:
    Compare the ORF amino acid results with the ORF amino acid results from the NCBI site

'''