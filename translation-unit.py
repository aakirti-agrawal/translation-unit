dna=raw_input("Enter your DNA Sequence\n")

print "\n"

flag=0

#convert all alphabets into upper case
dna=dna.upper()

if(dna.count("U")>0):
	dna=dna.replace("U", "T")
	flag=1

#a is a list storing the lengths of all translation units     
a=[]

#dictionary for starting index of translation units corresponding to their lengths
b={}

#dictionary for all amino acid sequences
gcode={
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
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }


#list of start codons & stop codons
start_codon={'ATG', 'GTG'}
stop_codon={'TAG', 'TAA', 'TGA'}


t=0

#function to convert the corresponding codon sequence into amino acid sequence
def codon2protein(dna, gcode):
	for r in range(0,3):
		protein=''

		t=(len(dna)-r)%3

		for i in range(r, len(dna)-t, 3):
			codon=dna[i:i+3]
			protein+=gcode[codon]   # accessing dictionary
   		print "amino acid sequence in frame",r," ",protein
		print "\n"

codon2protein(dna, gcode)

print "\n"

#finding the translation units and their positions
def translation(dna, gcode):
	
	for k in range(0,3,1):
		i=k
		print "frame:",i
		t=(len(dna)-k)%3
		while(i<=len(dna)-t):
			j=i
			codon=dna[i:i+3]
			if codon in start_codon: # condition for start codon
				trans=''
				start=i+3
				j+=3
				codon=dna[start:start+3] 

				while codon not in stop_codon: #condition for stop codon
					trans+=codon
					j+=3
					codon=dna[j:j+3]

				j=j-1

				i=j+4
	
				length=j-start+1

				key=length
				b.setdefault(key, [])
				b[key].append(start)

				a.append(length)
				
				if(trans.count("T")>0 and flag==1):
					trans=trans.replace("T", "U")
			
				print "translation unit:",trans
				print "starting position:",start
				print "ending position:",j
				print "length of translation unit",length
				print "\n"

			else:
				i+=3
		
translation(dna, gcode)

print "Dictionary for starting index of translation units corresponding to their length in original sequence"

print b

print "\n"

#reversing the dna sequence
dna_reverse=dna[::-1]

#finding complementary of dna sequence
def compliment(dna_reverse):
    comp = []
    for i in dna_reverse:
        if i == "T":
            comp.append("A")
        if i == "A":
            comp.append("T")
        if i == "G":
            comp.append("C")
        if i == "C":
            comp.append("G")

    return ''.join(comp)

dna1=compliment(dna_reverse)

print "\n"

b={}

if(dna1.count("T")>0 and flag==1):
	dna1=dna1.replace("T", "U")

print "Reverse complementary dna strand:", dna1
print "\n"

print "For reverse complementary dna strand:"

print "\n"

if(dna1.count("U")>0 and flag==1):
	dna1=dna1.replace("U", "T")

codon2protein(dna1, gcode)

print "\n"

translation(dna1, gcode)

print "Dictionary for starting index of translation units corresponding to their length in reverse complementary sequence"

print b

print "\n"

print "Length of longest translation unit", max(a)

print "\n"

print "Position of longest translation unit: ",b[max(a)], " assuming first element of sequence to be at index 0"




		
