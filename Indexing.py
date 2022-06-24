import sys, getopt
import random
import os

def Find_SA_for_Genome(L=100,l=100,G_file='example1.fasta'):    
#Input:  G_file: Reference Genome File,  L: Divide genome to substrings of length L,   l : Overlap length.
#Output:   suffix:Suffix Array

    genome = load_genome(G_file)
    first_ele_suffix=len(genome)
    genome=genome.upper()
    genome=genome.replace('N',random.choice("ACTG"))  

    genome=genome+'$'
    SA_s=[]
    print("Genome length:", len(genome))
    for i in range(0,len(genome),L):
        S1=SACA_DNA_Alphabets(genome,i,i+L,l)
        SA_s.append(S1)
    print("Number of substrings:", int(len(genome)/L)+1)
    i=0
    while len(SA_s)!=1:
            if i==len(SA_s) or i==len(SA_s)-1:i=0
            SA1=SA_s[i]
            SA2=SA_s[i+1]
            SA_s.remove(SA1)      
            SA_s[i]=Sorted_2_SA(SA1,SA2,genome,k2=50)
            i=i+1
    SA_s[0].insert(0,first_ele_suffix) 
    head_tail = os.path.split(G_file)
    x=head_tail[1]+".SA.txt"
    file_sufix_array = open(x, 'w')
    file_sufix_array.write(str(SA_s[0]))
    file_sufix_array.close()



def SACA_DNA_Alphabets(t,S,E,l=100):
#Input: DNA sequence(t) (substring of genome t[S:E] )
#Output: Suffix Array for the substring
    suffix=[]
    E2=E
    E=E+l
    suffix.append(len(t)-1)
    
    if E>len(t)-1:
        E=len(t)-1

    F_t=Find_index(t, 'A',S,E)
    Tree=[]
    if (len(F_t)!=0):Tree.append(F_t)
    suffix=SA_char(suffix,t,Tree)


    F_t=Find_index(t, 'C',S,E)
    Tree=[]
    if (len(F_t)!=0):Tree.append(F_t)
    suffix=SA_char(suffix,t,Tree)


    F_t=Find_index(t, 'G',S,E)
    Tree=[]
    if (len(F_t)!=0):Tree.append(F_t)
    suffix=SA_char(suffix,t,Tree)


    F_t=Find_index(t, 'T',S,E)
    Tree=[]
    if (len(F_t)!=0):Tree.append(F_t)
    suffix=SA_char(suffix,t,Tree)


    suffix=suffix[1:]
    y=0
    while y<=E2-S:
        if y<len(suffix)and suffix[y]>=E2 and suffix[y]<=E2+l:
            suffix.remove(suffix[y])
        else:
            y=y+1
    return suffix




def load_genome(genome_file):
    #load genome all bases together and return string
    genome=""
    with open(genome_file) as content:
        for line in content:
            line = line.strip()
            if line!="" and line[0]!='>':
                genome+= line
    return genome


def Find_index(G, k,START,END):
#return all indexes of k alphabet in Substring of G returned as a pair of values -1 and y 
    x=[]
    for i in range(START,END):
        if G[i] == k:
            x.append((-1,(i)))
    return x



def SA_char(suffix,t,F_index):
#Input: suffix: Initial values of Suffix Array, t: DNA sequence , 	F_index: all indexes of k alphabet in S returned as a pair of values x and y.	
#Output: Final values of Suffix Array

    while(len(F_index)!=0):
        F_t=F_index[0]
        
        AA=[]    
        CC=[]
        GG=[]
        TT=[]
         
        c=F_t[0][0]*-1
        F_index.remove(F_t)
        j=len(F_t)

        if j==1:
            suffix.append(F_t[0][1])
            F_t=[]
            j=len(F_t)

        else:
            x=True
            while(j!=0 and x == True):
                for i in range(0,j): 
                    if (suffix.count(F_t[i][1]+c)!=0):
                       F_t[i]=((suffix.index(F_t[i][1]+c)),(F_t[i][1]))
                F_t.sort()
                if F_t[-1][0]>=0: 
                    i=0
                    while F_t[i][0]<0:
                        i+=1
                    j=i
                    while i< len(F_t):
                        suffix.append(F_t[i][1])
                        i+=1
                    F_t=F_t[0:j]
                else:
                    x=False
                    
                    j=len(F_t)
                    for i in range(0,j):
                        if t[F_t[i][1]+c]=='$':
                            suffix.append(F_t[i][1])
                        if t[F_t[i][1]+c]=='A': 
                           AA.append((-(c+1),(F_t[i][1])))
                        elif t[F_t[i][1]+c]=='C': 
                           CC.append((-(c+1),(F_t[i][1])))
                        elif t[F_t[i][1]+c]=='G':
                            GG.append((-(c+1),(F_t[i][1])))
                        else:
                           TT.append((-(c+1),(F_t[i][1])))
                
                    if len(TT)!=0:F_index.insert(0,TT)
                    if len(GG)!=0:F_index.insert(0,GG)
                    if len(CC)!=0:F_index.insert(0,CC)
                    if len(AA)!=0:F_index.insert(0,AA)
                    F_T=[]
    return suffix


def Sorted_2_SA(X,Y,t,k2=50):
#Input:  X: First Suffix Array 1, Y: Second Suffix Array 2 , 	t : Reference Genome , K2 :constant
#Output: Merged two Suffix Array

    kemer=k2
    kemer2=k2  
    z=[]
    l=0
    j=0
    while l<len(X) and j<len(Y):
            while t[X[l]:X[l]+kemer]==t[Y[j]:Y[j]+kemer]:
                kemer=kemer+kemer2
            if t[X[l]:X[l]+kemer]<t[Y[j]:Y[j]+kemer]:
                z.append(X[l])
                l=l+1
            else:
                z.append(Y[j])
                j=j+1
    z=z+X[l:]
    z=z+Y[j:]
    return z

if __name__ == "__main__":
    argv = sys.argv[1:]
    Length =100;
    Overlap=100;
    file='example1.fasta';
    try:
        opts, args = getopt.getopt(argv,"hL:G:",["Length=","Genome="])
    except getopt.GetoptError:
        print('Indexing.py -L <Substrings_Length> -G <Genome_file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Indexing.py -L <Substrings_Length> -G <Genome_file>')
            sys.exit()
        elif opt in ("-L", "--Length"):
            Length = int(arg)
        elif opt in ("-G", "--Genome"):
            file = arg
            
    Find_SA_for_Genome(L=Length,l=Overlap,G_file=file)
    """
    G_file: Reference Genome File default=example1.fasta
    L: Divide genome to substrings of length L default=100
    l : Overlap length =100
    """
     







    

