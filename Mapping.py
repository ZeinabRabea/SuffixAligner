import numpy as np
import pandas as pd
import sys, getopt
import os
import random


def load_sam(sam_file):
# load sam file 
    header=[]
    read=[]
    with open(sam_file) as content:
        for line in content:
            line = line.strip()
            if line[0]!='@':
                read.append(line)
                for line in content:
                    line = line.strip()
                    if line[0]!='@':
                        read.append(line)
            else:
                header.append(line)
    return header,read

def load_Suffix_Array(Suffix_file):
#load Suffix_Array, output array of SA
    with open(Suffix_file) as content:
        s = content.readline()
    s=s[1:-1]
    x=s.split(", ")
    for i in range (0,len(x)):
      x[i]=int(x[i])
    return x
  
def load_genome(genome_file): 
#load genome all bases together return genome as string and ref_ADDRESS
    genome=""
    with open(genome_file) as content:
        line = content.readline().strip()
        genome_name=line[1:]
        genome_name=genome_name.rsplit()[0]
        for line in content:
            line = line.strip()
            if line[0]!='>':
                genome+= line          
    return genome,genome_name
  
def load_read(read_file):
#return array of read, read_name, read_qual
    read=[]
    read_name=[]
    read_qual=[]
    read1=""
    with open(read_file) as file:
        while True:
            n=file.readline().strip()
            if not n:
                break
            seq=file.readline().strip()
            file.readline()
            qual=file.readline().strip()
            read_name.append(n[1:].rsplit()[0])
            read_qual.append(qual)
            read.append(seq)
    return read,read_name,read_qual


def bwtViaSa(t,sss):
#input string(t) and SA return BWT(t)
    bw=[]
    for si in sss:
        if si==0: bw.append('$')
        else: bw.append(t[si-1])
    return ''.join(bw)

def rankBwt(bw):
# return the repeat of each DNA element and the rank of each element
    repeat=dict()
    ranks=[]
    for c in bw:
        if c not in repeat: repeat[c]=0
        ranks.append(repeat[c])
        repeat[c]+=1
    return ranks,repeat
  
def cc3(bw):
#return C[c]
    bw=sorted(bw)
    cc=dict()
    x=0
    for c in bw:
        if c not in cc: cc[c]=x
        x+=1
    return cc
  
def firstCol(tots):
#input repeat of element, output first colume of (BWM)
    first={}
    totc=0
    for c in sorted(tots):
        count=tots[c]
        first[c]=(totc,totc+count)
        totc+=count
    return first
  
def LF(bw,rank,C):
# return last to first array
    lf=[]
    for x in range(0,len(bw)):
        c2=bw[x]
        lf.append(rank[x]+C[c2])
    return lf

def reverse_read(x):
# A<>T   C<>G
  y=""
  for i in range(0,len(x)):
        if x[-1-i]=='A':
            y=y+'T'
        if x[-1-i]=='C':
            y=y+ 'G'
        if x[-1-i]=='G':
            y=y+ 'C'
        if x[-1-i]=='T':
            y=y+ 'A'
  return y
def remove_N(r):
    r=r.upper()
    r=r.replace('N',random.choice("ACTG"))
    return r
 
def search_read(i):
        g_read[i]=remove_N(g_read[i])
        x,y=Find_matched_seeds(g_bw,g_ranks,g_repeat,g_C,g_f,g_LF_Array,g_S_A,g_t,g_read[i])

        if y==True:
            y=x[0]
            flag=0
            find_in_pos=y[0]
            print(y[0])
            score=y[0]
            CIGAR_XX=str(len(g_read[i]))+"M"
            print_sam_element(g_read_name[i],flag,g_ref_title,find_in_pos,255,CIGAR_XX,g_read[i],"q")
            
        elif x==[]:
            search_read_reverse(i)
                       
        elif x[-1][0]<0:                   
            search_read_reverse(i)
            
        else:
            flag=0
            p_seed=det_window(g_t,g_read[i],x)
            score,find_in_pos,t_final,p_final=align_between_seeds(g_t,g_read[i],x,p_seed)
            CIGAR_XX=find_cigar(t_final,p_final)
            if score==0 and find_in_pos==0 and t_final=="4":
              search_read_reverse(i)
            else:
              print_sam_element(g_read_name[i],flag,g_ref_title,find_in_pos,score,CIGAR_XX,g_read[i],g_read_quality[i])
            
def search_read_reverse(i):
        g_read[i]=reverse_read(g_read[i])
        x,y=Find_matched_seeds(g_bw,g_ranks,g_repeat,g_C,g_f,g_LF_Array,g_S_A,g_t,g_read[i])

        if y==True:
            y=x[0]
            flag=16
            find_in_pos=y[0]
            print(y[0])
            score=y[0]
            CIGAR_XX=str(len(g_read[i]))+"M"
            print_sam_element(g_read_name[i],flag,g_ref_title,find_in_pos,255,CIGAR_XX,g_read[i],"q")
            
        elif x==[]:
            flag=4
            print_sam_element(g_read_name[i],flag,"*","0","0","*",g_read[i],g_read_quality[i])
            
        elif x[-1][0]<0:                   
            flag=4
            print_sam_element(g_read_name[i],flag,"*","0","0","*",g_read[i],g_read_quality[i])
                    
        else:
            flag=16
            p_seed=det_window(g_t,g_read[i],x)
            score,find_in_pos,t_final,p_final=align_between_seeds(g_t,g_read[i],x,p_seed)
            CIGAR_XX=find_cigar(t_final,p_final)
            if score==0 and find_in_pos==0 and t_final=="4":
              print_sam_element(g_read_name[i],4,"*","0","0","*",g_read[i],g_read_quality[i])
            else:
              print_sam_element(g_read_name[i],flag,g_ref_title,find_in_pos,score,CIGAR_XX,g_read[i],g_read_quality[i])
            

def Find_matched_seeds(bw,ranks,repeat,C,f,LF_Array,sss,t,p):
# return list of exact match(len > k1_mer) in read 
    k1_mer=10
    index_first=len(p) 
    p2=p
    x=[]
    y=0
    while index_first>k1_mer and index_first!=0:
        p_in_t,index_first,index_last=Modify_fm_index2 (bw,ranks,repeat,C,f,LF_Array,sss,t,p2)
        
        match=index_last-index_first
        p2=p[0:index_first]
        if match>k1_mer:
            for s in p_in_t:
                x.append([s-index_first,s,index_first,index_last,match])
                if len(p)==match:
                    y=True    #exact match 
    x.sort()
    return x,y

def Modify_fm_index(bw,ranks,repeat,C,f,LF_Array,t,p):
# return the position of p in LF_array  or return the start index of subset p that exact in t 
    i=len(p)-1
    c=p[i]
    A=f[c]
    tob=A[0]
    bot=A[1]
    index_true=-1
    Posation_p=[]
    while i>0 and bot>=tob:
        Posation_p=[]
        for m in range(tob,bot):
            try:
              if bw[m]==p[i-1]:
                Posation_p.append(LF_Array[m])
            except:
                nothingg=0 #bw[m] is the last
        if len(Posation_p)!=0:
            tob=Posation_p[0]
            bot=Posation_p[len(Posation_p)-1]+1
        else: 
            index_true=i
            i=0
        i=i-1
    return Posation_p,index_true
  
def Modify_fm_index2 (bw,ranks,repeat,C,f,LF_Array,sss,t,p):
    ss=[]
    index_first=0
    index_last=len(p)
    ll,index_true=Modify_fm_index(bw,ranks,repeat,C,f,LF_Array,t,p)
    if index_true==-1:
        for m in range(0,len(ll)):
          try:
            y=ll[m]
            ss.append(sss[y])
          except:  
            nothingg=0     
    else:
        index_first=index_true
        ll,index_true=Modify_fm_index(bw,ranks,repeat,C,f,LF_Array,t,p[index_true:index_last])
        for m in range(0,len(ll)):
          try:
            y=ll[m]
            ss.append(sss[y])
          except:
            nothingg=0
    return ss,index_first,index_last

def det_window(t,p,x):
  #select the expected mapping window 
    p_seed=[]
    if len(x)==1:
        y=x[0]
        p_seed.append([1,y[0],0,y[4]])
    else:
        k1_mer=100             
        col=list(zip(*x))
        match_N=col[4]
        v=col[0]
        
        v_before=v[0]  
        j=0
        count=0     
        total_N_of_match=0
        for i in range (0,len(v)):
            if v[i]>=-k1_mer:     
                if v_before+k1_mer<v[i] :
                    if i==len(v)-1:
                        total_N_of_match+=match_N[i]
                    if  v[j]>=-k1_mer and count>0:     
                       p_seed.append([count,v[j],j,total_N_of_match])
                    v_before=v[i]
                    j=i
                    count=1
                    total_N_of_match=match_N[i]
                elif  i==len(v)-1:
                    if i==len(v)-1:
                        total_N_of_match+=match_N[i]
                    if  v[j]>=-k1_mer and count>0:
                        count+=1
                        p_seed.append([count,v[j],j,total_N_of_match])
                    v_before=v[i]
                    j=i
                    count=1
                    total_N_of_match=match_N[i]
                    
                else:
                    count+=1
                    total_N_of_match+=match_N[i]
    p3_seed=[]
    for seed in p_seed:
        extend=pd.DataFrame(p_seed)
        extend=extend.sort_values(by=[3])
        extend=np.array(extend,'int')
    better_seed=extend[-1]
    p3_seed.append(better_seed)
    return p3_seed
  

def align_between_seeds(t,p,x,p_seed):
# Needle_Man is used to align the regions between the matched seeds 
    t_position=p_seed[0]
    if t_position[0]==1 :
          return 0,0,"4",""
    else:
        t_final=""
        p_final=""      
        extend=pd.DataFrame(x[t_position[2]:t_position[0]+t_position[2]])
        extend=extend.sort_values(by=[2])
        extend=np.array(extend,'int')
        i=extend[0]
        p3=p[0:i[2]]
        t3=t[i[1]-i[2]:i[1]]
        t3,p3=Needle_Man(t3,p3,1,-1,-1)
        t3,p3,countt=filter_remve_dash_first(t3,p3)    
        first=i[1]-i[2]+countt
        
        t_final=t_final+t3
        p_final=p_final+p3
        
        for j in range(0,len(extend)-1):
            i=extend[j]
            t_final=t_final+p[i[2]:i[3]]
            p_final=p_final+p[i[2]:i[3]]
            i2=extend[j+1]
            p3=p[i[3]:i2[2]]
            t3=t[i[1]+i[4]:i2[1]]
            t3,p3=Needle_Man(t3,p3,1,-1,-1)
            t_final=t_final+t3
            p_final=p_final+p3
        i=extend[len(extend)-1]
        t_final=t_final+p[i[2]:i[3]]
        p_final=p_final+p[i[2]:i[3]]
       
        p3=p[i[3]:len(p)]
        t3=t[i[1]+i[4]:i[1]+i[4]+len(p)-i[3]]  
        t3,p3=Needle_Man(t3,p3,1,-1,-1)
        t3,p3=filter_remve_dash_end(t3,p3)
        t_final=t_final+t3
        p_final=p_final+p3
        return Find_score(t_final,p_final),first,t_final,p_final   
    return

def filter_remve_dash_first(referen,reads) :
  #follow align_between_seeds
  mi=min(len(referen),len(reads))
  i=0
  count=0
  while i <mi and reads[i]=='-':
            i+=1
            count+=1   
  reads=reads[i:]
  referen=referen[i:]
  return referen,reads,count

def filter_remve_dash_end(referen,reads) :
  #follow align_between_seeds
  mi=min(len(referen),len(reads))         
  i=mi-1            
  while  i>=0 and i <mi and reads[i]=='-' :
            i-=1
  reads=reads[:i+1]
  referen=referen[:i+1]
  return referen,reads

def Needle_Man(Sub_text,Sub_pat,M,UM,G):
  #Needleman–Wunsch Algorithm 
    Sub_text="$"+Sub_text
    Sub_pat="$"+Sub_pat
    V = np.zeros((len(Sub_text), len(Sub_pat)), dtype=int)
    for i in range(0, len(Sub_text)):
        V[i,0]=UM*i
    for j in range(0, len(Sub_pat)):
        V[0,j]=UM*j
    for i in range(1, len(Sub_text)):
        for j in range(1, len(Sub_pat)):
            V[i, j] = max(V[i-1, j-1] + s_Needle_Man(Sub_text[i], Sub_pat[j],M,UM,G), 
                          V[i-1, j  ] + UM,    
                          V[i  , j-1] + UM)
    i, j = len(Sub_text)-1,len(Sub_pat)-1
    align1, align2 = '', ''
    while (i > 0 or j > 0):
        if i > 0 and j > 0 and V[i,j] == V[i-1,j-1] + s_Needle_Man(Sub_text[i], Sub_pat[j],M,UM,G) :
            a1=Sub_text[i]
            a2=Sub_pat[j]
            i-=1
            j-=1
        elif i > 0 and V[i,j] == V[i-1,j] + UM:
            a1=Sub_text[i]
            a2='-'
            i-=1
        else :
            a1='-'
            a2=Sub_pat[j]
            j-=1
        align1 = a1 +align1
        align2 = a2 +align2
    return align1,align2
def s_Needle_Man(t_c, p_c,M,UM,G):
  #retutn M, UM, G for Needleman–Wunsch Algorithm 
    if t_c == p_c: return M 
    if t_c == '-' or p_c == '-': return G
    return UM
  

def Mapping(Type="r",G_file="Acinetobacter_ref.fasta", R_file=["ERR776852.fastq"],
            SA_file="bacteria150overlap100.txt",Start=0,End=10,Sam_file=["map_minimap2_52.sam"],Output_file="x_output.sam"):
               
    global g_sam_file  
    global g_bw,g_ranks,g_repeat,g_C,g_f,g_LF_Array,g_S_A,g_t
    global g_read,g_read_name,g_ref_title,g_read_quality

    genome,g_ref_title = load_genome(G_file)
    genome=genome.replace('N','')                                                                  # chr of 64 less than A
    genome=genome.replace('n','')                                                                  # chr of 64 less than A
    g_t=genome+"$"
    print("Genome length",len(genome))
    
    g_S_A=load_Suffix_Array(SA_file)
    g_bw=bwtViaSa(g_t,g_S_A)
    g_ranks,g_repeat=rankBwt(g_bw)
    g_C=cc3(g_t)
    g_f=firstCol(g_repeat)
    g_LF_Array=LF(g_bw,g_ranks,g_C)

    
    g_sam_file = open(Output_file, 'w')
    #cl="command" # need modify
    print_sam_title(len(genome),cl)
    
    if Type=="r":
     for i in range(0,len(R_file)):
        g_read,g_read_name,g_read_quality=load_read(R_file[i])
        print("Number of read in file",R_file[i],"=  ",len(g_read))   
        if Start==-1:
            Start=0
        if End==-1:
            End=len(g_read)
        for i in range(Start,End,1)  :
          search_read(i)
          #print("**")
    if Type=="s":
      for i in range(0,len(Sam_file)):
        index=0
        g_read=["ACGT"]
        g_read_name=["n"]
        g_read_quality=["7h8b"]
        
        x,y=load_sam(Sam_file[i])
        print("Number of read in sam file",Sam_file[i],"=   ",len(y))
        if Start==-1:
            Start=0
        if End==-1:
            End=len(y)
        counter=0
        for i in range(Start,End,1):   
             d=y[i]
             r=d.rsplit()
             flag=r[1]
             if flag=="4" :
                counter+=1 
                g_read_name[index]=r[0]
                g_read[index]=r[9]
                g_read_quality[index]=r[10]
                search_read(index)
                #print("**")
        print("Number of unmapped read in sam file",counter)        
    g_sam_file.close()
def print_sam_title(n_of_ref_char,cl):
  t1="@HD"+"  "+"VN:1"+"  "+"SO:unsorted"+"\n"
  g_sam_file.write(t1)
  t2="@SQ"+"  "+"SN:"+g_ref_title+"  "+"LN:"+str(n_of_ref_char)+"\n"
  g_sam_file.write(t2)
  t3="@PG"+"  "+"ID:SuffixAligner"+"  "+"PN:SuffixAligner"+"  "+"VN:1"+"  "+"CL:"+cl+"\n"
  g_sam_file.write(t3)

  
def print_sam_element(r_n,flag_XX,ref_named,find_in_pos,score,CIGAR_XX,readd,quality):
  #MAPQ replace by score   
  sam_read=r_n+"  "+str(flag_XX)+"  "+ref_named+"  "+str(find_in_pos)+"  "+str(score)+"  "+CIGAR_XX+"  "+"*"+"  "+"0"+"  "+"0"+"  "+readd+"  "+quality+"\n"
  g_sam_file.write(sam_read)


def Find_score(t_final,p_final):
#find score for sam file Match =1, unmatch=-1, gap=-1 and return the mapq
    score=0
    mi=min(len(p_final),len(t_final))
    for i in range(0,mi):
        if t_final[i]==p_final[i]:
            score+=1
        #elif t_final[i]=='-' or p_final[i]=='-':
        #    score-=1
        else :
            score-=1
    if int(score)/len(p_final)< 0.15:
            return 60
    else:
            x=round( 60 +5* (0.2 - (-1*int(score)/len(p_final))) / 0.2)
            return min (x,255)
    return 0

def find_cigar(t_final,p_final):
#find cigar for sam file
    M_UM_c=0
    R_gap_delete_c=0
    G_gap_insert_c=0
    cigar=""
    mi=min(len(p_final),len(t_final))
    for i in range(0,mi):
        if t_final[i]!='-' and  p_final[i]!='-':
             M_UM_c+=1
             if R_gap_delete_c!=0:
               cigar=cigar+str(R_gap_delete_c)+"D"
               R_gap_delete_c=0
             if G_gap_insert_c!=0:
               cigar=cigar+str(G_gap_insert_c)+"I"
               G_gap_insert_c=0
        else :
          if M_UM_c!=0:
            cigar=cigar+str(M_UM_c)+"M"
            M_UM_c=0
          if t_final[i]=='-' :
            if R_gap_delete_c==0:
                G_gap_insert_c+=1
            else:
              cigar=cigar+str(R_gap_delete_c)+"D"
              R_gap_delete_c=0
              G_gap_insert_c+=1
          if p_final[i]=='-' :
            if G_gap_insert_c==0:
                R_gap_delete_c+=1
            else:
              cigar=cigar+str(G_gap_insert_c)+"I"
              G_gap_insert_c=0
              R_gap_delete_c+=1
    if R_gap_delete_c!=0:
               cigar=cigar+str(R_gap_delete_c)+"D"
    if M_UM_c!=0:
            cigar=cigar+str(M_UM_c)+"M"
    if G_gap_insert_c!=0:
               cigar=cigar+str(G_gap_insert_c)+"I"
    return cigar
  









 
if __name__ == "__main__":
    Default_Type="r"
    Default_G_file="example1.fasta"
    Default_R_file=[]
    Default_Start=-1
    Default_End=-1
    Default_sam_file=[] 
    Default_Output=""
    argv = sys.argv[1:]
    global cl
    cl=""
    
    
    
    try:
        cl="python Mapping.py"
        opts, args = getopt.getopt(argv,"hT:G:R:F:B:E:O:",["Type=","Genome=","Reads=","Sam_File=","Begin=","End=","Output="])
        for z in opts :
            cl+=str(z[0])+" "+str(z[1])+" "
    except getopt.GetoptError:
        print('Mapping.py -T r -G <Genome_file>  -R <Read_file>   -B <begining_read>  -E <Ending_read> -O <Output_file>')
        print("or")
        print('Mapping.py -T s -G <Genome_file>  -F <sam_file>  -B <begining_read>  -E <Ending_read> -O <Output_file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("for Mapping read in one file")
            print('Mapping.py -T r -G <Genome_file>  -R <Read_file>  -B <begining_read>  -E <Ending_read> -O <Output_file>')
            print("for Mapping read in more than one file")
            print('Mapping.py -T r -G <Genome_file>  -R <Read_file>  -R <Read_file2>  -R <Read_file3>  -B <begining_read>  -E <Ending_read> -O <Output_file>')
            print("for Mapping sam file")
            print('Mapping.py -T s -G <Genome_file>  -F <sam_file> -B <begining_read>  -E <Ending_read> -O <Output_file>')
            print("for Mapping more than one file")
            print('Mapping.py -T s -G <Genome_file>  -F <sam_file> -F <sam_file2> -B <begining_read>  -E <Ending_read> -O <Output_file>')
        
            sys.exit()
        elif opt in ("-T", "--Type"):
            Default_Type=arg
        elif opt in ("-G", "--Genome"):
            Default_G_file = arg
        elif opt in ("-R", "--Reads"):
            Default_R_file.append(arg)
        elif opt in ("-F", "--Sam_File"):
            Default_sam_file.append(arg)
        elif opt in ("-B", "--Begin"):
            Default_Start = int(arg)
        elif opt in ("-E", "--End"):
            Default_End = int(arg)
        elif opt in ("-O", "--Output"):
            Default_Output = arg
    
         
    if Default_Output=="":
        head_tail = os.path.split(Default_G_file)
        Default_Output=head_tail[1]+".sam"
        

    if Default_R_file==[]:
        Default_R_file=["Read_example1.fastq"]
        
    Default_SA_file=Default_G_file+".SA.txt"
    if cl=="":
        cl="Mapping(Type="+Default_Type+", G_file="+Default_G_file+",SA_file="+Default_SA_file+",Start="+Default_Start+",End="+Default_End+",Sam_file="+Default_sam_file+",Output_file="+Default_Output+")"
      
    Mapping(Type=Default_Type, G_file=Default_G_file,
            R_file=Default_R_file,
            SA_file=Default_SA_file,
            Start=Default_Start,End=Default_End,
            Sam_file=Default_sam_file,Output_file=Default_Output)

    

    
    
