import psutil 
import os
import numpy as np
import pandas as pd
import datetime
from psutil import virtual_memory
import time
import multiprocessing.dummy as mp
import multiprocessing as mp2
import threading
from numba import jit

def load_sam(sam_file):
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
            read_name.append(n[1:])
            read_qual.append(qual)
            read.append(seq)
    return read,read_name,read_qual


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

def bwtViaSa(t,sss):
#input string(t) and SA return BWT(t)
    bw=[]
    for si in sss:
        if si==0: bw.append('$')
        else: bw.append(t[si-1])
    return ''.join(bw)


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


def search_read(i):
        x,y=devide_pattern(g_bw,g_ranks,g_repeat,g_C,g_f,g_LF_Array,g_S_A,g_t,g_read[i])

        if y==True:
            y=x[0]
            flag=0
            find_in_pos=y[0]
            print(y[0])
            score=y[0]
            CIGAR_XX=str(len(g_read[i]))+"M"
            print_sam_element(g_read_name[i],flag,g_ref_title,find_in_pos,score,CIGAR_XX,g_read[i],"q")
            
        elif x==[]:
            search_read_reverse(i)
                       
        elif x[-1][0]<0:                   
            search_read_reverse(i)
            
        else:
            flag=0
            p_seed=det_window(g_t,g_read[i],x)
            score,find_in_pos,t_final,p_final=xxxM_plus(g_t,g_read[i],x,p_seed)
            CIGAR_XX=find_cigar(t_final,p_final)
            if score==0 and find_in_pos==0 and t_final=="4":
              search_read_reverse(i)
            else:
              print_sam_element(g_read_name[i],flag,g_ref_title,find_in_pos,score,CIGAR_XX,g_read[i],g_read_quality[i])
            
def search_read_reverse(i):
        this_read=g_read[i]
        g_read[i]=reverse_read(g_read[i])
        x,y=devide_pattern(g_bw,g_ranks,g_repeat,g_C,g_f,g_LF_Array,g_S_A,g_t,this_read)

        if y==True:
            y=x[0]
            flag=16
            find_in_pos=y[0]
            print(y[0])
            score=y[0]
            CIGAR_XX=str(len(g_read[i]))+"M"
            print_sam_element(g_read_name[i],flag,g_ref_title,find_in_pos,score,CIGAR_XX,this_read,"q")
            
        elif x==[]:
            flag=4
            print_sam_element(g_read_name[i],flag,"*","0","0","*",this_read,g_read_quality[i])
            
        elif x[-1][0]<0:                   
            flag=4
            print_sam_element(g_read_name[i],flag,"*","0","0","*",this_read,g_read_quality[i])
                    
        else:
            flag=16
            p_seed=det_window(g_t,g_read[i],x)
            score,find_in_pos,t_final,p_final=xxxM_plus(g_t,g_read[i],x,p_seed)
            CIGAR_XX=find_cigar(t_final,p_final)
            if score==0 and find_in_pos==0 and t_final=="4":
              print_sam_element(g_read_name[i],4,"*","0","0","*",this_read,g_read_quality[i])
            else:
              print_sam_element(g_read_name[i],flag,g_ref_title,find_in_pos,score,CIGAR_XX,this_read,g_read_quality[i])
            



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
            if bw[m]==p[i-1]:
                Posation_p.append(LF_Array[m])
        if len(Posation_p)!=0:
            tob=Posation_p[0]
            bot=Posation_p[len(Posation_p)-1]+1
        else: 
            index_true=i
            i=0
        i=i-1
    return Posation_p,index_true
  
def Used_Modify_fm_index (bw,ranks,repeat,C,f,LF_Array,sss,t,p):
    ss=[]
    index_first=0
    index_last=len(p)
    ll,index_true=Modify_fm_index(bw,ranks,repeat,C,f,LF_Array,t,p)
    if index_true==-1:
        for m in range(0,len(ll)):
            y=ll[m]
            ss.append(sss[y])      
    else:
        index_first=index_true
        ll,index_true=Modify_fm_index(bw,ranks,repeat,C,f,LF_Array,t,p[index_true:index_last])
        for m in range(0,len(ll)):
            y=ll[m]
            ss.append(sss[y])
    return ss,index_first,index_last


def devide_pattern(bw,ranks,repeat,C,f,LF_Array,sss,t,p):
# return list of exact match(len > k1_mer) in read 
    k1_mer=10
    index_first=len(p) 
    p2=p
    x=[]
    y=0
    while index_first>k1_mer and index_first!=0:
        p_in_t,index_first,index_last=Used_Modify_fm_index (bw,ranks,repeat,C,f,LF_Array,sss,t,p2)
        
        match=index_last-index_first
        p2=p[0:index_first]
        if match>k1_mer:
            for s in p_in_t:
                x.append([s-index_first,s,index_first,index_last,match])
                if len(p)==match:
                    y=True    #exact match 
    x.sort()
    return x,y

def det_window(t,p,x):
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
  



def main(NN="all_me_ERR776851_me.sam",SS=5000,EE=5500):
    
    global g_sam_file  
    global g_bw,g_ranks,g_repeat,g_C,g_f,g_LF_Array,g_S_A,g_t
    global g_read,g_read_name,g_ref_title,g_read_quality
    g_sam_file = open(NN, 'a')
    genome,g_ref_title = load_genome("Acinetobacter_ref.fasta")      #ALL
    
    genome=genome.replace('N','')                                                                  # chr of 64 less than A
    genome=genome.replace('n','')                                                                  # chr of 64 less than A
    
    print_sam_title(len(genome))

    g_t=genome+"$"
    g_S_A=load_Suffix_Array("bacteria150overlap100.txt")
    g_bw=bwtViaSa(g_t,g_S_A)

    g_read,g_read_name,g_read_quality=load_read("ERR776852.fastq")

     
    g_ranks,g_repeat=rankBwt(g_bw)
    g_C=cc3(g_t)
    g_f=firstCol(g_repeat)
    g_LF_Array=LF(g_bw,g_ranks,g_C)

    for i in range(SS,EE,1)  :
      search_read(i)
      print("**")
    g_sam_file.close()
    
def read_sam(NN="all_me_ERR776851_me.sam",SS=5000,EE=5500):
    
    global g_sam_file  
    global g_bw,g_ranks,g_repeat,g_C,g_f,g_LF_Array,g_S_A,g_t
    global g_read,g_read_name,g_ref_title,g_read_quality
    
    g_sam_file = open(NN, 'a')
    genome,g_ref_title = load_genome("Acinetobacter_ref.fasta")      #ALL
    
    genome=genome.replace('N','')                                                                  # chr of 64 less than A
    genome=genome.replace('n','')                                                                  # chr of 64 less than A
    
    print_sam_title(len(genome))

    g_t=genome+"$"
    g_S_A=load_Suffix_Array("bacteria150overlap100.txt")
    g_bw=bwtViaSa(g_t,g_S_A)

    g_read,g_read_name,g_read_quality=load_read("ERR776852.fastq")

     
    g_ranks,g_repeat=rankBwt(g_bw)
    g_C=cc3(g_t)
    g_f=firstCol(g_repeat)
    g_LF_Array=LF(g_bw,g_ranks,g_C)




    
    x,y=load_sam("all_me_ERR776851_meee.sam")
    count_flag_4=0
    not_found=0
    last_name="last"
    countother=0
    for i in range(SS,EE,1):    #0,len(y),1
      
     d=y[i]
     r=d.split("  ")
     flag=r[1]
     #print(r)
     if flag=="4" and r[0]!=last_name:        # must be chage to ==
      name_r=r[0]
      #print(name_r)
      name_r2=name_r.rsplit()
      name_r2=name_r2[0].split(".")
      #print(name_r2)
      n_r_i_maf=int(name_r2[1])-1
      #print(n_r_i_maf)
      count_flag_4+=1
      flag_me=0
      not_found+=len(r[9])
      
      g_read[n_r_i_maf]=reverse_read(r[9])
        
        #print("reverse ",r[0],g_read_name[n_r_i_maf])

      search_read(n_r_i_maf)
      #print(r[0])
      #print("not reverse ",r[0],g_read_name[n_r_i_maf])
      last_name=r[0]
      
    print(count_flag_4,not_found)
    print(err0r)
    
    
    
    
    """
    x,y=load_sam("map_bwa.sam")
    count_flag_4=0
    not_found=0
    for i in range(0,len(y),2):
     d=y[i]
     r=d.rsplit()
     flag=r[1]
     if flag=="4":

      
      name_r=r[0]
      name_r2=name_r.split("_")
      n_r_i_maf=int(name_r2[1])-1
      
      count_flag_4+=1
      flag_me=0
      #not_found+=len(r[9])
      r_read=reverse_read(r[9])
      score_me,find_in_position,flag_me=read_multiprocessor_read2(r_read)
      if flag_me=="4":
         print(r[0],"  ",find_in_position,"  ",score_me,"  ",len(r[9]),flag_me)
         not_found+=1
      else:
        print(r[0],"  ",find_in_position,"  ",score_me,"  ",len(r[9]))
      
    

        
    print(count_flag_4,not_found)
    print(err0r)
    """
   
    
    
    g_sam_file.close()
def print_sam_title(n_of_ref_char):
  
  t1="@HD"+"  "+"VN:try"+"  "+"SO:unsorted"+"\n"
  g_sam_file.write(t1)
  t2="@SQ"+"  "+"SN:"+g_ref_title+"  "+"LN:"+str(n_of_ref_char)+"\n"
  g_sam_file.write(t2)
  t3="@PG"+"  "+"ID:Zeinab Rabea"+"  "+"PN:Zeinab"+"  "+"VN:try"+"  "+"CL:A Fast Algorithm for Constructing Suffix Arrays for DNA Alphabets"+"\n"
  g_sam_file.write(t3)

  
def print_sam_element(r_n,flag_XX,ref_named,find_in_pos,score,CIGAR_XX,readd,quality):
  sam_read=r_n+"  "+str(flag_XX)+"  "+ref_named+"  "+str(find_in_pos)+"  "+str(score)+"  "+CIGAR_XX+"  "+"*"+"  "+"0"+"  "+"0"+"  "+readd+"  "+quality+"\n"
  #flag   0 map                        4 unmapped            16     reverse
  #CIGAR_XX
  g_sam_file.write(sam_read)



def score(t_final,p_final):
#Match =1, unmatch=-1, gap=-1
    score=0
    mi=min(len(p_final),len(t_final))
    for i in range(0,mi):
        if t_final[i]==p_final[i]:
            score+=1
        elif t_final[i]=='-' or p_final[i]=='-':
            score-=1
        else :
            score-=1
    return score





def find_cigar(t_final,p_final):
    score=0
    #M,UM,G    1   -1   -1
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
  








def xxxM_plus(t,p,x,p_seed):
# alignment between the exact match
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
        return score(t_final,p_final),first,t_final,p_final   
    return

def filter_remve_dash_first(referen,reads) :
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
  mi=min(len(referen),len(reads))         
  i=mi-1            
  
  while  i>=0 and i <mi and reads[i]=='-' :
            i-=1
  reads=reads[:i+1]
  referen=referen[:i+1]
  return referen,reads

  
def s(t_c, p_c,M,UM,G):
    if t_c == p_c: return M 
    if t_c == '-' or p_c == '-': return G
    return UM
def Needle_Man(Sub_text,Sub_pat,M,UM,G):
    Sub_text="$"+Sub_text
    Sub_pat="$"+Sub_pat
    V = np.zeros((len(Sub_text), len(Sub_pat)), dtype=int)
    for i in range(0, len(Sub_text)):
        V[i,0]=UM*i
    for j in range(0, len(Sub_pat)):
        V[0,j]=UM*j
    for i in range(1, len(Sub_text)):
        for j in range(1, len(Sub_pat)):
            V[i, j] = max(V[i-1, j-1] + s(Sub_text[i], Sub_pat[j],M,UM,G), 
                          V[i-1, j  ] + UM,    
                          V[i  , j-1] + UM)
    i, j = len(Sub_text)-1,len(Sub_pat)-1
    align1, align2 = '', ''
    while (i > 0 or j > 0):
        if i > 0 and j > 0 and V[i,j] == V[i-1,j-1] + s(Sub_text[i], Sub_pat[j],M,UM,G) :
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

 
if __name__ == "__main__":
    main(NN="test.sam",SS=0,EE=700)
    
