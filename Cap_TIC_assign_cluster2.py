import re
import pickle
import bisect
import numpy as np

## use the profile to judge whether the pre-miRNA is capped or not
## 2014.10.16

miRNA_cluster='/mnt/home/chenjiao/research/Project-miRNA/Data/miRNA/cel.cluster.1000.txt'
cluster_file='/mnt/home/chenjiao/research/Project-miRNA/Codes/python/cluster/RNA_short_cap_cluster05_whole.bed'

chrom_num={'I':0,'II':1,'III':2,'IV':3,'V':4,'X':5}

def find_TSS(miRNA,cluster_file,distance,cluster_edge):
    (mi_name,mi_chrom,mi_strand,mi_start,mi_end)=(miRNA[0],miRNA[1],miRNA[4],int(miRNA[2]),int(miRNA[3]))
    pri_TSS=[]
    pre_cap=[]
    with open(cluster_file,'r') as f:
        for c_line in f:
            c_line=c_line.strip()
            #I   35537   35631   Cluster1    14  +   35549   0.142857142857
            map=c_line.split('\t')
            (c_chrom,c_start,c_end,c_name,c_reads_num,c_strand,c_mode,c_depth,c_score)=(map[0],int(map[1]),int(map[2]),map[3],int(map[4]),map[5],int(map[6]),int(map[7]),float(map[8]))
            if c_chrom==mi_chrom and c_strand==mi_strand:
                if c_strand=='+':
                    if (c_mode>mi_start and c_mode<mi_end):
                        # mode reads_num depth score distance
                        pre_cap_cluster=[c_mode,c_reads_num,c_depth,c_score]
                        pre_mode_distance=abs(c_mode-mi_start)
                        pre_cap_cluster.append(pre_mode_distance)
                        pre_cap.append(pre_cap_cluster)
                    elif c_mode>cluster_edge[0]-distance and c_mode<=cluster_edge[0]:  # consider the cluster start site
                        TIC=[c_mode,c_reads_num,c_depth,c_score]
                        pre_mode_distance=abs(c_mode-mi_start)
                        TIC.append(pre_mode_distance)
                        pri_TSS.append(TIC)
                else: #minus strand
                    if c_mode>mi_start and c_mode<mi_end:
                        pre_cap_cluster=[c_mode,c_reads_num,c_depth,c_score]
                        pre_mode_distance=abs(c_mode-mi_end)
                        pre_cap_cluster.append(pre_mode_distance)
                        pre_cap.append(pre_cap_cluster)
                    elif c_mode>=cluster_edge[1] and c_mode<cluster_edge[1]+distance:
                        TIC=[c_mode,c_reads_num,c_depth,c_score]
                        pre_mode_distance=abs(c_mode-mi_end)
                        TIC.append(pre_mode_distance)
                        pri_TSS.append(TIC)
    #print len(pri_TSS)
    return mi_name,mi_chrom,mi_strand,mi_start,mi_end,pre_cap,pri_TSS

# find the data number(profile_array[i]) within the region [start,end]
def reads_count(profile_array,start,end,strand='+'):
    index_left=bisect.bisect_left(profile_array,start)
    index_right=bisect.bisect_right(profile_array,end)
    reads_num=sum(profile_array[index_left:index_right])

    if reads_num>0:
        cluster=profile_array[index_left:index_right]
        cluster_unique=np.unique(cluster)
        cluster_hist=[]
        for x in cluster_unique:
            cluster_hist.append(cluster.count(x))
        max_index=[i for i,j in enumerate(cluster_hist) if j==max(cluster_hist)]
        print max_index
        if strand=='+':
            return reads_num,cluster_unique[max_index[0]]
        else:
            return reads_num,cluster_unique[max_index[-1]]
    else:
        return reads_num,'.'

# detect whether the pre-miRNA is capped based on the reads profile
def cap_detect(miRNA,starts_plus,starts_minus,depth):
    ## output format: is_capped/ maximum depth position/ distance to the pre-start/ reads number within the pre-miRNA
    (mi_name,mi_chrom,mi_strand,mi_start,mi_end)=(miRNA[0],miRNA[1],miRNA[4],int(miRNA[2]),int(miRNA[3]))
    pre_cap=[]
    chr_num=chrom_num[mi_chrom]

    if mi_strand=='+':
        profile_plus=starts_plus[chr_num]
        reads_num,max_index=reads_count(profile_plus,mi_start,mi_end,mi_strand)
        if reads_num>=depth:
            pre_cap.append(1)
            pre_cap.append(max_index)
            dis_cap=abs(mi_start-max_index)
            pre_cap.append(dis_cap)
            pre_cap.append(reads_num)
        else:
            pre_cap=[0,'.','.',0]
    else:
        profile_minus=starts_minus[chr_num]
        profile_minus=np.add(profile_minus,35)
        reads_num,max_index=reads_count(profile_minus,mi_start,mi_end,mi_strand)
        if reads_num>=depth:
            pre_cap.append(1)
            pre_cap.append(max_index)
            dis_cap=abs(mi_end-max_index)
            pre_cap.append(dis_cap)
            pre_cap.append(reads_num)
        else:
            pre_cap=[0,'.','.',0]

    return mi_name,mi_chrom,mi_strand,mi_start,mi_end,pre_cap
    

def output_cap_TSS(pre_cap,pri_TSS):
## format the pre-cursor cap and primary miRNA TSS information
    pre_cap_out=[0,'','','','','']
    pri_TSS_out=[0,'','','','','']
    if len(pre_cap)==0:
        pre_cap_out=[0,'.','.','.','.','.']
    else:
        pre_cap=sorted(pre_cap,key=lambda cap:cap[0])
        pre_cap_out[0]=len(pre_cap)
        pre_cap_mode,pre_cap_distance='',''
        for pre_cap_cluster in pre_cap:
            temp=str(pre_cap_cluster[0])+','
            pre_cap_out[1]=pre_cap_out[1]+str(pre_cap_cluster[0])+',' ## TIC mode
            pre_cap_out[2]=pre_cap_out[2]+str(pre_cap_cluster[4])+',' ## TIC to 5' end distance
            pre_cap_out[3]=pre_cap_out[3]+str(pre_cap_cluster[1])+',' ## TIC reads number
            pre_cap_out[4]=pre_cap_out[4]+str(pre_cap_cluster[2])+',' ## TIC depth
            pre_cap_out[5]=pre_cap_out[5]+str(pre_cap_cluster[3])+',' ## TIC score


    if len(pri_TSS)==0:
        pri_TSS_out=[0,'.','.','.','.','.']
    else:
        pri_TSS=sorted(pri_TSS,key=lambda TSS:TSS[0])
        pri_TSS_out[0]=len(pri_TSS)
        for TIC in pri_TSS:
            # mode reads_num depth score distance
            pri_TSS_out[1]=pri_TSS_out[1]+str(TIC[0])+',' ## TIC mode
            pri_TSS_out[2]=pri_TSS_out[2]+str(TIC[4])+',' ## TIC to 5' end distance
            pri_TSS_out[3]=pri_TSS_out[3]+str(TIC[1])+',' ## TIC reads number
            pri_TSS_out[4]=pri_TSS_out[4]+str(TIC[2])+',' ## TIC depth
            pri_TSS_out[5]=pri_TSS_out[5]+str(TIC[3])+',' ## TIC score
    return pre_cap_out,pri_TSS_out

def parse_miRNA_cluster(miRNA_file):
    #0       cel-mir-35      MI0000006       7.91e+04        chrII   11537608        11537704        +
    miRNA_cluster={}
    with open(miRNA_file,'r') as f:
        for line in f:
            line=line.strip()
            map=line.split('\t')
            (mc_num,mc_ID,mc_chrom,mc_start,mc_end,mi_strand)=(int(map[0]),map[1],map[4],int(map[5]),int(map[6]),map[7])
            m=re.match(r'chr(\w+)',mc_chrom)
            if m:
                mc_chrom=m.group(1)
            else:
                print "Cluster chromosome error!"
            # ID chrom start end strand
            miRNA=[mc_ID,mc_chrom,mc_start,mc_end,mi_strand]
            if mc_num in miRNA_cluster:
                miRNA_cluster[mc_num].append(miRNA)
            else:
                miRNA_cluster[mc_num]=[]
                miRNA_cluster[mc_num].append(miRNA)

    miRNA_new_cluster={}
    cluster_edge={}
    for i in miRNA_cluster.keys():
        mi_cluster=miRNA_cluster[i]
        #print mi_cluster
        mi_new_cluster=[]
        miRNA_strand={'+':0,'-':0}
        for miRNA_gene in mi_cluster:
            #print str(i)+'\t'+miRNA_gene[0]+'\t'+miRNA_gene[1]+'\t'+str(miRNA_gene[2])+'\t'+str(miRNA_gene[3])+'\t'+miRNA_gene[4]+'\n'
            if miRNA_gene[4] in miRNA_strand:
                miRNA_strand[miRNA_gene[4]]=miRNA_strand[miRNA_gene[4]]+1
            else:
                print 'miRNA strand error!'

        ## judge whether there is miRNA from different strands, if so ,delete the minor genes from the same strand
        if miRNA_strand['+']>0 and miRNA_strand['-']>0:
            if miRNA_strand['+']>miRNA_strand['-'] and miRNA_strand['+']>1:
                for miRNA_gene in mi_cluster:
                    # ID chrom start end strand
                    if miRNA_gene[4]=='+':
                        print miRNA_strand['+']
                        mi_new_cluster.append(miRNA_gene)
            elif miRNA_strand['-']>miRNA_strand['+'] and miRNA_strand['-']>1:
                for miRNA_gene in mi_cluster:
                    if miRNA_gene[4]=='-':
                        print miRNA_strand['-']
                        mi_new_cluster.append(miRNA_gene)
        else:
            mi_new_cluster=mi_cluster
        if len(mi_new_cluster)>0:
            miRNA_new_cluster[i]=mi_new_cluster
            cluster_edge[i]=[mi_new_cluster[0][2],mi_new_cluster[-1][3]]
    #for j in miRNA_new_cluster.keys():
    #    mi_new_cluster=miRNA_new_cluster[j]
    #    print mi_new_cluster
    #    for miRNA_gene in mi_new_cluster:
    #        print str(j)+'\t'+miRNA_gene[0]+'\t'+miRNA_gene[1]+'\t'+str(miRNA_gene[2])+'\t'+str(miRNA_gene[3])+'\t'+miRNA_gene[4]+'\n'
    return miRNA_new_cluster,cluster_edge
    


####################################################################################################
plus_start_file='/mnt/home/chenjiao/research/Project-miRNA/Combine/miRNA/TIC/starts_plus.pkl'
minus_start_file='/mnt/home/chenjiao/research/Project-miRNA/Combine/miRNA/TIC/starts_minus.pkl'
f_plus_start=open(plus_start_file,'rb')
starts_plus=pickle.load(f_plus_start)
f_plus_start.close()
f_minus_start=open(minus_start_file,'rb')
starts_minus=pickle.load(f_minus_start)
f_minus_start.close()

miRNA_clusters,miRNA_cluster_edge=parse_miRNA_cluster(miRNA_cluster)
f_out=open('miRNA_cluster_TIC_assign_format02.txt','w')
f_out.write('#cluster_number'+'\t'+'mi_name'+'\t'+'mi_chrom'+'\t'+'mi_strand'+'\t'+'mi_start'+'\t'+'mi_end'+'\t'+'pre-miRNA_capped'+'\t'+'pre_cap_pisition'+'\t'+'pre_cap_to_5end_distance'+'\t'+'reads_num_within_pre-miRNA'+'\n')
for i in miRNA_clusters.keys():
    mi_cluster=miRNA_clusters[i]
    #print mi_cluster
    for miRNA_gene in mi_cluster:
        #print str(i)+'\t'+miRNA_gene[0]+'\t'+miRNA_gene[1]+'\t'+str(miRNA_gene[2])+'\t'+str(miRNA_gene[3])+'\t'+miRNA_gene[4]+'\n'
        mi_name,mi_chrom,mi_strand,mi_start,mi_end,pre_cap=cap_detect(miRNA_gene,starts_plus,starts_minus,10)
        if mi_strand=='+':
            f_out.write(str(i)+'\t'+mi_name+'\t'+mi_chrom+'\t'+mi_strand+'\t'+str(mi_start)+'\t'+str(mi_end)+'\t'+str(pre_cap[0])+'\t'+str(pre_cap[1])+'\t'+str(pre_cap[2])+'\t'+str(pre_cap[3])+'\n')
        elif mi_strand=='-':
            f_out.write(str(i)+'\t'+mi_name+'\t'+mi_chrom+'\t'+mi_strand+'\t'+str(mi_start)+'\t'+str(mi_end)+'\t'+str(pre_cap[0])+'\t'+str(pre_cap[1])+'\t'+str(pre_cap[2])+'\t'+str(pre_cap[3])+'\n')
        else:
            print "Output error!"
    #break;
f_out.close()

