#include "main.h"
static pthread_mutex_t ReadLock,UpdateLock;

struct job_seed
{
    struct m_opt *opt;
    FILE *in_file_1,*in_file_2;
    FILE *out_file_1;
    FILE *unmap_file_1;
};
struct job_snp
{
    struct m_opt *opt;
    FILE *in_file[MAP_READ_BUF_LENGTH];
    struct hash_temp *temp;
    int *temp_fn;
};
uint16_t char2base_16(char *readseq,int length)
{
    uint16_t baseseq = 0;
    int i, base = 0;

    for (i = 0; i < length; i++)
    {
        base = nst_nt4_table[(int)readseq[i]];
        if(base == -1)
            return -1;
        baseseq = baseseq << 2 | base;
    }
    return baseseq;
}
int hash_c_find(unsigned int **hash_front,unsigned int *hash_num,uint16_t base,unsigned int site)
{
    int left = 0, right = hash_num[base]-1, middle;

    if (right == -1)
        return 0;
    while (left <= right)
    {
        middle = (left + right)/2;

        if (hash_front[base][middle] == site)
        {
            while ((middle-1>=0)&&(hash_front[base][middle-1] == site))
                middle--;
            return middle;
        }
        else if (hash_front[base][middle] > site)
            right = middle -1;
        else
            left = middle + 1;
    }
    return left;
}
void insert_c_hash(unsigned int **hash_front,unsigned int *hash_num,uint16_t base,FILE *hash_file)
{
    unsigned int *temp;
    temp = (unsigned int *)calloc(hash_num[base], sizeof(unsigned int));
    if(temp == NULL)
        fprintf(stdout,"There have no enough space...\n");

    int i;
    char f_line[MAX_STRING_LENGTH];
    for (i = 0;i<hash_num[base];i++)
    {
        if(fgets(f_line,MAX_STRING_LENGTH,hash_file)!=NULL)
            sscanf(f_line,"%u",temp+i);
    }

    hash_front[base] = temp;
}
void build_c_hash(struct m_opt *opt,int chr,FILE *hash_file)
{
    uint16_t i = 0;
    unsigned int num = 0;
    char f_line[MAX_STRING_LENGTH];
    while (fgets(f_line,MAX_STRING_LENGTH,hash_file)!=NULL)
    {
        if(f_line[0]=='>')
        {
            sscanf(f_line,">%hu %u",&i,&num);
            opt->chr->list[chr].c_num[i]=num;
        }
        insert_c_hash(opt->chr->list[chr].c_hash,opt->chr->list[chr].c_num,i,hash_file);
    }
}
void load_hash(struct m_opt *opt,int chr)
{
    char hash_file_name[MAX_NAME_LENGTH];
    FILE *Hash_file;

    strcpy(hash_file_name,opt->Hash_path);
    strcat(hash_file_name,"/\0");
    strcat(hash_file_name,opt->chr->list[chr].name);
    strcat(hash_file_name,".hash");
    Hash_file = fopen(hash_file_name,"r");

    opt->chr->list[chr].c_hash = (unsigned int **)calloc(pow(4,KMER_FRONT),sizeof(unsigned int *));
	opt->chr->list[chr].c_num= (unsigned int *)calloc(pow(4,KMER_FRONT),sizeof(unsigned int));

    build_c_hash(opt,chr,Hash_file);
    fclose(Hash_file);
}

int result_cmp(const void *a,const void *b)
{
    struct read_map_inf *EA,*EB;
    EA = (struct read_map_inf *)a;
    EB = (struct read_map_inf *)b;

    if (strcmp(EA->name,EB->name)==0)
    {
        return (EB->score-EA->score);
    }
    else
        return strcmp(EA->name,EB->name);
}
void HeapAdjust_Read(unsigned int *heap,unsigned int *heap_num,unsigned int *order,unsigned int heap_length,unsigned int i)
{
    int LC=2*i;
    int RC=2*i+1;
    int min=i; 
    if(i<=heap_length/2)
    {
        if((LC<=heap_length)
        &&((read_buf[heap[LC]].ref_site<read_buf[heap[min]].ref_site)
           ||((read_buf[heap[LC]].ref_site==read_buf[heap[min]].ref_site)&&(strcmp(read_buf[heap[LC]].name,read_buf[heap[min]].name)<0))))
            min=LC;
        if((RC<=heap_length)
        &&((read_buf[heap[RC]].ref_site<read_buf[heap[min]].ref_site)
           ||((read_buf[heap[RC]].ref_site==read_buf[heap[min]].ref_site)&&(strcmp(read_buf[heap[RC]].name,read_buf[heap[min]].name)<0))))
            min=RC;
        if(min!=i)
        {
            SwapHeap(i,min,heap,heap_num,order);
            HeapAdjust_Read(heap,heap_num,order,heap_length,min);
        }
    }
}
void ReHeap_Read(unsigned int *heap,unsigned int *heap_num,unsigned int *order,unsigned int *heap_length,FILE *InPut_file[])
{
    heap_num[1]--;
    if(heap_num[1]>0)
    {
        if(fread(&read_buf[order[1]],sizeof(struct read_map_inf),1,InPut_file[order[1]])!=1)
        {
            heap_num[1] = 0;
            SwapHeap(1,*heap_length,heap,heap_num,order);
            (*heap_length)--;
        }
    }
    else{
        SwapHeap(1,*heap_length,heap,heap_num,order);
        (*heap_length)--;
    }
    HeapAdjust_Read(heap,heap_num,order,*heap_length,1);
}
void HeapAdjust_Res(unsigned int *heap,unsigned int *heap_num,unsigned int *order,unsigned int heap_length,unsigned int i)
{
    int LC=2*i; 
    int RC=2*i+1;
    int min=i; 
    if(i<=heap_length/2)
    {
        if((LC<=heap_length)&&(result_cmp(read_buf[heap[LC]].name,read_buf[heap[min]].name)<0))
            min=LC;
        if((RC<=heap_length)&&(result_cmp(read_buf[heap[RC]].name,read_buf[heap[min]].name)<0))
            min=RC;
        if(min!=i)
        {
            SwapHeap(i,min,heap,heap_num,order);
            HeapAdjust_Res(heap,heap_num,order,heap_length,min);
        }
    }
}
void ReHeap_Res(unsigned int *heap,unsigned int *heap_num,unsigned int *order,unsigned int *heap_length,FILE *InPut_file[])
{
    heap_num[1]--;
    if(heap_num[1]>0)
    {
        if(fread(&read_buf[order[1]],sizeof(struct read_map_inf),1,InPut_file[order[1]])!=1)
        {
            heap_num[1] = 0;
            SwapHeap(1,*heap_length,heap,heap_num,order);
            (*heap_length)--;
        }
    }
    else{
        SwapHeap(1,*heap_length,heap,heap_num,order);
        (*heap_length)--;
    }
    HeapAdjust_Res(heap,heap_num,order,*heap_length,1);
}



#define MAX_EDGE_NUM 100
#define MAX_NODE_NUM 500000
struct node_table
{
    char base;
    unsigned int site;
    int in_edge[MAX_EDGE_NUM];
    uint8_t inum;
    int out_edge[MAX_EDGE_NUM];
    uint8_t onum;
    uint8_t change_flag;
    uint16_t rnum;
    float weigth;
    float score;
};
struct map_t
{
    struct node_table *node;
    int node_num;
};
struct stack_element
{
    int node;
    int edge;
};
struct tree_stack
{
    struct stack_element *data;
    int top;
};
void push_stack(struct m_opt *opt,struct tree_stack *t_stack,int node,int edge)
{
    t_stack->top++;
    t_stack->data[t_stack->top].node = node;
    t_stack->data[t_stack->top].edge = edge;
}
void pop_stack(struct tree_stack *t_stack)
{
    if(t_stack->top > -1)
        t_stack->top--;
}
int new_node(struct m_opt *opt,int front,int back,struct map_t *map,char base,unsigned int site)
{
    if(map->node_num<MAX_NODE_NUM)
    {
    map->node[map->node_num].base = base;
    map->node[map->node_num].site = site;
    if(front == -1)
        map->node[map->node_num].inum = 0;
    else
    {
        map->node[map->node_num].inum = 1;
        map->node[map->node_num].in_edge[0] = front;
        map->node[front].out_edge[map->node[front].onum] = map->node_num;
        map->node[front].onum++;
    }

    if(back == -1 )
        map->node[map->node_num].onum = 0;
    else
    {
        map->node[map->node_num].onum = 1;
        map->node[map->node_num].out_edge[0] = back;
        map->node[back].in_edge[map->node[back].inum] = map->node_num;
        map->node[back].inum++;
    }
    map->node[map->node_num].rnum = 1;
    map->node[map->node_num].weigth = 0;
    map->node[map->node_num].weigth += 1*1;//(qual);
    map->node[map->node_num].score = 0;

    map->node_num++;
    return 0;
    }
    else return 1;
}
int build_map(struct m_opt *opt,struct map_t *map,unsigned int start,unsigned int end,int chr,int branch_length,int r_root)
{
    int i = 0;
    unsigned int j = 0,k = 0;
    int base = 0;
    int front = 0,back = 0;
    unsigned int ref_pos = start-opt->chr->list[chr].start_site;

    int exon_order = find_exon(opt->exon,start);

    for(k = 0;k<(end-start+1);k++)
    {
        front = k-1;

        back = -1;
        new_node(opt,front,back,map,opt->chr->list[chr].seq[ref_pos+k],ref_pos+k);
        while((exon_order<opt->exon->total)&&(start+k>opt->exon->exon[exon_order].end)) exon_order++;
        if((exon_order!=-1)&&(opt->exon->exon[exon_order].start<=start+k)&&(opt->exon->exon[exon_order].end>=start+k))
            map->node[map->node_num-1].rnum = opt->exon->exon[exon_order].dep[start+k-opt->exon->exon[exon_order].start];
        else
            map->node[map->node_num-1].rnum = 0;
    }
    //insert SNP
    unsigned int LB;
    if(start>opt->area) LB= max(opt->chr->list[chr].start_site,start-opt->area);
    else LB = opt->chr->list[chr].start_site;

    int order = find_s(opt->snp,LB);
    int sorder = 0;
    unsigned int sstart = 0;
    unsigned int send = 0;
    int snode;

    for(i = order;i<opt->snp->total;i++)
    {
        if((opt->snp->snp[i].type=='I')&&(opt->snp->snp[i].start>=start)&&(opt->snp->snp[i].start<=end))
        {
            ref_pos = opt->snp->snp[i].start;
            for(j = 0;j<opt->snp->snp[i].length;j++)
            {
                base = opt->snp->snp[i].seq>>((opt->snp->snp[i].length-1-j)<<1)&3;

                if(j == 0)front = ref_pos-start;
                else front = map->node_num-1;

                if(j == opt->snp->snp[i].length-1)
                {
                    back = ref_pos-start+1;
                    if(back>r_root) back = -1;
                }
                else back = -1;

                if((front!=-1)&&(map->node[front].onum>=MAX_EDGE_NUM)) return 1;
                if((back!=-1)&&(map->node[back].inum>=MAX_EDGE_NUM)) return 1;

                new_node(opt,front,back,map,base2char[base],ref_pos-opt->chr->list[chr].start_site);
                map->node[map->node_num-1].rnum=opt->snp->snp[i].num;
            }
        }
        else if((opt->snp->snp[i].type=='X')&&(opt->snp->snp[i].start>=start)&&(opt->snp->snp[i].start<=end)
        &&(opt->snp->snp[i].end>=start)&&(opt->snp->snp[i].end<=end))
        {
            ref_pos = opt->snp->snp[i].start;
            for(j = 0;j<opt->snp->snp[i].length;j++)
            {
                base = opt->snp->snp[i].seq>>((opt->snp->snp[i].length-1-j)<<1)&3;

                if((j == 0)&&(ref_pos>=start)&&(ref_pos<=end))front = ref_pos-start-1;
                else
                {
                    if(j==0) front = -1;
                    else front = map->node_num-1;
                }

                if((j == opt->snp->snp[i].length-1)&&(ref_pos+1>=start)&&(ref_pos+1<=end))
                {
                    back = ref_pos+j-start+1;
                    if(back>r_root) back = -1;
                }
                else back = -1;

                if((front!=-1)&&(map->node[front].onum>=MAX_EDGE_NUM)) return 1;
                if((back!=-1)&&(map->node[back].inum>=MAX_EDGE_NUM)) return 1;

                new_node(opt,front,back,map,base2char[base],ref_pos+j-opt->chr->list[chr].start_site);
                map->node[map->node_num-1].rnum=opt->snp->snp[i].num;
            }
        }
        else if((opt->snp->snp[i].type=='D')&&(opt->snp->snp[i].start-1>=start)&&(opt->snp->snp[i].start-1<=end)
        &&(opt->snp->snp[i].end+1>=start)&&(opt->snp->snp[i].end+1<=end))
        {
            front = opt->snp->snp[i].start-1-start;
            back = opt->snp->snp[i].end+1-start;

            if((front!=-1)&&(map->node[front].onum>=MAX_EDGE_NUM)) return 1;
            if((back!=-1)&&(map->node[back].inum>=MAX_EDGE_NUM)) return 1;

            map->node[back].in_edge[map->node[back].inum] = front;
            map->node[back].inum++;

            map->node[front].out_edge[map->node[front].onum] = back;
            map->node[front].onum++;
        }
        else if((opt->snp->snp[i].type=='N')&&(((opt->snp->snp[i].start-1>=start)&&(opt->snp->snp[i].start-1<=end))
            ||((opt->snp->snp[i].end+1>=start)&&(opt->snp->snp[i].end+1<=end))))
        {
            if((opt->snp->snp[i].start-1>=start)&&(opt->snp->snp[i].start-1<=end)&&(opt->snp->snp[i].end+1>=start)&&(opt->snp->snp[i].end+1<=end))
            {
                front = opt->snp->snp[i].start-1-start;
                back = opt->snp->snp[i].end+1-start;

                if((front!=-1)&&(map->node[front].onum>=MAX_EDGE_NUM)) return 1;
                if((back!=-1)&&(map->node[back].inum>=MAX_EDGE_NUM)) return 1;

                map->node[back].in_edge[map->node[back].inum] = front;
                map->node[back].inum++;

                map->node[front].out_edge[map->node[front].onum] = back;
                map->node[front].onum++;
                continue;
            }
            if((opt->snp->snp[i].end+1>=start)&&(opt->snp->snp[i].end+1<=end))
            {
                if(opt->snp->snp[i].start-branch_length<opt->chr->list[chr].start_site) sstart = opt->chr->list[chr].start_site;
                else sstart = opt->snp->snp[i].start-branch_length;send = opt->snp->snp[i].start-1;ref_pos = opt->snp->snp[i].end+1;
            }
            else
                {sstart = opt->snp->snp[i].end+1;send = opt->snp->snp[i].end+branch_length;ref_pos = opt->snp->snp[i].start-1;}

            for(j = 0;j<branch_length;j++)
            {
                if((j == 0)&&(opt->snp->snp[i].start-1>=start)&&(opt->snp->snp[i].start-1<=end)) front = opt->snp->snp[i].start-1-start;
                else
                {
                    if(j==0) front = -1;
                    else front = map->node_num-1;
                }

                if((j == branch_length-1)&&(opt->snp->snp[i].end+1>=start)&&(opt->snp->snp[i].end+1<=end)) back = opt->snp->snp[i].end+1-start;
                else back = -1;

                if((front!=-1)&&(map->node[front].onum>=MAX_EDGE_NUM)) return 1;
                if((back!=-1)&&(map->node[back].inum>=MAX_EDGE_NUM)) return 1;

                new_node(opt,front,back,map,opt->chr->list[chr].seq[sstart+j-opt->chr->list[chr].start_site],sstart+j-opt->chr->list[chr].start_site);
                map->node[map->node_num-1].rnum=opt->snp->snp[i].num;
                if(j==0) snode = map->node_num-1;
            }
            sorder = find_s(opt->snp,sstart);
            while((sorder<opt->snp->total)&&(sorder!=-1))
            {
                if(opt->snp->snp[sorder].start>send) break;

                if((opt->snp->snp[sorder].start>=sstart)&&(opt->snp->snp[sorder].start<=send)&&(opt->snp->snp[sorder].end>=sstart)&&(opt->snp->snp[sorder].end<=send))
                {
                    if((opt->snp->snp[sorder].type=='I')&&(opt->snp->snp[sorder].start>=sstart)&&(opt->snp->snp[sorder].start<=send))
                    {
                        ref_pos = opt->snp->snp[sorder].start;
                        for(j = 0;j<opt->snp->snp[sorder].length;j++)
                        {
                            base = opt->snp->snp[sorder].seq>>((opt->snp->snp[sorder].length-1-j)<<1)&3;

                            if(j == 0)front = ref_pos-sstart+snode;
                            else front = map->node_num-1;

                            if(j == opt->snp->snp[sorder].length-1)
                            {
                                if(ref_pos+1<=send)back = ref_pos-sstart+1+snode;
                                else if((map->node[ref_pos-sstart+snode].onum>=1)
                                    &&(map->node[map->node[ref_pos-sstart+snode].out_edge[0]].site>=start)
                                    &&(map->node[map->node[ref_pos-sstart+snode].out_edge[0]].site<=end))
                                    back = map->node[ref_pos-sstart+snode].out_edge[0];
                                else back = -1;
                            }
                            else back = -1;

                            if((front!=-1)&&(map->node[front].onum>=MAX_EDGE_NUM)) return 1;
                            if((back!=-1)&&(map->node[back].inum>=MAX_EDGE_NUM)) return 1;

                            new_node(opt,front,back,map,base2char[base],ref_pos-opt->chr->list[chr].start_site);
                            map->node[map->node_num-1].rnum=opt->snp->snp[sorder].num;
                        }
                    }
                    else if((opt->snp->snp[sorder].type=='X')&&(opt->snp->snp[sorder].start>=sstart)&&(opt->snp->snp[sorder].start<=send)
                    &&(opt->snp->snp[sorder].end>=sstart)&&(opt->snp->snp[sorder].end<=send))
                    {
                        ref_pos = opt->snp->snp[sorder].start;
                        for(j = 0;j<opt->snp->snp[sorder].length;j++)
                        {
                            base = opt->snp->snp[sorder].seq>>((opt->snp->snp[sorder].length-1-j)<<1)&3;

                            if(j == 0)
                            {
                                if((ref_pos>sstart)&&(ref_pos<=send))front = ref_pos-sstart-1+snode;
                                else if((ref_pos==sstart)&&(map->node[snode].inum>0)) front = map->node[snode].in_edge[0];
                                else front = -1;
                            }
                            else front = map->node_num-1;


                            if(j == opt->snp->snp[sorder].length-1)
                            {
                                if((opt->snp->snp[sorder].end+1>=sstart)&&(opt->snp->snp[sorder].end+1<=send)) back = opt->snp->snp[sorder].end-sstart+1+snode;
                                else if((opt->snp->snp[sorder].end==send)&&(map->node[snode+branch_length-1].onum>0)) back = map->node[snode+branch_length-1].out_edge[0];
                                else back = -1;
                            }
                            else back = -1;

                            if((front!=-1)&&(map->node[front].onum>=MAX_EDGE_NUM)) return 1;
                            if((back!=-1)&&(map->node[back].inum>=MAX_EDGE_NUM)) return 1;

                            new_node(opt,front,back,map,base2char[base],ref_pos+j-opt->chr->list[chr].start_site);
                            map->node[map->node_num-1].rnum=opt->snp->snp[sorder].num;
                        }
                    }
                    else if((opt->snp->snp[sorder].type=='D')&&(opt->snp->snp[sorder].start-1>=sstart)&&(opt->snp->snp[sorder].start-1<=send)
                    &&(opt->snp->snp[sorder].end+1>=sstart)&&(opt->snp->snp[sorder].end+1<=send))
                    {
                        front = opt->snp->snp[sorder].start-1-sstart+snode;
                        back = opt->snp->snp[sorder].end+1-sstart+snode;

                        if((front!=-1)&&(map->node[front].onum>=MAX_EDGE_NUM)) return 1;
                        if((back!=-1)&&(map->node[back].inum>=MAX_EDGE_NUM)) return 1;

                        map->node[back].in_edge[map->node[back].inum] = front;
                        map->node[back].inum++;

                        map->node[front].out_edge[map->node[front].onum] = back;
                        map->node[front].onum++;
                    }
                }
                sorder++;
            }
        }
        else if(opt->snp->snp[i].start>end) break;
    }
    return 0;
}
void extend_seed_tail_forward(struct m_opt *opt,int chr_order,char *seq,int length,int mode,unsigned int start,unsigned int end,struct seed_t *seed,int *seed_num,unsigned int pos,int Sstart)
{
    int j;

    struct cigar_t t_cigar[21];
    int t_num = 0;

    struct seed_t seed_r;

    uint64_t ref_pos;
    int read_pos;

    int init_num = (*seed_num);

    read_pos = 0;
    ref_pos = pos;

    if(Sstart>=5)
        tail_align_seed(opt,chr_order,0,pos+opt->chr->list[chr_order].start_site,seq,Sstart,t_cigar,&t_num,20,&ref_pos,0);
    if(t_num!=0)
    {
        ref_pos-=opt->chr->list[chr_order].start_site;
        for(j = 0;j<t_num;j++)
        {
            if((t_cigar[j].c=='M')&&(t_cigar[j].l>=5))
            {
                seed_r.pos = ref_pos;
                seed_r.start = read_pos;
                seed_r.length = t_cigar[j].l;
                seed_r.abs = seed_r.pos-seed_r.start;
                seed_r.lnum = 0;
                seed_r.score = 0;
                if((mode!=1)&&(seed_r.pos<start))
                {
                    seed_r.start += start-seed_r.pos+1;
                    seed_r.length -= start-seed_r.pos+1;
                    seed_r.pos += start-seed_r.pos+1;
                }
                if(seed_r.length>0)
                insert_seed(seed,seed_num,SEED_BUF_LENGTH,seed_r);

                if((*seed_num)>init_num)
                    extend_seed_tail_forward(opt,chr_order,seq,seed_r.start,mode,start,seed_r.pos,seed,seed_num,seed_r.pos,seed_r.start);
            }
            if((t_cigar[j].c=='M')||(t_cigar[j].c=='X')||(t_cigar[j].c=='D')||(t_cigar[j].c=='N')) ref_pos+=t_cigar[j].l;
            if((t_cigar[j].c=='M')||(t_cigar[j].c=='X')||(t_cigar[j].c=='I')||(t_cigar[j].c=='S')) read_pos+=t_cigar[j].l;
        }
    }
}
void extend_seed_tail_backward(struct m_opt *opt,int chr_order,char *seq,int length,int mode,unsigned int start,unsigned int end,struct seed_t *seed,int *seed_num,unsigned int pos,int Sstart)
{
    int j;
    //int length = strlen(seq);

    struct cigar_t t_cigar[21];
    int t_num = 0;

    struct seed_t seed_r;

    uint64_t ref_pos,post;
    int read_pos;

    int init_num = (*seed_num);

    ref_pos =pos;
    read_pos =Sstart;

    if(length-read_pos>=5)
        tail_align_seed(opt,chr_order,0,ref_pos-1+opt->chr->list[chr_order].start_site,&seq[read_pos],length-read_pos,t_cigar,&t_num,20,&post,1);
    if(t_num!=0)
    {
        for(j = 0;j<t_num;j++)
        {
            if((t_cigar[j].c=='M')&&(t_cigar[j].l>=5))
            {
                seed_r.pos = ref_pos;
                seed_r.start = read_pos;
                seed_r.length = t_cigar[j].l;
                seed_r.abs = seed_r.pos-seed_r.start;
                seed_r.lnum = 0;
                seed_r.score = 0;
                if((mode!=0)&&(seed_r.pos+seed_r.length>end))
                    seed_r.length-=end-seed_r.pos-seed_r.length;
                if(seed_r.length>0)
                insert_seed(seed,seed_num,SEED_BUF_LENGTH,seed_r);

                if((*seed_num)>init_num)
                    extend_seed_tail_backward(opt,chr_order,seq,length,mode,seed_r.pos+seed_r.length-1,end,seed,seed_num,seed_r.pos+seed_r.length,seed_r.start+seed_r.length);
            }
            if((t_cigar[j].c=='M')||(t_cigar[j].c=='X')||(t_cigar[j].c=='D')||(t_cigar[j].c=='N')) ref_pos+=t_cigar[j].l;
            if((t_cigar[j].c=='M')||(t_cigar[j].c=='X')||(t_cigar[j].c=='I')||(t_cigar[j].c=='S')) read_pos+=t_cigar[j].l;
        }
    }
}
void extend_seed_tail(struct m_opt *opt,char *chr,int chr_order,char *seq,int mode,unsigned int start,unsigned int end,struct seed_t *seed,int *seed_num)
{
    int i;
    int init_num = (*seed_num);
    int length = strlen(seq);

    if(mode!=0)
        extend_seed_tail_forward(opt,chr_order,seq,length,mode,start,end,seed,seed_num,end,length);
    if(mode!=1)
        extend_seed_tail_backward(opt,chr_order,seq,length,mode,start,end,seed,seed_num,start+1,0);

    for(i = 0;i<init_num;i++)
    {
        extend_seed_tail_forward(opt,chr_order,seq,length,mode,start,seed[i].pos,seed,seed_num,seed[i].pos,seed[i].start);
        extend_seed_tail_backward(opt,chr_order,seq,length,mode,seed[i].pos+seed[i].length,end,seed,seed_num,seed[i].pos+seed[i].length,seed[i].start+seed[i].length);
    }
}
void cigar2site(unsigned int pos,unsigned int *site,struct cigar_t *cigar,int cigar_num)
{
    int read_site = 0;
    unsigned int ref_site = pos;
    int i,j;

    for(i = 0;i<cigar_num;i++)
    {
        switch(cigar[i].c)
        {
        case 'M':
            for(j = 0;j<cigar[i].l;j++)
                site[read_site+j] = ref_site+j;

            ref_site+=cigar[i].l;
            read_site+=cigar[i].l;
            break;
        case 'X':
            for(j = 0;j<cigar[i].l;j++)
                site[read_site+j] = ref_site+j;

            ref_site+=cigar[i].l;
            read_site+=cigar[i].l;
            break;
        case 'I':
            for(j = 0;j<cigar[i].l;j++)
                site[read_site+j] = ref_site-1;

            read_site+=cigar[i].l;
            break;
        case 'D':
            ref_site+=cigar[i].l;
            break;
        case 'N':
            ref_site+=cigar[i].l;
            break;
        case 'U':
            ref_site+=cigar[i].l;
            break;
        case 'S':
            read_site+=cigar[i].l;
            break;
        }
    }
}
#define AREA_SEED_CAND_NUM 30

int find6(struct m_opt *opt,char *chr,char *seq,int mode,unsigned int start,unsigned int end,
               unsigned int **hash,unsigned int *hash_num,struct seed_t *seed,unsigned int *site,int align_flag,int chr_order)
{
    int i,j,k;
    int length = strlen(seq);
    char ref[MAX_READ_LENGTH+20];
    char text[MAX_READ_LENGTH];
    int length_r,length_t;
    char cigarBuf[2*MAX_READ_LENGTH];

    char f_seq[KMER_FRONT+1];
    uint16_t f_base = 0;
    uint32_t fmax = pow(4,KMER_FRONT);
    int hash_site = 0;

    struct seed_t seed_r;
    int seed_num = 0;

    int temp_start = 0;
    int temp_end = 0;
    int temp_mode = 0;
    unsigned int start_site = 0;
    unsigned int end_site = 0;
    int flag = 0;

    int step = 0;

    if(align_flag==1)
    {
        struct cigar_t t_cigar[21];
        int t_num;

        if(mode==0)
        {
            length_r = min(end-start-1,length+20);
            strncpy(ref,&chr[start+1],length_r);
            ref[length_r] = '\0';
            cigarBuf[0] = '\0';
            memset(cigarBuf,0,length_r);
            if(computeEditDistanceWithCigar(ref,length_r,seq,length,MAX_K-1,cigarBuf,MAX_READ_LENGTH,0,1)!=-1)
                landau2site(opt,length_r,length,start+1,cigarBuf,site,0,0);
        }
        else if(mode==1)
        {
            length_r = min(end-start-1,length+20);
            strncpy(ref,&chr[end-length_r],length_r);
            ref[length_r] = '\0';

            char r_ref[4*MAX_READ_LENGTH];
            char r_text[4*MAX_READ_LENGTH];
            for(k = 0; k<length_r; k++)
                r_ref[k] = toupper(ref[length_r-1-k]);
            r_ref[length_r] = '\0';
            for(k = 0; k<length; k++)
                r_text[k] = toupper(seq[length-1-k]);
            r_text[length] = '\0';

            cigarBuf[0]='\0';
            memset(cigarBuf,0,length_r);
            if(computeEditDistanceWithCigar(r_ref,length_r,r_text,length,MAX_K-1,cigarBuf,MAX_READ_LENGTH,0,1)!=-1)
                landau2site(opt,length_r,length,end-1,cigarBuf,site,1,0);
        }
        else if(mode==2)
        {
            unsigned int chr_start = opt->chr->list[chr_order].start_site;
            t_num = 0;
            middle_align_seed(opt,chr_order,0,start+chr_start,end+chr_start,seq,length,t_cigar,&t_num,20);
            if(t_num!=0)
                cigar2site(start+1,site,t_cigar,t_num);

            else return 1;
        }
        else return 1;
    }
    else
    {
        step = KMER_FRONT/2;
        i = 0;
        while ((i+KMER_FRONT)<=length)
        {
            strncpy(f_seq,&seq[i],KMER_FRONT);
            f_seq[KMER_FRONT] = '\0';
            f_base = char2base_16(f_seq,KMER_FRONT);

            if(f_base>=fmax) return -1;

            hash_site = hash_c_find(hash,hash_num,f_base,start);
            if((hash_site+AREA_SEED_CAND_NUM<hash_num[f_base])&&(hash[f_base][hash_site+AREA_SEED_CAND_NUM]>=start)&&(hash[f_base][hash_site+AREA_SEED_CAND_NUM]<=end)) {i+=step;continue;}
            while ((hash_site < hash_num[f_base])&&(hash[f_base][hash_site]<=end))
            {
                if((hash[f_base][hash_site]>=start)&&(hash[f_base][hash_site]<=end))
                {
                    seed_r.pos = hash[f_base][hash_site];
                    seed_r.start = i;
                    seed_r.length = KMER_FRONT;
                    seed_r.abs = hash[f_base][hash_site]-i;
                    seed_r.lnum = 0;
                    seed_r.score = 0;
                    insert_seed(seed,&seed_num,SEED_BUF_LENGTH,seed_r);
                }
                hash_site++;
            }
            i+=step;
        }
        if(seed_num != 0)
        {
        int s = 0;
        int l = 0;

        for(i = 0;i<seed_num;i++)
        {
            k = 0;
            s = seed[i].start;
            for(j = 1;j<=min(step,s);j++)
            {
                if(nst_nt4_table[(int)chr[seed[i].pos-j]]!=nst_nt4_table[(int)seq[s-j]])
                {
                    k = j-1;
                    break;
                }
                k = j;
            }
            seed[i].pos-=k;
            seed[i].start-=k;
            seed[i].length+=k;

            s = seed[i].start;
            l = seed[i].length;
            for(j = 1;j<=min(step-1,length - l - s);j++)
            {
                if(nst_nt4_table[(int)chr[seed[i].pos+j+l-1]]!=nst_nt4_table[(int)seq[s+l+j-1]])
                {
                    k = j-1;
                    break;
                }
                k = j;
            }
            seed[i].length+=k;
        }
        }
        else
        {
        step = 1;
        i = 0;
        while ((i+KMER_FRONT)<=length)
        {
            strncpy(f_seq,&seq[i],KMER_FRONT);
            f_seq[KMER_FRONT] = '\0';
            f_base = char2base_16(f_seq,KMER_FRONT);

            if(f_base>=fmax) return -1;

            hash_site = 0;
            hash_site = hash_c_find(hash,hash_num,f_base,start);
            if((hash_site+AREA_SEED_CAND_NUM<hash_num[f_base])&&(hash[f_base][hash_site+AREA_SEED_CAND_NUM]>=start)&&(hash[f_base][hash_site+AREA_SEED_CAND_NUM]<=end)) {i+=step;continue;}
            while ((hash_site < hash_num[f_base])&&(hash[f_base][hash_site]<=end))
            {
                if((hash[f_base][hash_site]>=start)&&(hash[f_base][hash_site]<=end))
                {
                    seed_r.pos = hash[f_base][hash_site];
                    seed_r.start = i;
                    seed_r.length = KMER_FRONT;
                    seed_r.abs = hash[f_base][hash_site]-i;
                    seed_r.lnum = 0;
                    seed_r.score = 0;
                    insert_seed(seed,&seed_num,SEED_BUF_LENGTH,seed_r);
                }
                hash_site++;
            }
            i+=step;;
        }
        }

        extend_seed_tail(opt,chr,chr_order,seq,mode,start,end,seed,&seed_num);

        int max_o = -1;
        int max_s = -1000;

        if(mode!=1) qsort(seed,seed_num,sizeof(struct seed_t),seed_cmp);
        else qsort(seed,seed_num,sizeof(struct seed_t),seed_cmp_r);
        for(i = 0; i<seed_num; i++)
        {
            max_o = -1;
            max_s = -1000;
            seed[i].score = seed[i].length;

            if(mode!=1)
            {
                if((seed[i].abs>=start+1)&&(seed[i].abs-start>opt->change_length)) seed[i].score-=opt->splice;

                else if ((seed[i].abs>=start+1)&&(seed[i].abs-start+1<=opt->change_length)) seed[i].score-=opt->gap;
                else if ((seed[i].abs<start+1)) seed[i].score-=(start+1-seed[i].abs)*opt->match+opt->gap;
            }

            for(j = 0;j<i;j++)
            {
                if(seed[j].length<opt->min_exon)
                {
                    if(abs(seed[i].abs-seed[j].abs)>opt->change_length)
                    {
                        if(((seed[j].lnum==0)&&(seed[i].abs-start>opt->change_length))||((seed[j].lnum!=0)&&(abs(seed[j].abs-seed[seed[j].last[0]].abs)>opt->change_length)))continue;
                    }
                }

                if((seed[i].start > seed[j].start)&&((seed[i].start+seed[i].length)>(seed[j].start+seed[j].length)&&(seed[i].pos > seed[j].pos)))
                {
                    if(seed[i].start<=seed[j].start+seed[j].length)
                        max_s = seed[j].score + (seed[i].start+seed[i].length-seed[j].start-seed[j].length);
                    else
                        max_s = seed[j].score + seed[i].length;

                    if(seed[i].abs<seed[j].abs)
                        max_s -=(seed[j].abs-seed[i].abs)+opt->gap;
                    else if((seed[i].abs>seed[j].abs)&&(seed[i].abs-seed[j].abs<opt->change_length))
                        max_s -=opt->gap;
                    else if((seed[i].abs>seed[j].abs)&&(seed[i].abs-seed[j].abs>opt->change_length))
                        max_s -=opt->splice;

                    if(seed[i].score < max_s)
                    {
                        max_o = j;
                        seed[i].score = max_s;
                    }
                }
            }
            if(max_o == -1)
            {
                {
                    seed[i].score = seed[i].length;
                    if((mode!=1)&&(seed[i].abs-start>opt->change_length)) seed[i].score-=opt->splice;
                }
            }
            else
            {
                seed[i].last[0] = max_o;
                seed[i].lnum = 1;
            }
        }

        max_o = -1;
        max_s = 0;
        int best_num = 0;
        for(i = seed_num-1; i>=0; i--)
        {


            if((mode!=0)&&(end-length>seed[i].abs)&&(end-length-seed[i].abs>opt->change_length)) seed[i].score-=opt->splice;

            else if((mode!=0)&&(end-length>seed[i].abs)&&(end-length-seed[i].abs<=opt->change_length)) seed[i].score-=opt->gap;
            else if((mode!=0)&&(end-length<seed[i].abs)) seed[i].score-=(seed[i].abs-(end-length))+opt->gap;
            if(seed[i].score>=max_s)
            {
                max_o = i;
                max_s = seed[i].score;
                best_num = 1;
            }
        }
        if(best_num>1) max_o = -1;
        while(max_o!=-1)
        {
            for(i = 0;i<seed[max_o].length;i++)
                site[seed[max_o].start+i] = seed[max_o].pos + i;
            if(seed[max_o].lnum>0)
                max_o = seed[max_o].last[0];
            else max_o = -1;
        }
        for(j = 0;j<length;j++)
        {
            if(site[j]>=end) site[j] = end-1;
            else if((site[j]!=0)&&(site[j]<start)) site[j] = start;

            if((j!=0)&&(site[j]!=0)&&(site[j]<site[j-1]))
                site[j]=site[j-1];

            if ((flag==0)&&(site[j]==0))
            {
                temp_start = j;
                flag = 1;
            }
            if ((flag==1)&&((site[j]!=0)||(j==(length-1))))
            {
                if(site[j]==0)
                    temp_end = j;
                else
                    temp_end = j-1;

                temp_mode = mode;
                if(temp_start == 0)
                    start_site = start;
                else
                {
                    start_site = site[temp_start-1];
                    if(temp_mode == 1)
                        temp_mode = 2;
                }

                if(j == length-1)
                    end_site = end;
                else
                {
                    if(site[temp_end+1]>end)
                        end_site = end;
                    else
                        end_site = site[temp_end+1];
                    if(temp_mode == 0)
                        temp_mode = 2;
                }

                length_t = temp_end-temp_start+1;
                strncpy(text,&seq[temp_start],length_t);
                text[length_t] = '\0';
                for(k = 0; k<length_t; k++)
                    text[k] = toupper(text[k]);

                if(end_site<=start_site+1)
                {
                    if((mode==1)&&(temp_end==(length-1)))
                    {
                        if(end_site<start_site+1)
                        {
                            i = 1;
                            while(((temp_start-i)>=0)&&(site[temp_start-i] >= end_site))
                            {
                                site[temp_start-i] = end_site-1;
                                i++;
                            }
                        }
                        for(i = temp_start;i<=temp_end;i++)
                            site[i] = end_site-1;
                    }
                    else
                    {
                        if(end_site<start_site+1)
                        {
                            i = 1;
                            while(((temp_end+i)<length)&&(site[temp_end+i] < start_site))
                            {
                                site[temp_end+i] = start_site;
                                i++;
                            }
                        }
                        for(i = temp_start;i<=temp_end;i++)
                            site[i] = start_site;
                    }

                    flag = 0;
                    continue;
                }
                if (find6(opt,chr,text,temp_mode,start_site,end_site,hash,hash_num,seed,&(site[temp_start]),1,chr_order)==1) return 1;
                flag = 0;
            }
        }
    }
    return 0;
};
int area_align(struct m_opt *opt,char *chr,char *seq,int mode,unsigned int start,unsigned int end,
               unsigned int **hash,unsigned int *hash_num,struct seed_t *seed,unsigned int *site,int chr_order)//mode 0 front 1 back 2 both
{
    if(end<=start) return 0;

    int length = strlen(seq);
    int i = 0,k = 0;

    int flag = 0;

    unsigned int lsite[MAX_READ_LENGTH];
    int score = 0,lscore = -4000;
    char ref[MAX_READ_LENGTH+20];
    int length_r ;
    char cigarBuf[2*MAX_READ_LENGTH];

    if(mode==0)
    {
            memset(lsite, 0, length*sizeof(unsigned int)/sizeof(char));
            length_r = min(end-start-1,length+20);
            strncpy(ref,&chr[start+1],length_r);
            ref[length_r] = '\0';
            cigarBuf[0] = '\0';
            memset(cigarBuf,0,length_r);
            if(computeEditDistanceWithCigar(ref,length_r,seq,length,MAX_K-1,cigarBuf,MAX_READ_LENGTH,0,1)!=-1)
                lscore = landau2site(opt,length_r,length,start+1,cigarBuf,lsite,0,0);

            if(lsite[0]-start>opt->change_length)
                lscore -= opt->splice;
            else if ((lsite[0]-start<opt->change_length)&&(lsite[0]-start>1))
                lscore -= opt->gap;
    }
    else if(mode==1)
    {
            lscore = 0;
            memset(lsite, 0, length*sizeof(unsigned int)/sizeof(char));
            length_r = min(end-start-1,length+20);
            strncpy(ref,&chr[end-length_r],length_r);
            ref[length_r] = '\0';

            char r_ref[4*MAX_READ_LENGTH];
            char r_text[4*MAX_READ_LENGTH];
            for(k = 0; k<length_r; k++)
                r_ref[k] = toupper(ref[length_r-1-k]);
            r_ref[length_r] = '\0';
            for(k = 0; k<length; k++)
                r_text[k] = toupper(seq[length-1-k]);
            r_text[length] = '\0';

            cigarBuf[0]='\0';
            memset(cigarBuf,0,length_r);
            if(computeEditDistanceWithCigar(r_ref,length_r,r_text,length,MAX_K-1,cigarBuf,MAX_READ_LENGTH,0,1)!=-1)
                lscore = landau2site(opt,length_r,length,end-1,cigarBuf,lsite,1,0);

            if(end-lsite[length-1]>opt->change_length)
                lscore -= opt->splice;
            else if ((end-lsite[length-1]<opt->change_length)&&(end-lsite[length-1]>1))
                lscore -= opt->gap;
    }

    if(find6(opt,chr,seq,mode,start,end,hash,hash_num,seed,site,0,chr_order)==1) return -1000;

    int out_flag = 0;
    score = 0;
    flag = 1;
    if(mode!=1)
    {
        if((site[0]-start>1)&&(site[0]-start<opt->change_length))score-=opt->gap;
        else if (site[0]-start>opt->change_length) score-=opt->splice;

        if(start==site[0]) {score-=opt->miss;flag = 0;}
        else
        {
            if(nst_nt4_table[(int)chr[site[0]]]==nst_nt4_table[(int)seq[0]]) score+=opt->match;
            else score-=opt->miss;
        }
    }
    else
    {
        if(nst_nt4_table[(int)chr[site[0]]]==nst_nt4_table[(int)seq[0]]) score+=opt->match;
        else score-=opt->miss;
    }
    if (site[0]==0) {score=0;out_flag = 1;}
    for(i = 1;i<length;i++)
    {
        if (site[i]==0) {score=0;out_flag = 1;break;}

        if(site[i]-site[i-1]==0)
        {
            if(flag==1) {score-=opt->miss;flag = 0;}
        }
        else
        {
            if(nst_nt4_table[(int)chr[site[i]]]==nst_nt4_table[(int)seq[i]]) score+=opt->match;
            else score-=opt->miss;
            flag = 1;
        }

        if((site[i]-site[i-1]>1)&&(site[i]-site[i-1]<=opt->change_length))
        {
            flag = 1;
            score-=opt->gap;
        }
        else if(site[i]-site[i-1]>opt->change_length)
        {
            flag = 1;
            score-=opt->splice;
        }
    }
    if((mode!=0)&&(out_flag==0))
    {
        if((end-site[length-1]>1)&&(end-site[length-1]<opt->change_length))score-=opt->gap;
        else if (end-site[length-1]>opt->change_length) score-=opt->splice;
        else if((end==site[length-1])&&(flag)) score-=opt->miss;
    }
    if((mode!=2)&&((lscore>score)||(out_flag==1))) {memcpy(site,lsite,sizeof(unsigned int)*length);score = lscore;}
    if((mode==2)&&(out_flag==1)) return -1000;
    return score;
}
int site2cigar(char *chr,char *seq,unsigned int *site,int length,struct cigar_t *cigar,int *cigar_num,int max)
{
    int change_length;
    int snp_num = 0;
    int i;

    if (site[0]==0) {(*cigar_num) = 0;return -1000;}
    {
        if(nst_nt4_table[(int)chr[site[0]]]==nst_nt4_table[(int)seq[0]])
        {
            cigar[(*cigar_num)].c = 'M';
            cigar[(*cigar_num)].l = 1;
        }
        else
        {
            cigar[(*cigar_num)].c = 'X';
            cigar[(*cigar_num)].l = 1;
            snp_num++;
        }
    }

    for(i = 1;i<length;i++)
    {
        if (site[i]==0) {(*cigar_num)=0;return -1000;}

        change_length = site[i] - site[i-1]-1;
        if((change_length > 0)&&(change_length<=opt->change_length))
        {
            (*cigar_num)++;
            if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}
            cigar[(*cigar_num)].c = 'D';
            cigar[(*cigar_num)].l = change_length;
            snp_num++;
        }
        else if(change_length>opt->change_length)
        {
            (*cigar_num)++;
            cigar[(*cigar_num)].c = 'N';
            cigar[(*cigar_num)].l = change_length;
            snp_num++;
        }

        if(site[i] == site[i-1])
        {
            if(cigar[(*cigar_num)].c!='I')
            {
                (*cigar_num)++;
                if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}
                cigar[(*cigar_num)].c = 'I';
                cigar[(*cigar_num)].l = 1;
                snp_num++;
            }
            else cigar[(*cigar_num)].l++;
        }
        else
        {
            if(nst_nt4_table[(int)seq[i]]==nst_nt4_table[(int)chr[site[i]]])
            {
                if(cigar[(*cigar_num)].c!='M')
                {
                    (*cigar_num)++;
                    if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}

                    cigar[(*cigar_num)].c = 'M';
                    cigar[(*cigar_num)].l = 1;
                }
                else cigar[(*cigar_num)].l++;
            }
            else
            {
                if(cigar[(*cigar_num)].c!='X')
                {
                    (*cigar_num)++;
                    if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}

                    cigar[(*cigar_num)].c = 'X';
                    cigar[(*cigar_num)].l = 1;
                }
                else cigar[(*cigar_num)].l++;
                snp_num++;
            }
        }
    }
    (*cigar_num)++;
    if((*cigar_num)>=max) {(*cigar_num)=0;return -1000;}
    return snp_num;
}
int insert_front_tail(struct m_opt *opt,struct map_t *map,struct read_map_inf *read,char *chr,unsigned int l,unsigned int *site)
{
    int i,j;
    int front,back,node;
    int flag;

                    back = read->ref_site-map->node[0].site;
                    for(i = 0;i<read->front;i++)
                    {
                        flag = 1;
                        for(j = 0;j<map->node[back].inum;j++)
                        {
                            node = map->node[back].in_edge[j];
                            if((map->node[node].base==read->seq[read->front-1-i])&&(map->node[node].site==site[read->front-1-i]))
                                {map->node[node].rnum++;back = node;flag = 0;read->node[read->front-1-i] = node;break;}
                        }
                        if(flag)
                        {
                            for(j = i;j<read->front;j++)
                            {
                                front  = -1;
                                if(j != i) back = map->node_num-1;

                                if((front!=-1)&&(map->node[front].onum>=MAX_EDGE_NUM)) return 1;
                                if((back!=-1)&&(map->node[back].inum>=MAX_EDGE_NUM)) return 1;

                                if(new_node(opt,front,back,map,read->seq[read->front-j-1],site[read->front-1-j])==1) return 2;
                                read->node[read->front-1-j] = map->node_num-1;
                            }
                            break;
                        }
                    }
                    return 0;
}
int insert_read2map_front(struct m_opt *opt,struct map_t *map,struct read_map_inf *read,struct tree_stack *t_stack)
{
    int i = 0,j = 0;

    int front = 0;
    int back = 0;
    int flag = 0,i_flag;

    int node;

    if(read->front>0)
    {
        {
            int branch = 0;
            branch = back = read->ref_site-map->node[0].site;

            for(i = 0;i<read->front;i++)
            {
                flag = 1;
                for(j = 0;j<map->node[back].inum;j++)
                {
                    node = map->node[back].in_edge[j];
                    if((map->node[node].base==read->seq[read->front-1-i])&&(map->node[node].site==read->site[read->front-1-i]))
                        {map->node[node].rnum++;branch = back = node;flag = 0;read->node[read->front-1-i] = node;break;}
                }
                if(flag)
                {
                    i_flag = 0;
                    t_stack->top = -1;
                    push_stack(opt,t_stack,branch,0);

                    while((t_stack->top >= 0)&&(read->site[read->front-1-i]!=read->site[read->front-i])&&(read->site[read->front-1-i]!=0))//xian deep search
                    {
                        while((t_stack->data[t_stack->top].edge < map->node[t_stack->data[t_stack->top].node].inum)&&(t_stack->top<=opt->change_length))
                        {
                            node = map->node[t_stack->data[t_stack->top].node].in_edge[t_stack->data[t_stack->top].edge];
                            if((map->node[node].base==read->seq[read->front-1-i])&&(map->node[node].site==read->site[read->front-1-i]))
                            {
                                if((map->node[node].onum>=MAX_EDGE_NUM)) {i_flag = 0;break;}
                                if((back!=-1)&&(map->node[back].inum>=MAX_EDGE_NUM)) return 1;
                                read->node[read->front-1-i] = node;

                                map->node[back].in_edge[map->node[back].inum] = node;
                                map->node[back].inum++;

                                map->node[node].out_edge[map->node[node].onum] = back;
                                map->node[node].onum++;

                                map->node[node].rnum++;
                                back = node;
                                i_flag = 1;
                                break;
                            }
                            else push_stack(opt,t_stack,node,0);
                        }
                        if(i_flag)break;

                        while(t_stack->top >= 0)
                        {
                            t_stack->data[t_stack->top].edge++;
                            if (t_stack->data[t_stack->top].edge >= map->node[t_stack->data[t_stack->top].node].inum){pop_stack(t_stack);}
                            else break;
                        }
                    }
                    if(!i_flag)
                    {
                        front  = -1;

                        if((front!=-1)&&(map->node[front].onum>=MAX_EDGE_NUM)) return 1;
                        if((back!=-1)&&(map->node[back].inum>=MAX_EDGE_NUM)) return 1;

                        if(new_node(opt,front,back,map,read->seq[read->front-i-1],read->site[read->front-1-i])==1) return 2;
                        read->node[read->front-1-i] = map->node_num-1;
                        back = map->node_num-1;
                    }
                }
            }
        }
    }
    return 0;
}
int insert_back_tail(struct m_opt *opt,struct map_t *map,struct read_map_inf *read,char *chr,unsigned int l,unsigned int *site)
{
    int i,j;
    int front,back,node;
    int flag;

                for(i = read->front;i<=read->back;i++)
                    read->site[i] = read->ref_site-map->node[0].site+i-read->front;

                front = read->ref_site+read->back-read->front-1-map->node[0].site;
                for(i = read->back;i<read->length;i++)
                {
                flag = 1;
                for(j = 0;j<map->node[front].onum;j++)
                {
                    node = map->node[front].out_edge[j];
                    if((map->node[node].base==read->seq[i])&&(map->node[node].site==site[i]))
                        {map->node[node].rnum++;front = node;flag = 0;read->node[i] = node;break;}
                }
                if(flag)
                {
                    for(j = i;j<read->length;j++)
                    {
                        back  = -1;

                        if((front!=-1)&&(map->node[front].onum>=MAX_EDGE_NUM)) return 1;
                        if((back!=-1)&&(map->node[back].inum>=MAX_EDGE_NUM))return 1;

                        if(new_node(opt,front,back,map,read->seq[j],site[j])==1) return 2;
                        read->node[j] = map->node_num-1;
                        front = map->node_num-1;
                    }
                    break;
                }
                }
    return 0;
}
int insert_read2map_back(struct m_opt *opt,struct map_t *map,struct read_map_inf *read,struct tree_stack *t_stack)
{
    int i = 0,j = 0;

    int front = 0;
    int back = 0;
    int flag = 0,i_flag;
    int node;

    if(read->back<read->length)
    {
        {
            int branch;
            branch = front = read->ref_site+read->back-read->front-map->node[0].site;

            for(i = read->front;i<=read->back;i++)
                read->node[i] = read->ref_site-map->node[0].site+i-read->front;

            for(i = read->back+1;i<read->length;i++)
            {
                flag = 1;
                for(j = 0;j<map->node[front].onum;j++)
                {
                    node = map->node[front].out_edge[j];
                    if((map->node[node].base==read->seq[i])&&(map->node[node].site==read->site[i]))
                        {map->node[node].rnum++;branch=front = node;flag = 0;read->node[i] = node;break;}
                }
                if(flag)
                {
                    i_flag = 0;
                    t_stack->top = -1;
                    push_stack(opt,t_stack,branch,0);
                    while((t_stack->top >= 0)&&(read->site[i]!=read->site[i-1])&&(read->site[i]!=0))
                    {
                        while((t_stack->data[t_stack->top].edge < map->node[t_stack->data[t_stack->top].node].onum)&&(t_stack->top<=opt->change_length))
                        {
                            node = map->node[t_stack->data[t_stack->top].node].out_edge[t_stack->data[t_stack->top].edge];
                            if((map->node[node].base==read->seq[i])&&(map->node[node].site==read->site[i]))
                            {
                                if((front!=-1)&&(map->node[front].onum>=MAX_EDGE_NUM)) return 1;
                                if((map->node[node].inum>=MAX_EDGE_NUM)) {i_flag = 0;break;}

                                read->node[i] = node;

                                map->node[node].in_edge[map->node[node].inum] = front;
                                map->node[node].inum++;

                                map->node[front].out_edge[map->node[front].onum] = node;
                                map->node[front].onum++;

                                map->node[node].rnum++;
                                front = node;
                                i_flag = 1;
                                break;
                            }
                            else push_stack(opt,t_stack,node,0);
                        }
                        if(i_flag)break;

                        while(t_stack->top >= 0)
                        {
                            t_stack->data[t_stack->top].edge++;
                            if (t_stack->data[t_stack->top].edge >= map->node[t_stack->data[t_stack->top].node].onum){pop_stack(t_stack);}
                            else break;
                        }
                    }
                    if(!i_flag)
                    {
                        back  = -1;

                        if((front!=-1)&&(map->node[front].onum>=MAX_EDGE_NUM)) return 1;
                        if((back!=-1)&&(map->node[back].inum>=MAX_EDGE_NUM)) return 1;

                        if(new_node(opt,front,back,map,read->seq[i],read->site[i])==1) return 2;
                        read->node[i] = map->node_num-1;
                        front = map->node_num-1;
                    }
                }
            }
        }
    }
    return 0;
}

void insert_read2map(struct m_opt *opt,struct map_t *map,struct read_map_inf *read,char *chr,unsigned int l)
{
    //int i;

    struct cigar_t cigar[MAX_CIGAR_BUF];
    int cigar_total;

    unsigned int end_pos;

    unsigned int site[MAX_READ_LENGTH];

    process_cigar(read->cigar,cigar,&cigar_total,&end_pos);
    cigar2site(read->start_site,site,cigar,cigar_total);

    insert_front_tail(opt,map,read,chr,l,site);
    insert_back_tail(opt,map,read,chr,l,site);
}
void read2site(struct read_map_inf *read,char *chr,unsigned int l)
{
    struct cigar_t cigar[MAX_CIGAR_BUF];
    int cigar_total;

    unsigned int end_pos;

    process_cigar(read->cigar,cigar,&cigar_total,&end_pos);
    cigar2site(read->start_site,read->site,cigar,cigar_total);

    read->sstart = 0;
    read->send = 0;
    int read_pos = 0;
    int i;
    for(i = 0;i<cigar_total;i++)
    {
        if(cigar[i].c=='M')
        {
            if(cigar[i].l>6)
            {
                if(read->sstart==0)read->sstart=read_pos+5;
            }
            read_pos += cigar[i].l;
            if(cigar[i].l>6) read->send = read_pos-6;
        }
        else if((cigar[i].c=='I')||(cigar[i].c=='S')||(cigar[i].c=='X')) read_pos += cigar[i].l;
    }
    if(read->sstart==0)read->sstart=read->front;
    if(read->send==0)read->send=read->back;
}
void stat_map(struct m_opt *opt,struct map_t *map,int r_root,
              char *chr,unsigned int l,
            struct tree_stack *t_stack,struct seed_t *seed)
{
    int i = 0,x;
    int total = 0;
    int node;

    int dep_site = 0;

    for (i = 0;i<map->node_num;i++)
        map->node[i].change_flag = 0;

    for(x = 0;x<=r_root;x++)
    {
        map->node[x].change_flag = 1;
        map->node[x].score = 1-pow((1-(float)map->node[x].rnum/(map->node[x].rnum+SCORE_BALANCE)),map->node[x].rnum);

        if(i+1<=r_root) total = map->node[x].rnum-map->node[x+1].rnum;
        else total = map->node[x].rnum;
        if(map->node[x].onum>1)
        {
            t_stack->top = -1;
            push_stack(opt,t_stack,x,1);
            while(t_stack->top>=0)
            {
                while(t_stack->data[t_stack->top].edge<map->node[t_stack->data[t_stack->top].node].onum)
                {
                    if(t_stack->top>=1999)
                        break;
                    node = map->node[t_stack->data[t_stack->top].node].out_edge[t_stack->data[t_stack->top].edge];
                    if((node<=r_root)||(map->node[node].change_flag==1))
                        push_stack(opt,t_stack,node,map->node[node].onum);
                    else
                        push_stack(opt,t_stack,node,0);
                }
                if(opt->hard==1)
                {
                    dep_site = 0;
                if((t_stack->top<8)&&(map->node[t_stack->data[t_stack->top].node].onum==0))
                {
                        for(i = 1;i<=t_stack->top;i++)
                        {
                            if (map->node[t_stack->data[i].node].site-map->node[t_stack->data[i-1].node].site>opt->change_length)
                                map->node[t_stack->data[i].node].site = 0;
                        }
                }
                else if((t_stack->top>=8)&&(map->node[t_stack->data[t_stack->top].node].onum==0))
                {
                    for(i = t_stack->top;i>=0;i--)
                    {
                        if(map->node[t_stack->data[i].node].rnum >= 2) {dep_site = i;break;}
                    }
                    if(dep_site<8)
                    {
                        for(i = 1;i<=t_stack->top;i++)
                        {
                            if (map->node[t_stack->data[i].node].site-map->node[t_stack->data[i-1].node].site>opt->change_length)
                                map->node[t_stack->data[i].node].site = 0;
                        }
                    }
                }
                }

                for(i = 1;i<=t_stack->top;i++)
                {
                    if(map->node[t_stack->data[i].node].change_flag==0)
                    {
                        map->node[t_stack->data[i].node].change_flag=1;
                        if(map->node[t_stack->data[i].node].rnum>total) map->node[t_stack->data[i].node].score = 1;
                        else map->node[t_stack->data[i].node].score =1-pow((1-(float)map->node[t_stack->data[i].node].rnum/(total+SCORE_BALANCE)),map->node[t_stack->data[i].node].rnum);
                       
                    }
                }
                while(t_stack->top >= 0)
                {
                    t_stack->data[t_stack->top].edge++;
                    if (t_stack->data[t_stack->top].edge >= map->node[t_stack->data[t_stack->top].node].onum) pop_stack(t_stack);
                    else break;
                }
            }
        }
    }
    for(x = 0;x<=r_root;x++)
    {
        if(x==0) total = map->node[x].rnum;
        else total = map->node[x].rnum-map->node[x-1].rnum;
        if(map->node[x].inum>1)
        {
            t_stack->top = -1;
            push_stack(opt,t_stack,x,1);
            while(t_stack->top>=0)
            {
                while(t_stack->data[t_stack->top].edge<map->node[t_stack->data[t_stack->top].node].inum)
                {
                    if(t_stack->top>=1999)
                    break;
                    node = map->node[t_stack->data[t_stack->top].node].in_edge[t_stack->data[t_stack->top].edge];
                    if((node<=r_root)||(map->node[node].change_flag==1))
                        push_stack(opt,t_stack,node,map->node[node].inum);
                    else
                        push_stack(opt,t_stack,node,0);
                }
                if(opt->hard==1)
                {
                    dep_site = 0;
                if((t_stack->top<8)&&(map->node[t_stack->data[t_stack->top].node].inum==0))
                {
                    for(i = 1;i<=t_stack->top;i++)
                    {
                        if (map->node[t_stack->data[i-1].node].site-map->node[t_stack->data[i].node].site>opt->change_length)
                            map->node[t_stack->data[i].node].site = 0;
                    }
                }
                else if((t_stack->top>=8)&&(map->node[t_stack->data[t_stack->top].node].inum==0))
                {
                    for(i = t_stack->top;i>=0;i--)
                    {
                        if(map->node[t_stack->data[i].node].rnum >= 2) {dep_site = i;break;}
                    }
                    if(dep_site<8)
                    {
                        for(i = 1;i<=t_stack->top;i++)
                        {
                            if (map->node[t_stack->data[i-1].node].site-map->node[t_stack->data[i].node].site>opt->change_length)
                            map->node[t_stack->data[i].node].site = 0;
                        }
                    }
                }
                }

                for(i = 1;i<=t_stack->top;i++)
                {
                    if(map->node[t_stack->data[i].node].change_flag==0)
                    {
                        map->node[t_stack->data[i].node].change_flag=1;
                        if(map->node[t_stack->data[i].node].rnum>total) map->node[t_stack->data[i].node].score = 1;
                        else map->node[t_stack->data[i].node].score =1-pow((1-(float)map->node[t_stack->data[i].node].rnum/(total+SCORE_BALANCE)),map->node[t_stack->data[i].node].rnum);
                    }
                }
                while(t_stack->top >= 0)
                {
                    t_stack->data[t_stack->top].edge++;
                    if (t_stack->data[t_stack->top].edge >= map->node[t_stack->data[t_stack->top].node].inum) pop_stack(t_stack);
                    else break;
                }
            }
        }
    }
}
int output_map_align(struct m_opt *opt,struct read_map_inf *read,int read_num,struct map_t *map,char *chr,unsigned int chr_end,struct tree_stack *t_stack,int r_root)
{
    int i = 0,j = 0,k = 0,x = 0;
    int flag = 0;
    int node;

    int change_length = 0,change_type;//change_type 0 I,1 D,2 X,3 N,4 M,5 S
    float change_score = 0;

    float score = 0;
    int cover = 0;
    float gap_score = 0;
    char cigarS[2*MAX_READ_LENGTH];
    char temp_cigar[2*MAX_READ_LENGTH];

    float best_score = 0;
    float best_gap = 0;
    char best_cigar[2*MAX_READ_LENGTH];
    unsigned int best_ref_site = 0;
    unsigned int best_end_site = 0;
    unsigned int ref_site = 0;

    int best_fscore = 0;
    //int best_start = 0;

    unsigned int start_site = 0;
    int snp_num;
    int change_num;
    int gap_num;

    struct cigar_t cigar[MAX_CIGAR_BUF];
    int cigar_num = 0;
    int read_site = 0;
    int length = 0;

    int temp_node[MAX_READ_LENGTH];
    int start_flag = 0;

    for(i = 0; i<read_num; i++)
    {
        if(read[i].output_flag==1) continue;

        best_cigar[0] = '*';
        best_cigar[1] = '\0';
        best_score = 0;
        best_gap = 0;
        best_ref_site = 0;

        {
        best_fscore = -100;
        //best_start = 0;
        best_gap = 0;
        start_flag = 0;

        node = read[i].ref_site-map->node[0].site;
        //if(read[i].front == 0) best_start = node;
        if((read[i].front != 0)&&(read[i].node[0]==0))
        {
            t_stack->top = -1;
            push_stack(opt,t_stack,node,0);
            j = read[i].front-1;

            while(t_stack->top >= 0)
            {
                while((t_stack->data[t_stack->top].edge<map->node[t_stack->data[t_stack->top].node].inum)&&(j>=0))
                {
                    node = map->node[t_stack->data[t_stack->top].node].in_edge[t_stack->data[t_stack->top].edge];
                    if(nst_nt4_table[(int)map->node[node].base]== nst_nt4_table[(int)read[i].seq[j]])
                    {
                        if ((map->node[node].site<0)||(map->node[node].site>chr_end)) break;
                        if(j>=read[i].sstart)//if(read[i].site[j]!=0)
                        {
                            if(map->node[node].site==read[i].site[j])
                            {
                                push_stack(opt,t_stack,node,0);
                                j--;
                            }
                            else t_stack->data[t_stack->top].edge++;
                        }
                        else
                        {
                            push_stack(opt,t_stack,node,0);
                            j--;
                        }
                    }
                    else t_stack->data[t_stack->top].edge++;
                }
                if(j < 0)
                {
                    score = 0;
                    change_num = 0;
                    gap_score = 0;
                    gap_num = 0;

                    for(k = t_stack->top;k>0;k--)
                    {
                        if(k == t_stack->top)
                        {
                            if (map->node[t_stack->data[k].node].site == 0)
                            {
                                change_type = 5;
                                change_length = 1;
                            }
                            else if(nst_nt4_table[(int)map->node[t_stack->data[k].node].base]
                            ==nst_nt4_table[(int)chr[map->node[t_stack->data[k].node].site]])
                            {
                                change_type = 4;
                                change_length = 1;
                                score+=opt->match;
                            }
                            else
                            {
                                change_num++;
                                change_type = 2;
                                change_length = 1;
                                if(map->node[t_stack->data[k].node].score * map->node[t_stack->data[k].node].weigth>(0.9))
                                    score-=(opt->miss-1);
                                else {score-=opt->miss;snp_num++;}
                            }
                        }
                        else
                        {
                            if (map->node[t_stack->data[k].node].site == 0)
                            {
                                if(change_type!=5)
                                {
                                    if(change_type==0)
                                    {
                                        change_num++;
                                        if(change_score > (0.9)) score -= (opt->miss-1);
                                        else {score -= opt->miss;snp_num++;}
                                    }
                                    change_type = 5;
                                    change_length = 1;
                                }
                                else
                                    change_length++;
                            }
                            else if(map->node[t_stack->data[k].node].site == (map->node[t_stack->data[k-1].node].site-1))//MX
                            {
                                if(nst_nt4_table[(int)map->node[t_stack->data[k].node].base]
                                ==nst_nt4_table[(int)chr[map->node[t_stack->data[k].node].site]])//M
                                {
                                    if(change_type!=4)//IX
                                    {
                                        if(change_type==0)
                                        {
                                            change_num++;
                                            if(change_score > (0.9)) score -= (opt->miss-1);
                                            else {score -= opt->miss;snp_num++;}
                                        }
                                        change_type = 4;
                                        change_length = 0;
                                    }
                                    score+=opt->match;
                                    change_length++;
                                }
                                else//X
                                {
                                    if(change_type!=2)//IX
                                    {
                                        if(change_type==0)
                                        {
                                            if(change_score > (0.9)) score -= (opt->miss-1);
                                            else {score -= opt->miss;snp_num++;}
                                        }
                                        change_type = 2;
                                        change_length = 0;
                                    }
                                    if(map->node[t_stack->data[k].node].score * map->node[t_stack->data[k].node].weigth>(0.9))
                                        score-=(opt->miss-1);
                                    else {score-=opt->miss;snp_num++;}
                                    change_length++;
                                    change_num++;
                                }
                            }
                            else if (map->node[t_stack->data[k].node].site == (map->node[t_stack->data[k-1].node].site))//I
                            {
                                if(change_type!=0)//IX
                                {
                                    change_num++;
                                    change_type = 0;
                                    change_length = 0;
                                    change_score = (float)1;
                                }
                                if(map->node[t_stack->data[k].node].score * map->node[t_stack->data[k].node].weigth<(0.9))
                                    change_score=((float)change_score)*map->node[t_stack->data[k].node].score * ((float)map->node[t_stack->data[k].node].weigth);
                                change_length++;
                            }
                            else if (map->node[t_stack->data[k].node].site <= (map->node[t_stack->data[k-1].node].site-1))//DN
                            {
                                change_num++;
                                if(change_type==5)
                                {
                                    if(nst_nt4_table[(int)map->node[t_stack->data[k].node].base]
                                    ==nst_nt4_table[(int)chr[map->node[t_stack->data[k].node].site]])
                                    {
                                        change_type = 4;
                                        change_length = 1;
                                        score+=opt->match;
                                    }
                                    else
                                    {
                                        change_type = 2;
                                        change_length = 1;
                                        if(map->node[t_stack->data[k].node].score * map->node[t_stack->data[k].node].weigth>(0.9))
                                            score-=(opt->miss-1);
                                        else {score-=opt->miss;snp_num++;}
                                    }
                                    continue;
                                }
                                if(change_type==0)
                                {
                                    change_num++;
                                    if(change_score > (0.9)) score -= (opt->miss-1);
                                    else {score -= opt->miss;snp_num++;}
                                }

                                change_length = map->node[t_stack->data[k].node].site - map->node[t_stack->data[k-1].node].site-1;
                                if(change_length<=opt->change_length)
                                {
                                    if((map->node[t_stack->data[k].node].score * map->node[t_stack->data[k].node].weigth>(0.9))
                                    &&(map->node[t_stack->data[k-1].node].score * map->node[t_stack->data[k-1].node].weigth>(0.9)))
                                        score-=(opt->gap-1);
                                    else
                                        {score-=opt->gap;snp_num++;}
                                }
                                else
                                {
                                    if((map->node[t_stack->data[read_site-1].node].rnum>=2)&&(map->node[t_stack->data[read_site].node].rnum>=2))
                                        score-=(opt->splice-1.5);
                                    else
                                        {score-=opt->splice;snp_num++;}
                                    gap_score+=map->node[t_stack->data[read_site-1].node].score+map->node[t_stack->data[read_site].node].score;
                                    gap_num++;
                                }

                                if(nst_nt4_table[(int)map->node[t_stack->data[k].node].base]
                                ==nst_nt4_table[(int)chr[map->node[t_stack->data[k].node].site]])
                                {
                                    change_type = 4;
                                    change_length = 1;
                                    score+=opt->match;
                                }
                                else
                                {
                                    change_num++;
                                    change_type = 2;
                                    change_length = 1;
                                    if(map->node[t_stack->data[k].node].score * map->node[t_stack->data[k].node].weigth>(0.9))
                                        score-=(opt->miss-1);
                                    else {score-=opt->miss;snp_num++;}
                                }
                            }
                            else if (map->node[t_stack->data[k].node].site > (map->node[t_stack->data[k-1].node].site))
                            {
                                score = -100;break;
                            }
                        }
                    }
                    gap_score = gap_score/gap_num;
                    if(score>best_fscore)
                    {
                        best_fscore = score;
                       // best_start = t_stack->data[t_stack->top].node;
                        best_gap = gap_score;
                        start_flag = 1;
                        for(k = 0;k<=t_stack->top;k++)
                            temp_node[k] = t_stack->data[t_stack->top-k].node;
                    }
                    else if(score==best_fscore)
                    {
                        if(gap_score>best_score)
                        {
                        best_fscore = score;
                        //best_start = t_stack->data[t_stack->top].node;
                        best_gap = gap_score;
                        start_flag = 1;
                        for(k = 0;k<=t_stack->top;k++)
                            temp_node[k] = t_stack->data[t_stack->top-k].node;
                        }
                    }
                }
                while(t_stack->top >= 0)
                {
                    t_stack->data[t_stack->top].edge++;
                    if (t_stack->data[t_stack->top].edge >= map->node[t_stack->data[t_stack->top].node].inum) {pop_stack(t_stack);j++;}
                    else break;
                }
            }

        if(start_flag == 0)
            continue;
        memcpy(&(read[i].node[0]),temp_node,sizeof(unsigned int)*(read[i].front));
        }
        }


        {
            if(read[i].node[read[i].length-1]==0) start_flag = 0;
            else {start_flag = 1;j = read[i].length;}

            start_site = 0;
            t_stack->top = -1;
            node = read[i].node[0];
            if((map->node[node].site>chr_end)||(map->node[node].site<0)) continue;
            //if(read[i].front != 0)
            {
                push_stack(opt,t_stack,node,0);
                for(k = 1;k<read[i].length;k++)
                {
                    if(read[i].node[k] == 0)
                    {
                        j = k;
                        break;
                    }
                    for(x = 0;x<map->node[read[i].node[k-1]].onum;x++)
                    {
                        if(map->node[read[i].node[k-1]].out_edge[x]==read[i].node[k])
                            t_stack->data[t_stack->top].edge = x;
                    }
                    push_stack(opt,t_stack,read[i].node[k],0);
                }
            }

            while(t_stack->top >= read[i].back)
            {
                while((t_stack->data[t_stack->top].edge<map->node[t_stack->data[t_stack->top].node].onum)&&(j<read[i].length))
                {
                    node = map->node[t_stack->data[t_stack->top].node].out_edge[t_stack->data[t_stack->top].edge];
                    if(nst_nt4_table[(int)map->node[node].base]
                           == nst_nt4_table[(int)read[i].seq[j]])
                    {
                        if ((map->node[node].site<0)||(map->node[node].site>chr_end)) break;
                        if(j<=read[i].send)
                        {
                            if(map->node[map->node[t_stack->data[t_stack->top].node].out_edge[t_stack->data[t_stack->top].edge]].site==read[i].site[j])
                            {
                                push_stack(opt,t_stack,map->node[t_stack->data[t_stack->top].node].out_edge[t_stack->data[t_stack->top].edge],0);
                                j++;
                            }
                            else t_stack->data[t_stack->top].edge++;
                        }
                        else
                        {
                            push_stack(opt,t_stack,map->node[t_stack->data[t_stack->top].node].out_edge[t_stack->data[t_stack->top].edge],0);
                            j++;
                        }
                    }
                    else t_stack->data[t_stack->top].edge++;
                }
                if(j == read[i].length)
                {
                    score = 0;
                    gap_score = 0;
                    cigarS[0] = '\0';

                    flag = 1;
                    snp_num = 0;
                    change_num = 0;

                    cigar_num = 0;
                    if(map->node[t_stack->data[0].node].site==0)
                    {
                        cigar[cigar_num].c = 'S';
                        cigar[cigar_num].l = 1;
                    }
                    else
                    {
                        start_site = map->node[t_stack->data[0].node].site;
                        if(nst_nt4_table[(int)read[i].seq[0]]==nst_nt4_table[(int)chr[map->node[t_stack->data[0].node].site]])
                        {
                            cigar[cigar_num].c = 'M';
                            cigar[cigar_num].l = 1;
                        }
                        else
                        {
                            cigar[cigar_num].c = 'X';
                            cigar[cigar_num].l = 1;
                        }
                    }

                    for(k = 1;k<= t_stack->top;k++)
                    {
                        if((start_site==0)&&(map->node[t_stack->data[k].node].site!=0)) start_site = map->node[t_stack->data[k].node].site;

                        if (map->node[t_stack->data[k].node].site < map->node[t_stack->data[k-1].node].site) {flag = 0;break;}
                        if(map->node[t_stack->data[k].node].site==0)
                        {
                            if(cigar[cigar_num].c!='S')
                            {
                                cigar_num++;
                                if(cigar_num>=MAX_CIGAR_BUF) {flag = 0;break;}

                                cigar[cigar_num].c = 'S';
                                cigar[cigar_num].l = 1;
                            }
                            else cigar[cigar_num].l++;
                        }
                        else
                        {
                            if(cigar[cigar_num].c=='S')
                            {
                                cigar_num++;
                                if(cigar_num>=MAX_CIGAR_BUF) {flag = 0;break;}
                                if(nst_nt4_table[(int)read[i].seq[k]]==nst_nt4_table[(int)chr[map->node[t_stack->data[k].node].site]])
                                {
                                    cigar[cigar_num].c = 'M';
                                    cigar[cigar_num].l = 1;
                                }
                                else
                                {
                                    cigar[cigar_num].c = 'X';
                                    cigar[cigar_num].l = 1;
                                }
                            }
                            else
                            {
                                change_length = map->node[t_stack->data[k].node].site - map->node[t_stack->data[k-1].node].site-1;
                                if((change_length > 0)&&(change_length<=opt->change_length))
                                {
                                    cigar_num++;
                                    if(cigar_num>=MAX_CIGAR_BUF) {flag = 0;break;}
                                    cigar[cigar_num].c = 'D';
                                    cigar[cigar_num].l = change_length;
                                }
                                else if(change_length>opt->change_length)
                                {
                                    cigar_num++;
                                    if(cigar_num>=MAX_CIGAR_BUF) {flag = 0;break;}
                                    cigar[cigar_num].c = 'N';
                                    cigar[cigar_num].l = change_length;
                                }

                                if(map->node[t_stack->data[k].node].site == (map->node[t_stack->data[k-1].node].site))
                                {
                                    if(cigar[cigar_num].c!='I')
                                    {
                                        cigar_num++;
                                        if(cigar_num>=MAX_CIGAR_BUF) {flag = 0;break;}

                                        cigar[cigar_num].c = 'I';
                                        cigar[cigar_num].l = 1;
                                    }
                                    else cigar[cigar_num].l++;
                                }
                                else
                                {
                                    if(nst_nt4_table[(int)read[i].seq[k]]==nst_nt4_table[(int)chr[map->node[t_stack->data[k].node].site]])
                                    {
                                        if(cigar[cigar_num].c!='M')
                                        {
                                            cigar_num++;
                                            if(cigar_num>=MAX_CIGAR_BUF) {flag = 0;break;}

                                            cigar[cigar_num].c = 'M';
                                            cigar[cigar_num].l = 1;
                                        }
                                        else cigar[cigar_num].l++;
                                    }
                                    else
                                    {
                                        if(cigar[cigar_num].c!='X')
                                        {
                                            cigar_num++;
                                            if(cigar_num>=MAX_CIGAR_BUF) {flag = 0;break;}

                                            cigar[cigar_num].c = 'X';
                                            cigar[cigar_num].l = 1;
                                        }
                                        else cigar[cigar_num].l++;
                                    }
                                }
                            }
                        }
                    }
                    cigar_num++;
                    if(flag)
                    {
                        cigarS[0] = '\0';
                        read_site = 0;
                        ref_site = start_site;

                        score = 0;
                        cover = 0;
                        snp_num = 0;
                        change_num = 0;
                        length = 0;
                        gap_num = 0;
                        gap_score = 0;

                        for(k = 0;k<cigar_num;k++)
                        {
                            if((cigar[k].c=='N')&&(read_site<=15))
                            {
                                if((map->node[t_stack->data[read_site-1].node].rnum<2)||(map->node[t_stack->data[read_site].node].rnum<2))
                                {
                                cigarS[0] = '\0';
                                sprintf(temp_cigar,"%dS",read_site);
                                strcat(cigarS,temp_cigar);
                                k++;
                                score = 0;
                                cover = 0;
                                snp_num = 0;
                                change_num = 0;
                                length = 0;
                                gap_num = 0;
                                gap_score = 0;
                                start_site = map->node[t_stack->data[read_site].node].site;
                                ref_site=start_site;
                                }
                            }
                            else if((cigar[k].c=='N')&&(read[i].length-read_site<=15))
                            {
                                if((map->node[t_stack->data[read_site-1].node].rnum<2)||(map->node[t_stack->data[read_site].node].rnum<2))
                                {
                                sprintf(temp_cigar,"%dS",read[i].length-read_site);
                                strcat(cigarS,temp_cigar);
                                break;
                                }
                            }

                            if(cigar[k].c=='M') {score+=cigar[k].l*opt->match;cover+=cigar[k].l*opt->match;ref_site+=cigar[k].l;read_site+=cigar[k].l;length+=cigar[k].l;}
                            else if(cigar[k].c=='X')
                            {
                                length+=cigar[k].l;
                                change_num+=cigar[k].l;

                                for(x = 0;x<cigar[k].l;x++)
                                {
                                    if(map->node[t_stack->data[x+read_site].node].score * map->node[t_stack->data[x+read_site].node].weigth>(0.9))
                                        score-=(opt->miss-1);
                                    else {score-=opt->miss;cover-=opt->miss;snp_num++;}
                                }
                                read_site+=cigar[k].l;
                                ref_site+=cigar[k].l;
                            }
                            else if(cigar[k].c=='I')
                            {
                                length+=cigar[k].l;
                                change_num++;

                                change_score = (float)1;
                                for(x = 0;x<cigar[k].l;x++)
                                {
                                    if(map->node[t_stack->data[x+read_site].node].score * map->node[t_stack->data[x+read_site].node].weigth<(0.9))
                                        change_score=((float)change_score)*map->node[t_stack->data[x+read_site].node].score * ((float)map->node[t_stack->data[x+read_site].node].weigth);
                                }
                                if(change_score > (0.9)) score -= (opt->miss-1);
                                else {score -= opt->miss;cover-=opt->miss;snp_num++;}

                                read_site+=cigar[k].l;
                            }
                            else if(cigar[k].c=='D')
                            {
                                ref_site+=cigar[k].l;
                                change_num++;
                                if((map->node[t_stack->data[read_site].node].score * map->node[t_stack->data[read_site].node].weigth<(0.9))
                                    ||(map->node[t_stack->data[read_site-1].node].score * map->node[t_stack->data[read_site-1].node].weigth<(0.9)))
                                   {score-=opt->gap;}
                                else
                                    score-=(opt->gap-1);
                            }
                            else if(cigar[k].c=='N')
                            {
                                ref_site+=cigar[k].l;
                                change_num++;
                                if((map->node[t_stack->data[read_site-1].node].rnum<2)||(map->node[t_stack->data[read_site].node].rnum<2))
                                    {score-=opt->splice;}
                                else
                                    score-=(opt->splice-1.5);
                                gap_score += map->node[t_stack->data[read_site-1].node].score+map->node[t_stack->data[read_site].node].score;
                                gap_num++;
                            }
                            else if(cigar[k].c=='S') read_site+=cigar[k].l;

                            sprintf(temp_cigar,"%d%c",cigar[k].l,cigar[k].c);
                            strcat(cigarS,temp_cigar);
                        }
                        gap_score = gap_score/gap_num;
                        if(score>best_score)
                        {
                            if((change_num<=(float)(0.12*read[i].length))&&(score>=(float)(0.8*read[i].length)))
                            {
                                strcpy(best_cigar,cigarS);
                                best_score = score;
                                best_ref_site = start_site;
                                best_end_site = ref_site-1;
                                best_gap = gap_score;
                            }
                        }
                        else if(score==best_score)
                        {
                            if(best_gap<gap_score)
                            {
                                strcpy(best_cigar,cigarS);
                                best_score = score;
                                best_ref_site = start_site;
                                best_end_site = ref_site-1;
                                best_gap = gap_score;
                            }
                        }
                    }
                }
                if(start_flag==1) break;
                pop_stack(t_stack);
                j--;
                while(t_stack->top >= read[i].back)
                {
                    t_stack->data[t_stack->top].edge++;
                    if (t_stack->data[t_stack->top].edge >= map->node[t_stack->data[t_stack->top].node].onum) {pop_stack(t_stack);j--;}
                    else break;
                }
            }
        }
        if((best_cigar[0]!='*')&&(best_score>read[i].score)&&(strlen(best_cigar)<MAX_READ_LENGTH))
        {
            strcpy(read[i].cigar,best_cigar);
            read[i].score = best_score;
            read[i].start_site = best_ref_site;
            read[i].end_site = best_end_site;
        }
        if(read[i].psite !=0) read[i].score +=2;
    }
    return 0;
}
int front_cmp(const void *a,const void *b)
{
    struct read_map_inf *EA,*EB;
    EA = (struct read_map_inf *)a;
    EB = (struct read_map_inf *)b;

    if(EA->ref_site>EB->ref_site) return -1;
    else if (EA->ref_site<EB->ref_site) return 1;
    else return (EB->front-EA->front);
}
int back_cmp(const void *a,const void *b)
{
    struct read_map_inf *EA,*EB;
    EA = (struct read_map_inf *)a;
    EB = (struct read_map_inf *)b;

    if(EA->ref_site+EA->back-EA->front >=EB->ref_site+EB->back-EB->front) return 1;
    else if (EA->ref_site+EA->back-EA->front <EB->ref_site+EB->back-EB->front) return -1;
    else return (EA->back-EB->back);
}
#define DIF_LENGTH 10//12//15
struct contig_t
{
    int base[4][2*MAX_READ_LENGTH];
    unsigned int site[2*MAX_READ_LENGTH];
    unsigned int start;
    int contig_start;
};
char BASE[5] = {'A','C','G','T'};
int align_read(struct m_opt *opt,char *chr,unsigned int l,int chr_order,struct read_map_inf *read,int read_num,struct contig_t *contig,unsigned int **hash,unsigned int *hash_num,struct seed_t *seed)
{
    int i = 0,j,k;
    int flag;
    int NW,TH;
    int istart,jstart;
    int abs;
    int length;

    int contig_num = 0;
    //int read_start = 0;

    unsigned int site[2*MAX_READ_LENGTH];
    char text[2*MAX_READ_LENGTH];
    char c;
    int r;

    struct cigar_t t_cigar[MAX_CIGAR_BUF];
    int t_num;

    int N;

    qsort(read,read_num,sizeof(struct read_map_inf),front_cmp);
    for(i = 0;i<read_num;i++)
    {
        read2site(read+i,chr,l);
        flag = 1;
        //while(read[read_start].ref_site>read[i].ref_site+read[i].back-read[i].front) read_start++;
        N = 0;
        for(k = 0;k<read[i].front;k++)
        {
            if(nst_nt4_table[(int)(read[i].seq[k])]>=4)
                {N=1;break;}
        }
        if(N)
        {
            read[i].contig_order=-1;
            //memset(read[i].site, 0, (read[i].front)*sizeof(unsigned int)/sizeof(char));
            continue;
        }
        for(j = 0;j<i;j++)
        {
            if(read[j].contig_order==-1)continue;
            if(read[j].ref_site>read[i].ref_site+read[i].back-read[i].front) continue;
            if(read[i].ref_site>read[j].ref_site+read[j].back-read[j].front) continue;

            abs = (read[j].ref_site-read[i].ref_site);
            jstart = read[j].front-abs;
            istart = read[i].front;
            length = min(jstart,istart);

            if(istart+length-5<=read[i].front) continue;

            NW = 0;
            TH = 1+length/DIF_LENGTH;
            for(k = 0;k<length;k++)
            {
                if(nst_nt4_table[(int)read[i].seq[istart-k]]!=nst_nt4_table[(int)read[j].seq[jstart-k]]) NW++;
                if(NW>TH) break;
            }
            if(NW<=TH)
            {
                if(read[i].ref_site>contig[read[j].contig_order].start) continue;

                read[i].contig_order = read[j].contig_order;
                read[i].contig_start = read[j].contig_start+abs;
                if (read[i].contig_start>contig[read[i].contig_order].contig_start)
                    contig[read[i].contig_order].contig_start = read[i].contig_start;

                for(k = read[i].front-1+abs;k>=0;k--)
                    contig[read[i].contig_order].base[nst_nt4_table[(int)(read[i].seq[k])]][read[j].contig_start+abs+read[i].front-1-k]++;

                flag = 0;
                break;
            }
        }
        if((flag)&&(read[i].front!=0))
        {
            read[i].contig_order = contig_num;
            read[i].contig_start = 0;
            contig[contig_num].start = read[i].ref_site;
            contig[contig_num].contig_start = read[i].contig_start;
            for(k = read[i].front-1;k>=0;k--)
                contig[contig_num].base[nst_nt4_table[(int)(read[i].seq[k])]][read[i].front-1-k]++;
            contig_num++;
        }
    }
    int max = 0;
    for(i = 0;i<contig_num;i++)
    {
        length = 0;
        for(j = 0;j<2*MAX_READ_LENGTH;j++)
        {
            max = 0;
            for(k = 1;k<4;k++)
            {
                if(contig[i].base[k][j]>contig[i].base[max][j]) max = k;
            }
            if(contig[i].base[max][j]==0) {length = j;break;}
            text[j] = BASE[max];
        }
        text[j] = '\0';
        for(j = 0;j<length;j++)
        {
            for(k = 0;k<4;k++)
                contig[i].base[k][j] = 0;
        }
        for(j = 0;j<length/2;j++)
        {
            c = text[j];
            text[j]=text[length-1-j];
            text[length-1-j] = c;
        }
        for(j = 0;j<contig[i].contig_start;j++)
            contig[i].site[j] = contig[i].start-1-j;
        if((length-contig[i].contig_start)<=0) {memset(site, 0, MAX_READ_LENGTH*sizeof(unsigned int)/sizeof(char));continue;}
        memset(site, 0, (length-contig[i].contig_start)*sizeof(unsigned int)/sizeof(char));
        text[length-contig[i].contig_start] = '\0';
        r=area_align(opt,chr,text,1,(contig[i].start-contig[i].contig_start>AREA)?(contig[i].start-contig[i].contig_start-AREA):0,contig[i].start-contig[i].contig_start,hash,hash_num,seed,site,chr_order);

        if(r!=-1000)
        {
            t_num = 0;
            if(site2cigar(chr,text,site,length-contig[i].contig_start,t_cigar,&t_num,MAX_CIGAR_BUF)>=0)
            {
                check_read_splice(text,text,t_cigar,&t_num,site[0]+opt->chr->list[chr_order].start_site,chr_order,0);
                cigar2site(site[0],site,t_cigar,t_num);
            }

            for(j = 0;j<length-contig[i].contig_start;j++)
                contig[i].site[j+contig[i].contig_start] = site[length-contig[i].contig_start-j-1];
        }
        else
        {
            for(j = 0;j<length-contig[i].contig_start;j++)
                contig[i].site[j+contig[i].contig_start] = 0;
        }
    }
    for(i = 0;i<read_num;i++)
    {
        if(read[i].contig_order==-1) continue;
        if(contig[read[i].contig_order].site[read[i].contig_start+read[i].front-1]==0) continue;
        for(k = 0;k<read[i].front;k++)
           read[i].site[read[i].front-1-k] = contig[read[i].contig_order].site[read[i].contig_start+k];
    }

    qsort(read,read_num,sizeof(struct read_map_inf),back_cmp);
    //read_start = 0;
    contig_num = 0;
    for(i = 0;i<read_num;i++)
    {
        flag = 1;
        N = 0;
        for(k = read[i].back;k<read[i].length;k++)
        {
            if(nst_nt4_table[(int)(read[i].seq[k])]>=4)
                {N=1;break;}
        }
        if(N)
        {
            read[i].contig_order=-1;
            continue;
        }

        for(j = 0;j<i;j++)
        {
            if(read[j].contig_order==-1)continue;
            if(read[j].ref_site>read[i].ref_site+read[i].back-read[i].front) continue;
            if(read[i].ref_site>read[j].ref_site+read[j].back-read[j].front) continue;

            abs = (read[i].ref_site+read[i].back-read[i].front)-(read[j].ref_site+read[j].back-read[j].front);
            jstart = read[j].back+1;
            istart = read[i].back+1-abs;
            length = min(read[j].length-jstart,read[i].length-istart);

            if(istart+length-5<=read[i].back) continue;

            NW = 0;
            TH = 1+length/DIF_LENGTH;
            for(k = 0;k<length;k++)
            {
                if(nst_nt4_table[(int)read[i].seq[istart+k]]!=nst_nt4_table[(int)read[j].seq[jstart+k]]) NW++;
                if(NW>TH) break;
            }
            if(NW<=TH)
            {
                if(read[i].ref_site+read[i].back-read[i].front<contig[read[j].contig_order].start) continue;

                read[i].contig_order = read[j].contig_order;
                read[i].contig_start = read[j].contig_start+abs;
                if (read[i].contig_start>contig[read[i].contig_order].contig_start)
                    contig[read[i].contig_order].contig_start = read[i].contig_start;

                for(k = 0;k<read[i].length-read[i].back-1+abs;k++)
                    contig[read[i].contig_order].base[nst_nt4_table[(int)(read[i].seq[read[i].back-abs+1+k])]][read[j].contig_start+k]++;

                flag = 0;
                break;
            }
        }
        if((flag)&&(read[i].back+1!=read[i].length))
        {
            read[i].contig_order = contig_num;
            read[i].contig_start = 0;
            contig[contig_num].start = read[i].ref_site+read[i].back-read[i].front;
            contig[contig_num].contig_start = read[i].contig_start;
            for(k = 0;k<read[i].length-read[i].back-1;k++)
                contig[contig_num].base[nst_nt4_table[(int)(read[i].seq[k+read[i].back+1])]][k]++;
            contig_num++;
        }
    }
    for(i = 0;i<contig_num;i++)
    {
        length = 0;
        for(j = contig[i].contig_start;j<2*MAX_READ_LENGTH;j++)
        {
            max = 0;
            for(k = 1;k<4;k++)
            {
                if(contig[i].base[k][j]>contig[i].base[max][j]) max = k;
            }
            if(contig[i].base[max][j]==0) {length = j;break;}
            text[j] = BASE[max];
        }
        text[j] = '\0';
        for(j = 0;j<length;j++)
        {
            for(k = 0;k<4;k++)
                contig[i].base[k][j] = 0;
        }
        for(j = 0;j<contig[i].contig_start;j++)
            contig[i].site[j] = contig[i].start+1+j;
        if((length-contig[i].contig_start)<=0) {memset(site, 0, MAX_READ_LENGTH*sizeof(unsigned int)/sizeof(char));continue;}
        memset(&(contig[i].site[contig[i].contig_start]), 0, (length-contig[i].contig_start)*sizeof(unsigned int)/sizeof(char));
        r=area_align(opt,chr,&(text[contig[i].contig_start]),0,contig[i].start+contig[i].contig_start,contig[i].start+contig[i].contig_start+AREA,hash,hash_num,seed,&(contig[i].site[contig[i].contig_start]),chr_order);

        if(r!=-1000)
        {
            t_num = 0;
            if(site2cigar(chr,&(text[contig[i].contig_start]),&(contig[i].site[contig[i].contig_start]),(length-contig[i].contig_start),t_cigar,&t_num,MAX_CIGAR_BUF)>=0)
            {
            check_read_splice(&(text[contig[i].contig_start]),&(text[contig[i].contig_start]),t_cigar,&t_num,(contig[i].site[contig[i].contig_start])+opt->chr->list[chr_order].start_site,chr_order,0);
            cigar2site((contig[i].site[contig[i].contig_start]),&(contig[i].site[contig[i].contig_start]),t_cigar,t_num);
            }
        }
        else memset(&(contig[i].site[contig[i].contig_start]), 0, (length-contig[i].contig_start)*sizeof(unsigned int)/sizeof(char));
    }
    for(i = 0;i<read_num;i++)
    {
        if(read[i].contig_order==-1) continue;
        if(contig[read[i].contig_order].site[read[i].contig_start+(read[i].length-read[i].back-1)-1]==0) continue;
        if(read[i].back+1>=read[i].length) continue;
         memcpy(&(read[i].site[read[i].back+1]),&(contig[read[i].contig_order].site[read[i].contig_start]),sizeof(unsigned int)*(read[i].length-read[i].back-1));
    }

    return 0;
}
#define MAX_COVERAGE (5000)
void *map_align_core(void* arg)
{
    struct job_snp *job = (struct job_snp *)arg;

    struct map_t *map = (struct map_t *)calloc(1,sizeof(struct map_t));
    map->node = (struct node_table *)calloc(MAX_NODE_NUM+1,sizeof(struct node_table));
    map->node_num = 0;

    struct tree_stack *t_stack = (struct tree_stack *)calloc(1,sizeof(struct tree_stack));
    t_stack->data = (struct stack_element *)calloc(2000,sizeof(struct stack_element));


	struct seed_t *seed = (struct seed_t *)calloc(SEED_BUF_LENGTH+1,sizeof(struct seed_t));
	struct read_map_inf *read = (struct read_map_inf *)calloc(MAX_COVERAGE,sizeof(struct read_map_inf));
	struct read_map_inf *out = (struct read_map_inf *)calloc(MAP_READ_BUF_LENGTH,sizeof(struct read_map_inf));
	int out_num = 0;
	int read_num = 0;

	struct contig_t *contig = (struct contig_t *)calloc(SEED_BUF_LENGTH+1,sizeof(struct contig_t));

	int i = 0;

	for(i = 0;i<MAX_COVERAGE;i++)
    {
        read[i].site = (unsigned int *)calloc(MAX_READ_LENGTH,sizeof(unsigned int));
        read[i].node = (unsigned int *)calloc(MAX_READ_LENGTH,sizeof(unsigned int));
    }

	int r = 0;
	unsigned int ref_pos = 0,end_pos = 0;
	int r_root = 0;
	int chr = -1;
	int branch_length = 0;
	int flag = 0;

	while (1)
	{
        pthread_mutex_lock(&ReadLock);
        read_num = 0;
        branch_length = 0;
        r = heap->heap[1];
        ref_pos = 0;
        if(heap->heap_num[1]>0) ref_pos = read_buf[r].ref_site;
        if(ref_pos!=0)
        {
            if((chr==-1)||(ref_pos>=opt->chr->list[chr].start_site+opt->chr->list[chr].length))
            {
                if(chr!=-1)
                {
                    opt->chr->list[chr].thread_num--;
                    if(opt->chr->list[chr].thread_num==0)
                    {
                        for (i=0; i<pow(4,KMER_FRONT); i++)
                        {
                            if(opt->chr->list[chr].c_num[i]!=0) free(opt->chr->list[chr].c_hash[i]);
                        }
                        free(opt->chr->list[chr].c_num);
                    }
                }
                for(i = chr+1;i<opt->chr->total;i++)
                {
                    if((ref_pos>=opt->chr->list[i].start_site)&&(ref_pos<opt->chr->list[i].start_site+opt->chr->list[i].length))break;
                }
                chr = i;
                if(opt->chr->list[chr].c_hash==NULL) load_hash(opt,chr);
                opt->chr->list[chr].thread_num++;
                //fprintf(stdout,"process chr %d\n",chr);
            }
            read_buf[r].site = read[0].site;
            read_buf[r].node = read[0].node;
            memcpy(read,&(read_buf[r]),sizeof(struct read_map_inf));
            read_num++;
            branch_length = max(read[0].front,read[0].back);
            ReHeap_Read(heap->heap,heap->heap_num,heap->order,&heap->heap_length,job->in_file);
        }
        while (heap->heap_num[1]>0)
        {
            r = heap->heap[1];
            if((read_num<MAX_COVERAGE)
            &&(read_buf[r].ref_site<read[read_num-1].ref_site+100)
            &&(read_buf[r].ref_site<ref_pos+400)
            &&(read_buf[r].ref_site<opt->chr->list[chr].start_site+opt->chr->list[chr].length))
            {
                if(strcmp(read[read_num-1].name,read_buf[r].name)!=0)
                {
                    read_buf[r].site = read[read_num].site;
                    read_buf[r].node = read[read_num].node;
                    memcpy(read+read_num,&(read_buf[r]),sizeof(struct read_map_inf));
                    read_num++;
                    branch_length = max(read[read_num-1].front,branch_length);
                    branch_length = max(read[read_num-1].back,branch_length);
                }
                ReHeap_Read(heap->heap,heap->heap_num,heap->order,&heap->heap_length,job->in_file);
            }
            else break;
        }
		pthread_mutex_unlock(&ReadLock);

		if(read_num==0) break;

		memset(map->node,0,map->node_num*sizeof(struct node_table));
		map->node_num = 0;
		if(ref_pos>=99) ref_pos=max(opt->chr->list[chr].start_site,ref_pos-99);
		else ref_pos = 0;
		end_pos = min(read[read_num-1].ref_site+99,opt->chr->list[chr].start_site+opt->chr->list[chr].length);
		r_root = end_pos-ref_pos;
		flag = 0;

        if(build_map(job->opt,map,ref_pos,end_pos,chr,branch_length+opt->change_length,r_root)==1) flag = 1;

        if(flag == 0)
        {
            for(i = 0; i<read_num; i++)
            {
            read[i].chr_order = chr;
            read[i].ref_site-=opt->chr->list[chr].start_site;
            read[i].start_site-=opt->chr->list[chr].start_site;
            memset(read[i].site, 0, read[i].length*sizeof(unsigned int)/sizeof(char));
            memset(read[i].node, 0, read[i].length*sizeof(unsigned int)/sizeof(char));
            }
            align_read(job->opt,opt->chr->list[chr].seq,opt->chr->list[chr].length,chr,read,read_num,contig,opt->chr->list[chr].c_hash,opt->chr->list[chr].c_num,seed);

            for(i = 0; i<read_num; i++)
            {
                if(insert_read2map_front(job->opt,map,read+i,t_stack)==2){flag = 1;break;}
                if(insert_read2map_back(job->opt,map,read+i,t_stack)==2){flag = 1;break;}
            }

            if(flag==0) stat_map(opt,map,r_root,opt->chr->list[chr].seq,opt->chr->list[chr].length,t_stack,seed);
            if(flag==0) output_map_align(opt,read,read_num,map,opt->chr->list[chr].seq,opt->chr->list[chr].length-1,t_stack,r_root);
        }
        else
        {
            for(i = 0; i<read_num; i++)
            {
            read[i].chr_order = chr;
            read[i].ref_site-=opt->chr->list[chr].start_site;
            read[i].start_site-=opt->chr->list[chr].start_site;
            }
        }
        for(i = 0;i<read_num;i++)
        {
            if(out_num>=MAP_READ_BUF_LENGTH)
            {
                pthread_mutex_lock(&UpdateLock);
                qsort(out,MAP_READ_BUF_LENGTH,sizeof(struct read_map_inf),result_cmp);

                FILE *result_file;
                char result_file_name[MAX_NAME_LENGTH];
                sprintf(result_file_name,"%s/%d.result",opt->Temp_path,opt->result_fn);
                result_file = fopen(result_file_name,"wb");

                if(result_file == NULL)
                {
                    fprintf(stdout,"cannot open file %s...\n",result_file_name);
                }

                if(fwrite(out,sizeof(struct read_map_inf),MAP_READ_BUF_LENGTH,result_file)!=MAP_READ_BUF_LENGTH)
                    fprintf(stdout,"[DEEP]cannot write file...\n");

                fclose(result_file);
                fflush(result_file);
                opt->result_fn++;
                out_num = 0;
                pthread_mutex_unlock(&UpdateLock);
            }
            memcpy(&(out[out_num]),&(read[i]),sizeof(struct read_map_inf));
            out_num++;
        }
	}
	if(out_num>0)
    {
        pthread_mutex_lock(&UpdateLock);
        qsort(out,out_num,sizeof(struct read_map_inf),result_cmp);

        FILE *result_file;
        char result_file_name[MAX_NAME_LENGTH];
        sprintf(result_file_name,"%s/%d.result",opt->Temp_path,opt->result_fn);
        result_file = fopen(result_file_name,"wb");

        if(result_file == NULL)
        {
            fprintf(stdout,"cannot open file %s...\n",result_file_name);
        }
        if(fwrite(out,sizeof(struct read_map_inf),out_num,result_file)!=out_num)
            fprintf(stdout,"[DEEP]cannot write file...\n");
        fclose(result_file);
        fflush(result_file);
        opt->result_fn++;

        pthread_mutex_unlock(&UpdateLock);
    }

    free(map->node);
    free(map);

    free(t_stack->data);
    free(t_stack);

    for(i = 0;i<MAX_COVERAGE;i++)
    {
        free(read[i].site);
        free(read[i].node);
    }

	free(seed);
	free(read);
	free(out);

    return (void*)(0);
}
struct job_result
{
    struct m_opt *opt;
    FILE *in_file[MAP_READ_BUF_LENGTH];
    FILE *out;
    FILE *un;
};
void find_pair_map(struct read_map_inf *read,int start,int length,int pstart,int plength)
{
    int i = 0,j = 0;
    int best,best_1;

    if((opt->output_mode==OUTPUT_BEST))
    {
        for(j = pstart;j<pstart+plength;j++)
        {
            read[j].porder = -1;
            read[j].psite = 0;
            if(read[j].output_flag==0)continue;
        }
        for(i = start;i<start+length;i++)
        {
            read[i].porder = -1;
            read[i].psite = 0;
            if(read[i].output_flag==0)continue;

            for(j = pstart;j<pstart+plength;j++)
            {
                if(read[i].chr_order!=read[j].chr_order) continue;
                if(((read[i].end_site>=read[j].start_site-1000)&&(read[i].end_site<read[j].start_site-100))
                    ||((read[i].start_site>read[j].end_site+100)&&(read[i].start_site<=read[j].end_site+1000)))
                {
                    if((read[i].porder==-1)||(read[j].score>read[read[i].porder].score))
                    {
                        read[i].porder = j;
                        read[i].psite = read[j].start_site;
                    }
                    if((read[j].porder==-1)||(read[i].score>read[read[j].porder].score))
                    {
                        read[j].porder = i;
                        read[j].psite = read[i].start_site;
                    }
                }
            }
        }
        for(i = start;i<start+length;i++)
        {
            if(read[i].output_flag==0)continue;
            if(read[i].porder!=-1) continue;

            for(j = pstart;j<pstart+plength;j++)
            {
                if(read[i].chr_order!=read[j].chr_order) continue;

                if((read[i].start_site>=read[j].start_site-opt->area)&&(read[i].start_site<=read[j].start_site+opt->area))
                {
                    if((read[i].porder==-1)||(read[j].score>read[read[i].porder].score))
                    {
                        read[i].porder = j;
                        read[i].psite = read[j].start_site;
                    }
                    if((read[j].porder==-1))
                    {
                        read[j].porder = i;
                        read[j].psite = read[i].start_site;
                    }
                }
            }
        }
        {
            best = -1;
            best_1 = -1;
            for(i = start;i<start+length;i++)
            {
                if(read[i].output_flag==0)continue;
                if(read[i].porder==-1) continue;
                if((read[read[i].porder].score+read[i].score)>best)
                {
                    best_1= i;
                    best = read[read[i].porder].score+read[i].score;
                }
            }
            if(best_1!=-1)
            {

                for(i = start;i<start+length;i++)
                {
                    if(read[i].output_flag==0)continue;
                    if(read[i].porder==-1) read[i].output_flag = 0;
                    else if((read[read[i].porder].score+read[i].score)<best-6)
                    {
                        read[i].output_flag = 0;
                    }
                }
                for(i = pstart;i<pstart+plength;i++)
                {
                    if(read[i].output_flag==0)continue;
                    if(read[i].porder==-1) read[i].output_flag = 0;
                    else if((read[read[i].porder].score+read[i].score)<best-6)
                    {
                        read[i].output_flag = 0;
                    }
                }
            }
            else
            {
                for(i = start;i<start+length;i++)
                {
                    if(read[i].score!=read[start].score) read[i].output_flag=0;
                }
                for(i = pstart;i<pstart+plength;i++)
                {
                    if(read[i].score!=read[pstart].score) read[i].output_flag=0;
                }
            }
        }
    }
    else
    {
        for(j = pstart;j<pstart+plength;j++)
        {
            read[j].porder = -1;
            read[j].psite = 0;
            if(read[j].output_flag==0)continue;
        }
        for(i = start;i<start+length;i++)
        {
            read[i].porder = -1;
            read[i].psite = 0;
            if(read[i].output_flag==0)continue;

            for(j = pstart;j<pstart+plength;j++)
            {
                if(read[i].chr_order!=read[j].chr_order) continue;

                if(((read[i].end_site>=read[j].start_site-1000)&&(read[i].end_site<read[j].start_site))
                    ||((read[i].start_site>read[j].end_site)&&(read[i].start_site<=read[j].end_site+1000)))
                {
                    if((read[i].porder==-1)||(read[j].score>read[read[i].porder].score))
                    {
                        read[i].porder = j;
                        read[i].psite = read[j].start_site;
                    }
                    if((read[j].porder==-1)||(read[i].score>read[read[j].porder].score))
                    {
                        read[j].porder = i;
                        read[j].psite = read[i].start_site;
                    }
                }
            }
        }
    }
}
void *result_core(void* arg)
{
    struct job_result *job = (struct job_result *)arg;

	struct read_map_inf *read = (struct read_map_inf *)calloc(MAP_READ_BUF_LENGTH+100,sizeof(struct read_map_inf));
	int read_num = 0;

	char TAG[MAX_STRING_LENGTH];

	int i,j,r;
	int length,plength;
	int start,pstart;
	int flag;

	while (1)
	{
        pthread_mutex_lock(&ReadLock);
        read_num = 0;
        while (heap->heap_num[1]>0)
        {
            r = heap->heap[1];
            if(read_num<MAP_READ_BUF_LENGTH)
            {
                pstart = read_num;
                memcpy(read+read_num,&read_buf[r],sizeof(struct read_map_inf));
                read_num++;
                ReHeap_Res(heap->heap,heap->heap_num,heap->order,&heap->heap_length,job->in_file);
            }
            else if(strcmp(read_buf[r].name,read[read_num-1].name)==0)
            {
                memcpy(read+read_num,&read_buf[r],sizeof(struct read_map_inf));
                read_num++;
                ReHeap_Res(heap->heap,heap->heap_num,heap->order,&heap->heap_length,job->in_file);
            }
            else if((opt->pair==1)&&(strcmp(read_buf[r].name,read[pstart].pname)==0))
            {
                memcpy(read+read_num,&read_buf[r],sizeof(struct read_map_inf));
                read_num++;
                ReHeap_Res(heap->heap,heap->heap_num,heap->order,&heap->heap_length,job->in_file);
            }
            else break;
        }
		pthread_mutex_unlock(&ReadLock);

		if(read_num==0) break;

		length = 0;
		start = 0;
		pstart = 0;
		plength = 0;
		for(i = 0;i<read_num;i++)
        {
            if(opt->pair==1)
            {
                if(strcmp(read[i].name,read[start].name)==0)
                {
                    if(read[i].cigar[0]=='*') {read[i].output_flag = 0;continue;}
                    {
                        flag = 1;
                        for(j = start;j<i;j++)
                        {
                            if((read[j].output_flag==1)&&(read[i].chr_order==read[j].chr_order)&&(read[i].start_site==read[j].start_site))
                            {
                                if(strcmp(read[i].cigar,read[j].cigar)==0)
                                    read[i].output_flag = 0;
                                else
                                {
                                    if(read[i].score>read[j].score) {read[i].output_flag = 1;read[j].output_flag = 0;}
                                    else if(read[i].score==read[j].score) read[i].output_flag = 1;
                                    else read[i].output_flag = 0;
                                }
                                flag = 0;break;
                            }
                        }
                        if(flag) {read[i].output_flag = 1;length++;}
                    }
                }
                else if(strcmp(read[i].name,read[start].pname)==0)
                {
                    if(pstart==0) pstart = i;
                    if(read[i].cigar[0]=='*') {read[i].output_flag = 0;continue;}
                    {
                        flag = 1;
                        for(j = pstart;j<i;j++)
                        {
                            if((read[j].output_flag==1)&&(read[i].chr_order==read[j].chr_order)&&(read[i].start_site==read[j].start_site))
                            {
                                if(strcmp(read[i].cigar,read[j].cigar)==0)
                                    read[i].output_flag = 0;
                                else
                                {
                                    if(read[i].score>read[j].score) {read[i].output_flag = 1;read[j].output_flag = 0;}
                                    else if(read[i].score==read[j].score) read[i].output_flag = 1;
                                    else read[i].output_flag = 0;
                                }
                                flag = 0;break;
                            }
                        }
                        if(flag) {read[i].output_flag = 1;plength++;}
                    }
                }
                else
                {
                    if(length==0) read[start].output_flag = 2;
                    else if((plength==0)&&(pstart!=0)) {read[pstart].output_flag = 2;}
                    else
                    {
                        if((length!=0)&&(plength!=0))find_pair_map(read,start,pstart-start,pstart,i-pstart);
                        else if(opt->output_mode==OUTPUT_BEST)
                        {
                            if(length!=0)
                            {
                                if (pstart!=0)
                                {
                                    for(j = start+1;j<pstart;j++)
                                    {
                                        if(read[j].score!=read[start].score) read[j].output_flag = 0;
                                    }
                                }
                                else
                                {
                                    for(j = start+1;j<i;j++)
                                    {
                                        if(read[j].score!=read[start].score) read[j].output_flag = 0;
                                    }
                                }
                            }
                            else if(plength!=0)
                            {
                                for(j = pstart+1;j<i;j++)
                                {
                                    if(read[j].score!=read[pstart].score) read[j].output_flag = 0;
                                }
                            }
                        }
                        for(j = start;j<i;j++)
                        {
                            if(read[j].output_flag>0)
                            {
                                write_TAG(read+j);
                                write_MapQ(read+j,length);
                            }
                        }
                    }
                    pstart = 0;plength = 0;
                    start = i;
                    length = 0;
                    if(read[i].cigar[0]!='*'){read[i].output_flag = 1;length++;}
                }
            }
            else
            {
            if(strcmp(read[i].name,read[start].name)==0)
            {
                if(read[i].cigar[0]=='*') {read[i].output_flag = 0;continue;}
                if(opt->output_mode==OUTPUT_BEST)
                {
                    if(read[i].score==read[start].score)
                    {
                        flag = 1;
                        for(j = start;j<i;j++)
                        {
                            if((read[i].chr_order==read[j].chr_order)&&(read[i].start_site==read[j].start_site)&&(strcmp(read[i].cigar,read[j].cigar)==0))
                                {read[i].output_flag = 0;flag = 0;break;}
                        }
                        if(flag) {read[i].output_flag = 1;length++;}
                    }
                    else read[i].output_flag = 0;
                }
                else
                {
                    flag = 1;
                    for(j = start;j<i;j++)
                    {
                        if((read[i].chr_order==read[j].chr_order)&&(read[i].start_site==read[j].start_site))
                        {
                            if(strcmp(read[i].cigar,read[j].cigar)==0)
                                read[i].output_flag = 0;
                            else
                            {
                                if(read[i].score>read[j].score) {read[i].output_flag = 1;read[j].output_flag = 0;}
                                else read[i].output_flag = 0;
                            }
                            flag = 0;break;
                        }
                    }
                    if(flag) {read[i].output_flag = 1;length++;}
                }
            }
            else
            {
                if(length==0) read[start].output_flag = 2;
                else
                {
                    for(j = start;j<i;j++)
                    {
                        if(read[j].output_flag>0)
                        {
                            write_TAG(read+j);
                            write_MapQ(read+j,length);
                        }
                    }
                }
                start = i;
                length = 0;
                if(read[i].cigar[0]!='*'){read[i].output_flag = 1;length++;}
            }
            }
        }
        if(length==0) read[start].output_flag = 2;
        else if((plength==0)&&(pstart!=0)) read[pstart].output_flag = 2;
        else
        {
            if((length!=0)&&(plength!=0)&&(opt->pair==1)) find_pair_map(read,start,pstart-start,pstart,i-pstart);
            else if(opt->pair==1)
            {
                if(length!=0)
                {
                    if (pstart!=0)
                    {
                        for(j = start+1;j<pstart;j++)
                        read[j].output_flag = 0;
                    }
                    else
                    {
                        for(j = start+1;j<i;j++)
                        read[j].output_flag = 0;
                    }

                }
                else if(plength!=0)
                {
                    for(j = pstart+1;j<i;j++)
                        read[j].output_flag = 0;
                }
            }
            for(j = start;j<i;j++)
            {
                if(read[j].output_flag>0)
                {
                    write_TAG(read+j);
                    write_MapQ(read+j,length);
                }
            }
        }

		pthread_mutex_lock(&UpdateLock);
        for(i = 0;i<read_num;i++)
        {
            length = strlen(read[i].name);
            read[i].name[length-1] = '\0';
            if(read[i].pname[0]!='*')
            {
                length = strlen(read[i].pname);
                read[i].pname[length-1] = '\0';
            }
            if(read[i].output_flag==1)
            {
                sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read[i].score,read[i].dis,read[i].TAG);
                fprintf(job->out,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%s\n",
                    read[i].name,read[i].flag,opt->chr->list[read[i].chr_order].name,read[i].start_site+1,read[i].MapQ,read[i].cigar,
                    read[i].pname,(read[i].psite==0)?read[i].psite:read[i].psite,read[i].seq,read[i].qual,TAG);
            }
            else if(read[i].output_flag==2)
            {
                read[i].flag |= 0x0008;
                fprintf(job->un,"%s\t%d\t*\t0\t0\t*\t%s\t%u\t0\t%s\t%s\n",
                    read[i].name,read[i].flag,read[i].pname,(read[i].psite==0)?read[i].psite:read[i].psite,read[i].seq,read[i].qual);
            }
        }
        pthread_mutex_unlock(&UpdateLock);
	}
	free(read);

    return (void*)(0);
}
int map_align(struct m_opt *opt)
{
    int i = 0;

    FILE *InPut_file[MAP_READ_BUF_LENGTH];
    char in_file[MAX_NAME_LENGTH];

    for(i = 0;i<opt->find_fn;i++)
    {
        sprintf(in_file,"%s/%d.read",opt->Temp_path,i);
        InPut_file[i] = fopen(in_file,"rb");
        if (InPut_file[i] == NULL)
        {
            fprintf(stdout,"[DEEP]cannot open file %s...\n",in_file);
            return 1;
        }
        if(fread(&read_buf[i],sizeof(struct read_map_inf),1,InPut_file[i])!=1)
        {
            fprintf(stdout,"[DEEP]cannot read seq...\n");
            return 1;
        }

        heap->heap[i+1] = i;
        heap->heap_num[i+1] = MAP_READ_BUF_LENGTH;
        heap->order[i+1] = i;
    }
    heap->heap_length = opt->find_fn;

    for(i=heap->heap_length/2;i>=1;i--)
        HeapAdjust_Read(heap->heap,heap->heap_num,heap->order,heap->heap_length,i);

    FILE *out;
    char out_file_name_1[MAX_NAME_LENGTH];


    sprintf(out_file_name_1,"%s/map_out.sam",opt->Output_path);
    out = fopen(out_file_name_1,"a+");
    if (out == NULL)
        return 1;

    int r;
    while (heap->heap_num[1]>0)
    {
        r = heap->heap[1];
        if((read_buf[r].ref_site==0)&&(read_buf[r].cigar[0]=='*'))
        {
            fprintf(out,"%s\t%d\t*\t0\t0\t*\t%s\t0\t0\t%s\t%s\n",read_buf[r].name,read_buf[r].flag,read_buf[r].pname,read_buf[r].seq,read_buf[r].qual);
            ReHeap_Read(heap->heap,heap->heap_num,heap->order,&heap->heap_length,InPut_file);
        }
        else break;
    }

    opt->result_fn = 0;

    struct job_snp *job = (struct job_snp *)calloc(1,sizeof(struct job_snp));
    for(i = 0;i<opt->find_fn;i++) job->in_file[i] = InPut_file[i];
    job->opt = opt;

	pthread_t *pthreads = malloc(sizeof(pthread_t) * opt->thread_num);

    for (i = 0; i < opt->thread_num; i++) pthread_create(&pthreads[i], NULL, map_align_core, job);
	for (i = 0; i < opt->thread_num; i++) pthread_join(pthreads[i], NULL);

	for(i = 0;i<opt->find_fn;i++)
    {
        if(InPut_file[i]!= NULL)
        fclose(InPut_file[i]);
        InPut_file[i] = NULL;
        sprintf(in_file,"%s/%d.read",opt->Temp_path,i);
        remove(in_file);
    }
    free(job);
    free(pthreads);

    if(opt->result_fn==0)
    {
        fclose(out);
        return 0;
    }

    for(i = 0;i<opt->result_fn;i++)
    {
        sprintf(in_file,"%s/%d.result",opt->Temp_path,i);
        InPut_file[i] = fopen(in_file,"rb");
        if (InPut_file[i] == NULL)
        {
            fprintf(stdout,"[DEEP]cannot open file %s...\n",in_file);
            return 1;
        }
        if(fread(&read_buf[i],sizeof(struct read_map_inf),1,InPut_file[i])!=1)
        {
            fprintf(stdout,"[DEEP]cannot read seq...\n");
            return 1;
        }

        heap->heap[i+1] = i;
        heap->heap_num[i+1] = MAP_READ_BUF_LENGTH;
        heap->order[i+1] = i;
    }
    heap->heap_length = opt->result_fn;

    for(i=heap->heap_length/2;i>=1;i--)
        HeapAdjust_Res(heap->heap,heap->heap_num,heap->order,heap->heap_length,i);

    struct job_result *jobr = (struct job_result *)calloc(1,sizeof(struct job_result));
    for(i = 0;i<opt->result_fn;i++) jobr->in_file[i] = InPut_file[i];
    jobr->opt = opt;
    jobr->out = out;
    jobr->un = out;

    pthread_t *pthreadsr = malloc(sizeof(pthread_t) * opt->thread_num);
    for (i = 0; i < opt->thread_num; i++) pthread_create(&pthreadsr[i], NULL, result_core, jobr);
	for (i = 0; i < opt->thread_num; i++) pthread_join(pthreadsr[i], NULL);

    for(i = 0;i<opt->result_fn;i++)
    {
        if(InPut_file[i]!= NULL)
        fclose(InPut_file[i]);
        sprintf(in_file,"%s/%d.result",opt->Temp_path,i);
        remove(in_file);
    }
    free(jobr);
    free(pthreadsr);
    fclose(out);
    return 0;
}

