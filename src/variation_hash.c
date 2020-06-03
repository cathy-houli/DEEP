#include "main.h"
static pthread_mutex_t ReadLock,UpdateLock;
struct snp_t *snp_buf;
struct snp_seq_t *snpb_buf;
void HeapAdjust_SNP(unsigned int *heap,unsigned int *heap_num,unsigned int *order,unsigned int heap_length,unsigned int i)
{
    int LC=2*i;
    int RC=2*i+1;
    int min=i;
    if(i<=heap_length/2)
    {
        if((LC<=heap_length)
        &&((snp_buf[heap[LC]].start<snp_buf[heap[min]].start)
           ||((snp_buf[heap[LC]].start==snp_buf[heap[min]].start)&&(snp_buf[heap[LC]].end<snp_buf[heap[min]].end))
           ||((snp_buf[heap[LC]].start==snp_buf[heap[min]].start)&&(snp_buf[heap[LC]].end==snp_buf[heap[min]].end)&&(snp_buf[heap[LC]].type<snp_buf[heap[min]].type))
           ||((snp_buf[heap[LC]].start==snp_buf[heap[min]].start)&&(snp_buf[heap[LC]].end==snp_buf[heap[min]].end)&&(snp_buf[heap[LC]].type==snp_buf[heap[min]].type)&&(snp_buf[heap[LC]].length<snp_buf[heap[min]].length))))
            min=LC;
        if((RC<=heap_length)
        &&((snp_buf[heap[RC]].start<snp_buf[heap[min]].start)
           ||((snp_buf[heap[RC]].start==snp_buf[heap[min]].start)&&(snp_buf[heap[RC]].end<snp_buf[heap[min]].end))
           ||((snp_buf[heap[RC]].start==snp_buf[heap[min]].start)&&(snp_buf[heap[RC]].end==snp_buf[heap[min]].end)&&(snp_buf[heap[RC]].type<snp_buf[heap[min]].type))
           ||((snp_buf[heap[RC]].start==snp_buf[heap[min]].start)&&(snp_buf[heap[RC]].end==snp_buf[heap[min]].end)&&(snp_buf[heap[RC]].type==snp_buf[heap[min]].type)&&(snp_buf[heap[RC]].length<snp_buf[heap[min]].length))))
            min=RC;
        if(min!=i)
        {
            SwapHeap(i,min,heap,heap_num,order);
            HeapAdjust_SNP(heap,heap_num,order,heap_length,min);
        }
    }
}
void ReHeap_SNP(unsigned int *heap,unsigned int *heap_num,unsigned int *order,unsigned int *heap_length,FILE *InPut_file[])
{
    heap_num[1]--;
    if(heap_num[1]>0)
    {
        if(fread(&snp_buf[order[1]],sizeof(struct snp_t),1,InPut_file[order[1]])!=1)
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
    HeapAdjust_SNP(heap,heap_num,order,*heap_length,1);
}
int exon_hash_cmp(const void *a,const void *b)
{
    struct hash_temp *EA,*EB;
    EA = (struct hash_temp *)a;
    EB = (struct hash_temp *)b;
    int i;

    if(EA->front==EB->front)
    {
        if (EA->back==EB->back)
        {
            if(EA->start==EB->start)
            {
                if(EA->end==EB->end)
                {
                    for(i = 0;i<3;i++)
                    {
                        if(EA->c[i] < EB->c[i])return -1;
                        else if(EA->c[i] > EB->c[i]) return 1;

                        if(EA->l[i] < EB->l[i])return -1;
                        else if(EA->l[i] > EB->l[i]) return 1;
                    }
                    return 0;
                }
                else if(EA->end<EB->end) return -1;
                else return 1;
            }
            else if(EA->start<EB->start) return -1;
            else return 1;
        }
        else if (EA->back<EB->back) return -1;
        else  return 1;
    }
    else
    {
        if(EA->front < EB->front) return -1;
        else return 1;
    }
}
int hash_e_find(struct exon_hash **hash_front,unsigned int *hash_num,uint32_t f_base,uint32_t b_base)
{
    int left = 0, right = hash_num[f_base]-1, middle;

    if (right == -1)
        return -1;
    while (left <= right)
    {
        middle = (left + right)/2;

        if ((hash_front[f_base][middle].back) == b_base)
        {
            while((middle-1>=0)&&((hash_front[f_base][middle-1].back) == b_base))
                middle--;
            return middle;
        }
        else if ((hash_front[f_base][middle].back) > b_base)
        {
            right = middle -1;
        }
        else
            left = middle + 1;
    }
    return -1;
}
void build_SNP(struct m_opt *opt)
{
    opt->total = 10*MAX_SNP_NUM;

    snp_buf = (struct snp_t *)calloc(opt->snp_fn,sizeof(struct snp_t));
    opt->snp = (struct snp_list_t *)calloc(1,sizeof(struct snp_list_t));
    opt->snp->snp = (struct snp_t *)calloc(opt->total,sizeof(struct snp_t));
    opt->snp->total = 0;

    int i = 0;
    FILE *InPut_file[MAP_READ_BUF_LENGTH];
    char in_file[MAX_NAME_LENGTH];

    for(i = 0;i<opt->snp_fn;i++)
    {
        sprintf(in_file,"%s/%d.snp",opt->Temp_path,i);
        InPut_file[i] = fopen(in_file,"rb");
        if (InPut_file[i] == NULL)
        {
            fprintf(stdout,"[DEEP]cannot open file %s...\n",in_file);
            continue;
        }
        if(fread(&snp_buf[i],sizeof(struct snp_t),1,InPut_file[i])!=1)
        {
            fprintf(stdout,"[DEEP]cannot read seq...\n");
            continue;
        }

        heap->heap[i+1] = i;
        heap->heap_num[i+1] = MAX_SNP_NUM;
        heap->order[i+1] = i;
    }
    heap->heap_length = opt->snp_fn;
    for(i=heap->heap_length/2;i>=1;i--)
        HeapAdjust_SNP(heap->heap,heap->heap_num,heap->order,heap->heap_length,i);

    int r = heap->heap[1];

    unsigned int exon_order = 0;
    int dep = 0;
    int score = 0;

    memcpy(opt->snp->snp+opt->snp->total,&snp_buf[r],sizeof(struct snp_t));
    opt->snp->total++;
    ReHeap_SNP(heap->heap,heap->heap_num,heap->order,&heap->heap_length,InPut_file);

    while (heap->heap_num[1]>0)
    {
        r = heap->heap[1];
        if((snp_buf[r].start==opt->snp->snp[opt->snp->total-1].start)
           &&(snp_buf[r].end==opt->snp->snp[opt->snp->total-1].end)
           &&(snp_buf[r].type==opt->snp->snp[opt->snp->total-1].type)
           &&(snp_buf[r].length==opt->snp->snp[opt->snp->total-1].length))
            opt->snp->snp[opt->snp->total-1].num+=snp_buf[r].num;
        else
        {
            while(opt->snp->snp[opt->snp->total-1].start>opt->exon->exon[exon_order].end) exon_order++;
            if(opt->snp->snp[opt->snp->total-1].start<opt->exon->exon[exon_order].start)
            {
                while(opt->snp->snp[opt->snp->total-1].start-1<opt->exon->exon[exon_order].start) exon_order--;
                dep = opt->exon->exon[exon_order].dep[opt->snp->snp[opt->snp->total-1].start-1-opt->exon->exon[exon_order].start];
            }
            else
                dep = opt->exon->exon[exon_order].dep[opt->snp->snp[opt->snp->total-1].start-opt->exon->exon[exon_order].start];
            opt->snp->snp[opt->snp->total-1].dep = dep;

            if(opt->snp->snp[opt->snp->total-1].type!='N')
            {
                //if(opt->snp->snp[opt->snp->total-1].num<4) opt->snp->total--;
                score =1-pow((1-(float)opt->snp->snp[opt->snp->total-1].num/(dep+SCORE_BALANCE)),opt->snp->snp[opt->snp->total-1].num);
                if((score<0.9)||(opt->snp->snp[opt->snp->total-1].num<4))opt->snp->total--;
            }

            memcpy(opt->snp->snp+opt->snp->total,&snp_buf[r],sizeof(struct snp_t));
            opt->snp->total++;
            if(opt->snp->total>=opt->total)
            {
                opt->total+=MAX_SNP_NUM;
                opt->snp->snp = (struct snp_t *)realloc(opt->snp->snp,opt->total*sizeof(struct snp_t));
            }
        }
        ReHeap_SNP(heap->heap,heap->heap_num,heap->order,&heap->heap_length,InPut_file);
    }
    for(i = 0;i<opt->snp_fn;i++)
    {
        fclose(InPut_file[i]);
        sprintf(in_file,"%s/%d.snp",opt->Temp_path,i);
        remove(in_file);
    }

    free(snp_buf);
}
void insert_temp_hash(struct m_opt *opt,struct hash_temp *temp,uint64_t *temp_num,int *temp_fn)
{
    qsort(temp,*temp_num,sizeof(struct hash_temp),exon_hash_cmp);

    FILE *hash_file;
    char hash_file_name[MAX_NAME_LENGTH];
    sprintf(hash_file_name,"%s/%d.hash",opt->Temp_path,*temp_fn);
    hash_file = fopen(hash_file_name,"wb");

    if(hash_file == NULL)
    {
        fprintf(stdout,"[DEEP]cannot open file %s...\n",hash_file_name);
    }
    unsigned int order = 0;
    unsigned int i = 0;
    for(i = 0;i<*temp_num;i++)
    {
        if(exon_hash_cmp(&(temp[order]),&(temp[i]))!=0)
        {
            if(fwrite(&temp[order],sizeof(struct hash_temp),1,hash_file)!=1)
            fprintf(stdout,"[DEEP]cannot write file...\n");
            order = i;
        }
    }
    if(fwrite(&temp[order],sizeof(struct hash_temp),1,hash_file)!=1)
        fprintf(stdout,"[DEEP]cannot write file...\n");

    fclose(hash_file);
    (*temp_fn)++;
    *temp_num = 0;
}
struct job_snp
{
    struct m_opt *opt;
    FILE *in_file[MAP_READ_BUF_LENGTH];
    struct hash_temp *temp;
    int *temp_fn;
};
uint64_t TOTAL;
#define MAX_SNP_HASH_NUM 5000
#define MAX_SNP_TEMP 10000000
void *snp_hash(void* arg)
{
    struct job_snp *job = (struct job_snp *)arg;

    unsigned int snp_order = 0;
    unsigned int snp_length = 0;

    struct hash_temp *snpb_temp = (struct hash_temp *)calloc(MAX_SNP_TEMP,sizeof(struct hash_temp));
    unsigned int temp_num = 0;

	int i,j,k,x;
	int r;

	char seq[MAX_READ_LENGTH];
	unsigned int site[MAX_READ_LENGTH];
	char flag[MAX_READ_LENGTH];

	uint32_t front = 0,back = 0;
	unsigned int Fmark = (1<<(opt->hash_front<<1))-1;
    unsigned int Bmark = (1<<(opt->hash_back<<1))-1;
	unsigned int ref_pos;
	int read_pos;
	int end_pos;
	int length;

	int order;

	while (1)
	{
        //read_file;
        pthread_mutex_lock(&ReadLock);
        snp_order = opt->snp_num;
        snp_length = min(opt->snp->total,snp_order+MAX_SNP_HASH_NUM)-snp_order;
        opt->snp_num+=snp_length;
		pthread_mutex_unlock(&ReadLock);

		if(snp_length==0)break;

        for (i = snp_order;i<snp_order+snp_length;i++)
        {
            {
                read_pos = 0;
                length = opt->hash_front+opt->hash_back-1;
                ref_pos = (opt->snp->snp[i].start>length)?(opt->snp->snp[i].start-length):0;
                if(opt->snp->snp[i].type=='I') length++;

                for(k = 0;k<length;k++)
                {
                    seq[read_pos+k] = (opt->idx->pac[(ref_pos+k) >> 2] >> ((~(ref_pos+k) & 3) << 1) & 3);
                    site[read_pos+k] = ref_pos+k;
                    flag[read_pos+k] = length-k;
                }
                ref_pos += length;
                read_pos += length;

                if(opt->snp->snp[i].type=='I')
                {
                    length = opt->snp->snp[i].length;
                    for(k = 0;k<length;k++)
                    {
                        seq[read_pos+k] = opt->snp->snp[i].seq>>((length-k-1)<<1)&3;
                        site[read_pos+k] = ref_pos-1;
                        flag[read_pos+k] = 'I';
                    }
                    read_pos+=length;
                }
                else if (opt->snp->snp[i].type=='X')
                {
                    length = opt->snp->snp[i].length;
                    for(k = 0;k<length;k++)
                    {
                        seq[read_pos+k] = opt->snp->snp[i].seq>>((length-k-1)<<1)&3;
                        site[read_pos+k] = ref_pos+k;
                        flag[read_pos+k] = 'X';
                    }
                    ref_pos+=length;
                    read_pos+=length;
                }
                else ref_pos+=opt->snp->snp[i].length;
            }
            end_pos = read_pos;
            length = opt->hash_front+opt->hash_back-1;
            for(k = 0;k<length;k++)
            {
                seq[read_pos+k] = (opt->idx->pac[(ref_pos+k) >> 2] >> ((~(ref_pos+k) & 3) << 1) & 3);
                site[read_pos+k] = ref_pos+k;
                flag[read_pos+k] = length-k;
            }
            read_pos+=length;

            front = 0;back = 0;
            for(k = 0;k<opt->hash_front-1;k++)
                front = (front<<2) | seq[k];
            for(k = 0;k<opt->hash_back-1;k++)
                back = (back<<2) | seq[k+opt->hash_front];

            for(k = 0;k<read_pos-opt->hash_front-opt->hash_back;k++)
            {
                front = (front<<2) | seq[k+opt->hash_front-1];
                front = front&Fmark;
                back = (back<<2) | seq[k+opt->hash_front+opt->hash_back-1];
                back = back&Bmark;

                snpb_temp[temp_num].front = front;
                snpb_temp[temp_num].back = back;
                snpb_temp[temp_num].start = site[k];
                snpb_temp[temp_num].end = site[k+opt->hash_front+opt->hash_back-1];
                snpb_temp[temp_num].c[0] = 0;
                snpb_temp[temp_num].c[1] = 0;
                snpb_temp[temp_num].l[0] = 0;
                snpb_temp[temp_num].l[1] = 0;

                x = 0;
                length = 0;
                for(j = k;j<k+opt->hash_front+opt->hash_back;j++)
                {
                    if(flag[j] == 'I')
                    {
                        snpb_temp[temp_num].c[x] = 0;
                        snpb_temp[temp_num].c[x] = snpb_temp[temp_num].c[x]<<6 | length;
                        snpb_temp[temp_num].l[x] = 0;
                        while(flag[j]=='I') {snpb_temp[temp_num].l[x]++;j++;}
                        j--;
                        x++;
                    }
                    else if(flag[j] == 'X')
                    {
                        snpb_temp[temp_num].c[x] = 1;
                        snpb_temp[temp_num].c[x] = snpb_temp[temp_num].c[x]<<6 | length;
                        snpb_temp[temp_num].l[x] = 0;
                        while(flag[j]=='X') {snpb_temp[temp_num].l[x]++;j++;}
                        j--;
                        x++;
                    }
                    else
                    {
                        if((j>0)&&(flag[j-1] == 1))
                        {
                            if(site[j]-site[j-1]<opt->change_length)
                            {
                                snpb_temp[temp_num].c[x] = 2;
                                snpb_temp[temp_num].c[x] = snpb_temp[temp_num].c[x]<<6 | length;
                                snpb_temp[temp_num].l[x] = site[j]-site[j-1]-1;
                                x++;
                            }
                            else
                            {
                                snpb_temp[temp_num].c[x] = 3;
                                snpb_temp[temp_num].c[x] = snpb_temp[temp_num].c[x]<<6 | length;
                                snpb_temp[temp_num].l[x] = 0;
                                snpb_temp[temp_num].end = site[j];
                                x++;
                            }
                        }
                        length=flag[j];
                        j+=flag[j];
                        j--;
                    }
                }
                temp_num++;

                if(temp_num>=MAX_SNP_TEMP)
                {
                    pthread_mutex_lock(&UpdateLock);
                    for(j = 0;j<MAX_SNP_TEMP;j++)
                    {
                        if(TOTAL>=104857600)
                            insert_temp_hash(job->opt,job->temp,&TOTAL,job->temp_fn);
                        memcpy(job->temp+TOTAL,&(snpb_temp[j]),sizeof(struct hash_temp));
                        TOTAL++;
                    }
                    pthread_mutex_unlock(&UpdateLock);
                    temp_num = 0;
                }
                k++;
                front = (front<<2) | seq[k+opt->hash_front-1];
                front = front&Fmark;
                back = (back<<2) | seq[k+opt->hash_front+opt->hash_back-1];
                back = back&Bmark;
                k++;
                front = (front<<2) | seq[k+opt->hash_front-1];
                front = front&Fmark;
                back = (back<<2) | seq[k+opt->hash_front+opt->hash_back-1];
                back = back&Bmark;
            }

        //double snp
            order = find_s(opt->snp,opt->snp->snp[i].end+1);
            length = opt->hash_front+opt->hash_back-1;
            for(r = order;r<opt->snp->total;r++)
            {
                if((opt->snp->snp[i].type=='N')&&(opt->snp->snp[r].type=='N')) continue;
                if(opt->snp->snp[i].start==opt->snp->snp[r].start) continue;

                if((opt->snp->snp[r].start-opt->snp->snp[i].end)>length) break;
                else
                {
                    read_pos = end_pos;
            {
                length = opt->snp->snp[r].start-opt->snp->snp[i].end-1;
                ref_pos = opt->snp->snp[r].start-length;
                if(opt->snp->snp[r].type=='I') length++;

                for(k = 0;k<length;k++)
                {
                    seq[read_pos+k] = (opt->idx->pac[(ref_pos+k) >> 2] >> ((~(ref_pos+k) & 3) << 1) & 3);
                    site[read_pos+k] = ref_pos+k;
                    flag[read_pos+k] = length-k;
                }
                ref_pos += length;
                read_pos += length;

                if(opt->snp->snp[r].type=='I')
                {
                    length = opt->snp->snp[r].length;
                    for(k = 0;k<length;k++)
                    {
                        seq[read_pos+k] = opt->snp->snp[r].seq>>((length-k-1)<<1)&3;
                        site[read_pos+k] = ref_pos-1;
                        flag[read_pos+k] = 'I';
                    }
                    read_pos+=length;
                }
                else if (opt->snp->snp[r].type=='X')
                {
                    length = opt->snp->snp[r].length;
                    for(k = 0;k<length;k++)
                    {
                        seq[read_pos+k] = opt->snp->snp[r].seq>>((length-k-1)<<1)&3;
                        site[read_pos+k] = ref_pos+k;
                        flag[read_pos+k] = 'X';
                    }
                    ref_pos+=length;
                    read_pos+=length;
                }
                else ref_pos+=opt->snp->snp[r].length;
            }
            length = opt->hash_front+opt->hash_back-1;
            for(k = 0;k<length;k++)
            {
                seq[read_pos+k] = (opt->idx->pac[(ref_pos+k) >> 2] >> ((~(ref_pos+k) & 3) << 1) & 3);
                site[read_pos+k] = ref_pos+k;
                flag[read_pos+k] = length-k;
            }
            read_pos+=length;

            front = 0;back = 0;
            for(k = 0;k<opt->hash_front-1;k++)
                front = (front<<2) | seq[k];
            for(k = 0;k<opt->hash_back-1;k++)
                back = (back<<2) | seq[k+opt->hash_front];

            for(k = 0;k<read_pos-opt->hash_front-opt->hash_back;k++)
            {
                front = (front<<2) | seq[k+opt->hash_front-1];
                front = front&Fmark;
                back = (back<<2) | seq[k+opt->hash_front+opt->hash_back-1];
                back = back&Bmark;

                if(k<opt->snp->snp[r].start-site[0]-opt->hash_front-opt->hash_back+1) continue;

                snpb_temp[temp_num].front = front;
                snpb_temp[temp_num].back = back;
                snpb_temp[temp_num].start = site[k];
                snpb_temp[temp_num].end = site[k+opt->hash_front+opt->hash_back-1];
                snpb_temp[temp_num].c[0] = 0;
                snpb_temp[temp_num].c[1] = 0;
                snpb_temp[temp_num].l[0] = 0;
                snpb_temp[temp_num].l[1] = 0;

                x = 0;
                length = 0;
                for(j = k;j<k+opt->hash_front+opt->hash_back;j++)
                {
                    if(flag[j] == 'I')
                    {
                        snpb_temp[temp_num].c[x] = 0;
                        snpb_temp[temp_num].c[x] = snpb_temp[temp_num].c[x]<<6 | length;
                        snpb_temp[temp_num].l[x] = 0;
                        while(flag[j]=='I') {snpb_temp[temp_num].l[x]++;j++;}
                        j--;
                        x++;
                        length = 0;
                    }
                    else if(flag[j] == 'X')
                    {
                        snpb_temp[temp_num].c[x] = 1;
                        snpb_temp[temp_num].c[x] = snpb_temp[temp_num].c[x]<<6 | length;
                        snpb_temp[temp_num].l[x] = 0;
                        while(flag[j]=='X') {snpb_temp[temp_num].l[x]++;j++;}
                        j--;
                        x++;
                        length = 0;
                    }
                    else
                    {
                        if((j>0)&&(flag[j-1] == 1))
                        {
                            if(site[j]-site[j-1]<opt->change_length)
                            {
                                snpb_temp[temp_num].c[x] = 2;
                                snpb_temp[temp_num].c[x] = snpb_temp[temp_num].c[x]<<6 | length;
                                snpb_temp[temp_num].l[x] = site[j]-site[j-1]-1;
                                x++;
                                length = 0;
                            }
                            else
                            {
                                snpb_temp[temp_num].c[x] = 3;
                                snpb_temp[temp_num].c[x] = snpb_temp[temp_num].c[x]<<6 | length;
                                snpb_temp[temp_num].l[x] = 0;
                                snpb_temp[temp_num].end = site[j];
                                x++;
                                length = 0;
                            }
                        }
                        length=flag[j];
                        j+=flag[j];
                        j--;
                    }
                }
                temp_num++;

                if(temp_num>=MAX_SNP_TEMP)
                {
                    pthread_mutex_lock(&UpdateLock);
                    for(j = 0;j<MAX_SNP_TEMP;j++)
                    {
                        if(TOTAL>=104857600)
                            insert_temp_hash(job->opt,job->temp,&TOTAL,job->temp_fn);
                        memcpy(job->temp+TOTAL,&(snpb_temp[j]),sizeof(struct hash_temp));
                        TOTAL++;
                    }
                    pthread_mutex_unlock(&UpdateLock);
                    temp_num = 0;
                }
                k++;
                front = (front<<2) | seq[k+opt->hash_front-1];
                front = front&Fmark;
                back = (back<<2) | seq[k+opt->hash_front+opt->hash_back-1];
                back = back&Bmark;
                k++;
                front = (front<<2) | seq[k+opt->hash_front-1];
                front = front&Fmark;
                back = (back<<2) | seq[k+opt->hash_front+opt->hash_back-1];
                back = back&Bmark;
            }
                }
            }
        }
	}
    free(snpb_temp);

    return (void*)(0);
}
#define MAX_HASH_SAME 100
int build_hash(struct m_opt *opt)
{
    unsigned int i = 0,j = 0;
    int k = 0,x = 0;
    uint32_t front = 0,back = 0;
    unsigned int LB,RB;
    unsigned int site_t = 0;
    int chr_order = 0;

    struct hash_temp *temp = (struct hash_temp *)calloc((100*1024*1024ll),sizeof(struct hash_temp));
    if(temp==NULL) return 1;
    uint64_t temp_num = 0;
    int temp_fn = 0;

    unsigned int start = 0;

    int hash_seed_length = opt->hash_front+opt->hash_back;

    unsigned int Fmark = (1<<(opt->hash_front<<1))-1;
    unsigned int Bmark = (1<<(opt->hash_back<<1))-1;

    opt->exon = (struct exon_array *)calloc(1,sizeof(struct exon_array));
    opt->total = MAX_EXON_NUM;
    opt->exon->exon = (struct exon_inf_t *)calloc(MAX_EXON_NUM,sizeof(struct exon_inf_t));
    opt->exon->total = 0;

    start = 0;
    for(i = 0;i<opt->idx->bns->l_pac;i++)
    {
        if(opt->dep[i]!=0)
        {
            if(start == 0) start = i;
        }
        else
        {
            if(start!=0)
            {
                opt->exon->exon[opt->exon->total].start = start;
                opt->exon->exon[opt->exon->total].end = i-1;
                opt->exon->exon[opt->exon->total].dep = (uint16_t *)calloc(i-start,sizeof(uint16_t));

                memcpy(opt->exon->exon[opt->exon->total].dep,&(opt->dep[start]),(i-start)*sizeof(uint16_t));
                opt->exon->total++;
                if(opt->exon->total>=opt->total)
                {
                    opt->total+=MAX_EXON_NUM;
                    opt->exon->exon = (struct exon_inf_t *)realloc(opt->exon->exon,opt->total*sizeof(struct exon_inf_t));
                }

                while(start >= opt->chr->list[chr_order].start_site) chr_order++;
                chr_order--;
                LB = max(opt->chr->list[chr_order].start_site,(start>20)?start-6:0);
                RB = min(opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length-1,i+6);

                site_t = LB;
                for(k = 0;k<opt->hash_front-1;k++)
                    front = (front<<2) | (opt->idx->pac[(site_t+k) >> 2] >> ((~(site_t+k) & 3) << 1) & 3);
                for(k = 0;k<opt->hash_back-1;k++)
                    back = (back<<2) | (opt->idx->pac[(site_t+opt->hash_front+k) >> 2] >> ((~(site_t+opt->hash_front+k) & 3) << 1) & 3);

                for(j = LB;j<RB-hash_seed_length;j++)
                {
                    front = (front<<2) | (opt->idx->pac[(j+opt->hash_front-1) >> 2] >> ((~(j+opt->hash_front-1) & 3) << 1) & 3);
                    front = front&Fmark;
                    back = (back<<2) | (opt->idx->pac[(j+hash_seed_length-1) >> 2] >> ((~(j+hash_seed_length-1) & 3) << 1) & 3);
                    back = back&Bmark;

                    temp[temp_num].front = front;
                    temp[temp_num].back = back;
                    temp[temp_num].start = j;
                    temp[temp_num].end = j+hash_seed_length-1;
                    temp[temp_num].c[0] = 0;
                    temp[temp_num].l[0] = 0;
                    temp[temp_num].c[1] = 0;
                    temp[temp_num].l[1] = 0;
                    temp_num++;
                    if(temp_num>=104857600)
                        insert_temp_hash(opt,temp,&temp_num,&temp_fn);

                    j++;
                    front = (front<<2) | (opt->idx->pac[(j+opt->hash_front-1) >> 2] >> ((~(j+opt->hash_front-1) & 3) << 1) & 3);
                    front = front&Fmark;
                    back = (back<<2) | (opt->idx->pac[(j+hash_seed_length-1) >> 2] >> ((~(j+hash_seed_length-1) & 3) << 1) & 3);
                    back = back&Bmark;

                    j++;
                    front = (front<<2) | (opt->idx->pac[(j+opt->hash_front-1) >> 2] >> ((~(j+opt->hash_front-1) & 3) << 1) & 3);
                    front = front&Fmark;
                    back = (back<<2) | (opt->idx->pac[(j+hash_seed_length-1) >> 2] >> ((~(j+hash_seed_length-1) & 3) << 1) & 3);
                    back = back&Bmark;
                }
                start = 0;
            }
        }
    }
    free(opt->dep);

    build_SNP(opt);

    struct job_snp *job = (struct job_snp *)calloc(1,sizeof(struct job_snp));
    job->opt = opt;
    job->temp = temp;
    job->temp_fn = &temp_fn;
    TOTAL = temp_num;

	pthread_t *pthreads = malloc(sizeof(pthread_t) * opt->thread_num);

    for (i = 0; i < opt->thread_num; i++) pthread_create(&pthreads[i], NULL, snp_hash, job);
	for (i = 0; i < opt->thread_num; i++) pthread_join(pthreads[i], NULL);

	temp_num = TOTAL;

    insert_temp_hash(opt,temp,&temp_num,&temp_fn);
    struct hash_temp *temp_buf = (struct hash_temp *)calloc(temp_fn,sizeof(struct hash_temp));
    char buf_flag[temp_fn];

    FILE *InPut_file[MAP_READ_BUF_LENGTH];
    char in_file[MAX_NAME_LENGTH];

    for(i = 0;i<temp_fn;i++)
    {
        sprintf(in_file,"%s/%d.hash",opt->Temp_path,i);
        InPut_file[i] = fopen(in_file,"rb");
        if (InPut_file[i] == NULL)
        {
            fprintf(stdout,"[DEEP]cannot open file %s...\n",in_file);
            continue;
        }
        if(fread(&temp_buf[i],sizeof(struct hash_temp),1,InPut_file[i])!=1)
        {
            fprintf(stdout,"[DEEP]cannot read seq...\n");
            continue;
        }
        buf_flag[i] = 1;
    }

    int temp_order = 0;
    x = 0;
    front = -1;
    back = -1;
    temp_num = 0;
    int buf_length = temp_fn;
    int back_num;
    int flag = 1;
    while(buf_length>0)
    {
        for(x = 0;x<temp_fn;x++)
        {
            if ((buf_flag[x]==1)&&(exon_hash_cmp(&(temp_buf[x]),&(temp_buf[temp_order]))<0))
                temp_order = x;
        }
        if(temp_buf[temp_order].front ==front)
        {
            if(back!=temp_buf[temp_order].back)
            {
                back = temp_buf[temp_order].back;
                back_num = 0;
                flag = 1;
            }
            else if((flag)&&(back_num>=MAX_HASH_SAME)) {flag = 0;temp_num-=back_num;}

            if(flag)
            {
                if((temp_num==0)||(exon_hash_cmp(&(temp[temp_num-1]),&(temp_buf[temp_order]))!=0))
                {
                    memcpy(temp+temp_num,temp_buf+temp_order,sizeof(struct hash_temp));
                    temp_num++;
                    back_num++;
                }
            }
        }
        else
        {
            if(temp_num!=0)
            {
                opt->e_num[front]=temp_num;
                opt->e_hash[front] = (struct exon_hash *)calloc(opt->e_num[front],sizeof(struct exon_hash));
                if(opt->e_hash[front]==NULL)
                    return 1;
                x = 0;
                while(x<temp_num)
                {
                    opt->e_hash[front][x].back = temp[x].back;
                    opt->e_hash[front][x].start = temp[x].start;
                    opt->e_hash[front][x].end = temp[x].end;
                    opt->e_hash[front][x].c[0] = temp[x].c[0];
                    opt->e_hash[front][x].l[0] = temp[x].l[0];
                    opt->e_hash[front][x].c[1] = temp[x].c[1];
                    opt->e_hash[front][x].l[1] = temp[x].l[1];
                    x++;
                }
                temp_num = 0;
            }

            memcpy(temp,temp_buf+temp_order,sizeof(struct hash_temp));
            temp_num++;
            front = temp_buf[temp_order].front;
            back = temp_buf[temp_order].back;
            back_num = 1;
            flag = 1;
        }

        if(fread(&temp_buf[temp_order],sizeof(struct hash_temp),1,InPut_file[temp_order])!=1)
        {
            buf_length--;
            buf_flag[temp_order] = 0;
            for(i = 0;i<temp_fn;i++)
            {
                if(buf_flag[i]==1){temp_order = i;break;}
            }
        }
    }
    opt->e_num[front]=temp_num;
    opt->e_hash[front] = (struct exon_hash *)calloc(opt->e_num[front],sizeof(struct exon_hash));
    if(opt->e_hash[front]==NULL) return 1;
    x = 0;
    while(x<temp_num)
    {
        opt->e_hash[front][x].back = temp[x].back;
        opt->e_hash[front][x].start = temp[x].start;
        opt->e_hash[front][x].end = temp[x].end;
        opt->e_hash[front][x].c[0] = temp[x].c[0];
        opt->e_hash[front][x].l[0] = temp[x].l[0];
        opt->e_hash[front][x].c[1] = temp[x].c[1];
        opt->e_hash[front][x].l[1] = temp[x].l[1];
        x++;
    }

    free(temp);
    free(temp_buf);
    for(i = 0;i<temp_fn;i++)
    {
        fclose(InPut_file[i]);
        sprintf(in_file,"%s/%d.hash",opt->Temp_path,i);
        remove(in_file);
    }
    return 0;
}

