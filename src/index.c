#include "main.h"

#define MAX_CHR_HASH_SAME 100000
#define BUCKET_LENGTH (1<<20)
#define MAX_INDEX_LENGTH 300*1024*1024
#define CHR_MARK ((1LL<<((KMER_FRONT)<<1))-1)

struct sort_chr
{
    uint16_t mark;
    unsigned int site;
};
struct chr_sort_hash
{
    struct sort_chr *buf;
    unsigned int Enum;
    unsigned int Bnum;
};
struct chr_sort_hash *sort_array;

void SwapHeap(unsigned int a,unsigned int b,unsigned int *heap,unsigned int *heap_num,unsigned int *order)
{
    unsigned int temp = 0;
    temp = heap[a];
    heap[a] = heap[b];
    heap[b] = temp;

    temp = heap_num[a];
    heap_num[a] = heap_num[b];
    heap_num[b] = temp;

    temp = order[a];
    order[a] = order[b];
    order[b] = temp;
}
void HeapAdjust_Chr(unsigned int *heap,unsigned int *heap_num,unsigned int *order,unsigned int heap_length,unsigned int i)
{
    int LC=2*i;
    int RC=2*i+1; 
    int min=i;
    if(i<=heap_length/2)
    {
        if((LC<=heap_length)
        &&((sort_array->buf[heap[LC]].mark<sort_array->buf[heap[min]].mark)||(sort_array->buf[heap[LC]].mark==sort_array->buf[heap[min]].mark&&sort_array->buf[heap[LC]].site<sort_array->buf[heap[min]].site)))
            min=LC;
        if((RC<=heap_length)
        &&((sort_array->buf[heap[RC]].mark<sort_array->buf[heap[min]].mark)||(sort_array->buf[heap[RC]].mark==sort_array->buf[heap[min]].mark&&sort_array->buf[heap[RC]].site<sort_array->buf[heap[min]].site)))
            min=RC;
        if(min!=i)
        {
            SwapHeap(i,min,heap,heap_num,order);
            HeapAdjust_Chr(heap,heap_num,order,heap_length,min);
        }
    }
}
void BuildHeap_Chr(unsigned int *sort_buf,unsigned int *heap,unsigned int *heap_num,unsigned int *order,unsigned int *heap_length,unsigned int Bnum,unsigned int Enum)
{
    unsigned int i;
    for(i=1;i<=Bnum;i++)
    {
        heap[i] = sort_buf[i*BUCKET_LENGTH-1];
        heap_num[i] = BUCKET_LENGTH;
        order[i] = i-1;
    }
    if(Enum!=0)
    {
        (*heap_length) = Bnum+1;
        heap[i] = sort_buf[Bnum*BUCKET_LENGTH+Enum-1];
        heap_num[i] = Enum;
        order[i] = Bnum;
    }
    else (*heap_length) = Bnum;

    for(i=(*heap_length)/2;i>=1;i--)
    {
        HeapAdjust_Chr(heap,heap_num,order,(*heap_length),i);
    }
}
unsigned int OutputHeap_Chr(unsigned int *sort_buf,unsigned int *heap,unsigned int *heap_num,unsigned int *order,unsigned int *heap_length)
{
    unsigned int result = 0;
    if((heap_num[1]>0)&&(*heap_length>0))
        result = heap[1];
    else
        return -1;

    heap_num[1]--;
    if(heap_num[1]>0)
        heap[1] = sort_buf[order[1]*BUCKET_LENGTH+heap_num[1]-1];
    else{
        SwapHeap(1,*heap_length,heap,heap_num,order);
        (*heap_length)--;
    }
    HeapAdjust_Chr(heap,heap_num,order,*heap_length,1);

    return result;
}

int chr_hash_cmp(const void *a,const void *b)
{
    struct sort_chr EA,EB;
    EA = sort_array->buf[(*(unsigned int *)a)];
    EB = sort_array->buf[(*(unsigned int *)b)];

    if(EA.mark==EB.mark)
    {
        if(EA.site==EB.site) return 0;
        if(EA.site>EB.site) return -1;
        if(EA.site<EB.site) return 1;
    }
    else if(EA.mark>EB.mark) return -1;
    else if(EA.mark<EB.mark) return 1;
    return 0;
}
KSEQ_INIT(gzFile, gzread)

int index_main(int argc, char *argv[])
{
    if(argc<3)return usage();

    time_t present_time;
    struct tm *present_tm;

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"build hash index...%s\n",asctime(present_tm));

    char system_order[MAX_STRING_LENGTH];

    if((access(argv[1],0))== -1)
    {
        fprintf(stdout,"cannot access %s, please check your input.\n",argv[1]);
        return 1;
    }
    if((access(argv[2],0))== -1)
    {
        fprintf(stdout,"cannot access %s, create output files folder %s...\n",argv[2],argv[2]);
        sprintf(system_order,"mkdir %s",argv[2]);
        system(system_order);

        if((access(argv[2],0))== -1)
        {
            fprintf(stdout,"error!!! cannot create output files folder...\n");
            return 1;
        }
    }

    strcpy(system_order,"ls ");
    strcat(system_order,argv[1]);
    strcat(system_order,"/*.fa >");
    strcat(system_order,argv[2]);
    strcat(system_order,"/ref_name.txt");
    system(system_order);

    FILE *Chr_info;
    FILE *Name_info;
    FILE *Output_hash;

    gzFile ref;
    kseq_t *q;

    char chr_file[MAX_NAME_LENGTH];
    char name_file[MAX_NAME_LENGTH];
    char in_file[MAX_NAME_LENGTH];
    char out_file[MAX_NAME_LENGTH];

    strcpy(chr_file,argv[2]);
    strcat(chr_file,"/chr_list.txt");

    strcpy(name_file,argv[2]);
    strcat(name_file,"/ref_name.txt");

    Chr_info = fopen(chr_file,"w");
    if (Chr_info == NULL)
        return 1;

    Name_info = fopen(name_file,"r");
    if (Name_info == NULL)
        return 1;

    strcpy(system_order,"sort -k1,1 ");
    strcat(system_order,name_file);
    strcat(system_order," -o ");
    strcat(system_order,name_file);
    system(system_order);

    char f_line[MAX_STRING_LENGTH];
    uint16_t base = 0;
    uint32_t temp_base = -1;
    unsigned int i = 0;
    unsigned int rel_site = 0;
    unsigned int read_site = 0;
    int j = 0,move = 0;
    int N = 0;

    unsigned int *buffer=(unsigned int *)calloc(MAX_CHR_HASH_SAME,sizeof(unsigned int));
    int buffer_num = 0;
    int buffer_flag = 0;

    int *num;
    num= (int *)calloc(pow(4,KMER_FRONT),sizeof(int));

    sort_array = (struct chr_sort_hash *)calloc(1,sizeof(struct chr_sort_hash));
    sort_array->buf = (struct sort_chr *)calloc(MAX_INDEX_LENGTH,sizeof(struct sort_chr));
    sort_array->Enum = 0;
    sort_array->Bnum = 0;

    heap = (struct heap_array *)calloc(1,sizeof(struct heap_array));
    heap->sort_buf = (unsigned int *)calloc(MAX_INDEX_LENGTH,sizeof(unsigned int));
    heap->heap = (unsigned int *)calloc(MAX_INDEX_LENGTH/BUCKET_LENGTH+1,sizeof(unsigned int));
    heap->heap_num = (unsigned int *)calloc(MAX_INDEX_LENGTH/BUCKET_LENGTH+1,sizeof(unsigned int));
    heap->order = (unsigned int *)calloc(MAX_INDEX_LENGTH/BUCKET_LENGTH+1,sizeof(unsigned int));
    heap->heap_length = 0;

    while (fgets(f_line,MAX_STRING_LENGTH,Name_info)!=NULL)
    {
        time(&present_time);
        present_tm = localtime(&present_time);

        sscanf(f_line,"%s",in_file);
        ref = gzopen(in_file,"r");
        q = kseq_init(ref);
        if (kseq_read(q) < 0)
        {
            fprintf(stdout,"cann't open reference file %s...\n",in_file);
            return 1;
        }
        fprintf(stdout,"process %s...%s\n",q->name.s,asctime(present_tm));
        fprintf(Chr_info,"%s\t%u\t%u\n",q->name.s,(unsigned int)q->seq.l,rel_site);
        rel_site += q->seq.l;

        strcpy(out_file,argv[2]);
        strcat(out_file,"/");
        strcat(out_file,q->name.s);
        strcat(out_file,".hash");

        memset(num,0,sizeof(int)*pow(4,KMER_FRONT));

        sort_array->Enum = 0;
        sort_array->Bnum = 0;
        heap->heap_length = 0;

        move = 0;
        base = 0;

        for (i = 0; i < q->seq.l-KMER_FRONT; i++)
        {
            N = 0;
            if((sort_array->Enum+sort_array->Bnum*BUCKET_LENGTH)>0)
            {
                move =i-sort_array->buf[sort_array->Enum+sort_array->Bnum*BUCKET_LENGTH-1].site;
                base = sort_array->buf[sort_array->Enum+sort_array->Bnum*BUCKET_LENGTH-1].mark;
            }
            if((move>0)&&(move<(KMER_FRONT))) j = KMER_FRONT-move;
            else j = 0;
            for (; j <KMER_FRONT; j++)
            {
                if(nst_nt4_table[(int)q->seq.s[i+j]] == -1)
                {
                    N = 1;
                    break;
                }
                base = base << 2 | nst_nt4_table[(int)q->seq.s[i+j]];
            }
            if(N) continue;
            base = base & CHR_MARK;


            if(sort_array->Enum>=BUCKET_LENGTH)
            {
                qsort(&(heap->sort_buf)[sort_array->Bnum*BUCKET_LENGTH],BUCKET_LENGTH,sizeof(unsigned int),chr_hash_cmp);
                sort_array->Enum = 0;
                sort_array->Bnum++;
            }
            sort_array->buf[sort_array->Enum+sort_array->Bnum*BUCKET_LENGTH].mark = base;
            sort_array->buf[sort_array->Enum+sort_array->Bnum*BUCKET_LENGTH].site = i;
            heap->sort_buf[sort_array->Enum+sort_array->Bnum*BUCKET_LENGTH] = sort_array->Enum+sort_array->Bnum*BUCKET_LENGTH;
            sort_array->Enum++;
        }
        if(sort_array->Enum!=0)
            qsort(&(heap->sort_buf)[sort_array->Bnum*BUCKET_LENGTH],sort_array->Enum,sizeof(unsigned int),chr_hash_cmp);
        gzclose(ref);

        BuildHeap_Chr(heap->sort_buf,heap->heap,heap->heap_num,heap->order,&heap->heap_length,sort_array->Bnum,sort_array->Enum);

        Output_hash = fopen(out_file,"wb");
        i = OutputHeap_Chr(heap->sort_buf,heap->heap,heap->heap_num,heap->order,&heap->heap_length);
        int x = 0;
        temp_base = -1;
        while ((int)i!=-1)
        {
            base = sort_array->buf[i].mark;
            read_site = sort_array->buf[i].site;

            if (base!=temp_base)
            {
                if (temp_base == -1)
                {
                    temp_base = base;

                    buffer_num = 0;
                    buffer_flag = 0;

                    buffer[buffer_num] = read_site;
                    buffer_num++;
                }
                else
                {
                    if (buffer_flag == 0)
                    {
                        num[temp_base] = buffer_num;
                        fprintf(Output_hash,">%hu %d\n",temp_base,buffer_num);
                        for(x= 0;x<buffer_num;x++)
                            fprintf(Output_hash,"%u\n",buffer[x]);
                    }
                    buffer_num = 0;
                    buffer_flag = 0;

                    temp_base = base;

                    buffer[buffer_num] = read_site;
                    buffer_num++;
                }
            }
            else
            {
                if (buffer_flag == 1)
                {
                    i = OutputHeap_Chr(heap->sort_buf,heap->heap,heap->heap_num,heap->order,&heap->heap_length);
                    continue;
                }
                if ((buffer_num+1)>MAX_CHR_HASH_SAME)
                {
                    buffer_flag = 1;
                    buffer_num = 0;
                }
                else
                {
                    buffer[buffer_num] = read_site;
                    buffer_num++;
                }
            }
            i = OutputHeap_Chr(heap->sort_buf,heap->heap,heap->heap_num,heap->order,&heap->heap_length);
        }
        if (buffer_flag == 0)
        {
            num[temp_base] = buffer_num;
            fprintf(Output_hash,">%hu %d\n",temp_base,buffer_num);
            for(x= 0;x<buffer_num;x++)
                fprintf(Output_hash,"%u\n",buffer[x]);
        }
        fclose(Output_hash);
    }
    fclose(Chr_info);
    fclose(Name_info);

    remove(name_file);
    kseq_destroy(q);
    free(num);

    free(sort_array->buf);
    free(sort_array);

    free(heap->sort_buf);
    free(heap->heap);
    free(heap->heap_num);
    free(heap->order);
    free(heap);

    free(buffer);

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"complete build hash index...%s\n",asctime(present_tm));
    return 0;
}

