#include "main.h"
static pthread_mutex_t ReadLock,UpdateLock;

struct job_seed
{
    struct m_opt *opt;
    FILE *in_file_1,*in_file_2;
    FILE *out_file_1;
    FILE *unmap_file_1;
};
int read_file_sam(struct read_map_inf *read,FILE *file)
{
    int i = 0;
    char f_line[MAX_STRING_LENGTH];

    while((i<READ_BUF_LENGTH)&&fgets(f_line,MAX_STRING_LENGTH,file)!=NULL)
    {
        sscanf(f_line,"%s\t%d\t%*s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%d\n",
                read[i].name,&read[i].flag,&read[i].start_site,&read[i].MapQ,read[i].cigar,
                read[i].pname,&read[i].psite,read[i].seq,read[i].qual,&read[i].score);
        read[i].length = strlen(read[i].seq);
        read[i].output_flag = 0;
        read[i].TAG[0] = '\0';
        i++;
    }
    if(i<READ_BUF_LENGTH) return i;

    fpos_t pos;
    fgetpos(file,&pos);
    while(fgets(f_line,MAX_STRING_LENGTH,file)!=NULL)
    {
        sscanf(f_line,"%s\t%d\t%*s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%d\n",
                read[i].name,&read[i].flag,&read[i].start_site,&read[i].MapQ,read[i].cigar,
                read[i].pname,&read[i].psite,read[i].seq,read[i].qual,&read[i].score);
        read[i].length = strlen(read[i].seq);
        read[i].output_flag = 0;
        read[i].TAG[0] = '\0';

        if(strcmp(read[i].name,read[READ_BUF_LENGTH-1].name)!=0) {fsetpos(file,&pos);fflush(file);return i;}
        if(i>=READ_BUF_LENGTH+MAX_CAND_NUM) {fsetpos(file,&pos);fflush(file);return i;}
        i++;
        fgetpos(file,&pos);
    }
    return i;
}
uint32_t char2base_32(char *readseq,int length)
{
    uint32_t baseseq = 0;
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
int read_site_cmp(const void *a,const void *b)
{
    struct read_map_inf *EA,*EB;
    EA = (struct read_map_inf *)a;
    EB = (struct read_map_inf *)b;

    if(EA->ref_site==EB->ref_site)
    {
        if(strcmp(EA->name,EB->name)==0)
            return EB->score-EA->score;
        else return strcmp(EA->name,EB->name);
    }
    else if(EA->ref_site < EB->ref_site) return -1;
    else return 1;

}
void find_cigar(struct exon_hash *hash,struct seed_t *seed,int *seed_num,int read_site)
{
    struct seed_t seed_r;
    int i = 0,x = 0;
    unsigned int ref_site = hash->start;
    char L_MARK = (1<<6)-1;
    int length;

    for (i = 0;i<2;i++)
    {
        if((hash->c[i]==0)&&(hash->l[i]==0)) break;

        length = hash->c[i]&L_MARK;
        if(x+length>opt->hash_front+opt->hash_back)
            continue;
        if(length!=0)
        {
            seed_r.pos = ref_site;
            seed_r.start = read_site;
            seed_r.length = length;
            seed_r.score = length;
            seed_r.abs = ref_site-read_site;
            seed_r.lnum = 0;
            insert_seed(seed,seed_num,SEED_BUF_LENGTH,seed_r);
        }
        x+=length;
        ref_site += length;
        read_site += length;

        if(((hash->c[i])>>6&3)==0) {read_site += hash->l[i];x+=hash->l[i];}
        else if(((hash->c[i])>>6&3)==1) {read_site += hash->l[i];ref_site+=hash->l[i];x+=hash->l[i];}
        else if(((hash->c[i])>>6&3)==2) {ref_site+=hash->l[i];}
        else if(((hash->c[i])>>6&3)==3) {ref_site =hash->end;}
    }
    length = opt->hash_front+opt->hash_back-x;
    if(length!=0)
    {
        seed_r.pos = ref_site;
        seed_r.start = read_site;
        seed_r.length = length;
        seed_r.score = length;
        seed_r.abs = ref_site-read_site;
        seed_r.lnum = 0;
        insert_seed(seed,seed_num,SEED_BUF_LENGTH,seed_r);
    }
}

void extend_seed(struct seed_t *seed,int *seed_num,char *read,int length)
{
    int i ,j;
    int x = 0;
    int chr_order = 0;
    unsigned int chr_start = 0;

    for (i = 0;i<*seed_num;i++)
    {
        for(j = 0;j<opt->chr->total;j++)
        {
            if((seed[i].pos>=opt->chr->list[j].start_site)&&(seed[i].pos<opt->chr->list[j].start_site+opt->chr->list[j].length)) break;
        }
        chr_order = j;chr_start = opt->chr->list[j].start_site;
        int front = 0;
        front = seed[i].start;
        for(x = 0;x<front;x++)
        {
            if((seed[i].pos>chr_start)&&(nst_nt4_table[(int)opt->chr->list[chr_order].seq[seed[i].pos-1-chr_start]]==nst_nt4_table[(int)read[seed[i].start-1]])) {seed[i].pos--;seed[i].start--;seed[i].length++;}
            else break;
        }
        for(x = seed[i].start+seed[i].length;x<length;x++)
        {
            if((seed[i].pos+seed[i].length<opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length)
                &&(nst_nt4_table[(int)opt->chr->list[chr_order].seq[seed[i].pos+seed[i].length-chr_start]]==nst_nt4_table[(int)read[x]])) seed[i].length++;
            else break;
        }
    }
}
void insert_read(struct read_map_inf *read,struct read_map_inf *out,int *out_num)
{
    if((*out_num)>=MAP_READ_BUF_LENGTH)
    {
        pthread_mutex_lock(&UpdateLock);
        qsort(out,(*out_num),sizeof(struct read_map_inf),read_site_cmp);

        FILE *read_file;
        char read_file_name[MAX_NAME_LENGTH];
        sprintf(read_file_name,"%s/%d.read",opt->Temp_path,opt->find_fn);
        read_file = fopen(read_file_name,"wb");

        if(read_file == NULL)
        {
            fprintf(stdout,"[DEEP]cannot open file %s...\n",read_file_name);
        }
        if(fwrite(out,sizeof(struct read_map_inf),MAP_READ_BUF_LENGTH,read_file)!=MAP_READ_BUF_LENGTH)
            fprintf(stdout,"[DEEP]cannot write file...\n");

        fclose(read_file);
        opt->find_fn++;
        (*out_num) = 0;
        pthread_mutex_unlock(&UpdateLock);
    }
    memcpy(out+(*out_num),read,sizeof(struct read_map_inf));
    (*out_num)++;
}
void find_hash_seed(struct m_opt *opt,struct read_map_inf *read,struct read_inf_t *read_a,struct seed_t *seed,int *seed_num,struct seed_t *rseed,int *rseed_num,struct read_map_inf *out,int *out_num)
{
    int i = 0,j = 0;
    int seed_length = opt->hash_front+opt->hash_back;
    uint32_t front = 0,back = 0;
    uint32_t fmax = pow(4,opt->hash_front);
    uint32_t bmax = pow(4,opt->hash_back);
    int length;
    int hash_site = 0;

    struct seed_t seed_r;

    unsigned int SAREA = 0;
    unsigned int EAREA = 0;

    int chr_order = 0;

    if(read->cigar[0]!='*')
    {

        for(i = 0;i<opt->chr->total;i++)
        {
            if((read->start_site>=opt->chr->list[i].start_site)&&(read->start_site<opt->chr->list[i].start_site+opt->chr->list[i].length))
            {
                chr_order = i;
                break;
            }
        }

        struct cigar_t cigar[MAX_CIGAR_BUF];
        int cigar_total = 0;
        unsigned int ref_site = read->start_site;
        unsigned int end_site;
        int read_site = 0;

        process_cigar(read->cigar,cigar,&cigar_total,&end_site);

        for(i = 0;i<cigar_total;i++)
        {
            switch(cigar[i].c)
            {
                case 'M':
                    seed_r.pos = ref_site;
                    seed_r.start = read_site;
                    seed_r.length = cigar[i].l;
                    seed_r.score = cigar[i].l;
                    seed_r.abs = ref_site-read_site;
                    seed_r.lnum = 0;
                    insert_seed(seed,seed_num,SEED_BUF_LENGTH,seed_r);

                    ref_site+=cigar[i].l;
                    read_site+=cigar[i].l;
                    break;
                case 'X':
                    ref_site+=cigar[i].l;
                    read_site+=cigar[i].l;
                    break;
                case 'I':
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
        if(cigar[0].c=='S')
        {
            if(read->start_site>opt->area)
                SAREA = max(read->start_site-opt->area,opt->chr->list[chr_order].start_site);
            else SAREA = opt->chr->list[chr_order].start_site;
            EAREA = read->start_site;

            i = 0;
            while ((i<cigar[0].l)&&(i<=read->length-seed_length))
            {
                front = char2base_32(read->seq+i,opt->hash_front);
                back = char2base_32(read->seq+i+opt->hash_front,opt->hash_back);
                if((front<fmax)&&(back<bmax))
                    hash_site = hash_e_find(opt->e_hash,opt->e_num,front,back);
                else
                    hash_site = -1;

                if(hash_site!=-1)
                {
                    while((opt->e_num[front]>hash_site)&&(back == opt->e_hash[front][hash_site].back))
                    {
                        if((opt->e_hash[front][hash_site].start>=SAREA)&&(opt->e_hash[front][hash_site].start<=EAREA))
                            find_cigar(&(opt->e_hash[front][hash_site]),seed,seed_num,i);
                        hash_site++;
                    }
                }
                i++;
            }
        }
        if(cigar[cigar_total-1].c=='S')
        {
            SAREA = end_site;
            EAREA = min(end_site+opt->area,opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length-1);

            i = max(0,read->length-cigar[cigar_total-1].l-seed_length+1);
            while (i<=read->length-seed_length)
            {
                front = char2base_32(read->seq+i,opt->hash_front);
                back = char2base_32(read->seq+i+opt->hash_front,opt->hash_back);
                if((front<fmax)&&(back<bmax))
                    hash_site = hash_e_find(opt->e_hash,opt->e_num,front,back);
                else
                    hash_site = -1;

                if(hash_site!=-1)
                {
                    while((opt->e_num[front]>hash_site)&&(back == opt->e_hash[front][hash_site].back))
                    {
                        if((opt->e_hash[front][hash_site].start>=SAREA)&&(opt->e_hash[front][hash_site].start<=EAREA))
                            find_cigar(&(opt->e_hash[front][hash_site]),seed,seed_num,i);
                        hash_site++;
                    }
                }
                i++;
            }
        }
    }
    else
    {
        length = read->length;
    for(j = 0;j<read->length;j++)
    {
        if((read->seq[length-j-1]=='A')||(read->seq[length-j-1]=='a')) read_a->rseq[j] = 'T';
        if((read->seq[length-j-1]=='C')||(read->seq[length-j-1]=='c')) read_a->rseq[j] = 'G';
        if((read->seq[length-j-1]=='G')||(read->seq[length-j-1]=='g')) read_a->rseq[j] = 'C';
        if((read->seq[length-j-1]=='T')||(read->seq[length-j-1]=='t')) read_a->rseq[j] = 'A';
        if((read->seq[length-j-1]=='N')||(read->seq[length-j-1]=='n')) read_a->rseq[j] = 'N';
        read_a->rqual[j] = read_a->rqual[length-j-1];
    }
    if(read->qual[0]=='*') {read_a->rqual[0] = '*';read_a->rqual[1] = '\0';}
    else read_a->rqual[length] = '\0';
    read_a->rseq[length] = '\0';

    if(read->psite!=0)
    {
        int chr_order = 0;
        for(i = 0;i<opt->chr->total;i++)
        {
            if((read->psite>=opt->chr->list[i].start_site)&&(read->psite<opt->chr->list[i].start_site+opt->chr->list[i].length))
            {
                chr_order = i;
                break;
            }
        }
        if(read->psite>2*opt->area)
            SAREA = max(read->psite-2*opt->area,opt->chr->list[chr_order].start_site);
        else SAREA = opt->chr->list[chr_order].start_site;
        EAREA = min(read->psite+2*opt->area,opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length-1);
    }

    i = 0;
    while (i+seed_length<=length)
    {
        front = char2base_32(read->seq+i,opt->hash_front);
        back = char2base_32(read->seq+i+opt->hash_front,opt->hash_back);
        if((front<fmax)&&(back<bmax))
            hash_site = hash_e_find(opt->e_hash,opt->e_num,front,back);
        else
            hash_site = -1;

        if(hash_site!=-1)
        {
            while((opt->e_num[front]>hash_site)&&(back == opt->e_hash[front][hash_site].back))
            {
                if((SAREA!=0)&&(EAREA!=0)&&(opt->e_hash[front][hash_site].start>=SAREA)&&(opt->e_hash[front][hash_site].start<=EAREA))
                    find_cigar(&(opt->e_hash[front][hash_site]),seed,seed_num,i);
                else if((SAREA==0)&&(EAREA==0))
                    find_cigar(&(opt->e_hash[front][hash_site]),seed,seed_num,i);
                hash_site++;
            }
        }

        front = char2base_32(&(read_a->rseq[i]),opt->hash_front);
        back = char2base_32(&(read_a->rseq[i+opt->hash_front]),opt->hash_back);
        if((front<fmax)&&(back<bmax))
            hash_site = hash_e_find(opt->e_hash,opt->e_num,front,back);
        else
            hash_site = -1;

        if(hash_site!=-1)
        {
            while((opt->e_num[front]>hash_site)&&(back == opt->e_hash[front][hash_site].back))
            {
                if((SAREA!=0)&&(EAREA!=0)&&(opt->e_hash[front][hash_site].start>=SAREA)&&(opt->e_hash[front][hash_site].start<=EAREA))
                    find_cigar(&(opt->e_hash[front][hash_site]),rseed,rseed_num,i);
                else if((SAREA==0)&&(EAREA==0))
                    find_cigar(&(opt->e_hash[front][hash_site]),rseed,rseed_num,i);
                hash_site++;
            }
        }
        i+=seed_length;
    }
    if(i<length)
    {
        i = length-seed_length;
        front = char2base_32(read->seq+i,opt->hash_front);
        back = char2base_32(read->seq+i+opt->hash_front,opt->hash_back);
        if((front<fmax)&&(back<bmax))
            hash_site = hash_e_find(opt->e_hash,opt->e_num,front,back);
        else
            hash_site = -1;

        if(hash_site!=-1)
        {
            while((opt->e_num[front]>hash_site)&&(back == opt->e_hash[front][hash_site].back))
            {
                if((SAREA!=0)&&(EAREA!=0)&&(opt->e_hash[front][hash_site].start>=SAREA)&&(opt->e_hash[front][hash_site].start<=EAREA))
                    find_cigar(&(opt->e_hash[front][hash_site]),seed,seed_num,i);
                else if((SAREA==0)&&(EAREA==0))
                    find_cigar(&(opt->e_hash[front][hash_site]),seed,seed_num,i);
                hash_site++;
            }
        }

        front = char2base_32(&(read_a->rseq[i]),opt->hash_front);
        back = char2base_32(&(read_a->rseq[i+opt->hash_front]),opt->hash_back);
        if((front<fmax)&&(back<bmax))
            hash_site = hash_e_find(opt->e_hash,opt->e_num,front,back);
        else
            hash_site = -1;

        if(hash_site!=-1)
        {
            while((opt->e_num[front]>hash_site)&&(back == opt->e_hash[front][hash_site].back))
            {
                if((SAREA!=0)&&(EAREA!=0)&&(opt->e_hash[front][hash_site].start>=SAREA)&&(opt->e_hash[front][hash_site].start<=EAREA))
                    find_cigar(&(opt->e_hash[front][hash_site]),rseed,rseed_num,i);
                else if((SAREA==0)&&(EAREA==0))
                    find_cigar(&(opt->e_hash[front][hash_site]),rseed,rseed_num,i);
                hash_site++;
            }
        }
    }

    if((*seed_num==0)&&(*rseed_num==0))
    {
        i = 0;
    while (i+seed_length<=length)
    {
        front = char2base_32(read->seq+i,opt->hash_front);
        back = char2base_32(read->seq+i+opt->hash_front,opt->hash_back);
        if((front<fmax)&&(back<bmax))
            hash_site = hash_e_find(opt->e_hash,opt->e_num,front,back);
        else
            hash_site = -1;

        if(hash_site!=-1)
        {
            while((opt->e_num[front]>hash_site)&&(back == opt->e_hash[front][hash_site].back))
            {
                if((SAREA!=0)&&(EAREA!=0)&&(opt->e_hash[front][hash_site].start>=SAREA)&&(opt->e_hash[front][hash_site].start<=EAREA))
                    find_cigar(&(opt->e_hash[front][hash_site]),seed,seed_num,i);
                else if((SAREA==0)&&(EAREA==0))
                    find_cigar(&(opt->e_hash[front][hash_site]),seed,seed_num,i);
                hash_site++;
            }
        }

        front = char2base_32(&(read_a->rseq[i]),opt->hash_front);
        back = char2base_32(&(read_a->rseq[i+opt->hash_front]),opt->hash_back);
        if((front<fmax)&&(back<bmax))
            hash_site = hash_e_find(opt->e_hash,opt->e_num,front,back);
        else
            hash_site = -1;

        if(hash_site!=-1)
        {
            while((opt->e_num[front]>hash_site)&&(back == opt->e_hash[front][hash_site].back))
            {
                if((SAREA!=0)&&(EAREA!=0)&&(opt->e_hash[front][hash_site].start>=SAREA)&&(opt->e_hash[front][hash_site].start<=EAREA))
                    find_cigar(&(opt->e_hash[front][hash_site]),rseed,rseed_num,i);
                else if((SAREA==0)&&(EAREA==0))
                    find_cigar(&(opt->e_hash[front][hash_site]),rseed,rseed_num,i);
            hash_site++;
            }
        }
        i++;
    }
    }

        extend_seed(seed,seed_num,read->seq,read->length);
        extend_seed(rseed,rseed_num,read_a->rseq,read->length);

        for(i = 0;i<(*rseed_num);i++)
        {
            if((*seed_num)>=SEED_BUF_LENGTH) break;
            seed[(*seed_num)].pos = (opt->idx->bns->l_pac<<1)-rseed[i].length-rseed[i].pos;
            seed[(*seed_num)].start = read->length-rseed[i].length-rseed[i].start;
            seed[(*seed_num)].length = rseed[i].length;
            seed[(*seed_num)].score = rseed[i].length;
            seed[(*seed_num)].abs = seed[(*seed_num)].pos-seed[(*seed_num)].start;
            seed[(*seed_num)].lnum = 0;
            (*seed_num)++;
        }
        qsort(seed,*seed_num,sizeof(struct seed_t),seed_cmp);
    }
}
void write_TAG(struct read_map_inf *read);
void out_put_hash_read_deep(struct read_map_inf *read,struct read_map_inf *out,int *out_num,FILE *out_file)
{
    int i;
    struct cigar_t cigar[MAX_CIGAR_BUF];
    int cigar_total;
    unsigned int ref_site;
    unsigned int end_site;
    int read_site;
    int flag;

    ref_site = read->start_site;
    read_site = 0;
    flag = 1;

    /*int N = 0;
    for(i = 0;i<read->length;i++)
    {
        if(nst_nt4_table[(int)(read->seq[i])]>=4){N =1;break;}
    }
    if(N==1)
    {
        if(read->cigar[0] =='*')
        {
            fprintf(out_file,"%s\t%d\t*\t0\t0\t*\t%s\t0\t0\t%s\t%s\t0\n",read[i].name,read[i].flag,read[i].pname,read[i].seq,read[i].qual);
        }
        else
        {
            for(i = 0;i<opt->chr->total;i++)
            {
                if((read->start_site>=opt->chr->list[i].start_site)&&(read->start_site<opt->chr->list[i].start_site+opt->chr->list[i].length))
                {
                    read->chr_order = i;
                    break;
                }
            }

            read->start_site-=opt->chr->list[read->chr_order].start_site;
            if(read->TAG[0]=='\0') write_TAG(read);
            char TAG[MAX_STRING_LENGTH];
            sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read->score,read->dis,read->TAG);
                //+MD
            fprintf(out_file,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%s\n",
                read->name,read->flag,opt->chr->list[read->chr_order].name,read->start_site+1,read->MapQ,read->cigar,
                read->pname,read->psite,read->seq,read->qual,TAG);
        }
    }*/
    if(read->cigar[0]=='*')
    {
        read->ref_site = 0;
        read->front = 0;
        read->back = 0;
        insert_read(read,out,out_num);
    }
    else
    {
    process_cigar(read->cigar,cigar,&cigar_total,&end_site);
    for(i = 0;i<cigar_total;i++)
    {
        switch(cigar[i].c)
        {
        case 'M':
            if((flag == 1)&&(read_site+cigar[i].l>10)&&(cigar[i].l>=12))
            {

                read->ref_site = ref_site+opt->change_length;
                read->front = read_site+opt->change_length;
                read->back = read_site+cigar[i].l-opt->change_length;
                read->end_site = ref_site + cigar[i].l-1;

                flag = 0;
                insert_read(read,out,out_num);
            }
            ref_site+=cigar[i].l;
            read_site+=cigar[i].l;
            break;
        case 'X':
            ref_site+=cigar[i].l;
            read_site+=cigar[i].l;
            break;
        case 'I':
            read_site+=cigar[i].l;
            break;
        case 'D':
            ref_site+=cigar[i].l;
            break;
        case 'N':
            flag = 1;
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
}
void out_put_hash_read(struct read_map_inf *read,int *out_num,FILE *out_file)
{
    int i,x;
    char TAG[MAX_STRING_LENGTH];
    int length = 0;

    pthread_mutex_lock(&UpdateLock);
    for (i = 0;i<(*out_num);i++)
    {
        length = strlen(read[i].name);
        read[i].name[length-1] = '\0';
        if(read[i].pname[0]!='*')
        {
            length = strlen(read[i].pname);
            read[i].pname[length-1] = '\0';
        }
        if(read[i].cigar[0] =='*')
        {
            fprintf(out_file,"%s\t%d\t*\t0\t0\t*\t%s\t0\t0\t%s\t%s\t0\n",read[i].name,read[i].flag,read[i].pname,read[i].seq,read[i].qual);
        }
        else
        {
            for(x = 0;x<opt->chr->total;x++)
            {
                if((read[i].start_site>=opt->chr->list[x].start_site)&&(read[i].start_site<opt->chr->list[x].start_site+opt->chr->list[x].length))
                {
                    read[i].chr_order = x;
                    break;
                }
            }

            read[i].start_site-=opt->chr->list[read[i].chr_order].start_site;
            if(read[i].TAG[0]=='\0') write_TAG(read+i);
            sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read[i].score,read[i].dis,read[i].TAG);
                //+MD
            fprintf(out_file,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%s\n",
                read[i].name,read[i].flag,opt->chr->list[read[i].chr_order].name,read[i].start_site+1,read[i].MapQ,read[i].cigar,
                read[i].pname,read[i].psite,read[i].seq,read[i].qual,TAG);
        }
    }
    (*out_num) = 0;
    pthread_mutex_unlock(&UpdateLock);

}
void out_put_N_read(struct read_map_inf *read,FILE *out_file)
{
    int i;

    pthread_mutex_lock(&UpdateLock);
    if(read->cigar[0] =='*')
    {
            fprintf(out_file,"%s\t%d\t*\t0\t0\t*\t%s\t0\t0\t%s\t%s\t0\n",read->name,read->flag,read->pname,read->seq,read->qual);
    }
    else
    {
            for(i = 0;i<opt->chr->total;i++)
            {
                if((read->start_site>=opt->chr->list[i].start_site)&&(read->start_site<opt->chr->list[i].start_site+opt->chr->list[i].length))
                {
                    read->chr_order = i;
                    break;
                }
            }

            read->start_site-=opt->chr->list[read->chr_order].start_site;
            if(read->TAG[0]=='\0') write_TAG(read);
            char TAG[MAX_STRING_LENGTH];
            sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read->score,read->dis,read->TAG);
                //+MD
            fprintf(out_file,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%s\n",
                read->name,read->flag,opt->chr->list[read->chr_order].name,read->start_site+1,read->MapQ,read->cigar,
                read->pname,read->psite,read->seq,read->qual,TAG);
    }
    pthread_mutex_unlock(&UpdateLock);
}
int is_output(struct read_map_inf *read)
{
    int i;
    int l = strlen(read->cigar);
    for(i = 0;i<l;i++)
    {
        if(read->cigar[i]=='S') {read->output_flag = 0;return 0;}
    }
    {read->output_flag = 1;return 1;}
}
void *find_site_core(void* arg)
{
    struct job_seed *job = (struct job_seed *)arg;

	struct seed_t *seed = (struct seed_t *)calloc(SEED_BUF_LENGTH,sizeof(struct seed_t));
	int seed_num;
	struct seed_t *rseed = (struct seed_t *)calloc(SEED_BUF_LENGTH,sizeof(struct seed_t));
	int rseed_num;
	struct read_map_inf *read = (struct read_map_inf *)calloc(READ_BUF_LENGTH+MAX_CAND_NUM+1,sizeof(struct read_map_inf));
	int read_num = 0;
	struct read_map_inf *out = (struct read_map_inf *)calloc(MAP_READ_BUF_LENGTH+1,sizeof(struct read_map_inf));
	int out_num = 0;

	struct read_inf_t *read_a = (struct read_inf_t *)calloc(1,sizeof(struct read_inf_t));
	struct snp_list_t *snp = NULL;
    struct exon_array *exon = NULL;


	int i = 0,x = 0,j = 0;
	int start = 0;
	int max = 0;
	int max_n = 0;

	int flag = 0;

	while (1)
	{
        //read_file;
        pthread_mutex_lock(&ReadLock);
        read_num = read_file_sam(read,job->in_file_1);
		pthread_mutex_unlock(&ReadLock);

		if(read_num==0) break;
		start = 0;
		strcpy(read_a->name,read[0].name);
        strcpy(read_a->seq,read[0].seq);
        strcpy(read_a->qual,read[0].qual);
        read_a->length = read[0].length;
        read_a->cand_num = 0;
        read_a->out_flag = 0;
        //seed_align
        for (i = 0;i<=read_num;i++)
        {
            if((i!=read_num)&&(strcmp(read[start].name,read[i].name)==0))
            {
                read_a->cand_num = 0;
                seed_num = 0;
                rseed_num = 0;
                find_hash_seed(job->opt,read+i,read_a,seed,&seed_num,rseed,&rseed_num,out,&out_num);
                qsort(seed,seed_num,sizeof(struct seed_t),seed_cmp);
                find_cand(job->opt,read_a,seed,&seed_num,exon,snp,0.35*strlen(read->seq));

                if(read[i].cigar[0]!='*')
                {
                    if((read_a->cand_num==0)&&(opt->deep_mode==1))
                    {
                        for(x = 0;x<opt->chr->total;x++)
                        {
                            if((read[i].start_site>=opt->chr->list[x].start_site)&&(read[i].start_site<opt->chr->list[x].start_site+opt->chr->list[x].length))
                            {
                                read[i].chr_order = x;
                                break;
                            }
                        }
                    }
                    max = -100;
                    for(x = 0;x<read_a->cand_num;x++)
                    {
                        if(read_a->cand[x].score>max) max = read_a->cand[x].score;
                    }
                    for(x = 0;x<read_a->cand_num;x++)
                    {
                        if((read_a->cand[x].score==max)&&(read_a->cand[x].score>read[i].score))
                        {
                            read[i].start_site = read_a->cand[x].pos;
                            read[i].score = read_a->cand[x].score;
                            read[i].dis = read_a->cand[x].dis;
                            read[i].chr_order = read_a->cand[x].chr_order;
                            strcpy(read[i].cigar,read_a->cand[x].cigar);
                            strcpy(read[i].TAG,read_a->cand[x].TAG);
                            break;
                        }
                    }
                }
                else
                {
                    generate_MapQ(read_a);
                    read[i].MapQ = read_a->MapQ;
                }
            }
            else
            {
                if(opt->deep_mode==1)
                {
                    int N = 0;
                    for(x = 0;x<read[start].length;x++)
                    {
                        if(nst_nt4_table[(int)(read[start].seq[x])]>=4){N =1;break;}
                    }
                    if(N==1)
                    {
                        if(read[start].cigar[0]!='*')
                        {
                            max = -100;
                            max_n = 0;
                            for(x = start;x<i;x++)
                            {
                                if(read[x].score>max) {max = read[x].score;max_n = 1;}
                                else if(read[x].score==max) {max_n++;}
                            }
                            if (max_n >= 5) read_a->MapQ = 0;
                            else if (max_n == 4) read_a->MapQ = 1;
                            else if (max_n== 3) read_a->MapQ = 2;
                            else if (max_n == 2) read_a->MapQ = 3;
                            else read_a->MapQ = 50;

                            for(x = start;x<i;x++)
                            {
                                if(((opt->output_mode==OUTPUT_BEST)&&(read[x].score==max))||(opt->output_mode==OUTPUT_ALL))
                                {
                                    flag = 1;
                                    for(j = start;j<x;j++)
                                    {
                                        if((read[x].chr_order==read[j].chr_order)&&(read[x].start_site==read[j].start_site))
                                        {
                                            if(strcmp(read[x].cigar,read[j].cigar)==0){flag = 0;break;}
                                        }
                                    }
                                    if(flag)
                                    {
                                        read[x].MapQ = read_a->MapQ;
                                        out_put_N_read(read+x,job->out_file_1);
                                    }
                                }
                            }
                        }
                        else
                        {
                            if(read_a->cand_num==0)
                            {
                                read[start].flag|=0x0008;
                                out_put_N_read(read+start,job->out_file_1);
                            }
                            else
                            {
                                max = -100;
                                max_n = 0;
                                for(x = 0;x<read_a->cand_num;x++)
                                {
                                    if(read_a->cand[x].score>max) {max = read_a->cand[x].score;max_n = 1;}
                                    else if(read_a->cand[x].score==max) {max_n++;}
                                }
                                if (max_n >= 5) read_a->MapQ = 0;
                                else if (max_n == 4) read_a->MapQ = 1;
                                else if (max_n== 3) read_a->MapQ = 2;
                                else if (max_n == 2) read_a->MapQ = 3;
                                else read_a->MapQ = 50;

                                flag = read[start].flag;
                                for(x = 0;x<read_a->cand_num;x++)
                                {
                                    if(((opt->output_mode==OUTPUT_BEST)&&(read_a->cand[x].score==max))||(opt->output_mode==OUTPUT_ALL))
                                    {
                                    if(read_a->cand[x].strand)
                                    {
                                        read[start].flag = flag|0x0010;
                                        strcpy(read[start].seq,read_a->rseq);
                                        strcpy(read[start].qual,read_a->rqual);
                                    }
                                    else
                                    {
                                        read[start].flag = flag;
                                        strcpy(read[start].seq,read_a->seq);
                                        strcpy(read[start].qual,read_a->qual);
                                    }
                                    read[start].MapQ = read_a->MapQ;
                                    read[start].start_site = read_a->cand[x].pos;
                                    read[start].chr_order = read_a->cand[x].chr_order;
                                    read[start].score = read_a->cand[x].score;
                                    read[start].dis = read_a->cand[x].dis;
                                    strcpy(read[start].cigar,read_a->cand[x].cigar);
                                    strcpy(read[start].TAG,read_a->cand[x].TAG);

                                    out_put_N_read(read+start,job->out_file_1);
                                    }

                                }
                            }
                        }
                    }
                    else
                    {
                    if(read[start].cigar[0]!='*')
                    {
                        for(x = start;x<i;x++)
                        {
                            flag = 1;
                            for(j = start;j<x;j++)
                            {
                                if((read[x].chr_order==read[j].chr_order)&&(read[x].start_site==read[j].start_site))
                                {
                                    if(strcmp(read[x].cigar,read[j].cigar)==0){flag = 0;break;}
                                }
                            }
                            if(flag) {is_output(read+x);out_put_hash_read_deep(read+x,out,&out_num,job->out_file_1);}
                        }
                    }
                    else
                    {
                        if(read_a->cand_num==0)
                        {
                            read[start].flag|=0x0008;
                            out_put_hash_read_deep(read+start,out,&out_num,job->out_file_1);
                        }
                        else
                        {
                            max = -100;
                            for(x = 0;x<read_a->cand_num;x++)
                            {
                                if(read_a->cand[x].score>max) max = read_a->cand[x].score;
                            }

                            flag = read[start].flag;
                            for(x = 0;x<read_a->cand_num;x++)
                            {
                                if(read_a->cand[x].strand)
                                {
                                    read[start].flag = flag|0x0010;
                                    strcpy(read[start].seq,read_a->rseq);
                                    strcpy(read[start].qual,read_a->rqual);
                                }
                                else
                                {
                                    read[start].flag = flag;
                                    strcpy(read[start].seq,read_a->seq);
                                    strcpy(read[start].qual,read_a->qual);
                                }
                                read[start].start_site = read_a->cand[x].pos;
                                read[start].chr_order = read_a->cand[x].chr_order;
                                read[start].score = read_a->cand[x].score;
                                read[start].dis = read_a->cand[x].dis;
                                strcpy(read[start].cigar,read_a->cand[x].cigar);
                                strcpy(read[start].TAG,read_a->cand[x].TAG);

                                is_output(read+start);
                                out_put_hash_read_deep(read+start,out,&out_num,job->out_file_1);
                            }
                        }
                    }
                    }
                }
                else if(opt->output_mode==OUTPUT_BEST)
                {
                    if(read[start].cigar[0]!='*')
                    {
                        max = -100;
                        max_n = 0;
                        for(x = start;x<i;x++)
                        {
                            if(read[x].score>max) {max = read[x].score;max_n = 1;}
                            else if(read[x].score==max) {max_n++;}
                        }
                        if (max_n >= 5) read_a->MapQ = 0;
                        else if (max_n == 4) read_a->MapQ = 1;
                        else if (max_n== 3) read_a->MapQ = 2;
                        else if (max_n == 2) read_a->MapQ = 3;
                        else read_a->MapQ = 50;

                        for(x = start;x<i;x++)
                        {
                            if(read[x].score==max)
                            {
                                flag = 1;
                                for(j = start;j<x;j++)
                                {
                                    if((read[x].chr_order==read[j].chr_order)&&(read[x].start_site==read[j].start_site))
                                    {
                                        if(strcmp(read[x].cigar,read[j].cigar)==0){flag = 0;break;}
                                    }
                                }
                                if(flag)
                                {
                                    read[x].MapQ = read_a->MapQ;
                                    memcpy(out+out_num,read+x,sizeof(struct read_map_inf));
                                    out_num++;
                                    if(out_num>=MAP_READ_BUF_LENGTH)
                                        out_put_hash_read(out,&out_num,job->out_file_1);
                                }
                            }
                        }
                    }
                    else
                    {
                        if(read_a->cand_num==0)
                        {
                            read[start].flag|=0x0008;
                            memcpy(out+out_num,read+start,sizeof(struct read_map_inf));
                            out_num++;
                            if(out_num>=MAP_READ_BUF_LENGTH)
                                out_put_hash_read(out,&out_num,job->out_file_1);
                        }
                        else
                        {
                            flag = read[start].flag;
                        max = -100;
                        for(x = 0;x<read_a->cand_num;x++)
                        {
                            if(read_a->cand[x].score>max) max = read_a->cand[x].score;
                        }
                        for(x = 0;x<read_a->cand_num;x++)
                        {
                            if(read_a->cand[x].score==max)
                            {
                                if(read_a->cand[x].strand)
                                {
                                    read[start].flag = flag|0x0010;
                                    strcpy(read[start].seq,read_a->rseq);
                                    strcpy(read[start].qual,read_a->rqual);
                                }
                                else
                                {
                                    read[start].flag = flag;
                                    strcpy(read[start].seq,read_a->seq);
                                    strcpy(read[start].qual,read_a->qual);
                                }
                                read[start].start_site = read_a->cand[x].pos;
                                read[start].chr_order = read_a->cand[x].chr_order;
                                read[start].score = read_a->cand[x].score;
                                read[start].dis = read_a->cand[x].dis;
                                strcpy(read[start].cigar,read_a->cand[x].cigar);
                                strcpy(read[start].TAG,read_a->cand[x].TAG);

                                memcpy(out+out_num,read+start,sizeof(struct read_map_inf));
                                out_num++;
                                if(out_num>=MAP_READ_BUF_LENGTH)
                                    out_put_hash_read(out,&out_num,job->out_file_1);
                            }
                        }
                        }
                    }
                }
                else
                {
                    if(read[start].cigar[0]!='*')
                    {
                        max = -100;
                        max_n = 0;
                        for(x = start;x<i;x++)
                        {
                            if(read[x].score>max) {max = read[x].score;max_n = 1;}
                            else if(read[x].score==max) {max_n++;}
                        }
                        if (max_n >= 5) read_a->MapQ = 0;
                        else if (max_n == 4) read_a->MapQ = 1;
                        else if (max_n== 3) read_a->MapQ = 2;
                        else if (max_n == 2) read_a->MapQ = 3;
                        else read_a->MapQ = 50;

                        for(x = start;x<i;x++)
                        {
                            flag = 1;
                            for(j = start;j<x;j++)
                            {
                                if((read[x].chr_order==read[j].chr_order)&&(read[x].start_site==read[j].start_site))
                                {
                                    if(strcmp(read[x].cigar,read[j].cigar)==0){flag = 0;break;}
                                }
                            }
                            if(flag)
                            {
                                memcpy(out+out_num,read+x,sizeof(struct read_map_inf));
                                out_num++;
                                if(out_num>=MAP_READ_BUF_LENGTH)
                                    out_put_hash_read(out,&out_num,job->out_file_1);
                            }
                        }
                    }
                    else
                    {
                        if(read_a->cand_num==0)
                        {
                            read[start].flag|=0x0008;
                            memcpy(out+out_num,read+start,sizeof(struct read_map_inf));
                            out_num++;
                            if(out_num>=MAP_READ_BUF_LENGTH)
                                out_put_hash_read(out,&out_num,job->out_file_1);
                        }
                        else
                        {
                            flag = read[start].flag;
                            for(x = 0;x<read_a->cand_num;x++)
                            {
                                if(read_a->cand[x].strand)
                                {
                                    read[start].flag = flag|0x0010;
                                    strcpy(read[start].seq,read_a->rseq);
                                    strcpy(read[start].qual,read_a->rqual);
                                }
                                else
                                {
                                    read[start].flag = flag;
                                    strcpy(read[start].seq,read_a->seq);
                                    strcpy(read[start].qual,read_a->qual);
                                }
                                read[start].start_site = read_a->cand[x].pos;
                                read[start].chr_order = read_a->cand[x].chr_order;
                                read[start].score = read_a->cand[x].score;
                                read[start].dis = read_a->cand[x].dis;
                                strcpy(read[start].cigar,read_a->cand[x].cigar);
                                strcpy(read[start].TAG,read_a->cand[x].TAG);

                                memcpy(out+out_num,read+start,sizeof(struct read_map_inf));
                                out_num++;
                                if(out_num>=MAP_READ_BUF_LENGTH)
                                    out_put_hash_read(out,&out_num,job->out_file_1);
                            }
                        }
                    }
                }
                start = i;

                strcpy(read_a->name,read[i].name);
                strcpy(read_a->seq,read[i].seq);
                strcpy(read_a->qual,read[i].qual);
                read_a->length = read[i].length;
                read_a->cand_num = 0;
                read_a->out_flag = 0;

                if(i!=read_num)i--;
            }
        }
	}
	if((opt->deep_mode==1)&&(out_num>0))
    {
        pthread_mutex_lock(&UpdateLock);
        qsort(out,out_num,sizeof(struct read_map_inf),read_site_cmp);

        FILE *read_file;
        char read_file_name[MAX_NAME_LENGTH];
        sprintf(read_file_name,"%s/%d.read",opt->Temp_path,opt->find_fn);
        read_file = fopen(read_file_name,"wb");

        if(read_file == NULL)
        {
            fprintf(stdout,"[DEEP]cannot open file %s...\n",read_file_name);
        }
        if(fwrite(&out[0],sizeof(struct read_map_inf),out_num,read_file)!=out_num)
            fprintf(stdout,"[DEEP]cannot write file...\n");

        fclose(read_file);
        opt->find_fn++;
        pthread_mutex_unlock(&UpdateLock);
    }
    else if((opt->deep_mode==0)&&(out_num>0))
        out_put_hash_read(out,&out_num,job->out_file_1);

	free(seed);
	free(rseed);
	free(read);
	free(out);
	free(read_a);
    return (void*)(0);
}
int find_site(struct m_opt *opt)
{
    opt->step_flag = SNP_STEP;
    //open_file input output unmap
    FILE *in_file_1;
    FILE *out_file_1;
    char out_file_name_1[MAX_NAME_LENGTH+50];

    //char unmap[MAX_NAME_LENGTH+50];
    //sprintf(unmap,"%s/unmap_out.sam",opt->Temp_path);
    //in_file_1 = fopen(unmap,"r");

    in_file_1 = fopen(opt->input_file_1->file[0].name,"r");
    if (in_file_1 == NULL)
        return 1;

    sprintf(out_file_name_1,"%s/map_out.sam",opt->Output_path);
    //sprintf(out_file_name_1,"%s/map_out.sam",opt->Temp_path);
    out_file_1 = fopen(out_file_name_1,"a+");
    if (out_file_1 == NULL)
        return 1;

    struct job_seed *job = (struct job_seed *)calloc(1,sizeof(struct job_seed));
    job->opt = opt;
    job->in_file_1 = in_file_1;
    job->out_file_1 = out_file_1;

	pthread_t *pthreads = malloc(sizeof(pthread_t) * opt->thread_num);

	int i = 0;
    for (i = 0; i < opt->thread_num; i++) pthread_create(&pthreads[i], NULL, find_site_core, job);
	for (i = 0; i < opt->thread_num; i++) pthread_join(pthreads[i], NULL);

	free(job);
	free(pthreads);
    fclose(in_file_1);
    remove(opt->input_file_1->file[0].name);
    fclose(out_file_1);

    return 0;
}

