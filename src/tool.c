#include "main.h"
int generate_MapQ(struct read_inf_t *read)
{
    int i = 0;
    int num = 0;
    int max_score = 0;

    for (i = 0; i < read->cand_num; i++)
    {
        if (max_score >read->cand[i].score) max_score = read->cand[i].score;
    }
    for (i = 0; i < read->cand_num; i++)
    {
        if (max_score == read->cand[i].score) num++;
    }
    if (num >= 5) return 0;
    else if (num == 4) return 1;
    else if (num == 3) return 2;
    else if (num == 2) return 3;
    else return 50;
}
int generate_flag_single(struct m_opt *opt,struct read_inf_t *read,int cand_order)
{
    int flag = 0;

    if(read->cand[cand_order].cigar[0]=='\0') {flag = 0x0004;return flag;}
    else
    {
        if(read->cand[cand_order].strand==1) flag|=0x0010;
    }
    return flag;
}
int generate_flag_paired(struct m_opt *opt,struct read_inf_t *read,struct read_inf_t *read_p,int cand_order,int read_order)
{
    int flag = 0x0001;

    if(read->cand[cand_order].cigar[0]=='\0')
    {
        flag |= 0x0004;
        if(read_p->cand_num==0) flag |= 0x0008;
        return flag;
    }
    else
    {
        if(read->cand[cand_order].pair!=-1) flag|=0x0002;
        //if(read_p->cand_num==0) flag |= 0x0008;
        if(read->cand[cand_order].strand==1) flag|=0x0010;
        if((read_p->cand[read->cand[cand_order].pair].strand==1)&&(read->cand[cand_order].pair!=-1)) flag|=0x0020;
        if(read_order==1) flag|=0x0040;
        else flag|=0x0080;
    }
    return flag;
}
void write_TAG(struct read_map_inf *read)
{
    int i,j;
    struct cigar_t cigar[400];
    int cigar_num = 0;
    char temp[40];

    int cigar_l = 0;
    int cigar_n = 0;

    unsigned int ref_pos = read->start_site;

    if(read->cigar[0]!='*')
    {
        j = 0;
        cigar_num = 0;
        for(i = 0;i<strlen(read->cigar);i++)
        {
            if((read->cigar[i]>='0')&&(read->cigar[i]<='9'))
            {
                temp[j] = read->cigar[i];
                j++;
            }
            else
            {
                temp[j] = '\0';
                cigar[cigar_num].c = read->cigar[i];
                cigar[cigar_num].l = atoi(temp);
                cigar_num++;
                j = 0;
            }
        }
        read->cigar[0]='\0';
        for(i = 0;i<cigar_num;i++)
        {
            if((cigar[i].c=='M')||(cigar[i].c=='X'))
            {
                cigar_l += cigar[i].l;
                if (cigar[i].c=='M')
                {
                    if((i+1<cigar_num)&&((cigar[i+1].c=='N')||(cigar[i+1].c=='I')||(cigar[i+1].c=='M'))) {cigar_n += cigar[i].l;continue;}
                    if(cigar_n!=0) {sprintf(temp,"%d",cigar[i].l+cigar_n);cigar_n = 0;}
                    else sprintf(temp,"%d",cigar[i].l);
                    strcat(read->TAG,temp);
                }
                else
                {
                    read->dis+=cigar[i].l;
                    if(cigar_n!=0) {sprintf(temp,"%d",cigar_n);strcat(read->TAG,temp);cigar_n = 0;}
                    for(j = 0;j<cigar[i].l;j++)
                    {
                        sprintf(temp,"%c",opt->chr->list[read->chr_order].seq[j+ref_pos]);
                        strcat(read->TAG,temp);
                    }
                }
                ref_pos+=cigar[i].l;
            }
            else
            {
                if(cigar_l!=0)
                {
                    sprintf(temp,"%dM",cigar_l);
                    strcat(read->cigar,temp);
                }
                sprintf(temp,"%d%c",cigar[i].l,cigar[i].c);
                strcat(read->cigar,temp);
                cigar_l = 0;
                if(cigar[i].c=='I')
                {
                    read->dis+=cigar[i].l;
                }
                if(cigar[i].c=='D')
                {
                    if(cigar_n!=0) {sprintf(temp,"%d",cigar_n);strcat(read->TAG,temp);cigar_n = 0;}
                    read->dis+=cigar[i].l;
                    strcat(read->TAG,"^");
                    for(j = 0;j<cigar[i].l;j++)
                    {
                        sprintf(temp,"%c",opt->chr->list[read->chr_order].seq[j+ref_pos]);
                        strcat(read->TAG,temp);
                    }
                    ref_pos+=cigar[i].l;
                }
                else if (cigar[i].c=='N') ref_pos+=cigar[i].l;
            }
        }
        if(cigar_l!=0)
        {
            sprintf(temp,"%dM",cigar_l);
            strcat(read->cigar,temp);
        }
    }
}
void write_MapQ(struct read_map_inf *read,int num)
{
    if (num >= 5) read->MapQ = 0;
    else if (num == 4) read->MapQ = 0;
    else if (num == 3) read->MapQ = 0;
    else if (num == 2) read->MapQ = 0;
    else read->MapQ = 50;
}
void process_cigar(char *cigarBuf,struct cigar_t *cigar,int *cigar_total,unsigned int *end_pos)
{
    int i,j;
    char temp[20];
    {
        j = 0;
        (*cigar_total) = 0;
        for(i = 0;i<strlen(cigarBuf);i++)
        {
            if((cigarBuf[i]>='0')&&(cigarBuf[i]<='9'))
            {
                temp[j] = cigarBuf[i];
                j++;
            }
            else
            {
                temp[j] = '\0';
                cigar[(*cigar_total)].c = cigarBuf[i];
                cigar[(*cigar_total)].l = atoi(temp);

                if((cigar[(*cigar_total)].c=='M')||(cigar[(*cigar_total)].c=='X')||(cigar[(*cigar_total)].c=='D')||(cigar[(*cigar_total)].c=='N')) (*end_pos)+= cigar[(*cigar_total)].l;
                (*cigar_total)++;
                j = 0;
            }
        }
    }
}

