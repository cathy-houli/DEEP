#include "main.h"
static pthread_mutex_t ReadLock,OutputLock,UpdateLock;

int read_file(struct read_inf_t *read,struct m_opt *opt,FILE *file_1,FILE *file_2)
{
    int i = 0,j = 0;
    char f_line[MAX_STRING_LENGTH];

    if (opt->input_mode==FA_FILE)
    {
        while((i<READ_BUF_LENGTH)&&fgets(f_line,MAX_STRING_LENGTH,file_1)!=NULL)
        {
            if(f_line[0] == '>')
            {
                sscanf(f_line,">%s",read[i].name);
                strcat(read[i].name,"a\0");
                if(fgets(f_line,MAX_STRING_LENGTH,file_1)!=NULL)
                {
                    read[i].length = strlen(f_line)-1;
                    strncpy(read[i].seq,f_line,read[i].length);
                    read[i].seq[read[i].length] = '\0';
                    read[i].qual[0] = '*';
                    read[i].qual[1] = '\0';
                    read[i].rqual[0] = '*';
                    read[i].rqual[1] = '\0';
                    read[i].cand_num = 0;
                    read[i].out_flag = 0;

                    for(j = 0;j<read[i].length;j++)
                    {
                        if((read[i].seq[read[i].length-j-1]=='A')||(read[i].seq[read[i].length-j-1]=='a')) read[i].rseq[j] = 'T';
                        if((read[i].seq[read[i].length-j-1]=='C')||(read[i].seq[read[i].length-j-1]=='c')) read[i].rseq[j] = 'G';
                        if((read[i].seq[read[i].length-j-1]=='G')||(read[i].seq[read[i].length-j-1]=='g')) read[i].rseq[j] = 'C';
                        if((read[i].seq[read[i].length-j-1]=='T')||(read[i].seq[read[i].length-j-1]=='t')) read[i].rseq[j] = 'A';
                    }
                    i++;
                }
                else return 0;
                if(opt->pair==1)
                {
                if(fgets(f_line,MAX_STRING_LENGTH,file_2)!=NULL)
                {
                    if(f_line[0] == '>') {sscanf(f_line,">%s",read[i].name);strcat(read[i].name,"b\0");}
                    else return 0;
                    if(fgets(f_line,MAX_STRING_LENGTH,file_2)!=NULL)
                    {
                        read[i].length = strlen(f_line)-1;
                        strncpy(read[i].seq,f_line,read[i].length);
                        read[i].seq[read[i].length] = '\0';
                        read[i].qual[0] = '*';
                        read[i].qual[1] = '\0';
                        read[i].rqual[0] = '*';
                        read[i].rqual[1] = '\0';
                        read[i].cand_num = 0;
                        read[i].out_flag = 0;

                        for(j = 0;j<read[i].length;j++)
                        {
                            if((read[i].seq[read[i].length-j-1]=='A')||(read[i].seq[read[i].length-j-1]=='a')) read[i].rseq[j] = 'T';
                            if((read[i].seq[read[i].length-j-1]=='C')||(read[i].seq[read[i].length-j-1]=='c')) read[i].rseq[j] = 'G';
                            if((read[i].seq[read[i].length-j-1]=='G')||(read[i].seq[read[i].length-j-1]=='g')) read[i].rseq[j] = 'C';
                            if((read[i].seq[read[i].length-j-1]=='T')||(read[i].seq[read[i].length-j-1]=='t')) read[i].rseq[j] = 'A';
                        }
                        i++;
                    }
                    else return 0;
                }
                else return 0;
                }
            }
        }
    }
    else if (opt->input_mode==FQ_FILE)
    {
        while((i<READ_BUF_LENGTH)&&fgets(f_line,MAX_STRING_LENGTH,file_1)!=NULL)
        {
            if(f_line[0] == '@')
            {
                sscanf(f_line,"@%s",read[i].name);
                strcat(read[i].name,"a\0");
                if(fgets(f_line,MAX_STRING_LENGTH,file_1)!=NULL)
                {
                    read[i].length = strlen(f_line)-1;
                    strncpy(read[i].seq,f_line,read[i].length);
                    read[i].seq[read[i].length] = '\0';
                }else return 0;
                if(fgets(f_line,MAX_STRING_LENGTH,file_1)==NULL) return 0;
                if(fgets(f_line,MAX_STRING_LENGTH,file_1)!=NULL)
                {
                    strncpy(read[i].qual,f_line,read[i].length);
                    read[i].qual[read[i].length] = '\0';
                    read[i].cand_num = 0;
                    read[i].out_flag = 0;

                    for(j = 0;j<read[i].length;j++)
                    {
                        if((read[i].seq[read[i].length-j-1]=='A')||(read[i].seq[read[i].length-j-1]=='a')) read[i].rseq[j] = 'T';
                        if((read[i].seq[read[i].length-j-1]=='C')||(read[i].seq[read[i].length-j-1]=='c')) read[i].rseq[j] = 'G';
                        if((read[i].seq[read[i].length-j-1]=='G')||(read[i].seq[read[i].length-j-1]=='g')) read[i].rseq[j] = 'C';
                        if((read[i].seq[read[i].length-j-1]=='T')||(read[i].seq[read[i].length-j-1]=='t')) read[i].rseq[j] = 'A';
                        read[i].rqual[j] = read[i].qual[read[i].length-j-1];
                    }
                    read[i].rseq[read[i].length] = '\0';
                    read[i].rqual[read[i].length] = '\0';
                    i++;
                }
                else return 0;
                if(opt->pair==1)
                {
                if(fgets(f_line,MAX_STRING_LENGTH,file_2)!=NULL)
                {
                    if(f_line[0] == '@') {sscanf(f_line,"@%s",read[i].name);strcat(read[i].name,"b\0");}
                    else return 0;
                    if(fgets(f_line,MAX_STRING_LENGTH,file_2)!=NULL)
                    {
                        read[i].length = strlen(f_line)-1;
                        strncpy(read[i].seq,f_line,read[i].length);
                        read[i].seq[read[i].length] = '\0';
                    }else return 0;
                    if(fgets(f_line,MAX_STRING_LENGTH,file_2)==NULL) return 0;
                    if(fgets(f_line,MAX_STRING_LENGTH,file_2)!=NULL)
                    {
                        strncpy(read[i].qual,f_line,read[i].length);
                        read[i].qual[read[i].length] = '\0';
                        read[i].cand_num = 0;
                        read[i].out_flag = 0;

                        for(j = 0;j<read[i].length;j++)
                        {
                            if((read[i].seq[read[i].length-j-1]=='A')||(read[i].seq[read[i].length-j-1]=='a')) read[i].rseq[j] = 'T';
                            if((read[i].seq[read[i].length-j-1]=='C')||(read[i].seq[read[i].length-j-1]=='c')) read[i].rseq[j] = 'G';
                            if((read[i].seq[read[i].length-j-1]=='G')||(read[i].seq[read[i].length-j-1]=='g')) read[i].rseq[j] = 'C';
                            if((read[i].seq[read[i].length-j-1]=='T')||(read[i].seq[read[i].length-j-1]=='t')) read[i].rseq[j] = 'A';
                            read[i].rqual[j] = read[i].qual[read[i].length-j-1];
                        }
                        read[i].rseq[read[i].length] = '\0';
                        read[i].rqual[read[i].length] = '\0';
                        i++;
                    }
                    else return 0;
                }
                else return 0;
                }
            }
        }
    }
    else return 0;

    return i;
}

void out_put_read(struct m_opt *opt,struct read_inf_t *read,int read_num,FILE *out,FILE *un)
{
    int i = 0,j = 0,x = 0;
    char TAG[MAX_STRING_LENGTH];
    int flag = 0;
    int max_s = 0;//,max_p = 0;
    int pair_flag = 0;

    char name1[MAX_NAME_LENGTH];
    char name2[MAX_NAME_LENGTH];
    int length = 0;

    opt->total_read+=read_num;
    for (i = 0;i<read_num;i++)
    {
        if(opt->pair)
        {
            if(read[i].out_flag==0) opt->unmapped_read++;
            if(read[i+1].out_flag==0) opt->unmapped_read++;

            if((read[i].out_flag==0)&&(read[i+1].out_flag==0))
            {
                if (read[i].cand_num==0)
                    {
                        flag = generate_flag_paired(opt,read+i,read+i+1,0,1);
                        fprintf(un,"%s\t%d\t*\t0\t0\t*\t%s\t0\t0\t%s\t%s\t0\n",read[i].name,flag,read[i+1].name,read[i].seq,read[i].qual);
                    }
                    else
                    {
                        for(x = 0;x<read[i].cand_num;x++)
                        {
                            if((x>5)&&(read[i].cand[x].score<80)) continue;
                            flag = generate_flag_paired(opt,read+i,read+i+1,x,1);

                            fprintf(un,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%d\n",
                                read[i].name,flag,opt->chr->list[read[i].cand[x].chr_order].name,read[i].cand[x].pos,read[i].MapQ,read[i].cand[x].cigar,
                                read[i+1].name,(read[i].cand[x].pair==-1)?0:read[i+1].cand[read[i].cand[x].pair].pos-opt->chr->list[read[i].cand[x].chr_order].start_site+1,
                                (read[i].cand[x].strand==0)?read[i].seq:read[i].rseq,(read[i].cand[x].strand==0)?read[i].qual:read[i].rqual,read[i].cand[x].score);
                        }
                    }
                if (read[i+1].cand_num==0)
                    {
                        flag = generate_flag_paired(opt,read+i+1,read+i,0,2);
                        fprintf(un,"%s\t%d\t*\t0\t0\t*\t%s\t0\t0\t%s\t%s\t0\n",read[i+1].name,flag,read[i].name,read[i+1].seq,read[i+1].qual);
                    }
                    else
                    {
                        for(x = 0;x<read[i+1].cand_num;x++)
                        {
                            if((x>5)&&(read[i+1].cand[x].score<80)) continue;
                            flag = generate_flag_paired(opt,read+i+1,read+i,x,2);

                            fprintf(un,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%d\n",
                                read[i+1].name,flag,opt->chr->list[read[i+1].cand[x].chr_order].name,read[i+1].cand[x].pos,read[i+1].MapQ,read[i+1].cand[x].cigar,
                                read[i].name,(read[i+1].cand[x].pair==-1)?0:read[i].cand[read[i+1].cand[x].pair].pos-opt->chr->list[read[i+1].cand[x].chr_order].start_site+1,
                                (read[i+1].cand[x].strand==0)?read[i+1].seq:read[i+1].rseq,(read[i+1].cand[x].strand==0)?read[i+1].qual:read[i+1].rqual,read[i+1].cand[x].score);
                        }
                    }
            }
            else if((read[i].out_flag==0)||(read[i+1].out_flag==0))
            {
                if(read[i].out_flag!=0)
                {
                    read[i].MapQ = generate_MapQ(read+i);
                    length = strlen(read[i].name);
                    strncpy(name1,read[i].name,length-1);
                    name1[length-1] = '\0';
                    length = strlen(read[i+1].name);
                    strncpy(name2,read[i+1].name,length-1);
                    name2[length-1] = '\0';
                    {
                        pair_flag = 0;
                        max_s = 0;
                        for(x = 0;x< read[i].cand_num;x++)
                        {
                            if(read[i].cand[x].out_put_flag==0) continue;
                            if(read[i].cand[x].pair!=-1) pair_flag=1;

                            if((read[i].cand[x].pair!=-1)&&(read[i].cand[x].score+2>max_s))
                                max_s = read[i].cand[x].score+2;
                            else if(read[i].cand[x].score>max_s)
                                max_s = read[i].cand[x].score;
                        }
                        for(x = 0;x< read[i].cand_num;x++)
                        {
                            if(read[i].cand[x].out_put_flag==0) continue;
                            if(((read[i].cand[x].pair!=-1)&&(read[i].cand[x].score+2==max_s))||(read[i].cand[x].score==max_s))
                            //if(read[i].cand[x].score==max_s)
                            {
                                flag = generate_flag_paired(opt,read+i,read+i+1,x,1);
                                sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read[i].cand[x].score,read[i].cand[x].dis,read[i].cand[x].TAG);
                            //+MD
                                fprintf(out,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t0\t0\t%s\t%s\t%s\n",
                                name1,flag,opt->chr->list[read[i].cand[x].chr_order].name,read[i].cand[x].pos-opt->chr->list[read[i].cand[x].chr_order].start_site+1,read[i].MapQ,read[i].cand[x].cigar,
                                name2,(read[i].cand[x].strand==0)?read[i].seq:read[i].rseq,(read[i].cand[x].strand==0)?read[i].qual:read[i].rqual,TAG);
                            }
                            else read[i].cand[x].out_put_flag = 0;
                        }
                    }
                    if (read[i+1].cand_num==0)
                    {
                        flag = generate_flag_paired(opt,read+i+1,read+i,0,2);
                        fprintf(un,"%s\t%d\t*\t0\t0\t*\t%s\t0\t0\t%s\t%s\t0\n",read[i+1].name,flag,read[i].name,read[i+1].seq,read[i+1].qual);
                    }
                    else
                    {
                        for(x = 0;x<read[i+1].cand_num;x++)
                        {
                            if((pair_flag==1)&&((read[i+1].cand[x].pair==-1)||(read[i].cand[read[i+1].cand[x].pair].out_put_flag!=1))) continue;
                            if((x>10)&&(read[i+1].cand[x].score<50)) continue;
                            flag = generate_flag_paired(opt,read+i+1,read+i,x,2);

                            fprintf(un,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%d\n",
                                read[i+1].name,flag,opt->chr->list[read[i+1].cand[x].chr_order].name,read[i+1].cand[x].pos,read[i+1].MapQ,read[i+1].cand[x].cigar,
                                read[i].name,(read[i+1].cand[x].pair==-1)?0:read[i].cand[read[i+1].cand[x].pair].pos-opt->chr->list[read[i+1].cand[x].chr_order].start_site+1,
                                (read[i+1].cand[x].strand==0)?read[i+1].seq:read[i+1].rseq,(read[i+1].cand[x].strand==0)?read[i+1].qual:read[i+1].rqual,read[i+1].cand[x].score);
                        }
                    }
                }
                if(read[i+1].out_flag!=0)
                {
                    read[i+1].MapQ = generate_MapQ(read+i+1);
                    length = strlen(read[i].name);
                    strncpy(name1,read[i].name,length-1);
                    name1[length-1] = '\0';
                    length = strlen(read[i+1].name);
                    strncpy(name2,read[i+1].name,length-1);
                    name2[length-1] = '\0';
                    {
                        max_s = 0;
                        pair_flag = 0;
                        for(x = 0;x< read[i+1].cand_num;x++)
                        {
                            if(read[i+1].cand[x].out_put_flag==0) continue;
                            if(read[i+1].cand[x].pair!=-1) pair_flag = 1;

                            if((read[i+1].cand[x].pair!=-1)&&(read[i+1].cand[x].score+2>max_s))
                                max_s = read[i+1].cand[x].score+2;
                            else if(read[i+1].cand[x].score>max_s)
                                max_s = read[i+1].cand[x].score;
                        }
                        for(x = 0;x< read[i+1].cand_num;x++)
                        {
                            if(read[i+1].cand[x].out_put_flag==0) continue;
                            //if(read[i+1].cand[x].score==max_s)
                            if(((read[i+1].cand[x].pair!=-1)&&(read[i+1].cand[x].score+2==max_s))||(read[i+1].cand[x].score==max_s))
                            {
                                flag = generate_flag_paired(opt,read+i+1,read+i,x,2);
                                sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read[i+1].cand[x].score,read[i+1].cand[x].dis,read[i+1].cand[x].TAG);
                            //+MD
                                fprintf(out,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t0\t0\t%s\t%s\t%s\n",
                                name2,flag,opt->chr->list[read[i+1].cand[x].chr_order].name,read[i+1].cand[x].pos-opt->chr->list[read[i+1].cand[x].chr_order].start_site+1,read[i+1].MapQ,read[i+1].cand[x].cigar,
                                name1,(read[i+1].cand[x].strand==0)?read[i+1].seq:read[i+1].rseq,(read[i+1].cand[x].strand==0)?read[i+1].qual:read[i+1].rqual,TAG);
                            }
                            else read[i+1].cand[x].out_put_flag = 0;
                        }
                    }
                    if (read[i].cand_num==0)
                    {
                        flag = generate_flag_paired(opt,read+i,read+i+1,0,1);
                        fprintf(un,"%s\t%d\t*\t0\t0\t*\t%s\t0\t0\t%s\t%s\t0\n",read[i].name,flag,read[i+1].name,read[i].seq,read[i].qual);
                    }
                    else
                    {
                        for(x = 0;x<read[i].cand_num;x++)
                        {
                            if((pair_flag==1)&&((read[i].cand[x].pair==-1)||(read[i+1].cand[read[i].cand[x].pair].out_put_flag!=1))) continue;
                            if((x>10)&&(read[i].cand[x].score<50)) continue;
                            flag = generate_flag_paired(opt,read+i,read+i+1,x,1);

                            fprintf(un,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%d\n",
                                read[i].name,flag,opt->chr->list[read[i].cand[x].chr_order].name,read[i].cand[x].pos,read[i].MapQ,read[i].cand[x].cigar,
                                read[i+1].name,(read[i].cand[x].pair==-1)?0:read[i+1].cand[read[i].cand[x].pair].pos-opt->chr->list[read[i].cand[x].chr_order].start_site+1,
                                (read[i].cand[x].strand==0)?read[i].seq:read[i].rseq,(read[i].cand[x].strand==0)?read[i].qual:read[i].rqual,read[i].cand[x].score);
                        }
                    }
                }
            }
            else
            {
                read[i].MapQ = generate_MapQ(read+i);
                read[i+1].MapQ = generate_MapQ(read+i+1);

                length = strlen(read[i].name);
                strncpy(name1,read[i].name,length-1);
                name1[length-1] = '\0';
                length = strlen(read[i+1].name);
                strncpy(name2,read[i+1].name,length-1);
                name2[length-1] = '\0';

                flag = 0;
                max_s = 0;
                //max_p = 0;
                for(x = 0;x<read[i].cand_num;x++)
                {
                    if((read[i].cand[x].out_put_flag==1)&&(read[i].cand[x].pair!=-1))flag = 1;
                    if(read[i].cand[x].score>max_s) max_s=read[i].cand[x].score;
                    //if((read[i].cand[x].out_put_flag==1)&&(read[i].cand[x].pair!=-1))
                    //{
                        //if(read[i].cand[x].score+read[i].cand[read[i].cand[x].pair].score>max_p)
                            //max_p = read[i].cand[x].score+read[i].cand[read[i].cand[x].pair].score;
                    //}
                }
                for(x = 0;x<read[i+1].cand_num;x++)
                {
                    if((read[i+1].cand[x].out_put_flag==1)&&(read[i+1].cand[x].pair!=-1)){flag = 1;break;}
                }
                if(flag)
                {
                for(x = 0;x<read[i].cand_num;x++)
                {
                    if((read[i].cand[x].out_put_flag==1)&&(read[i].cand[x].pair!=-1))
                    {
                        if(read[i].cand[x].score<read[i].cand[read[i+1].cand[read[i].cand[x].pair].pair].score) continue;
                        //if(read[i].cand[x].score+read[i+1].cand[read[i].cand[x].pair].score<max_p-10) continue;
                        flag = generate_flag_paired(opt,read+i,read+i+1,x,1);
                        sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read[i].cand[x].score,read[i].cand[x].dis,read[i].cand[x].TAG);
                        //+MD
                        fprintf(out,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%d\t%s\t%s\t%s\n",
                        name1,flag,opt->chr->list[read[i].cand[x].chr_order].name,read[i].cand[x].pos-opt->chr->list[read[i].cand[x].chr_order].start_site+1,read[i].MapQ,read[i].cand[x].cigar,
                        name2,(read[i].cand[x].pair==-1)?0:read[i+1].cand[read[i].cand[x].pair].pos-opt->chr->list[read[i].cand[x].chr_order].start_site+1,
                        (read[i].cand[x].pos>read[i+1].cand[read[i].cand[x].pair].pos)?read[i].cand[x].pos-read[i+1].cand[read[i].cand[x].pair].pos:read[i+1].cand[read[i].cand[x].pair].pos-read[i].cand[x].pos,
                        (read[i].cand[x].strand==0)?read[i].seq:read[i].rseq,(read[i].cand[x].strand==0)?read[i].qual:read[i].rqual,TAG);
                    }
                    //else if((read[i].cand[x].out_put_flag==1)&&(read[i].cand[x].score==max_s))
                    //{
                        //flag = generate_flag_paired(opt,read+i,read+i+1,x,1);
                        //sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read[i].cand[x].score,read[i].cand[x].dis,read[i].cand[x].TAG);
                        //+MD
                        //fprintf(out,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%s\n",
                        //read[i].name,flag,opt->chr->list[read[i].cand[x].chr_order].name,read[i].cand[x].pos-opt->chr->list[read[i].cand[x].chr_order].start_site+1,read[i].MapQ,read[i].cand[x].cigar,
                        //read[i+1].name,(read[i].cand[x].pair==-1)?0:read[i+1].cand[read[i].cand[x].pair].pos-opt->chr->list[read[i].cand[x].chr_order].start_site+1,
                        //(read[i].cand[x].strand==0)?read[i].seq:read[i].rseq,(read[i].cand[x].strand==0)?read[i].qual:read[i].rqual,TAG);
                    //}
                    else if((read[i].cand[x].out_put_flag==0)&&(read[i].cand[x].pair!=-1))
                    {
                        if(read[i].cand[x].score<read[i].cand[read[i+1].cand[read[i].cand[x].pair].pair].score) continue;
                        fprintf(un,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%d\n",
                        read[i].name,flag,opt->chr->list[read[i].cand[x].chr_order].name,read[i].cand[x].pos,read[i].MapQ,read[i].cand[x].cigar,
                        read[i+1].name,read[i+1].cand[read[i].cand[x].pair].pos-opt->chr->list[read[i].cand[x].chr_order].start_site+1,
                        (read[i].cand[x].strand==0)?read[i].seq:read[i].rseq,(read[i].cand[x].strand==0)?read[i].qual:read[i].rqual,read[i].cand[x].score);
                    }
                }
                max_s = 0;
                for(x = 0;x<read[i+1].cand_num;x++)
                {
                    if(read[i+1].cand[x].score>max_s) max_s=read[i+1].cand[x].score;
                }
                for(x = 0;x<read[i+1].cand_num;x++)
                {
                    if((read[i+1].cand[x].out_put_flag==1)&&(read[i+1].cand[x].pair!=-1))
                    {
                        if(read[i+1].cand[x].score<read[i+1].cand[read[i].cand[read[i+1].cand[x].pair].pair].score) continue;
                        //if(read[i+1].cand[x].score+read[i].cand[read[i+1].cand[x].pair].score<max_p-10) continue;
                        flag = generate_flag_paired(opt,read+i+1,read+i,x,2);
                        sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read[i+1].cand[x].score,read[i+1].cand[x].dis,read[i+1].cand[x].TAG);
                        //+MD
                        fprintf(out,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%d\t%s\t%s\t%s\n",
                        name2,flag,opt->chr->list[read[i+1].cand[x].chr_order].name,read[i+1].cand[x].pos-opt->chr->list[read[i+1].cand[x].chr_order].start_site+1,read[i+1].MapQ,read[i+1].cand[x].cigar,
                        name1,(read[i+1].cand[x].pair==-1)?0:read[i].cand[read[i+1].cand[x].pair].pos-opt->chr->list[read[i+1].cand[x].chr_order].start_site+1,
                        (read[i].cand[read[i+1].cand[x].pair].pos>read[i+1].cand[x].pos)?read[i].cand[read[i+1].cand[x].pair].pos-read[i+1].cand[x].pos:read[i+1].cand[x].pos-read[i].cand[read[i+1].cand[x].pair].pos,
                        (read[i+1].cand[x].strand==0)?read[i+1].seq:read[i+1].rseq,(read[i+1].cand[x].strand==0)?read[i+1].qual:read[i+1].rqual,TAG);
                    }
                    //else if((read[i+1].cand[x].out_put_flag==1)&&(read[i+1].cand[x].score==max_s)&&(read[i+1].cand[x].pair==-1))
                    //{
                        //flag = generate_flag_paired(opt,read+i+1,read+i,x,2);
                        //sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read[i+1].cand[x].score,read[i+1].cand[x].dis,read[i+1].cand[x].TAG);
                        //+MD
                        //fprintf(out,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%s\n",
                        //read[i+1].name,flag,opt->chr->list[read[i+1].cand[x].chr_order].name,read[i+1].cand[x].pos-opt->chr->list[read[i+1].cand[x].chr_order].start_site+1,read[i+1].MapQ,read[i+1].cand[x].cigar,
                        //read[i].name,(read[i+1].cand[x].pair==-1)?0:read[i].cand[read[i+1].cand[x].pair].pos-opt->chr->list[read[i+1].cand[x].chr_order].start_site+1,
                        //(read[i+1].cand[x].strand==0)?read[i+1].seq:read[i+1].rseq,(read[i+1].cand[x].strand==0)?read[i+1].qual:read[i+1].rqual,TAG);
                    //}
                    else if((read[i+1].cand[x].out_put_flag==0)&&(read[i+1].cand[x].pair!=-1))
                    {
                        if(read[i+1].cand[x].score<read[i+1].cand[read[i].cand[read[i+1].cand[x].pair].pair].score) continue;
                        flag = generate_flag_paired(opt,read+i+1,read+i,x,2);

                        fprintf(un,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%d\n",
                        read[i+1].name,flag,opt->chr->list[read[i+1].cand[x].chr_order].name,read[i+1].cand[x].pos,read[i+1].MapQ,read[i+1].cand[x].cigar,
                        read[i].name,read[i].cand[read[i+1].cand[x].pair].pos-opt->chr->list[read[i+1].cand[x].chr_order].start_site+1,
                        (read[i+1].cand[x].strand==0)?read[i+1].seq:read[i+1].rseq,(read[i+1].cand[x].strand==0)?read[i+1].qual:read[i+1].rqual,read[i+1].cand[x].score);
                    }
                }
                }
                else
                {
                for(x = 0;x<read[i].cand_num;x++)
                {
                    if((read[i].cand[x].out_put_flag==1)&&(read[i].cand[x].score==max_s)&&(read[i].cand[x].pair!=-1))
                    {
                        flag = generate_flag_paired(opt,read+i,read+i+1,x,1);
                        sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read[i].cand[x].score,read[i].cand[x].dis,read[i].cand[x].TAG);
                        //+MD
                        fprintf(out,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%d\t%s\t%s\t%s\n",
                        name1,flag,opt->chr->list[read[i].cand[x].chr_order].name,read[i].cand[x].pos-opt->chr->list[read[i].cand[x].chr_order].start_site+1,read[i].MapQ,read[i].cand[x].cigar,
                        name2,(read[i].cand[x].pair==-1)?0:read[i+1].cand[read[i].cand[x].pair].pos-opt->chr->list[read[i].cand[x].chr_order].start_site+1,
                        (read[i].cand[x].pos>read[i+1].cand[read[i].cand[x].pair].pos)?read[i].cand[x].pos-read[i+1].cand[read[i].cand[x].pair].pos:read[i+1].cand[read[i].cand[x].pair].pos-read[i].cand[x].pos,
                        (read[i].cand[x].strand==0)?read[i].seq:read[i].rseq,(read[i].cand[x].strand==0)?read[i].qual:read[i].rqual,TAG);
                    }
                    else if((read[i].cand[x].out_put_flag==1)&&(read[i].cand[x].score==max_s)&&(read[i].cand[x].pair==-1))
                    {
                        flag = generate_flag_paired(opt,read+i,read+i+1,x,1);
                        sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read[i].cand[x].score,read[i].cand[x].dis,read[i].cand[x].TAG);
                        //+MD
                        fprintf(out,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%s\n",
                        name1,flag,opt->chr->list[read[i].cand[x].chr_order].name,read[i].cand[x].pos-opt->chr->list[read[i].cand[x].chr_order].start_site+1,read[i].MapQ,read[i].cand[x].cigar,
                        name2,(read[i].cand[x].pair==-1)?0:read[i+1].cand[read[i].cand[x].pair].pos-opt->chr->list[read[i].cand[x].chr_order].start_site+1,
                        (read[i].cand[x].strand==0)?read[i].seq:read[i].rseq,(read[i].cand[x].strand==0)?read[i].qual:read[i].rqual,TAG);
                    }
                }
                max_s = 0;
                for(x = 0;x<read[i+1].cand_num;x++)
                {
                    if(read[i+1].cand[x].score>max_s) max_s=read[i+1].cand[x].score;
                }
                for(x = 0;x<read[i+1].cand_num;x++)
                {
                    if((read[i+1].cand[x].out_put_flag==1)&&(read[i+1].cand[x].score==max_s)&&(read[i+1].cand[x].pair!=-1))
                    {
                        flag = generate_flag_paired(opt,read+i+1,read+i,x,2);
                        sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read[i+1].cand[x].score,read[i+1].cand[x].dis,read[i+1].cand[x].TAG);
                        //+MD
                        fprintf(out,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%d\t%s\t%s\t%s\n",
                        name2,flag,opt->chr->list[read[i+1].cand[x].chr_order].name,read[i+1].cand[x].pos-opt->chr->list[read[i+1].cand[x].chr_order].start_site+1,read[i+1].MapQ,read[i+1].cand[x].cigar,
                        name1,(read[i+1].cand[x].pair==-1)?0:read[i].cand[read[i+1].cand[x].pair].pos-opt->chr->list[read[i+1].cand[x].chr_order].start_site+1,
                        (read[i].cand[read[i+1].cand[x].pair].pos>read[i+1].cand[x].pos)?read[i].cand[read[i+1].cand[x].pair].pos-read[i+1].cand[x].pos:read[i+1].cand[x].pos-read[i].cand[read[i+1].cand[x].pair].pos,
                        (read[i+1].cand[x].strand==0)?read[i+1].seq:read[i+1].rseq,(read[i+1].cand[x].strand==0)?read[i+1].qual:read[i+1].rqual,TAG);
                    }
                    else if((read[i+1].cand[x].out_put_flag==1)&&(read[i+1].cand[x].score==max_s)&&(read[i+1].cand[x].pair==-1))
                    {
                        flag = generate_flag_paired(opt,read+i+1,read+i,x,2);
                        sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read[i+1].cand[x].score,read[i+1].cand[x].dis,read[i+1].cand[x].TAG);
                        //+MD
                        fprintf(out,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%s\n",
                        name2,flag,opt->chr->list[read[i+1].cand[x].chr_order].name,read[i+1].cand[x].pos-opt->chr->list[read[i+1].cand[x].chr_order].start_site+1,read[i+1].MapQ,read[i+1].cand[x].cigar,
                        name1,(read[i+1].cand[x].pair==-1)?0:read[i].cand[read[i+1].cand[x].pair].pos-opt->chr->list[read[i+1].cand[x].chr_order].start_site+1,
                        (read[i+1].cand[x].strand==0)?read[i+1].seq:read[i+1].rseq,(read[i+1].cand[x].strand==0)?read[i+1].qual:read[i+1].rqual,TAG);
                    }
                }
                }
            }
            i++;
        }
        else
        {
            if(read[i].out_flag==0) opt->unmapped_read++;
            if(read[i].out_flag==0)
            {
                if(read[i].cand_num==0)
                {
                    flag = generate_flag_single(opt,read+i,0);
                    fprintf(un,"%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\t0\n",read[i].name,flag,read[i].seq,read[i].seq);
                }
                else
                {
                    for(x = 0;x<read[i].cand_num;x++)
                    {
                        if((x>5)&&(read[i].cand[x].score<80)) continue;
                        flag = generate_flag_single(opt,read+i,x);
                        fprintf(un,"%s\t%d\t%s\t%u\t%d\t%s\t*\t0\t0\t%s\t%s\t%d\n",
                        read[i].name,flag,opt->chr->list[read[i].cand[x].chr_order].name,read[i].cand[x].pos,read[i].MapQ,read[i].cand[x].cigar,
                        (read[i].cand[x].strand==0)?read[i].seq:read[i].rseq,(read[i].cand[x].strand==0)?read[i].qual:read[i].rqual,read[i].cand[x].score);
                    }
                }

                continue;
            }
            read[i].MapQ = generate_MapQ(read+i);
            length = strlen(read[i].name);
            read[i].name[length-1] = '\0';
            if(opt->output_mode==OUTPUT_BEST)
            {
                max_s = 0;
                for (j = 0; j < read[i].cand_num; j++)
                {
                    if (max_s < read[i].cand[j].score) max_s = read[i].cand[j].score;
                }
                for (j = 0; j < read[i].cand_num; j++)
                {
                    if (max_s < read[i].cand[j].score)
                    {
                        flag = generate_flag_single(opt,read+i,x);
                sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read[i].cand[x].score,read[i].cand[x].dis,read[i].cand[x].TAG);
                //+MD
                fprintf(out,"%s\t%d\t%s\t%u\t%d\t%s\t*\t0\t0\t%s\t%s\t%s\n",
                read[i].name,flag,opt->chr->list[read[i].cand[x].chr_order].name,read[i].cand[x].pos-opt->chr->list[read[i].cand[x].chr_order].start_site+1,read[i].MapQ,read[i].cand[x].cigar,
                (read[i].cand[x].strand==0)?read[i].seq:read[i].rseq,(read[i].cand[x].strand==0)?read[i].qual:read[i].rqual,TAG);
                    }
                }
            }
            else
            {
            for(x = 0;x<read[i].cand_num;x++)
            {
                if(read[i].cand[x].out_put_flag==0) continue;

                flag = generate_flag_single(opt,read+i,x);
                sprintf(TAG,"AS:%d\tNM:%d",read[i].cand[x].score,read[i].cand[x].dis);
                sprintf(TAG,"AS:i:%d\tNM:i:%d\tMD:Z:%s",read[i].cand[x].score,read[i].cand[x].dis,read[i].cand[x].TAG);
                //+MD
                fprintf(out,"%s\t%d\t%s\t%u\t%d\t%s\t*\t0\t0\t%s\t%s\t%s\n",
                read[i].name,flag,opt->chr->list[read[i].cand[x].chr_order].name,read[i].cand[x].pos-opt->chr->list[read[i].cand[x].chr_order].start_site+1,read[i].MapQ,read[i].cand[x].cigar,
                (read[i].cand[x].strand==0)?read[i].seq:read[i].rseq,(read[i].cand[x].strand==0)?read[i].qual:read[i].rqual,TAG);
            }
            }
        }
    }
}


struct job_seed
{
    struct m_opt *opt;
    FILE *in_file_1,*in_file_2;
    FILE *out_file_1;
    FILE *unmap_file_1;
};
int snp_cmp(const void *a,const void *b)
{
    struct snp_t *EA,*EB;
    EA = (struct snp_t *)a;
    EB = (struct snp_t *)b;

    if(EA->start==EB->start)
    {
        if (EA->end==EB->end)
        {
            if(EA->type==EB->type)
            {
                if(EA->length==EB->length)
                {
                    if(EA->seq==EB->seq)return 0;
                    else if(EA->seq<EB->seq)return -1;
                    else return 1;
                }
                else if(EA->length<EB->length)return -1;
                else return 1;
            }
            else if(EA->type<EB->type)return -1;
            else return 1;
        }
        else if (EA->end < EB->end) return -1;
        else  return 1;
    }
    else if(EA->start < EB->start) return -1;
    else return 1;

}
void Update_Snp(struct m_opt *opt,struct snp_list_t *snp)
{
    if(snp->total>=0)
    {
        pthread_mutex_lock(&UpdateLock);
        qsort(snp->snp,snp->total,sizeof(struct snp_t),snp_cmp);

        FILE *snp_file;
        char snp_file_name[MAX_NAME_LENGTH];
        sprintf(snp_file_name,"%s/%d.snp",opt->Temp_path,opt->snp_fn);
        snp_file = fopen(snp_file_name,"wb");

        if(snp_file == NULL)
        {
            fprintf(stdout,"[DEEP]cannot open file %s...\n",snp_file_name);
        }
        unsigned int order = 0;
        unsigned int i = 0;
        snp->snp[0].num--;
        for(i = 0;i<snp->total;i++)
        {
            if(snp_cmp(&(snp->snp[order]),&(snp->snp[i]))==0)
                snp->snp[order].num++;
            else
            {
                if(fwrite(&snp->snp[order],sizeof(struct snp_t),1,snp_file)!=1)
                fprintf(stdout,"[DEEP]cannot write file...\n");
                order = i;
            }
        }
        if(fwrite(&snp->snp[order],sizeof(struct snp_t),1,snp_file)!=1)
            fprintf(stdout,"[DEEP]cannot write file...\n");

        fclose(snp_file);
        fflush(snp_file);
        opt->snp_fn++;
        pthread_mutex_unlock(&UpdateLock);
    }
}
void Update_Exon(struct m_opt *opt,struct exon_array *exon)
{
    unsigned int i = 0,j = 0;
    for (i = 0;i<exon->total;i++)
    {
        for (j = exon->exon[i].start;j<=exon->exon[i].end;j++)
            opt->dep[j]++;
    }
    exon->total = 0;
}
void insert_snp(struct snp_list_t *snp,struct snp_t snp_a)
{
    memcpy(&(snp->snp[snp->total]),&snp_a,sizeof(struct snp_t));
    snp->total++;
    if(snp->total>=MAX_SNP_NUM)
    {
        pthread_mutex_lock(&UpdateLock);
        qsort(snp->snp,MAX_SNP_NUM,sizeof(struct snp_t),snp_cmp);

        FILE *snp_file;
        char snp_file_name[MAX_NAME_LENGTH];
        sprintf(snp_file_name,"%s/%d.snp",opt->Temp_path,opt->snp_fn);
        snp_file = fopen(snp_file_name,"wb");

        if(snp_file == NULL)
        {
            fprintf(stdout,"[DEEP]cannot open file %s...\n",snp_file_name);
        }

        unsigned int order = 0;
        unsigned int i = 0;
        snp->snp[0].num--;
        for(i = 0;i<MAX_SNP_NUM;i++)
        {
            if(snp_cmp(&(snp->snp[order]),&(snp->snp[i]))==0)
                snp->snp[order].num++;
            else
            {
                if(fwrite(&snp->snp[order],sizeof(struct snp_t),1,snp_file)!=1)
                    fprintf(stdout,"[DEEP]cannot write file...\n");
                order = i;
            }
        }
        if(fwrite(&snp->snp[order],sizeof(struct snp_t),1,snp_file)!=1)
            fprintf(stdout,"[DEEP]cannot write file...\n");

        fclose(snp_file);
        opt->snp_fn++;
        snp->total = 0;
        pthread_mutex_unlock(&UpdateLock);
    }
}
void insert_exon(struct exon_array *exon,struct snp_list_t *snp,unsigned int start,unsigned int end)
{
    exon->exon[exon->total].end = end;
    exon->exon[exon->total].start = start;
    exon->total++;
    if(exon->total>=READ_BUF_LENGTH*3)
    {
        pthread_mutex_lock(&UpdateLock);
        Update_Exon(opt,exon);
        pthread_mutex_unlock(&UpdateLock);
        exon->total = 0;
    }
}
void find_snp2(struct m_opt *opt,struct read_inf_t *read,struct exon_array *exon,struct snp_list_t *snp,struct cigar_t *cigar,int cigar_total,unsigned int ref_site)
{
    struct snp_seq_t snp_r;
    struct snp_t snp_a;
    unsigned int exon_start = 0,exon_end = 0;

    int read_site = 0;

    unsigned int i = 0,j = 0;

    read_site = 0;
    snp_r.snp_num = 0;
    exon_start = ref_site;
    for(i = 0;i<cigar_total;i++)
    {
        if(cigar[i].l<0)
            fprintf(stdout,"%s\n",read->name);
        switch(cigar[i].c)
        {
        case 'M':
            ref_site+=cigar[i].l;
            read_site+=cigar[i].l;
            break;
        case 'X':
            snp_a.start = ref_site;
            snp_a.end = ref_site+cigar[i].l-1;
            snp_a.type = 'X';
            snp_a.seq = 0;
            snp_a.length = cigar[i].l;
            snp_a.num = 1;
            snp_a.dep = 0;

            if(read->cand[read->cand_num].strand)
            {
                for(j = 0;j<cigar[i].l;j++)
                snp_a.seq = (snp_a.seq<<2)|nst_nt4_table[(int)read->rseq[read_site+j]];
            }
            else
            {
                for(j = 0;j<cigar[i].l;j++)
                snp_a.seq = (snp_a.seq<<2)|nst_nt4_table[(int)read->seq[read_site+j]];
            }
            insert_snp(snp,snp_a);
            if(snp_r.snp_num==0)
            {memcpy(&snp_r.snp[snp_r.snp_num],&snp_a,sizeof(struct snp_t));
            //snp_r.snp_num++;
            }
            else
            {
                memcpy(&snp_r.snp[snp_r.snp_num],&snp_a,sizeof(struct snp_t));
                //snp_r.snp_num++;
            }

            ref_site+=cigar[i].l;
            read_site+=cigar[i].l;
            break;
        case 'I':
            snp_a.start = ref_site;
            snp_a.end = ref_site;
            snp_a.type = 'I';
            snp_a.seq = 0;
            snp_a.length = cigar[i].l;
            snp_a.num = 1;
            snp_a.dep = 0;
            if(read->cand[read->cand_num].strand)
            {
                for(j = 0;j<cigar[i].l;j++)
                snp_a.seq = (snp_a.seq<<2)|nst_nt4_table[(int)read->rseq[read_site+j]];
            }
            else
            {
                for(j = 0;j<cigar[i].l;j++)
                snp_a.seq = (snp_a.seq<<2)|nst_nt4_table[(int)read->seq[read_site+j]];
            }
            insert_snp(snp,snp_a);
            if(snp_r.snp_num==0)
            {memcpy(&snp_r.snp[snp_r.snp_num],&snp_a,sizeof(struct snp_t));
            //snp_r.snp_num++;
            }
            else
            {
                memcpy(&snp_r.snp[snp_r.snp_num],&snp_a,sizeof(struct snp_t));
                //snp_r.snp_num++;
            }

            read_site+=cigar[i].l;
            break;
        case 'D':
            snp_a.start = ref_site;
            snp_a.end = ref_site+cigar[i].l-1;
            snp_a.type = 'D';
            snp_a.seq = 0;
            snp_a.length = cigar[i].l;
            snp_a.num = 1;
            snp_a.dep = 0;
            insert_snp(snp,snp_a);
            if(snp_r.snp_num==0)
            {memcpy(&snp_r.snp[snp_r.snp_num],&snp_a,sizeof(struct snp_t));
            //snp_r.snp_num++;
            }
            else
            {
                memcpy(&snp_r.snp[snp_r.snp_num],&snp_a,sizeof(struct snp_t));
                //snp_r.snp_num++;
            }

            ref_site+=cigar[i].l;
            break;
        case 'N':
            snp_a.start = ref_site;
            snp_a.end = ref_site+cigar[i].l-1;
            snp_a.type = 'N';
            snp_a.seq = 0;
            snp_a.length = cigar[i].l;
            snp_a.num = 1;
            snp_a.dep = 0;
            insert_snp(snp,snp_a);
            if(snp_r.snp_num==0)
            {memcpy(&snp_r.snp[snp_r.snp_num],&snp_a,sizeof(struct snp_t));
            //snp_r.snp_num++;
            }
            else
            {
                memcpy(&snp_r.snp[snp_r.snp_num],&snp_a,sizeof(struct snp_t));
                //snp_r.snp_num++;
            }
            exon_end = ref_site-1;
            insert_exon(exon,snp,exon_start,exon_end);
            ref_site+=cigar[i].l;
            exon_start = ref_site;
            break;
        case 'U':
            ref_site+=cigar[i].l;
            break;
        case 'S':
            read_site+=cigar[i].l;
            break;
        }
    }

    exon_end = ref_site-1;
    insert_exon(exon,snp,exon_start,exon_end);
}
int seed_cmp(const void *a,const void *b)
{
    struct seed_t *EA,*EB;
    EA = (struct seed_t *)a;
    EB = (struct seed_t *)b;

    if(EA->start==EB->start)
    {
        if (EA->pos==EB->pos) return 0;
        else if (EA->pos<EB->pos) return -1;
        else  return 1;
    }
    else
    {
        if(EA->start < EB->start) return -1;
        else return 1;
    }
}
int seed_cmp_r(const void *a,const void *b)
{
    struct seed_t *EA,*EB;
    EA = (struct seed_t *)a;
    EB = (struct seed_t *)b;

    if(EA->start==EB->start)
    {
        if (EA->pos==EB->pos) return 0;
        else if (EA->pos>EB->pos) return -1;
        else  return 1;
    }
    else
    {
        if(EA->start < EB->start) return -1;
        else return 1;
    }
}
void insert_seed(struct seed_t *seed,int *seed_num,int max,struct seed_t seed_r)
{
    int i = 0;
    int flag = 1;
    for(i = 0;i<*seed_num;i++)
    {
        if((seed[i].abs==seed_r.abs)&&(seed_r.start>=seed[i].start)&&(seed_r.start<=seed[i].start+seed[i].length))
        {
            if((seed_r.start+seed_r.length)>(seed[i].start+seed[i].length))
                seed[i].length+=seed_r.start+seed_r.length-seed[i].start-seed[i].length;
            flag = 0;
            break;
        }
    }
    if((flag)&&((*seed_num)<max))
        {memcpy(seed+(*seed_num),&seed_r,sizeof(struct seed_t));(*seed_num)++;}
}
void find_seed(struct m_opt *opt,struct read_inf_t *read,struct seed_t *seed,int *seed_num)
{
    char flag[MAX_READ_LENGTH];
    int i = 0,j = 0,x = 0;//back_step = 0;;

    int c = 0;
    bwtint_t k, l, ok, ol;
    struct seed_a seed_w[MAX_READ_LENGTH];

    for (i = 0;i<read->length;i++) flag[i] = 0;
    //flag[50] = 1;flag[100] = 1;flag[150] = 1;flag[200] = 1;flag[250] = 1;flag[300] = 1;flag[350] = 1;
    flag[read->length-1] = 1;
    flag[read->length-5] = 1;

    //struct seed_t seed_a;

    for (i = read->length-1;i>15;i--)
    {
        if(flag[i])
        {
            //bwt_find_seed;
            k = 0; l = opt->idx->bwt->seq_len;
            for (j = i;j>=0;j--)
            {
                if(nst_nt4_table[(int)read->seq[j]]>=4) {flag[i-1]=1;break;}
                c = nst_nt4_table[(int)read->seq[j]];
                bwt_2occ(opt->idx->bwt, k - 1, l, c, &ok, &ol);
                k = opt->idx->bwt->L2[c] + ok + 1;
                l = opt->idx->bwt->L2[c] + ol;
                seed_w[j].start = k;
                seed_w[j].end = l;

                if (k > l)
                {
                    if(((seed_w[j+1].end-seed_w[j+1].start+1) >SEED_CAND_NUM)||(i-j<16)) {flag[i-1]=1;break;}

                    for (x = 0;x<=seed_w[j+1].end-seed_w[j+1].start;x++)
                    {
                        seed[*seed_num].start = j+1;
                        seed[*seed_num].length = i-j;
                        seed[*seed_num].pos = bwt_sa(opt->idx->bwt,seed_w[j+1].start + x);
                        seed[*seed_num].lnum = 0;
                        seed[*seed_num].score = 0;
                        seed[*seed_num].flag = 0;
                        seed[*seed_num].abs = seed[*seed_num].pos-seed[*seed_num].start;
                        if((*seed_num)<SEED_BUF_LENGTH) (*seed_num)++;
                        //insert_seed(seed,seed_num,SEED_BUF_LENGTH,seed_a);
                    }
                    flag[j]=1;
                    flag[j-4]=1;
                    flag[j+4]=1;
                    break;
                }
            }
            if(j == -1)
            {
                if (k <= l)
                {
                    if((i == read->length-1)&&((seed_w[j+1].end-seed_w[j+1].start+1) >SEED_CAND_NUM))
                        {read->out_flag = 1;break;}
                    if(((seed_w[j+1].end-seed_w[j+1].start+1) >SEED_CAND_NUM)||(i-j<16)) {flag[i-1]=1;continue;}

                    for (x = 0;x<=seed_w[j+1].end-seed_w[j+1].start;x++)
                    {
                        seed[*seed_num].start = j+1;
                        seed[*seed_num].length = i-j;
                        seed[*seed_num].pos = bwt_sa(opt->idx->bwt,seed_w[j+1].start + x);
                        seed[*seed_num].lnum = 0;
                        seed[*seed_num].score = 0;
                        seed[*seed_num].flag = 0;
                        seed[*seed_num].abs = seed[*seed_num].pos-seed[*seed_num].start;
                        if((*seed_num)<SEED_BUF_LENGTH)(*seed_num)++;
                        //insert_seed(seed,seed_num,SEED_BUF_LENGTH,seed_a);
                    }
                }
            }
        }
    }
    qsort(seed,*seed_num,sizeof(struct seed_t),seed_cmp);
}
void find_seed_r(struct m_opt *opt,struct read_inf_t *read,struct seed_t *seed,int *seed_num)
{
    char flag[MAX_READ_LENGTH];
    int i = 0,j = 0,x = 0;//back_step = 0;;

    int c = 0;
    bwtint_t k, l, ok, ol;
    struct seed_a seed_w[MAX_READ_LENGTH];


    for (i = 0;i<read->length;i++) flag[i] = 0;
    //flag[50] = 1;flag[100] = 1;flag[150] = 1;flag[200] = 1;flag[250] = 1;flag[300] = 1;flag[350] = 1;
    flag[read->length-1] = 1;
    flag[read->length-5] = 1;

    //struct seed_t seed_a;

    for (i = read->length-1;i>15;i--)
    {
        if(flag[i])
        {
            //bwt_find_seed;
            k = 0; l = opt->idx->bwt->seq_len;
            for (j = i;j>=0;j--)
            {
                if(nst_nt4_table[(int)read->rseq[j]]>=4) {flag[i-1]=1;break;}
                c = nst_nt4_table[(int)read->rseq[j]];
                bwt_2occ(opt->idx->bwt, k - 1, l, c, &ok, &ol);
                k = opt->idx->bwt->L2[c] + ok + 1;
                l = opt->idx->bwt->L2[c] + ol;
                seed_w[j].start = k;
                seed_w[j].end = l;

                if (k > l)
                {
                    if(((seed_w[j+1].end-seed_w[j+1].start+1) >SEED_CAND_NUM)||(i-j<16)) {flag[i-1]=1;break;}

                    for (x = 0;x<=seed_w[j+1].end-seed_w[j+1].start;x++)
                    {
                        seed[*seed_num].start = read->length-1-i;
                        seed[*seed_num].length = i-j;
                        seed[*seed_num].pos = (opt->idx->bns->l_pac<<1)-seed[*seed_num].length-bwt_sa(opt->idx->bwt,seed_w[j+1].start + x);
                        seed[*seed_num].lnum = 0;
                        seed[*seed_num].score = 0;
                        seed[*seed_num].flag = 0;
                        seed[*seed_num].abs = seed[*seed_num].pos-seed[*seed_num].start;
                        if((*seed_num)<SEED_BUF_LENGTH)(*seed_num)++;
                        //insert_seed(seed,seed_num,SEED_BUF_LENGTH,seed_a);
                    }
                    flag[j]=1;
                    flag[j-4]=1;
                    flag[j+4]=1;

                    break;
                }
            }
            if(j == -1)
            {
                if (k <= l)
                {
                    if((i == read->length-1)&&((seed_w[j+1].end-seed_w[j+1].start+1) >SEED_CAND_NUM))
                        {read->out_flag = 1;break;}
                    if(((seed_w[j+1].end-seed_w[j+1].start+1) >SEED_CAND_NUM)||(i-j<16)) {flag[i-1]=1;continue;}

                    for (x = 0;x<=seed_w[j+1].end-seed_w[j+1].start;x++)
                    {
                        seed[*seed_num].start = read->length-1-i;
                        seed[*seed_num].length = i-j;
                        seed[*seed_num].pos = (opt->idx->bns->l_pac<<1)-seed[*seed_num].length-bwt_sa(opt->idx->bwt,seed_w[j+1].start + x);
                        seed[*seed_num].lnum = 0;
                        seed[*seed_num].score = 0;
                        seed[*seed_num].flag = 0;
                        seed[*seed_num].abs = seed[*seed_num].pos-seed[*seed_num].start;
                        if((*seed_num)<SEED_BUF_LENGTH)(*seed_num)++;
                        //insert_seed(seed,seed_num,SEED_BUF_LENGTH,seed_a);
                    }
                }
            }
        }
    }
    qsort(seed,*seed_num,sizeof(struct seed_t),seed_cmp);
}
void find_seed_s(struct m_opt *opt,struct read_inf_t *read,struct seed_t *seed,int *seed_num)
{
    int i = 0,j = 0,x = 0;//back_step = 0;;

    int c = 0;
    bwtint_t k, l, ok, ol;
    struct seed_a seed_w[MAX_READ_LENGTH];

    int step = 5;
    int seed_length = 22;

    //struct seed_t seed_a;

    for (i = i*step;i*step+seed_length<read->length;i++)
    {
            //bwt_find_seed;
            k = 0; l = opt->idx->bwt->seq_len;
            for (j =i*step;j<i*step+seed_length;j++)
            {
                if(nst_nt4_table[(int)read->seq[j]]>=4) break;
                c = 3-nst_nt4_table[(int)read->seq[j]];
                bwt_2occ(opt->idx->bwt, k - 1, l, c, &ok, &ol);
                k = opt->idx->bwt->L2[c] + ok + 1;
                l = opt->idx->bwt->L2[c] + ol;
                seed_w[j].start = k;
                seed_w[j].end = l;

                if (k > l)
                    break;
            }
            if(j == i*step+seed_length)
            {
                if (k <= l)
                {
                    if((seed_w[j-1].end-seed_w[j-1].start+1) >SEED_CAND_NUM) continue;

                    for (x = 0;x<=seed_w[j-1].end-seed_w[j-1].start;x++)
                    {
                        seed[*seed_num].start = i*step;
                        seed[*seed_num].length = seed_length;
                        seed[*seed_num].pos = (opt->idx->bns->l_pac<<1)-seed[*seed_num].length-bwt_sa(opt->idx->bwt,seed_w[j-1].start + x);
                        seed[*seed_num].lnum = 0;
                        seed[*seed_num].score = 0;
                        seed[*seed_num].flag = 0;
                        seed[*seed_num].abs = seed[*seed_num].pos-seed[*seed_num].start;
                        if((*seed_num)<SEED_BUF_LENGTH)(*seed_num)++;
                        //insert_seed(seed,seed_num,SEED_BUF_LENGTH,seed_a);
                    }
                }
            }
    }
    qsort(seed,*seed_num,sizeof(struct seed_t),seed_cmp);
}
int tail_align(struct m_opt *opt,int chr,char strand,uint64_t ref_pos,char *seq,int length,int direction,struct cigar_t cigar[20],int *num,int *score)
{
    char ref[MAX_READ_LENGTH];
    char text[MAX_READ_LENGTH];
    int length_r;
    char cigarBuf[4*MAX_READ_LENGTH];

    int t = 0,r = 0,k = 0;
    char change_type;
    int change_num = 0;
    int change_site = 0;
    int flag = 1;

    if(direction^strand)//to back
    {
        length_r = min(opt->chr->list[chr].length+opt->chr->list[chr].start_site-ref_pos,length+opt->change_length);
        for(k = 0;k<length_r;k++)
            ref[k] = opt->chr->list[chr].seq[ref_pos+k];
    }
    else//to front
    {
        length_r = min(ref_pos,length+opt->change_length);
        for(k = 0;k<length_r;k++)
            ref[k] = opt->chr->list[chr].seq[ref_pos-k];
    }
    ref[k] = '\0';
    if(direction)//to back
    {
        for(k = 0;k<length;k++)
        {
            if(strand) text[k] = base2char[3-nst_nt4_table[(int)seq[k]]];
            else text[k] = seq[k];
        }
    }
    else//to front
    {
        for(k = 0;k<length;k++)
        {
            if(strand) text[k] = base2char[3-nst_nt4_table[(int)seq[length-k-1]]];
            else text[k] = seq[length-k-1];
        }
    }
    text[k] = '\0';

    (*score) = 0;

    cigarBuf[0] = '\0';
    if(computeEditDistanceWithCigar(ref,length_r,text,length,6,cigarBuf,4*MAX_READ_LENGTH,0,1)!=-1)
    {
        k = 0;
        while(t<length)//只允许有一个变异
        {
            if(cigarBuf[k]!= '=')
            {
                if (change_num ==0)
                {
                    change_type = cigarBuf[k];//IDX
                    change_site = k;
                    change_num++;

                    if(cigarBuf[k]=='I') t++;
                    else if(cigarBuf[k]=='X') {t++;r++;}
                    else if(cigarBuf[k]=='D') r++;
                }
                else if(((change_site+1) == k)&&(change_type == cigarBuf[k]))
                {
                    change_site++;
                    if(cigarBuf[k]=='I') t++;
                    else if(cigarBuf[k]=='X') {t++;r++;}
                    else if(cigarBuf[k]=='D') r++;
                    k++;
                    continue;
                }
                else{flag = 0; break;}
            }
            else {r++;t++;}
            k++;
        }
        if(flag)
        {
            t = 0;r = 0;k = 0;
            *num = 0;
            cigar[0].c = cigarBuf[0];
            cigar[0].l = 0;
            while((t<length)&&(r<length_r))//只允许有一个变异
            {
                if(cigarBuf[k]== cigar[*num].c) {cigar[*num].l++;}
                else
                {
                    (*num)++;
                    cigar[*num].c = cigarBuf[k];
                    cigar[*num].l = 1;
                }

                if(cigarBuf[k]=='I') t++;
                else if(cigarBuf[k]=='X') {t++;r++;}
                else if(cigarBuf[k]=='D') r++;
                else if(cigarBuf[k]=='=') {t++;r++;}
                k++;
            }
            if (cigar[(*num)].l!=0)(*num)++;
            if(t<length)
            {
                cigar[*num].c = 'I';
                cigar[*num].l = length-t;
                (*num)++;
            }
            change_num =0;
            for(k = 0;k<(*num);k++)
            {
                if(cigar[k].c=='=') {cigar[k].c='M';(*score)+=cigar[k].l*opt->match;}
                else if(cigar[k].c=='X') {(*score)-=cigar[k].l*opt->miss;change_num++;}
                else if(cigarBuf[k]=='I') {(*score)-=opt->miss;change_num++;}
                else if(cigarBuf[k]=='D') {(*score)-=opt->gap;change_num++;}
            }
            if(change_num>1) {(*num) = 0;return 0;}
            if(!direction)
            {
                for(k = 0;k<*num/2;k++)
                {
                    change_type = cigar[k].c;
                    change_num = cigar[k].l;

                    cigar[k].c = cigar[*num-k-1].c;
                    cigar[k].l = cigar[*num-k-1].l;

                    cigar[*num-k-1].c = change_type;
                    cigar[*num-k-1].l = change_num;

                }
            }
            for(k = 0;k<(*num);k++)
            {
                if((cigar[k].c=='X')||(cigar[k].c=='I'))
                {
                    if(cigar[k].l>opt->change_length){(*num) = 0;return 0;}
                }
                else if((cigar[k].c=='D')&&(cigar[k].l>opt->change_length)) cigar[k].c='N';
            }
            return r;
        }
        return 0;
    }
    return 0;
}
int middle_align(struct m_opt *opt,int chr,char strand,uint64_t pos1,uint64_t pos2,char *seq,int length,struct cigar_t cigar[20],int *num)
{
    if(strand==0)
    {
        if((pos2<pos1)||(pos2-pos1>MAX_READ_LENGTH))
        {(*num) = 0;return 0;}
    }
    if((strand==1)&&(pos2>pos1))
    {
        if((pos2>pos1)||(pos1-pos2>MAX_READ_LENGTH))
        {(*num) = 0;return 0;}
    }

    char ref[MAX_READ_LENGTH];
    char text[MAX_READ_LENGTH];
    int length_r = abs(pos2-pos1)-1;
    char cigarBuf[4*MAX_READ_LENGTH];

    int t = 0,r = 0,k = 0;

    if(strand)
    {
        for(k = 0;k<length_r;k++)
            ref[k] = base2char[3-nst_nt4_table[(int)opt->chr->list[chr].seq[pos1-1-k]]];
    }
    else
    {
        for(k = 0;k<length_r;k++)
            ref[k] = opt->chr->list[chr].seq[pos1+1+k];
    }
    ref[k] = '\0';
    for(k = 0;k<length;k++)
        text[k] = seq[k];
    text[k] = '\0';

    cigarBuf[0] = '\0';
    if(computeEditDistanceWithCigar(ref,length_r,text,length,MAX_K-1,cigarBuf,4*MAX_READ_LENGTH,0,1)!=-1)
    {

        t = 0;r = 0;k = 0;
        *num = 0;
        cigar[0].c = cigarBuf[0];
        cigar[0].l = 0;
        while((t<length)&&(r<length_r))
        {
            if(cigarBuf[k]== cigar[*num].c) {cigar[*num].l++;}
            else
            {
                if((*num)<20)(*num)++;
                else {(*num) = 0;return 0;}
                cigar[*num].c = cigarBuf[k];
                cigar[*num].l = 1;
            }

            if(cigarBuf[k]=='I') t++;
            else if(cigarBuf[k]=='X') {t++;r++;}
            else if(cigarBuf[k]=='D') r++;
            else if(cigarBuf[k]=='=') {t++;r++;}
            k++;
        }
        if((*num)<20)(*num)++;
        else {(*num) = 0;return 0;}
        if(t<length)
        {
            cigar[*num].c = 'I';
            cigar[*num].l = length-t;
            if((*num)<20)(*num)++;
            else {(*num) = 0;return 0;}
        }
        if(r<length_r)
        {
            cigar[*num].c = 'D';
            cigar[*num].l = length_r-r;
            if((*num)<20)(*num)++;
            else {(*num) = 0;return 0;}
        }
        for(k = 0;k<*num;k++)
        {
            if(cigar[k].c=='=') cigar[k].c='M';
        }
    }
    for(k = 0;k<(*num);k++)
    {
        if((cigar[k].c=='X')||(cigar[k].c=='I'))
        {
            if(cigar[k].l>opt->change_length){(*num) = 0;return 0;}
        }
        else if((cigar[k].c=='D')&&(cigar[k].l>opt->change_length)) cigar[k].c='N';
    }
    return 0;
}
int seed2cigar(char *seq,int length,struct seed_t *seed,struct cigar_t *cigar,int *cigar_num,int max,int seed_order[10],int seed_order_n,int chr_order,int strand,int mode);
int tail_align_seed(struct m_opt *opt,int chr,char strand,uint64_t pos,char *seq,int length,struct cigar_t *cigar,int *cigar_num,int max,uint64_t *start_pos,int mode)//0 front 1 back
{
    int i,j,k;
    struct seed_t seed_r;
    struct seed_t seed[101];
    int seed_num = 0;

    struct cigar_t t_cigar[21];
    int t_num;
    struct cigar_t b_cigar[21];
    int b_num;

    int seed_order[10];
    int seed_order_n = 0;

    unsigned int chr_start = opt->chr->list[chr].start_site;

    if(mode==0)
    {
    i = 0;
    while((i<=opt->change_length)&&(i+5<=length))
    {
        for(j = 0;j<=opt->change_length;j++)
        {
            if(strand==1)
            {
                k = 0;
                while((i+k<length)&&(3-nst_nt4_table[(int)seq[length-i-k-1]]==nst_nt4_table[(int)opt->chr->list[chr].seq[(opt->idx->bns->l_pac<<1)-(pos-j-k-1)-1-chr_start]]))
                    k++;
            }
            else
            {
                k = 0;
                while((i+k<length)&&(nst_nt4_table[(int)seq[length-i-k-1]]==nst_nt4_table[(int)opt->chr->list[chr].seq[pos-j-k-1-chr_start]]))
                    k++;
            }

            if(k>=5)
            {
                seed_r.pos = pos-j-k;
                seed_r.start = length-i-k;
                seed_r.length = k;
                seed_r.abs = seed_r.pos-seed_r.start;
                seed_r.lnum = 0;

                {
                    seed_r.score = seed_r.length;
                    insert_seed(seed,&seed_num,100,seed_r);
                }
            }
        }
        i++;
    }
    }
    else if(mode==1)
    {
    i = 0;
    while((i<=opt->change_length)&&(i+5<=length))
    {
        for(j = 0;j<=opt->change_length;j++)
        {
            if(strand==1)
            {
                k = 0;
                while((i+k<length)&&(3-nst_nt4_table[(int)seq[i+k]]==nst_nt4_table[(int)opt->chr->list[chr].seq[(opt->idx->bns->l_pac<<1)-(pos+1+j+k)-1-chr_start]]))
                    k++;
            }
            else
            {
                k = 0;
                while((i+k<length)&&(nst_nt4_table[(int)seq[i+k]]==nst_nt4_table[(int)opt->chr->list[chr].seq[pos+1+j+k-chr_start]]))
                    k++;

            }

            if(k>=5)
            {
                seed_r.pos = pos+1+j;
                seed_r.start = i;
                seed_r.length = k;
                seed_r.abs = pos+1+j-i;
                seed_r.lnum = 0;

                {
                    seed_r.score = seed_r.length;
                    insert_seed(seed,&seed_num,100,seed_r);
                }
            }
        }
        i++;
    }
    }

    int max_o,max_s;
    //回溯
    max_o = -1;
    max_s = -1000;
    for(i = seed_num-1; i>=0; i--)
    {
        if(seed[i].score>=max_s)
        {
            max_o = i;
            max_s = seed[i].score;
        }
    }
    if(seed_num!=0)
    {
        seed_order_n = 1;
        seed_order[0] = max_o;
    }
    else seed_order_n = 0;

    uint64_t ref_pos = 0;
    int read_pos = 0;
    int r;
    int score_t;

    t_num = 0;
    (*cigar_num) = 0;
    if ((seed_num==0)||(seed2cigar(seq,length,&seed[0],&t_cigar[0],&t_num,20,seed_order,seed_order_n,chr,strand,1)!=0))
    {
        if(length>15)
            {(*cigar_num)=0;return 1;}
        if(mode==1)
        {
            if((strand==0)&&(pos+1<chr_start)){(*cigar_num) = 0;return 1;}
            if((strand==1)&&((opt->idx->bns->l_pac<<1)<pos+1+1+chr_start)){(*cigar_num) = 0;return 1;}

            r = tail_align(opt,chr,strand,(strand==0)?pos+1-chr_start:(opt->idx->bns->l_pac<<1)-pos-1-1-chr_start,seq,length,1,t_cigar,&t_num,&score_t);
            if((t_num!=0))
            {
                memcpy(&(cigar[(*cigar_num)]),t_cigar,t_num*(sizeof(struct cigar_t)));
                (*cigar_num)+=t_num;

                (*start_pos) = pos+1;
                return 0;
            }
            else{(*cigar_num) = 0;return 1;}
        }
        if(mode==0)
        {
            if((strand==0)&&(pos<chr_start+1)){(*cigar_num) = 0;return 1;}
            if((strand==1)&&((opt->idx->bns->l_pac<<1)<pos+chr_start)){(*cigar_num) = 0;return 1;}
            r = tail_align(opt,chr,strand,(strand==0)?pos-1-chr_start:(opt->idx->bns->l_pac<<1)-pos+1-1-chr_start,seq,length,0,t_cigar,&t_num,&score_t);
            if((t_num!=0))
            {
                memcpy(&(cigar[(*cigar_num)]),t_cigar,t_num*(sizeof(struct cigar_t)));
                (*cigar_num)+=t_num;

                (*start_pos) = pos-r;
                return 0;
            }
            else{(*cigar_num) = 0;return 1;}
        }

    }
    else
    {
        (*start_pos) = ref_pos = seed[seed_order[seed_order_n-1]].pos;
        if((mode==1)&&(ref_pos>pos)&&(t_cigar[0].c!='S'))
        {
        if((ref_pos-pos>opt->change_length))
        {
            cigar[(*cigar_num)].c = 'N';
            cigar[(*cigar_num)].l = ref_pos-pos-1;
            (*cigar_num)++;
        }
        else if((ref_pos-pos<=opt->change_length)&&(ref_pos-pos>1))
        {
            cigar[(*cigar_num)].c = 'D';
            cigar[(*cigar_num)].l = ref_pos-pos-1;
            (*cigar_num)++;
        }
        }
        i = 0;
        if((mode==1)&&(t_cigar[0].c=='S'))
        {
            if(ref_pos-pos==1)
            {
                cigar[(*cigar_num)].c = 'I';
                cigar[(*cigar_num)].l = t_cigar[i].l;
                (*cigar_num)++;
                if((*cigar_num)>=max){(*cigar_num)=0;return 1;}
            }
            else
            {
            b_num = 0;
            middle_align(opt,chr,strand,(strand==0)?pos-chr_start:(opt->idx->bns->l_pac<<1)-pos-1-chr_start,
                             (strand==0)?ref_pos-chr_start:(opt->idx->bns->l_pac<<1)-ref_pos-1-chr_start,seq,t_cigar[i].l,b_cigar,&b_num);
            if(b_num==0){(*cigar_num)=0;return 1;}
            else
            {
                 memcpy(&(cigar[(*cigar_num)]),b_cigar,b_num*(sizeof(struct cigar_t)));
                (*cigar_num)+=b_num;
            }
            }
            i = 1;
        }
        for(;i<t_num;i++)
        {
            if((mode==0)&&(i==t_num-1)&&(t_cigar[t_num-1].c=='S'))
            {
                if(pos<=ref_pos)
                {
                    cigar[(*cigar_num)].c = 'I';
                    cigar[(*cigar_num)].l = t_cigar[i].l+ref_pos-pos;
                    cigar[(*cigar_num)-1].l -= ref_pos-pos;
                    (*cigar_num)++;
                    if((*cigar_num)>=max){(*cigar_num)=0;return 1;}
                    return 0;
                }
                b_num = 0;
                middle_align(opt,chr,strand,(strand==0)?ref_pos-1-chr_start:(opt->idx->bns->l_pac<<1)-ref_pos+1-1-chr_start,
                             (strand==0)?pos-chr_start:(opt->idx->bns->l_pac<<1)-pos-1-chr_start,&seq[read_pos],t_cigar[i].l,b_cigar,&b_num);
                if(b_num==0){(*cigar_num)=0;return 1;}
                else
                {
                    for(j = 0;j<b_num;j++)
                    {
                        cigar[(*cigar_num)].c = b_cigar[j].c;
                        cigar[(*cigar_num)].l = b_cigar[j].l;
                        (*cigar_num)++;
                        if((*cigar_num)>=max){(*cigar_num)=0;return 1;}
                    }
                }
                break;
            }
            if((t_cigar[i].c=='M')||(t_cigar[i].c=='X')){read_pos+=t_cigar[i].l;ref_pos+=t_cigar[i].l;}
            else if((t_cigar[i].c=='I')||(t_cigar[i].c=='S')){read_pos+=t_cigar[i].l;}
            else if((t_cigar[i].c=='D')||(t_cigar[i].c=='N')){ref_pos+=t_cigar[i].l;}

            cigar[(*cigar_num)].c = t_cigar[i].c;
            cigar[(*cigar_num)].l = t_cigar[i].l;
            (*cigar_num)++;
            if((*cigar_num)>=max){(*cigar_num)=0;return 1;}

            if((mode==0)&&(pos>ref_pos)&&(i == t_num-1)&&(t_cigar[t_num-1].c!='S'))
            {
            if((pos-ref_pos>opt->change_length))
            {
                cigar[(*cigar_num)].c = 'N';
                cigar[(*cigar_num)].l = pos-ref_pos;
                (*cigar_num)++;
            }
            else if(pos-ref_pos<=opt->change_length)
            {
                cigar[(*cigar_num)].c = 'D';
                cigar[(*cigar_num)].l = pos-ref_pos;
                (*cigar_num)++;
            }
            }
        }
        if((mode==0)&&(pos<ref_pos))
        {
            cigar[(*cigar_num)].c = 'I';
            cigar[(*cigar_num)].l = ref_pos-pos;
            cigar[(*cigar_num)-1].l -= ref_pos-pos;
            (*cigar_num)++;
            if((*cigar_num)>=max){(*cigar_num)=0;return 1;}
        }
    }
    //int snp_num = 0;
    //for(i = 0;i<(*cigar_num);i++)
    //{
        //if((cigar[i].c=='I')||(cigar[i].c=='X')||(cigar[i].c=='D')||(cigar[i].c=='N')) snp_num++;
    //}
    //if(snp_num>1){(*cigar_num)=0;return 1;}
    return 0;
}
struct splice_t
{
    unsigned start;
    unsigned end;
    int length;

    struct cigar_t cigar[21];
    int cigar_num;
    int score;
};
int landau2cigar(int length_r,int length_t,char *cigarBuf,struct cigar_t *t_cigar,int *cigar_num,int *score,int max)
{
        int t = 0,r = 0,k = 0;
        int t_num = 0;
        t_cigar[0].c = cigarBuf[0];
        t_cigar[0].l = 0;

        while((t<length_t)&&(r<length_r))
        {
            if(cigarBuf[k]== t_cigar[t_num].c) {t_cigar[t_num].l++;}
            else
            {
                t_num++;
                if(t_num>max){(*cigar_num) = 0;(*score) = 0;return 0;}
                t_cigar[t_num].c = cigarBuf[k];
                t_cigar[t_num].l = 1;
            }

            if(cigarBuf[k]=='I') t++;
            else if(cigarBuf[k]=='X') {t++;r++;}
            else if(cigarBuf[k]=='D') r++;
            else if(cigarBuf[k]=='=') {t++;r++;}
            k++;
        }
        t_num++;
        if(t_num>max){(*cigar_num) = 0;(*score) = 0;return 0;}
        if(t<length_t)
        {
            t_cigar[t_num].c = 'I';
            t_cigar[t_num].l = length_t-t;
            t_num++;
            if(t_num>max){(*cigar_num) = 0;(*score) = 0;return 0;}
        }
        if(r<length_r)
        {
            t_cigar[t_num].c = 'D';
            t_cigar[t_num].l = length_r-r;
            t_num++;
            if(t_num>max){(*cigar_num) = 0;(*score) = 0;return 0;}
        }

        (*score) = 0;
        for(k = 0;k<t_num;k++)
        {
            if(t_cigar[k].c=='=') {t_cigar[k].c='M';(*score)+=t_cigar[k].l*opt->match;}
            else if(t_cigar[k].c=='X') (*score)-=t_cigar[k].l*opt->miss;
            else if(t_cigar[k].c=='I') {(*score)-=opt->miss;if(t_cigar[k].l>opt->change_length) {(*cigar_num) = 0;(*score) = 0;return 0;}}
            else if(t_cigar[k].c=='D') {(*score)-=opt->gap;if(t_cigar[k].l>opt->change_length) {(*cigar_num) = 0;(*score) = 0;return 0;}}
        }
        (*cigar_num) = t_num;
        return 0;
}
int splice_pos(int chr_order,char strand,uint64_t start,uint64_t end,char *seq,int length,struct cigar_t *cigar,int *cigar_num,int *score,int max)
{
    unsigned int GT[20];
    unsigned int AG[20];
    unsigned int CT[20];
    unsigned int AC[20];

    int GTN = 0,AGN = 0,CTN = 0,ACN = 0;

    unsigned int i = 0,j;
    unsigned int Sstart,Send;

    int TH = 5;

    Sstart = start+length+TH;
    Send = end-length-TH;

    char *chr = opt->chr->list[chr_order].seq;

    for(i = start;i<Sstart;i++)
    {
        if((nst_nt4_table[(int)chr[i]]==2)&&(nst_nt4_table[(int)chr[i+1]]==3))
        {GT[GTN] = i;GTN++;if(GTN>20){(*cigar_num) = 0;(*score) = -1000;return 0;}}
        if((nst_nt4_table[(int)chr[i]]==1)&&(nst_nt4_table[(int)chr[i+1]]==3))
        {CT[CTN] = i;CTN++;if(CTN>20){(*cigar_num) = 0;(*score) = -1000;return 0;}}
    }
    for(i = Send;i<=end+1;i++)
    {
        if((nst_nt4_table[(int)chr[i-2]]==0)&&(nst_nt4_table[(int)chr[i-1]]==2))
        {AG[AGN] = i;AGN++;if(AGN>20){(*cigar_num) = 0;(*score) = -1000;return 0;}}
        if((nst_nt4_table[(int)chr[i-2]]==0)&&(nst_nt4_table[(int)chr[i-1]]==1))
        {AC[ACN] = i;ACN++;if(ACN>20){(*cigar_num) = 0;(*score) = -1000;return 0;}}
    }

    char ref[MAX_READ_LENGTH];
    char text[MAX_READ_LENGTH];
    for(j = 0;j<length;j++)
        text[j] = seq[j];
    text[length] = '\0';

    int istart;
    int length_r;
    char cigarBuf[4*MAX_READ_LENGTH];

    struct splice_t splice[100];
    int SN = 0;
    int k = 0;

    for(k = 0;k<GTN;k++)
    {
        for(j = 0;j<AGN;j++)
        {
            splice[SN].start = GT[k];
            splice[SN].end = AG[j];
            splice[SN].length = splice[SN].end-splice[SN].start;
            if(splice[SN].length<=0) continue;

            for(i = start;i<splice[SN].start;i++)
                ref[i-start] = chr[i];
            istart = i-start;
            for(i = splice[SN].end;i<end+1;i++)
                ref[istart+i-splice[SN].end] = chr[i];
            length_r = (splice[SN].start-start)+(end+1-splice[SN].end);
            ref[length_r] = '\0';

            if(length_r==0)
            {
                splice[SN].cigar_num = 1;
                splice[SN].cigar[0].c ='I';
                splice[SN].cigar[0].l =length;
                splice[SN].score = -1;
            }
            else if(computeEditDistanceWithCigar(ref,length_r,text,length,MAX_K-1,cigarBuf,4*MAX_READ_LENGTH,0,1)!=-1)
            {
                landau2cigar(length_r,length,cigarBuf,splice[SN].cigar,&(splice[SN].cigar_num),&(splice[SN].score),20);
                if(splice[SN].cigar_num!=0)
                {
                     SN++;
                    if(SN>100){(*cigar_num) = 0;(*score) = -1000;return 0;}
                }
            }

        }
    }
    for(k = 0;k<CTN;k++)
    {
        for(j = 0;j<ACN;j++)
        {
            splice[SN].start = CT[k];
            splice[SN].end = AC[j];
            splice[SN].length = splice[SN].end-splice[SN].start;
            if(splice[SN].length<=0) continue;

            for(i = start;i<splice[SN].start;i++)
                ref[i-start] = chr[i];
            istart = i-start;
            for(i = splice[SN].end;i<end+1;i++)
                ref[istart+i-splice[SN].end] = chr[i];
            length_r = (splice[SN].start-start)+(end+1-splice[SN].end);
            ref[length_r] = '\0';

            if(length_r==0)
            {
                splice[SN].cigar_num = 1;
                splice[SN].cigar[0].c ='I';
                splice[SN].cigar[0].l =length;
                splice[SN].score = -1;
            }
            else if(computeEditDistanceWithCigar(ref,length_r,text,length,MAX_K-1,cigarBuf,4*MAX_READ_LENGTH,0,1)!=-1)
            {
                landau2cigar(length_r,length,cigarBuf,splice[SN].cigar,&(splice[SN].cigar_num),&(splice[SN].score),20);
                if(splice[SN].cigar_num!=0)
                {
                     SN++;
                    if(SN>100){(*cigar_num) = 0;(*score) = -1000;return 0;}
                }
            }
        }
    }
    int best = 0;
    int best_score = splice[0].score;
    for(i = 1;i<SN;i++)
    {
        if(splice[i].score>best_score)
        {
            best_score = splice[i].score;
            best = i;
        }
    }
    if(SN!=0)
    {
        (*score) = splice[best].score;
        uint64_t ref_pos = start;
        int read_pos = 0;

        for(i = 0;i<splice[best].cigar_num;i++)
        {
            if((splice[best].cigar[i].c=='M')||(splice[best].cigar[i].c=='X')||(splice[best].cigar[i].c=='D'))
            {
                if(ref_pos+splice[best].cigar[i].l<=splice[best].start)
                {
                    cigar[(*cigar_num)].c = splice[best].cigar[i].c;
                    cigar[(*cigar_num)].l = splice[best].cigar[i].l;
                    (*cigar_num)++;
                    if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;return 0;}

                    ref_pos+=splice[best].cigar[i].l;
                    read_pos+=splice[best].cigar[i].l;
                }
                else
                {
                    if(splice[best].start-ref_pos>0)
                    {
                        cigar[(*cigar_num)].c = splice[best].cigar[i].c;
                        cigar[(*cigar_num)].l = splice[best].start-ref_pos;
                        (*cigar_num)++;
                        if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;return 0;}
                        splice[best].cigar[i].l-= splice[best].start-ref_pos;
                    }
                    break;
                }
            }
            else if(splice[best].cigar[i].c=='I')
            {
                cigar[(*cigar_num)].c = splice[best].cigar[i].c;
                cigar[(*cigar_num)].l = splice[best].cigar[i].l;
                (*cigar_num)++;
                if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;return 0;}
            }
            else {(*cigar_num) = 0;(*score) = -1000;return 0;}
        }
        cigar[(*cigar_num)].c = 'N';
        cigar[(*cigar_num)].l = splice[best].length;
        (*cigar_num)++;
        if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;return 0;}
        for(j = i;j<splice[best].cigar_num;j++)
        {
            cigar[(*cigar_num)].c = splice[best].cigar[j].c;
            cigar[(*cigar_num)].l = splice[best].cigar[j].l;
            (*cigar_num)++;
            if((*cigar_num)>max) {(*cigar_num) = 0;(*score) = -1000;return 0;}
        }
    }
    else {(*cigar_num) = 0;(*score) = -1000;return 0;}
    return 0;
}
int middle_pos(int chr,char strand,uint64_t start,uint64_t end,char *seq,int length,struct cigar_t *cigar,int *cigar_num,int *score_m,int maxN)
{
    if(start>end) {(*cigar_num) = 0;return 1;}

    int j;

    int score = 0;
    unsigned int chr_start = opt->chr->list[chr].start_site;

    for(j = 0;j<length;j++)
    {
        if(nst_nt4_table[(int)seq[j]]==nst_nt4_table[(int)opt->chr->list[chr].seq[start-chr_start+j]]) score+=opt->match;
        else score-=opt->miss;
    }
    int max = score;
    int max_l = length;

    for(j = 1;j<=length;j++)
    {
        if(nst_nt4_table[(int)seq[length-j]]==nst_nt4_table[(int)opt->chr->list[chr].seq[start-chr_start+(length-j)]]) score-=opt->match;
        else score+=opt->miss;

        if(nst_nt4_table[(int)seq[length-j]]==nst_nt4_table[(int)opt->chr->list[chr].seq[end-chr_start-j+1]]) score+=opt->match;
        else score-=opt->miss;

        if(((score>max)&&(strand==0))||((score>=max)&&(strand==1)))
        {
            max = score;
            max_l = length-j;
        }
    }

    int flag = 0;
    for(j = 0;j<length;j++)
    {
        if(j<max_l)
        {
            if((nst_nt4_table[(int)seq[j]]==nst_nt4_table[(int)opt->chr->list[chr].seq[start+j-chr_start]]))
            {
                if((cigar[(*cigar_num)].c=='M')&&(flag == 1))
                    cigar[(*cigar_num)].l++;
                else
                {
                    if(flag==1)
                    {
                        (*cigar_num)++;
                        if((*cigar_num)>maxN){(*cigar_num) = 0;return 1;}
                    }
                    cigar[(*cigar_num)].c='M';
                    cigar[(*cigar_num)].l = 1;
                    flag = 1;
                }
            }
            else
            {
                if((cigar[(*cigar_num)].c=='X')&&(flag == 1))
                    cigar[(*cigar_num)].l++;
                else
                {
                    if(flag==1)
                    {
                        (*cigar_num)++;
                        if((*cigar_num)>maxN){(*cigar_num) = 0;return 1;}
                    }
                    cigar[(*cigar_num)].c='X';
                    cigar[(*cigar_num)].l = 1;
                    flag = 1;
                }
            }
        }
        else
        {
            if(j==max_l)
            {
                if(flag==1)
                {
                    (*cigar_num)++;
                    if((*cigar_num)>maxN){(*cigar_num) = 0;return 1;}
                }
                cigar[(*cigar_num)].c='N';
                cigar[(*cigar_num)].l = end-start-length+1;
                flag = 1;
            }
            if(nst_nt4_table[(int)seq[j]]==nst_nt4_table[(int)opt->chr->list[chr].seq[end-length+j+1-chr_start]])
            {
                if((cigar[(*cigar_num)].c=='M')&&(flag == 1))
                    cigar[(*cigar_num)].l++;
                else
                {
                    if(flag==1)
                    {
                        (*cigar_num)++;
                        if((*cigar_num)>maxN){(*cigar_num) = 0;return 1;}
                    }
                    cigar[(*cigar_num)].c='M';
                    cigar[(*cigar_num)].l = 1;
                    flag = 1;
                }
            }
            else
            {
                if((cigar[(*cigar_num)].c=='X')&&(flag == 1))
                    cigar[(*cigar_num)].l++;
                else
                {
                    if(flag==1)
                    {
                        (*cigar_num)++;
                        if((*cigar_num)>maxN){(*cigar_num) = 0;return 1;}
                    }
                    cigar[(*cigar_num)].c='X';
                    cigar[(*cigar_num)].l = 1;
                    flag = 1;
                }
            }
        }
    }
    (*cigar_num)++;
    if(max_l==length)
    {
        cigar[(*cigar_num)].c='N';
        cigar[(*cigar_num)].l = end-start-length+1;
        (*cigar_num)++;
    }
    if((*cigar_num)>maxN){(*cigar_num) = 0;return 1;}
    *score_m = max;
    {
        struct cigar_t t_cigar[21];
        int t_num = 0;
        score  = -1000;

        splice_pos(chr,strand,start-opt->chr->list[chr].start_site,end-opt->chr->list[chr].start_site,seq,length,t_cigar,&t_num,&score,maxN);
        if(score>=max)
        {
            *score_m = score;
             memcpy(cigar,t_cigar,t_num*sizeof(struct cigar_t));
             (*cigar_num) = t_num;
        }
    }
    return 0;
}
int middle_align_seed(struct m_opt *opt,int chr,char strand,uint64_t start,uint64_t end,char *seq,int length,struct cigar_t *cigar,int *cigar_num,int max)
{
    int i,j,k;
    struct seed_t seed_r;
    struct seed_t seed[101];
    int seed_num = 0;

    int length_r = length+10;

    uint64_t ref1_start = start+1;
    uint64_t ref2_start = max((end<length_r+1)?0:end-1-length_r,start+1);

    struct cigar_t t_cigar[21];
    int t_num;
    struct cigar_t b_cigar[21];
    int b_num;

    int seed_order[30];
    int seed_order_n = 0;

    unsigned int chr_start = opt->chr->list[chr].start_site;

    i = 0;
    while(i<length-5)
    {
        for(j = 0;j<length_r-5;j++)
        {
            if(strand==1)
            {
                k = 0;
                while((ref1_start+j+k<end)&&(i+k<length)&&(3-nst_nt4_table[(int)seq[i+k]]==nst_nt4_table[(int)opt->chr->list[chr].seq[(opt->idx->bns->l_pac<<1)-(ref1_start+j+k)-1-chr_start]])) k++;
            }
            else
            {
                k = 0;
                while((ref1_start+j+k<end)&&(i+k<length)&&(nst_nt4_table[(int)seq[i+k]]==nst_nt4_table[(int)opt->chr->list[chr].seq[(ref1_start+j+k)-chr_start]])) k++;
            }

            if(k>=5)
            {
                seed_r.pos = ref1_start+j;
                seed_r.start = i;
                seed_r.length = k;
                seed_r.abs = ref1_start+j-i;
                seed_r.lnum = 0;
                seed_r.score = 0;
                insert_seed(seed,&seed_num,100,seed_r);
            }
        }
        for(j = 0;j<length_r-5;j++)
        {
            if(strand==1)
            {
                k = 0;
                while((ref2_start+j+k<end)&&(i+k<length)&&(3-nst_nt4_table[(int)seq[i+k]]==nst_nt4_table[(int)opt->chr->list[chr].seq[(opt->idx->bns->l_pac<<1)-(ref2_start+j+k)-1-chr_start]])) k++;
            }
            else
            {
                k = 0;
                while((ref2_start+j+k<end)&&(i+k<length)&&(nst_nt4_table[(int)seq[i+k]]==nst_nt4_table[(int)opt->chr->list[chr].seq[(ref2_start+j+k)-chr_start]])) k++;
            }

            if(k>=5)
            {
                seed_r.pos = ref2_start+j;
                seed_r.start = i;
                seed_r.length = k;
                seed_r.abs = ref2_start+j-i;
                seed_r.lnum = 0;
                seed_r.score = 0;
                insert_seed(seed,&seed_num,100,seed_r);
            }
        }
        i++;
    }
    int max_o = -1;
    int max_s = -1000;

    for(i = 0; i<seed_num; i++)
    {
        max_o = -1;
        max_s = -1000;
        for(j = 0;j<i;j++)
        {
            if((seed[i].start > seed[j].start)&&((seed[i].start+seed[i].length)>(seed[j].start+seed[j].length))&&(seed[i].pos > seed[j].pos))
            {
                if(seed[i].start<=seed[j].start+seed[j].length)
                    max_s = seed[j].score + (seed[i].start+seed[i].length-seed[j].start-seed[j].length);
                else
                    max_s = seed[j].score + seed[i].length;

                if(seed[i].abs<seed[j].abs)
                    max_s -=(seed[j].abs-seed[i].abs)+2;
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
            seed[i].score = seed[i].length;
            if((seed[i].abs-start>opt->change_length)) seed[i].score-=opt->splice;
        }
        else
        {
            seed[i].last[0] = max_o;
            seed[i].lnum = 1;
        }

    }
    //回溯
    max_o = -1;
    max_s = -1000;

    for(i = seed_num-1; i>=0; i--)
    {
        if(end-seed[i].pos-seed[i].length>opt->change_length) seed[i].score-=opt->splice;
        if(seed[i].score>=max_s)
        {
            max_o = i;
            max_s = seed[i].score;
        }
    }
    seed_order_n = 0;
    while(max_o!=-1)
    {
        seed_order[seed_order_n] = max_o;
        seed_order_n++;
        if(seed_order_n>=30) break;

        if(seed[max_o].lnum>0)
            max_o = seed[max_o].last[0];
        else max_o = -1;
    }
    uint64_t ref_pos = 0;
    int read_pos = 0;
    t_num = 0;
    (*cigar_num) = 0;
    if ((seed_num==0)||(seed2cigar(seq,length,&seed[0],&t_cigar[0],&t_num,20,seed_order,seed_order_n,chr,strand,1)!=0))
    {
        if(length+6>=end-start-1)
        {
            b_num = 0;
            middle_align(opt,chr,strand,(strand==0)?start-chr_start:(opt->idx->bns->l_pac<<1)-start-1-chr_start,
                             (strand==0)?end-chr_start:(opt->idx->bns->l_pac<<1)-end-1-chr_start,seq,length,cigar,cigar_num);
            if((*cigar_num)==0){(*cigar_num)=0;return 1;}
        }
        else
        {
            struct cigar_t m_cigar[21];
            int m_num;
            int score_t = -1000;
            int score_b = -1000;
            int score_m = -1000;
            int r = 0,b = 0;
            t_num = 0;
            b_num = 0;
            m_num = 0;

            if((strand==0)&&(start+1<chr_start)){(*cigar_num) = 0;return 1;}
            if((strand==1)&&((opt->idx->bns->l_pac<<1)<start+1+1+chr_start)){(*cigar_num) = 0;return 1;}

            if((strand==0)&&(end<chr_start+1)){(*cigar_num) = 0;return 1;}
            if((strand==1)&&((opt->idx->bns->l_pac<<1)<end+chr_start)){(*cigar_num) = 0;return 1;}

            r = tail_align(opt,chr,strand,(strand==0)?start+1-chr_start:(opt->idx->bns->l_pac<<1)-start-1-1-chr_start,seq,length,1,t_cigar,&t_num,&score_t);
            b = tail_align(opt,chr,strand,(strand==0)?end-1-chr_start:(opt->idx->bns->l_pac<<1)-end+1-1-chr_start,seq,length,0,b_cigar,&b_num,&score_b);
            if(strand==0)
                middle_pos(chr,strand,start+1,end-1,seq,length,m_cigar,&m_num,&score_m,20);
            else
            {
                char text[MAX_READ_LENGTH];
                for(j = 0;j<length;j++)
                {
                    if((seq[length-j-1]=='A')||(seq[length-j-1]=='a')) text[j]='T';
                    else if((seq[length-j-1]=='C')||(seq[length-j-1]=='c')) text[j]='G';
                    else if((seq[length-j-1]=='G')||(seq[length-j-1]=='g')) text[j]='C';
                    else if((seq[length-j-1]=='T')||(seq[length-j-1]=='t')) text[j]='A';
                    else text[j]='N';
                }
                text[length] = '\0';
                middle_pos(chr,strand,(opt->idx->bns->l_pac<<1)-end+1-1,(opt->idx->bns->l_pac<<1)-start-1-1,text,length,m_cigar,&m_num,&score_m,20);

                char c;
                int l;
                for(j = 0;j<m_num/2;j++)
                {
                    c = m_cigar[j].c;
                    l = m_cigar[j].l;

                    m_cigar[j].c = m_cigar[m_num-j-1].c;
                    m_cigar[j].l = m_cigar[m_num-j-1].l;

                    m_cigar[m_num-j-1].c = c;
                    m_cigar[m_num-j-1].l = l;
                }
            }

            if((score_m>=score_t)&&(score_m>=score_t)&&(m_num!=0))
            {
                memcpy(&(cigar[(*cigar_num)]),m_cigar,m_num*(sizeof(struct cigar_t)));
                (*cigar_num)+=m_num;
            }
            else if((score_t>=score_b)&&(score_t>score_m)&&(t_num!=0))
            {
                memcpy(&(cigar[(*cigar_num)]),t_cigar,t_num*(sizeof(struct cigar_t)));
                (*cigar_num)+=t_num;

                cigar[(*cigar_num)].c = 'N';
                cigar[(*cigar_num)].l = end-start-1-r;
                (*cigar_num)++;
            }
            else if((score_t<score_b)&&(score_b>score_m)&&(b_num!=0))
            {
                cigar[(*cigar_num)].c = 'N';
                cigar[(*cigar_num)].l = end-start-1-b;
                (*cigar_num)++;

                memcpy(&(cigar[(*cigar_num)]),b_cigar,b_num*(sizeof(struct cigar_t)));
                (*cigar_num)+=b_num;
            }
            else {(*cigar_num) = 0;return 1;}
        }
    }
    else
    {
        ref_pos = seed[seed_order[seed_order_n-1]].pos;
        if((ref_pos>start)&&(t_cigar[0].c!='S'))
        {
        if((ref_pos-start>opt->change_length))
        {
            cigar[(*cigar_num)].c = 'N';
            cigar[(*cigar_num)].l = ref_pos-start-1;
            (*cigar_num)++;
        }
        else if((ref_pos-start<=opt->change_length)&&(ref_pos-start>1))
        {
            cigar[(*cigar_num)].c = 'D';
            cigar[(*cigar_num)].l = ref_pos-start-1;
            (*cigar_num)++;
        }
        }

        i = 0;
        if(t_cigar[0].c=='S')
        {
            if(ref_pos-start==1)
            {
                cigar[(*cigar_num)].c = 'I';
                cigar[(*cigar_num)].l = t_cigar[i].l;
                (*cigar_num)++;
                if((*cigar_num)>=max){(*cigar_num)=0;return 1;}
            }
            else
            {
               b_num = 0;
               if(start>ref_pos)
                    {(*cigar_num)=0;return 1;}
            //middle_align(opt,chr,strand,(strand==0)?start-chr_start:(opt->idx->bns->l_pac<<1)-start-chr_start,
                             //(strand==0)?ref_pos-chr_start:(opt->idx->bns->l_pac<<1)-ref_pos-chr_start,seq,t_cigar[i].l,b_cigar,&b_num);
                middle_align_seed(opt,chr,strand,start,ref_pos,seq,t_cigar[i].l,b_cigar,&b_num,20);
                if(b_num==0){(*cigar_num)=0;return 1;}
                else
                {
                    memcpy(&(cigar[(*cigar_num)]),b_cigar,b_num*(sizeof(struct cigar_t)));
                    (*cigar_num)+=b_num;
                }
            }

            i = 1;
        }
        for(;i<t_num;i++)
        {
            if((i == t_num-1)&&(t_cigar[t_num-1].c=='S'))
            {
                if(end<=ref_pos)
                {
                    cigar[(*cigar_num)].c = 'I';
                    cigar[(*cigar_num)].l = t_cigar[i].l+ref_pos-end;
                    cigar[(*cigar_num)-1].l -= ref_pos-end;
                    (*cigar_num)++;
                    if((*cigar_num)>=max){(*cigar_num)=0;return 1;}
                    return 0;
                }
                b_num = 0;
                if(end<ref_pos)
                    {(*cigar_num)=0;return 1;}
                //middle_align(opt,chr,strand,(strand==0)?ref_pos-1-chr_start:(opt->idx->bns->l_pac<<1)-ref_pos+1-chr_start,
                             //(strand==0)?end-chr_start:(opt->idx->bns->l_pac<<1)-end-chr_start,&seq[read_pos],t_cigar[i].l,b_cigar,&b_num);
                middle_align_seed(opt,chr,strand,ref_pos-1,end,&seq[read_pos],t_cigar[i].l,b_cigar,&b_num,20);
                if(b_num==0){(*cigar_num)=0;return 1;}
                else
                {
                    for(j = 0;j<b_num;j++)
                    {
                        cigar[(*cigar_num)].c = b_cigar[j].c;
                        cigar[(*cigar_num)].l = b_cigar[j].l;
                        (*cigar_num)++;
                        if((*cigar_num)>=max){(*cigar_num)=0;return 1;}
                    }
                }
                break;
            }
            if((t_cigar[i].c=='M')||(t_cigar[i].c=='X')){read_pos+=t_cigar[i].l;ref_pos+=t_cigar[i].l;}
            else if((t_cigar[i].c=='I')||(t_cigar[i].c=='S')){read_pos+=t_cigar[i].l;}
            else if((t_cigar[i].c=='D')||(t_cigar[i].c=='N')){ref_pos+=t_cigar[i].l;}

            cigar[(*cigar_num)].c = t_cigar[i].c;
            cigar[(*cigar_num)].l = t_cigar[i].l;
            (*cigar_num)++;
            if((*cigar_num)>=max){(*cigar_num)=0;return 1;}

            if((end>ref_pos)&&(i == t_num-1)&&(t_cigar[t_num-1].c!='S'))
            {
            if((end-ref_pos>opt->change_length))
            {
                cigar[(*cigar_num)].c = 'N';
                cigar[(*cigar_num)].l = end-ref_pos;
                (*cigar_num)++;
            }
            else if(end-ref_pos<=opt->change_length)
            {
                cigar[(*cigar_num)].c = 'D';
                cigar[(*cigar_num)].l = end-ref_pos;
                (*cigar_num)++;
            }
            }
        }
    }
    return 0;
}
int check_splice_pos(struct m_opt *opt,char *chr,char strand,uint64_t splice_pos,int splice_length,int f,int b,char *seq,int pos,int flag)
{
    int i = 0;
    int start = f;
    int end = b;

    for(i = 0;i<start;i++)
    {
        if(strand)
        {
            if(3-nst_nt4_table[(int)chr[splice_pos-splice_length+i]]==nst_nt4_table[(int)seq[pos-i]])
            {
                if((nst_nt4_table[(int)chr[splice_pos+i]]==2)&&
                (nst_nt4_table[(int)chr[splice_pos+i-1]]==0)
                &&(nst_nt4_table[(int)chr[splice_pos-splice_length+i+2]]==3)
                &&(nst_nt4_table[(int)chr[splice_pos-splice_length+i+1]]==2))
                {return -i;}
            }
            else break;
        }
        else
        {
            if((nst_nt4_table[(int)chr[splice_pos+splice_length-i]]==nst_nt4_table[(int)seq[pos-i]])&&(splice_pos>=i))
            {
                if((nst_nt4_table[(int)chr[splice_pos-i]]==2)&&
                (nst_nt4_table[(int)chr[splice_pos-i+1]]==3)
                &&(nst_nt4_table[(int)chr[splice_pos+splice_length-i-2]]==0)
                &&(nst_nt4_table[(int)chr[splice_pos+splice_length-i-1]]==2))
                {return -i;}
            }
            else break;
        }
    }

    for(i = 0;i<end;i++)
    {
        if(strand)
        {
            if((3-nst_nt4_table[(int)chr[splice_pos-i]]==nst_nt4_table[(int)seq[pos+i]])&&(splice_pos-2>=i))
            {
                if((nst_nt4_table[(int)chr[splice_pos-i-1]]==2)&&
                (nst_nt4_table[(int)chr[splice_pos-i-2]]==0)
                &&(nst_nt4_table[(int)chr[splice_pos-splice_length-i+1]]==3)
                &&(nst_nt4_table[(int)chr[splice_pos-splice_length-i]]==2))
                {return i+1;}
            }
            else break;
        }
        else
        {
            if(nst_nt4_table[(int)chr[splice_pos+i]]==nst_nt4_table[(int)seq[pos+i]])
            {
                if((nst_nt4_table[(int)chr[splice_pos+i+1]]==2)&&
                (nst_nt4_table[(int)chr[splice_pos+i+2]]==3)
                &&(nst_nt4_table[(int)chr[splice_pos+splice_length+i-1]]==0)
                &&(nst_nt4_table[(int)chr[splice_pos+splice_length+i]]==2))
                {return i+1;}
            }
            else break;
        }
    }
    for(i = 0;i<start;i++)
    {
        if(strand)
        {
            if(3-nst_nt4_table[(int)chr[splice_pos-splice_length+i]]==nst_nt4_table[(int)seq[pos-i]])
            {
                if((nst_nt4_table[(int)chr[splice_pos+i]]==1)&&
                (nst_nt4_table[(int)chr[splice_pos+i-1]]==0)
                &&(nst_nt4_table[(int)chr[splice_pos-splice_length+i+2]]==3)
                &&(nst_nt4_table[(int)chr[splice_pos-splice_length+i+1]]==1))
                {return -i;}
            }
            else break;
        }
        else
        {
            if((nst_nt4_table[(int)chr[splice_pos+splice_length-i]]==nst_nt4_table[(int)seq[pos-i]])&&(splice_pos>=i))
            {
                if((nst_nt4_table[(int)chr[splice_pos-i]]==1)&&
                (nst_nt4_table[(int)chr[splice_pos-i+1]]==3)
                &&(nst_nt4_table[(int)chr[splice_pos+splice_length-i-2]]==0)
                &&(nst_nt4_table[(int)chr[splice_pos+splice_length-i-1]]==1))
                {return -i;}
            }
            else break;
        }
    }

    for(i = 0;i<end;i++)
    {
        if(strand)
        {
            if((3-nst_nt4_table[(int)chr[splice_pos-i]]==nst_nt4_table[(int)seq[pos+i]])&&(splice_pos-2>=i))
            {
                if((nst_nt4_table[(int)chr[splice_pos-i-1]]==1)&&
                (nst_nt4_table[(int)chr[splice_pos-i-2]]==0)
                &&(nst_nt4_table[(int)chr[splice_pos-splice_length-i+1]]==3)
                &&(nst_nt4_table[(int)chr[splice_pos-splice_length-i]]==1))
                {return i+1;}
            }
            else break;
        }
        else
        {
            if(nst_nt4_table[(int)chr[splice_pos+i]]==nst_nt4_table[(int)seq[pos+i]])
            {
                if((nst_nt4_table[(int)chr[splice_pos+i+1]]==1)&&
                (nst_nt4_table[(int)chr[splice_pos+i+2]]==3)
                &&(nst_nt4_table[(int)chr[splice_pos+splice_length+i-1]]==0)
                &&(nst_nt4_table[(int)chr[splice_pos+splice_length+i]]==1))
                {return i+1;}
            }
            else break;
        }
    }
    if(flag == 1)
    {
        if(strand)
        {
            for(i = 0;i<start;i++)
            {if(3-nst_nt4_table[(int)chr[splice_pos-splice_length+i+1]]!=nst_nt4_table[(int)seq[pos-i-1]]) break;}
            return -i;
        }
        else
        {
            for(i = 0;i<end;i++)
            {if(nst_nt4_table[(int)chr[splice_pos+i]]!=nst_nt4_table[(int)seq[pos+i]]) break;}
            return i;
        }
    }
    return 0;
}

int seed2cigar(char *seq,int length,struct seed_t *seed,struct cigar_t *cigar,int *cigar_num,int max,int seed_order[10],int seed_order_n,int chr_order,int strand,int mode)//mode 0 tail_align 1 tail_seed_align
{
    struct seed_t temp_seed;

    int i,j,x;
    int read_pos = 0;
    uint64_t ref_pos = seed[seed_order[seed_order_n-1]].pos;
    //(*start_pos) = seed[seed_order[seed_order_n-1]].pos;
    uint64_t chr_start = opt->chr->list[chr_order].start_site;

    struct cigar_t t_cigar[21];
    int t_num;
    int r = 0;

    int break_flag = 0;

    int flag;

    if(seed[seed_order[seed_order_n-1]].start!=0)
    {
        cigar[(*cigar_num)].c = 'S';
        cigar[(*cigar_num)].l = seed[seed_order[seed_order_n-1]].start;
        (*cigar_num)++;
    }
    cigar[(*cigar_num)].c = 'M';
    cigar[(*cigar_num)].l = seed[seed_order[seed_order_n-1]].length;
    read_pos = seed[seed_order[seed_order_n-1]].start+seed[seed_order[seed_order_n-1]].length;
    ref_pos = seed[seed_order[seed_order_n-1]].pos+seed[seed_order[seed_order_n-1]].length;

    for (i = seed_order_n-2;i>=0;i--)
    {
        if(seed[seed_order[i]].abs==seed[seed_order[i+1]].abs)
        {
            if(seed[seed_order[i]].start<=read_pos)//基本不可能
            {
                cigar[(*cigar_num)].l += seed[seed_order[i]].length-(read_pos-seed[seed_order[i]].start);

                ref_pos+=seed[seed_order[i]].length-(read_pos-seed[seed_order[i]].start);
                read_pos+=seed[seed_order[i]].length-(read_pos-seed[seed_order[i]].start);
            }
            else //XM *
            {
                x = seed[seed_order[i]].start-read_pos;
                for(j = 0;j<x;j++)
                {
                    if(((strand==0)&&(nst_nt4_table[(int)seq[read_pos+j]]==nst_nt4_table[(int)opt->chr->list[chr_order].seq[ref_pos+j-chr_start]]))
                        ||((strand==1)&&(nst_nt4_table[(int)seq[read_pos+j]]==3-nst_nt4_table[(int)opt->chr->list[chr_order].seq[(opt->idx->bns->l_pac<<1)-ref_pos-j-1-chr_start]])))
                    {
                        if(cigar[(*cigar_num)].c=='M')
                            cigar[(*cigar_num)].l++;
                        else
                        {
                            (*cigar_num)++;
                            if((*cigar_num)>max){break_flag = 1;break;}
                            cigar[(*cigar_num)].c='M';
                            cigar[(*cigar_num)].l = 1;
                        }
                    }
                    else
                    {
                        if(cigar[(*cigar_num)].c=='X')
                            cigar[(*cigar_num)].l++;
                        else
                        {
                            (*cigar_num)++;
                            if((*cigar_num)>max){break_flag = 1;break;}
                            cigar[(*cigar_num)].c='X';
                            cigar[(*cigar_num)].l = 1;
                        }
                    }
                }
                (*cigar_num)++;
                if((*cigar_num)>max){break_flag = 1;break;}
                cigar[(*cigar_num)].l = seed[seed_order[i]].length;
                cigar[(*cigar_num)].c = 'M';

                ref_pos = seed[seed_order[i]].pos+seed[seed_order[i]].length;
                read_pos = seed[seed_order[i]].start+seed[seed_order[i]].length;
            }
        }
        else if(seed[seed_order[i]].abs<=seed[seed_order[i+1]].abs)
        {
            (*cigar_num)++;
            if((*cigar_num)>max){break_flag = 1;break;}

            memcpy(&temp_seed,&(seed[seed_order[i]]),sizeof(struct seed_t));
            if(seed[seed_order[i]].start>read_pos) //IX *
            {
                if(seed[seed_order[i]].pos<=ref_pos)
                {
                    cigar[(*cigar_num)].l = seed[seed_order[i+1]].abs-seed[seed_order[i]].abs;
                    cigar[(*cigar_num)].c = 'I';

                    read_pos+=cigar[(*cigar_num)].l;
                    temp_seed.pos +=read_pos-seed[seed_order[i]].start;
                    temp_seed.length -=read_pos-seed[seed_order[i]].start;
                    temp_seed.start +=read_pos-seed[seed_order[i]].start;

                    if((temp_seed.length<=0)||(cigar[(*cigar_num)].l>opt->change_length)) {break_flag = 1;break;}

                    (*cigar_num)++;
                    if((*cigar_num)>max){break_flag = 1;break;}


                }
                else
                {
                    t_num = 0;
                    r = middle_align(opt,chr_order,strand,(strand==0)?ref_pos-chr_start-1:(opt->idx->bns->l_pac<<1)-ref_pos+1-chr_start,
                                 (strand==0)?seed[seed_order[i]].pos-chr_start:(opt->idx->bns->l_pac<<1)-seed[seed_order[i]].pos-chr_start,&(seq[read_pos]),seed[seed_order[i]].start-read_pos,t_cigar,&t_num);
                    if(t_num==0) {break_flag = 1;break;}
                    else
                    {
                        for(x = 0;x<t_num;x++)
                        {
                            cigar[(*cigar_num)].c = t_cigar[x].c;
                            cigar[(*cigar_num)].l = t_cigar[x].l;
                            (*cigar_num)++;if((*cigar_num)>=max){break_flag = 1;break;}
                        }
                    }
                }
            }
            else
            {
                cigar[(*cigar_num)].l = seed[seed_order[i+1]].abs-seed[seed_order[i]].abs;
                cigar[(*cigar_num)].c = 'I';

                read_pos+=cigar[(*cigar_num)].l;
                temp_seed.pos +=read_pos-seed[seed_order[i]].start;
                temp_seed.length -=read_pos-seed[seed_order[i]].start;
                temp_seed.start +=read_pos-seed[seed_order[i]].start;

                if((temp_seed.length<=0)||(cigar[(*cigar_num)].l>opt->change_length)) {break_flag = 1;break;}

                (*cigar_num)++;
                if((*cigar_num)>max){break_flag = 1;break;}
            }
            if(temp_seed.length>0)
            {
                cigar[(*cigar_num)].l = temp_seed.length;
                cigar[(*cigar_num)].c = 'M';

                ref_pos = seed[seed_order[i]].pos+seed[seed_order[i]].length;
                read_pos = seed[seed_order[i]].start+seed[seed_order[i]].length;
            }
        }
        else//DN
        {
            (*cigar_num)++;
            if((*cigar_num)>=max){break_flag = 1;break;}

            memcpy(&temp_seed,&(seed[seed_order[i]]),sizeof(struct seed_t));

            if(read_pos<seed[seed_order[i]].start)
            {
                middle_align_seed(opt,chr_order,strand,ref_pos-1,seed[seed_order[i]].pos,&seq[read_pos],seed[seed_order[i]].start-read_pos,t_cigar,&t_num,20);
                if(t_num==0) {break_flag = 1;break;}
                else
                {
                    for(x = 0;x<t_num;x++)
                    {
                        if(t_cigar[x].c=='N')
                        {
                            if(((strand==0)&&(cigar[(*cigar_num)-1].c=='M'))||((strand==1)&&((x+1>=t_num)||(t_cigar[x+1].c=='M')))) flag = 1;
                            else flag = 0;
                            r = 0;
                            r=check_splice_pos(opt,opt->chr->list[chr_order].seq,strand,(strand==0)?ref_pos-chr_start:(opt->idx->bns->l_pac<<1)-ref_pos-chr_start-1,t_cigar[x].l,cigar[(*cigar_num)-1].l,(x+1>=t_num)?seed[seed_order[i]].length:t_cigar[x+1].l,seq,read_pos,flag);

                            if(r>0)
                            {
                                if(cigar[(*cigar_num)-1].c=='M') cigar[(*cigar_num)-1].l+=r;
                                else
                                {
                                cigar[(*cigar_num)].c = 'M';
                                cigar[(*cigar_num)].l = r;
                                (*cigar_num)++;
                                if((*cigar_num)>max){break_flag = 1;break;}
                                }
                            }
                            else
                            {
                                cigar[(*cigar_num)-1].l+=r;
                                if(cigar[(*cigar_num)-1].l==0)
                                    (*cigar_num)--;
                            }
                            ref_pos+=r;
                            read_pos+=r;

                            cigar[(*cigar_num)].c = 'N';
                            cigar[(*cigar_num)].l = t_cigar[x].l;
                            (*cigar_num)++;
                            if((*cigar_num)>max){break_flag = 1;break;}

                            if(x+1>=t_num)
                            {
                                temp_seed.pos+=r;
                                temp_seed.length-=r;
                                temp_seed.start+=r;
                            }
                            else t_cigar[x+1].l-=r;
                            ref_pos+=t_cigar[x].l;
                            x++;
                            if(x>=t_num) break;
                        }

                        cigar[(*cigar_num)].l = t_cigar[x].l;cigar[(*cigar_num)].c = t_cigar[x].c;
                        (*cigar_num)++;
                        if((*cigar_num)>max){break_flag = 1;break;}

                        if((t_cigar[x].c=='M')||(t_cigar[x].c=='X')){read_pos+=t_cigar[x].l;ref_pos+=t_cigar[x].l;}
                        else if((t_cigar[x].c=='I')||(t_cigar[x].c=='S')){read_pos+=t_cigar[x].l;}
                        else if((t_cigar[x].c=='D')||(t_cigar[x].c=='N')){ref_pos+=t_cigar[x].l;}
                    }
                }
            }
            else if(read_pos>=seed[seed_order[i]].start)
            {
                temp_seed.pos+=read_pos-seed[seed_order[i]].start;
                temp_seed.length-=read_pos-seed[seed_order[i]].start;
                temp_seed.start+=read_pos-seed[seed_order[i]].start;

                cigar[(*cigar_num)].l = temp_seed.pos-ref_pos;
                if(cigar[(*cigar_num)].l<0)
                {
                    cigar[(*cigar_num)].l = 0-cigar[(*cigar_num)].l;
                    cigar[(*cigar_num)].c = 'I';
                    (*cigar_num)++;
                    if((*cigar_num)>max){break_flag = 1;break;}
                }
                else if(cigar[(*cigar_num)].l<=opt->change_length)
                {
                    if(cigar[(*cigar_num)-1].c == 'D')
                        {
                            cigar[(*cigar_num)-1].l+=cigar[(*cigar_num)].l;
                            (*cigar_num)--;
                        }
                    cigar[(*cigar_num)].c = 'D';
                    (*cigar_num)++;
                    if((*cigar_num)>max){break_flag = 1;break;}
                }
                else
                {
                    if(cigar[(*cigar_num)-1].c == 'N')
                        {
                            cigar[(*cigar_num)-1].l+=cigar[(*cigar_num)].l;
                            (*cigar_num)--;
                        }
                    cigar[(*cigar_num)].c = 'N';
                    if(((strand==0)&&(cigar[(*cigar_num)-1].c=='M'))||(strand==1)) flag = 1;
                    else flag = 0;
                    r = 0;
                    r=check_splice_pos(opt,opt->chr->list[chr_order].seq,strand,(strand==0)?ref_pos-chr_start:(opt->idx->bns->l_pac<<1)-ref_pos-chr_start-1,cigar[(*cigar_num)].l,cigar[(*cigar_num)-1].l,temp_seed.length,seq,read_pos,flag);

                    if(cigar[(*cigar_num)-1].c=='M') cigar[(*cigar_num)-1].l+=r;
                    else if(r>0)
                    {
                            (*cigar_num)++;
                            cigar[(*cigar_num)].c = 'N';
                            cigar[(*cigar_num)].l = cigar[(*cigar_num)-1].l;
                            cigar[(*cigar_num)-1].c = 'M';
                            cigar[(*cigar_num)-1].l = r;
                    }
                    (*cigar_num)++;
                    if((*cigar_num)>max){break_flag = 1;break;}
                    temp_seed.pos+=r;
                    temp_seed.length-=r;
                    temp_seed.start+=r;
                }
            }

            if(temp_seed.length>0)
            {
                cigar[(*cigar_num)].l = temp_seed.length;
                cigar[(*cigar_num)].c = 'M';

                read_pos=seed[seed_order[i]].start+seed[seed_order[i]].length;
                ref_pos=seed[seed_order[i]].pos+seed[seed_order[i]].length;
            }
            else {
                    if(i==0)
                    {
                        ref_pos -= cigar[(*cigar_num)-1].l;
                        (*cigar_num)--;
                    }
                    (*cigar_num)--;
            }
        }
    }
    (*cigar_num)++;
    if((read_pos!=length)&&(!break_flag))
    {
        if(read_pos>length)
            cigar[(*cigar_num)].l -= read_pos-length;
        else
        {
            cigar[(*cigar_num)].l = length-read_pos;
            cigar[(*cigar_num)].c = 'S';
            (*cigar_num)++;
        }
    }
    //if(strand) {(*start_pos) = ref_pos;(*start_pos) = (opt->idx->bns->l_pac<<1)-(*start_pos);}
    if(break_flag) {(*cigar_num) = 0;return 1;}
    return 0;
}

void insert_cand_seed(struct m_opt *opt,struct read_inf_t *read,struct seed_t *seed,int seed_order[10],int seed_order_n,int strand,int chr_order,struct exon_array *exon,struct snp_list_t *snp)
{
    int i,j;
    struct cigar_t cigar[MAX_CIGAR_BUF];
    int cigar_num = 0;
    char temp[10];
    int break_flag;
    int flag = 0;
    for (i = seed_order_n-1;i>=0;i--)
        {
            if(seed[seed_order[i]].length<12){flag = 0;continue;}
            if(seed[seed_order[i]].start+seed[seed_order[i]].length<10) {continue;}
            if(read->length-seed[seed_order[i]].start<10) {continue;}

            if((i!=seed_order_n-1)&&(flag==1))
            {
                if((seed[seed_order[i]].abs>=seed[seed_order[i+1]].abs)&&(seed[seed_order[i]].abs-seed[seed_order[i+1]].abs<opt->change_length))continue;
                else if((seed[seed_order[i]].abs<seed[seed_order[i+1]].abs)&&(seed[seed_order[i+1]].abs-seed[seed_order[i]].abs<opt->change_length))continue;
            }
            flag = 1;
            read->cand[read->cand_num].strand = strand;
            read->cand[read->cand_num].pos = seed[seed_order[i]].pos;
            read->cand[read->cand_num].chr_order = chr_order;
            read->cand[read->cand_num].dis = 0;
            read->cand[read->cand_num].cigar[0] = '\0';
            read->cand[read->cand_num].score = seed[seed_order[i]].length;
            read->cand[read->cand_num].out_put_flag = 0;

            cigar_num = 0;
            if(seed[seed_order[i]].start!=0)
            {
                cigar[cigar_num].c = 'S';
                cigar[cigar_num].l = seed[seed_order[i]].start;
                cigar_num++;
            }
            cigar[cigar_num].c = 'M';
            cigar[cigar_num].l = seed[seed_order[i]].length;
            cigar_num++;
            if((seed[seed_order[i]].start+seed[seed_order[i]].length)!=read->length)
            {
                cigar[cigar_num].c = 'S';
                cigar[cigar_num].l = read->length-seed[seed_order[i]].start-seed[seed_order[i]].length;
                cigar_num++;
            }
            if(strand)
            {
                for(j = cigar_num-1;j>=0;j--)
                {
                    sprintf(temp,"%d%c",cigar[j].l,cigar[j].c);
                    strcat(read->cand[read->cand_num].cigar,temp);
                }
                read->cand[read->cand_num].pos = (opt->idx->bns->l_pac<<1)-seed[seed_order[i]].pos-seed[seed_order[i]].length;
            }
            else
            {
                for(j = 0;j<cigar_num;j++)
                {
                    sprintf(temp,"%d%c",cigar[j].l,cigar[j].c);
                    strcat(read->cand[read->cand_num].cigar,temp);
                }
            }

            break_flag = 0;
            for(j = 0;j<read->cand_num;j++)
            {
                if((read->cand[read->cand_num].pos==read->cand[j].pos)&&(strcmp(read->cand[read->cand_num].cigar,read->cand[j].cigar)==0)){break_flag=1;break;}
            }
            if(break_flag) read->cand[read->cand_num].cigar[0] = '\0';
            else
            {
                if(opt->step_flag==SEED_STEP)
                    find_snp2(opt,read,exon,snp,cigar,cigar_num,read->cand[read->cand_num].pos);
                read->cand_num++;
            }
            if(read->cand_num>=MAX_CAND_NUM) break;
        }
}
void check_read_splice(char *read_seq,char *read_rseq,struct cigar_t *cigar,int *cigar_num,uint64_t start_pos,int chr_order,char strand)
{
    int i,j;
    unsigned int ref_pos = start_pos-opt->chr->list[chr_order].start_site;
    int read_pos = 0;

    char text[MAX_READ_LENGTH];

    struct cigar_t t_cigar[21];
    int t_num;

        unsigned int Start_ref_pos = ref_pos,End_ref_pos = ref_pos;
        int Start_read_pos = 0,End_read_pos = 0;
        int Cigar_start = 0,Cigar_end = 0;
        int score = 0;
        int Sscore = 0;
        int flag = 0;
        int length = 0;

        for(i = 0;i<(*cigar_num);i++)
        {
            if(cigar[i].c=='M')
            {
                ref_pos+=cigar[i].l;
                read_pos+=cigar[i].l;

                if(cigar[i].l>=5)
                {
                    Start_ref_pos = ref_pos;
                    Start_read_pos = read_pos;
                    Cigar_start = i;
                }
            }
            else if(cigar[i].c=='X'){ref_pos+=cigar[i].l;read_pos+=cigar[i].l;}
            else if((cigar[i].c=='I')||(cigar[i].c=='S')){read_pos+=cigar[i].l;}
            else if(cigar[i].c=='D'){ref_pos+=cigar[i].l;}
            else if(cigar[i].c=='N')
            {
                if((nst_nt4_table[(int)opt->chr->list[chr_order].seq[ref_pos]]==2)&&(nst_nt4_table[(int)opt->chr->list[chr_order].seq[ref_pos+1]]==3)
                   &&(nst_nt4_table[(int)opt->chr->list[chr_order].seq[ref_pos+cigar[i].l-2]]==0)&&(nst_nt4_table[(int)opt->chr->list[chr_order].seq[ref_pos+cigar[i].l-1]]==2))
                {ref_pos+=cigar[i].l;continue;}
                if((nst_nt4_table[(int)opt->chr->list[chr_order].seq[ref_pos]]==1)&&(nst_nt4_table[(int)opt->chr->list[chr_order].seq[ref_pos+1]]==3)
                   &&(nst_nt4_table[(int)opt->chr->list[chr_order].seq[ref_pos+cigar[i].l-2]]==0)&&(nst_nt4_table[(int)opt->chr->list[chr_order].seq[ref_pos+cigar[i].l-1]]==1))
                {ref_pos+=cigar[i].l;continue;}
                if(Start_ref_pos==start_pos-opt->chr->list[chr_order].start_site)
                    {ref_pos+=cigar[i].l;continue;}

                if(cigar[i-1].c=='S'){ref_pos+=cigar[i].l;}
                else
                {
                    ref_pos+=cigar[i].l;
                    j = i+1;
                    End_ref_pos = ref_pos;
                    End_read_pos = read_pos;
                    Cigar_end = i+1;
                    flag = 0;
                    score = 0;

                    while(j<(*cigar_num))
                    {
                        if(cigar[j].c=='M')
                        {
                            if(cigar[j].l>=5) break;
                            else
                            {
                                score+=cigar[j].l*opt->match;
                                End_ref_pos+=cigar[j].l;
                                End_read_pos+=cigar[j].l;
                            }
                        }
                        else if(cigar[j].c=='X'){End_ref_pos+=cigar[j].l;End_read_pos+=cigar[j].l;score-=opt->miss*cigar[j].l;}
                        else if(cigar[j].c=='I'){End_read_pos+=cigar[j].l;score-=opt->miss;}
                        else if(cigar[j].c=='D'){End_ref_pos+=cigar[j].l;score-=opt->miss;}
                        else
                            {flag = 1;break;}
                        Cigar_end = j+1;
                        j++;
                    }
                    if(!flag)
                    {
                        for(j = Start_read_pos-4;j<End_read_pos+4;j++)
                        {
                            if(strand) text[j-Start_read_pos+4] = read_rseq[j];
                            else text[j-Start_read_pos+4] = read_seq[j];
                        }
                        length = End_read_pos+8-Start_read_pos;
                        text[length] = '\0';

                        t_num = 0;
                        Sscore = -1000;
                        splice_pos(chr_order,strand,Start_ref_pos-4,End_ref_pos+3,text,length,t_cigar,&t_num,&Sscore,20);
                        if(Sscore>=score+4)
                        //if(Sscore!=-1000)
                        {
                            End_ref_pos += cigar[Cigar_end].l;
                            End_read_pos += cigar[Cigar_end].l;

                            cigar[Cigar_start].l-=4;
                            cigar[Cigar_end].l-=4;

                            int tstart = 0;
                            int tend  = t_num-1;

                            if(t_cigar[0].c==cigar[Cigar_start].c) {tstart = 1;cigar[Cigar_start].l+=t_cigar[0].l;}
                            if(t_cigar[tend].c==cigar[Cigar_end].c) {cigar[Cigar_end].l+=t_cigar[tend].l;tend --;}

                            int abs = (tend-tstart+1)-(Cigar_end-Cigar_start-1);

                            if(abs>0)
                            {
                                for(j = (*cigar_num)-1;j>=Cigar_end;j--)
                                {
                                    cigar[j+abs].c = cigar[j].c;
                                    cigar[j+abs].l = cigar[j].l;
                                }
                            }
                            else if(abs<0)
                            {
                                for(j = Cigar_end;j<(*cigar_num);j++)
                                {
                                    cigar[j+abs].c = cigar[j].c;
                                    cigar[j+abs].l = cigar[j].l;
                                }
                            }
                            (*cigar_num)+=abs;
                            Cigar_end+=abs;

                            for(j = tstart;j<=tend;j++)
                            {
                                cigar[Cigar_start+1+j-tstart].c = t_cigar[j].c;
                                cigar[Cigar_start+1+j-tstart].l = t_cigar[j].l;
                            }
                            i = Cigar_end;
                            ref_pos = End_ref_pos;
                            read_pos = End_read_pos;

                            if(cigar[i].l>=5)
                            {
                                Start_ref_pos = ref_pos;
                                Start_read_pos = read_pos;
                                Cigar_start = i;
                            }
                            if(cigar[Cigar_start].l<5)
                            {
                                read_pos = 0;
                                ref_pos = start_pos-opt->chr->list[chr_order].start_site;
                                Start_ref_pos = ref_pos;
                                Start_read_pos = 0;
                                Cigar_start = 0;
                                i = -1;
                            }
                        }
                    }
                }
            }
        }
}
void generate_cigar(struct m_opt *opt,struct read_inf_t *read,struct seed_t *seed,int seed_order[10],int seed_order_n,struct exon_array *exon,struct snp_list_t *snp)
{
    int i = 0,j = 0;
    //int read_pos = 0;
    uint64_t ref_pos = 0,start_pos = 0,end_pos = 0;
    int chr_order = 0;
    char strand = 0;

    int score = 0;
    struct cigar_t cigar[MAX_CIGAR_BUF];
    int cigar_num = 0;
    int distance = 0;

    int break_flag = 0;

    unsigned int chr_start = 0;

    struct cigar_t t_cigar[MAX_CIGAR_BUF];
    int t_num;

    struct cigar_t b_cigar[21];
    int b_num;

    char temp[10];

    if(seed[seed_order[0]].pos>=opt->idx->bns->l_pac)
    {
        for(j = 0;j<opt->chr->total;j++)
        {
            if(seed[seed_order[0]].pos <= (opt->idx->bns->l_pac<<1)-opt->chr->list[j].start_site-opt->chr->list[j].length){continue;}
            else break;
        }
        chr_order = j;
        strand = 1;
    }
    else{
        for(j = 0;j<opt->chr->total;j++)
        {
            if(seed[seed_order[0]].pos >= opt->chr->list[j].start_site+opt->chr->list[j].length){continue;}
            else break;
        }
        chr_order = j;
    }
    chr_start = opt->chr->list[chr_order].start_site;
    read->cand[read->cand_num].out_put_flag = 1;

    int snp_num = 0;
    t_num = 0;
    cigar_num = 0;

    ref_pos = 0;
    if((seed2cigar(read->seq,read->length,seed,t_cigar,&t_num,MAX_CIGAR_BUF-1,seed_order,seed_order_n,chr_order,strand,0)!=0))
    {
        if(read->out_flag==0)
            insert_cand_seed(opt,read,seed,seed_order,seed_order_n,strand,chr_order,exon,snp);
    }
    else
    {
        ref_pos = start_pos = seed[seed_order[seed_order_n-1]].pos;
        i = 0;
        if(t_cigar[0].c=='S')
        {
            b_num = 0;
            tail_align_seed(opt,chr_order,strand,start_pos,read->seq,t_cigar[i].l,b_cigar,&b_num,20,&start_pos,0);
            if(b_num==0)
            {
                ref_pos = start_pos = seed[seed_order[seed_order_n-1]].pos;
                cigar[cigar_num].c = t_cigar[0].c;
                cigar[cigar_num].l = t_cigar[0].l;
                cigar_num++;
                read->cand[read->cand_num].out_put_flag = 0;
            }
            else
            {
                for(j = 0;j<b_num;j++)
                {
                    cigar[cigar_num].c = b_cigar[j].c;
                    cigar[cigar_num].l = b_cigar[j].l;
                    cigar_num++;

                    if(b_cigar[j].c=='X') snp_num+=b_cigar[j].l;
                    else if((b_cigar[j].c=='I')||(b_cigar[j].c=='D')||(b_cigar[j].c=='N')) snp_num+=1;

                    if(b_cigar[j].c=='S') read->cand[read->cand_num].out_put_flag = 0;
                }
            }
            i = 1;
            if((t_cigar[0].c=='M')||(t_cigar[0].c=='X')||(t_cigar[0].c=='D')||(t_cigar[0].c=='N')) ref_pos+=t_cigar[0].l;
            if(t_cigar[0].c=='X') snp_num+=t_cigar[0].l;
            else if((t_cigar[0].c=='I')||(t_cigar[0].c=='D')||(t_cigar[0].c=='N')) snp_num+=1;
        }

        for(;i<t_num;i++)
        {
            if((i==t_num-1)&&(t_cigar[i].c=='S'))
            {
                b_num = 0;
                uint64_t pos = 0;
                tail_align_seed(opt,chr_order,strand,ref_pos-1,&read->seq[read->length-t_cigar[i].l],t_cigar[i].l,b_cigar,&b_num,20,&pos,1);

                if(b_num==0)
                {
                    cigar[cigar_num].c = t_cigar[i].c;
                    cigar[cigar_num].l = t_cigar[i].l;
                    cigar_num++;
                    read->cand[read->cand_num].out_put_flag = 0;
                }
                else
                {
                    for(j = 0;j<b_num;j++)
                    {
                        cigar[cigar_num].c = b_cigar[j].c;
                        cigar[cigar_num].l = b_cigar[j].l;
                        cigar_num++;

                        if((b_cigar[j].c=='M')||(b_cigar[j].c=='X')||(b_cigar[j].c=='D')||(b_cigar[j].c=='N')) ref_pos+=b_cigar[j].l;

                        if(b_cigar[j].c=='X') snp_num+=b_cigar[j].l;
                        else if((b_cigar[j].c=='I')||(b_cigar[j].c=='D')||(b_cigar[j].c=='N')) snp_num+=1;

                        if(b_cigar[j].c=='S') read->cand[read->cand_num].out_put_flag = 0;
                    }
                }
                break;
            }
            if((t_cigar[i].c=='M')||(t_cigar[i].c=='X')||(t_cigar[i].c=='D')||(t_cigar[i].c=='N')) ref_pos+=t_cigar[i].l;

            if(t_cigar[i].c=='X') snp_num+=t_cigar[i].l;
            else if((t_cigar[i].c=='I')||(t_cigar[i].c=='D')||(t_cigar[i].c=='N')) snp_num+=1;

            cigar[cigar_num].c = t_cigar[i].c;
            cigar[cigar_num].l = t_cigar[i].l;
            cigar_num++;
        }
    if(strand) {end_pos = (opt->idx->bns->l_pac<<1)-start_pos-1;start_pos = (opt->idx->bns->l_pac<<1)-ref_pos;}
    else end_pos = ref_pos-1;
    if(snp_num>4)
        read->cand[read->cand_num].out_put_flag = 0;

    //if(snp_num<7)
    {
        score = 0;
        for(i = 0;i<cigar_num;i++)
        {
            if(cigar[i].c=='M') score+=cigar[i].l;
            else if(cigar[i].c=='X') score-=cigar[i].l*opt->miss;
            else if(cigar[i].c=='D') score-=opt->gap;
            else if(cigar[i].c=='N') score-=opt->splice;
            else if(cigar[i].c=='I') score-=opt->miss;
        }
        if(score<(float)read->length*0.88) read->cand[read->cand_num].out_put_flag= 0;

        if(strand)
        {
            char c;
            int l;
            for(i = 0;i<cigar_num/2;i++)
            {
                c = cigar[i].c;
                l = cigar[i].l;

                cigar[i].c = cigar[cigar_num-i-1].c;
                cigar[i].l = cigar[cigar_num-i-1].l;

                cigar[cigar_num-i-1].c = c;
                cigar[cigar_num-i-1].l = l;
            }
        }
        check_read_splice(read->seq,read->rseq,cigar,&cigar_num,start_pos,chr_order,strand);

        int cigar_l = 0;
        int cigar_n = 0;
        char cigarS[MAX_STRING_LENGTH];
        char TAG[MAX_STRING_LENGTH];
        cigarS[0] = '\0';
        TAG[0] = '\0';
        if(read->cand[read->cand_num].out_put_flag==1)
        {
            ref_pos = start_pos;
            cigar_l = 0;
            cigar_n = 0;
            distance = 0;

            for(i = 0;i<cigar_num;i++)
            {
                if((cigar[i].c=='M')||(cigar[i].c=='X'))
                {
                    cigar_l += cigar[i].l;
                    if (cigar[i].c=='M')
                    {
                        if((i+1<cigar_num)&&((cigar[i+1].c=='N')||(cigar[i+1].c=='I'))) {cigar_n += cigar[i].l;continue;}
                        if(cigar_n!=0) {sprintf(temp,"%d",cigar[i].l+cigar_n);cigar_n = 0;}
                        else sprintf(temp,"%d",cigar[i].l);
                        strcat(TAG,temp);
                    }
                    else
                    {
                        distance+=cigar[i].l;
                        snp_num+=cigar[i].l;
                        if(cigar_n!=0) {sprintf(temp,"%d",cigar_n);strcat(TAG,temp);cigar_n = 0;}
                        for(j = 0;j<cigar[i].l;j++)
                        {
                            sprintf(temp,"%c",opt->chr->list[chr_order].seq[j+ref_pos-chr_start]);
                            strcat(TAG,temp);
                        }
                    }
                        ref_pos+=cigar[i].l;
                }
                else
                {
                    if(cigar_l!=0)
                    {
                        sprintf(temp,"%dM",cigar_l);
                        strcat(cigarS,temp);
                    }
                    if(cigar[i].c!='S') snp_num++;
                    sprintf(temp,"%d%c",cigar[i].l,cigar[i].c);
                    strcat(cigarS,temp);
                    cigar_l = 0;
                    if(cigar[i].c=='I')
                    {
                        distance+=cigar[i].l;
                    }
                    if(cigar[i].c=='D')
                    {
                        if(cigar_n!=0) {sprintf(temp,"%d",cigar_n);strcat(TAG,temp);cigar_n = 0;}
                        distance+=cigar[i].l;
                        strcat(TAG,"^");
                        for(j = 0;j<cigar[i].l;j++)
                        {
                            sprintf(temp,"%c",opt->chr->list[chr_order].seq[j+ref_pos-chr_start]);
                            strcat(TAG,temp);
                        }
                        ref_pos+=cigar[i].l;
                    }
                    else if (cigar[i].c=='N') {ref_pos+=cigar[i].l;}
                }
            }
            if(cigar_l!=0)
            {
                sprintf(temp,"%dM",cigar_l);
                strcat(cigarS,temp);
            }
        }
        else
        {
            for(i = 0;i<cigar_num;i++)
            {
                sprintf(temp,"%d%c",cigar[i].l,cigar[i].c);
                strcat(cigarS,temp);
            }
        }
        if(strlen(cigarS)<MAX_READ_LENGTH) strcpy(read->cand[read->cand_num].cigar,cigarS);
        else insert_cand_seed(opt,read,seed,seed_order,seed_order_n,strand,chr_order,exon,snp);
        if(strlen(TAG)<MAX_READ_LENGTH) strcpy(read->cand[read->cand_num].TAG,TAG);
        else insert_cand_seed(opt,read,seed,seed_order,seed_order_n,strand,chr_order,exon,snp);

        break_flag = 0;
        for(i = 0;i<read->cand_num;i++)
        {
            if((start_pos==read->cand[i].pos)&&(strcmp(read->cand[read->cand_num].cigar,read->cand[i].cigar)==0)){break_flag=1;break;}
        }
        if(break_flag) read->cand[read->cand_num].cigar[0] = '\0';
        else
        {
            if(read->cand[read->cand_num].out_put_flag==1)
                read->out_flag = 1;
            //if(((score>(float)read->length*0.9)&&(read->cand[read->cand_num].out_put_flag))||(!read->cand[read->cand_num].out_put_flag))
            {
                read->cand[read->cand_num].strand = strand;
                read->cand[read->cand_num].score = score;
                read->cand[read->cand_num].pos = start_pos;
                read->cand[read->cand_num].end = end_pos;
                read->cand[read->cand_num].chr_order = chr_order;
                read->cand[read->cand_num].dis = distance;
                if(opt->step_flag==SEED_STEP)
                {
                    if(((read->out_flag==0)&&(read->cand[read->cand_num].out_put_flag==0))||((read->out_flag==1)&&(read->cand[read->cand_num].out_put_flag==1)))
                    {if(snp_num<7)find_snp2(opt,read,exon,snp,cigar,cigar_num,start_pos);}
                }
                read->cand_num++;
            }
        }
    }
    //else if(read->out_flag==0)
        //insert_cand_seed(opt,read,seed,seed_order,seed_order_n,strand,chr_order,exon,snp);
    }
}

void unique_seed(struct seed_t *seed,int seed_num,int length)
{
    int i,j;
    int start = 0;

    for(i = 0;i<seed_num;i++)
    {
        for(j = start;j<i;j++)
        {
            if(seed[start].start==seed[i].start)

            if((seed[i].start==seed[j].start)&&(seed[i].length==seed[j].length)&&(seed[i].length!=length))
            {
                if((seed[i].pos+AREA>seed[j].pos)&&(((seed[i].pos>AREA)?(seed[i].pos-AREA):0)<seed[j].pos))
                {seed[j].flag = 1;seed[i].flag = 1;}
            }
        }
    }
    for(i = 0;i<seed_num;i++)
    {
        if(seed[i].flag == 1)
        {
            seed[i].start = 0;
            seed[i].pos = 0;
            seed[i].length = 0;
            seed[i].abs = 0;
            seed[i].score = 0;
            seed[i].lnum = 0;
        }
    }
}
void extend_seed_tail_forward_seed(struct m_opt *opt,int chr_order,char strand,char *seq,struct seed_t *seed,int *seed_num,uint64_t pos,int Sstart)
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
        tail_align_seed(opt,chr_order,strand,pos,seq,Sstart,t_cigar,&t_num,20,&ref_pos,0);
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
                seed_r.last[0] = chr_order;

                insert_seed(seed,seed_num,SEED_BUF_LENGTH,seed_r);

                if((*seed_num)>init_num)
                    extend_seed_tail_forward_seed(opt,chr_order,strand,seq,seed,seed_num,seed_r.pos,seed_r.start);
            }
            if((t_cigar[j].c=='M')||(t_cigar[j].c=='X')||(t_cigar[j].c=='D')||(t_cigar[j].c=='N')) ref_pos+=t_cigar[j].l;
            if((t_cigar[j].c=='M')||(t_cigar[j].c=='X')||(t_cigar[j].c=='I')||(t_cigar[j].c=='S')) read_pos+=t_cigar[j].l;
        }
    }
}
void extend_seed_tail_backward_seed(struct m_opt *opt,int chr_order,char strand,char *seq,int length,struct seed_t *seed,int *seed_num,uint64_t pos,int Sstart)
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
        tail_align_seed(opt,chr_order,strand,ref_pos-1,&seq[read_pos],length-read_pos,t_cigar,&t_num,20,&post,1);
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
                seed_r.last[0] = chr_order;

                insert_seed(seed,seed_num,SEED_BUF_LENGTH,seed_r);

                if((*seed_num)>init_num)
                    extend_seed_tail_backward_seed(opt,chr_order,strand,seq,length,seed,seed_num,seed_r.pos+seed_r.length,seed_r.start+seed_r.length);
            }
            if((t_cigar[j].c=='M')||(t_cigar[j].c=='X')||(t_cigar[j].c=='D')||(t_cigar[j].c=='N')) ref_pos+=t_cigar[j].l;
            if((t_cigar[j].c=='M')||(t_cigar[j].c=='X')||(t_cigar[j].c=='I')||(t_cigar[j].c=='S')) read_pos+=t_cigar[j].l;
        }
    }
}
int seed_pos_cmp(const void *a,const void *b)
{
    struct seed_t *EA,*EB;
    EA = (struct seed_t *)a;
    EB = (struct seed_t *)b;

    if(EA->pos==EB->pos)
    {
        if (EA->start==EB->start) return 0;
        else if (EA->start<EB->start) return -1;
        else  return 1;
    }
    else
    {
        if(EA->pos < EB->pos) return -1;
        else return 1;
    }
}
void find_cand(struct m_opt *opt,struct read_inf_t *read,struct seed_t *seed,int *seed_num,struct exon_array *exon,struct snp_list_t *snp,float score_T)
{
    int max_o = -1;
    int max_s = 0;
    int seed_order[30];
    int seed_order_n = 0;

    int i = 0,j = 0;
    int x = 0;
    int chr_order = 0;
    unsigned int chr_start = 0;

    uint64_t LB,RB;

    int init_num = *seed_num;
    char strand = 0;

    int start = 0;

    //unique_seed(seed,*seed_num,read->length);
    for(i = 0; i<(*seed_num); i++)
    {
        seed[i].flag = 0;
        if((seed[i].start == 0)&&(seed[i].pos == 0)&&(seed[i].length == 0)) continue;

        if(seed[i].pos>=opt->idx->bns->l_pac)
        {
            strand = 1;
            for(j = 0;j<opt->chr->total;j++)
            {
                if((seed[i].pos >= (opt->idx->bns->l_pac<<1)-opt->chr->list[j].start_site-opt->chr->list[j].length)
                    &&(seed[i].pos < (opt->idx->bns->l_pac<<1)-opt->chr->list[j].start_site))
                    break;
            }
        }
        else{
            strand = 0;
            for(j = 0;j<opt->chr->total;j++)
            {
                if((seed[i].pos >= opt->chr->list[j].start_site)
                    &&(seed[i].pos <opt->chr->list[j].start_site+opt->chr->list[j].length))
                    break;
            }
        }
        if((strand == 0)&&(seed[i].pos+seed[i].length>=opt->chr->list[j].start_site+opt->chr->list[j].length))
        {seed[i].start = 0;seed[i].pos = 0;seed[i].length = 0;seed[i].abs = 0;seed[i].score = 0;seed[i].lnum = 0;continue;}
        if((strand == 1)&&((opt->idx->bns->l_pac<<1)-seed[i].pos-seed[i].length<opt->chr->list[j].start_site))
        {seed[i].start = 0;seed[i].pos = 0;seed[i].length = 0;seed[i].abs = 0;seed[i].score = 0;seed[i].lnum = 0;continue;}
        //seed_sxtend
        seed[i].last[0] = j;
        chr_order = j;
        chr_start = opt->chr->list[j].start_site;
        int front = 0;
        front = seed[i].start;
        for(x = 0;x<front;x++)
        {
            if(seed[i].pos>=opt->idx->bns->l_pac)
            {
                if((seed[i].pos < (opt->idx->bns->l_pac<<1)-opt->chr->list[chr_order].start_site-opt->chr->list[chr_order].length)
                    ||(seed[i].pos >= (opt->idx->bns->l_pac<<1)-opt->chr->list[chr_order].start_site)) break;
                if ((3-nst_nt4_table[(int)opt->chr->list[chr_order].seq[(opt->idx->bns->l_pac<<1)-seed[i].pos-chr_start]])!=(nst_nt4_table[(int)read->seq[seed[i].start-1]])) break;
                else {seed[i].pos--;seed[i].start--;seed[i].length++;}
            }
            else
            {
                if((seed[i].pos < opt->chr->list[chr_order].start_site)
                    ||(seed[i].pos >=opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length))
                    break;
                if(nst_nt4_table[(int)opt->chr->list[chr_order].seq[seed[i].pos-1-chr_start]]!=nst_nt4_table[(int)read->seq[seed[i].start-1]])break;
                else {seed[i].pos--;seed[i].start--;seed[i].length++;}
            }
        }
        for(x = seed[i].start+seed[i].length;x<read->length;x++)
        {
            if(seed[i].pos>=opt->idx->bns->l_pac)
            {
                if((seed[i].pos+seed[i].length < (opt->idx->bns->l_pac<<1)-opt->chr->list[chr_order].start_site-opt->chr->list[chr_order].length)
                    ||(seed[i].pos+seed[i].length >= (opt->idx->bns->l_pac<<1)-opt->chr->list[chr_order].start_site)) break;
                if((3-nst_nt4_table[(int)opt->chr->list[chr_order].seq[(opt->idx->bns->l_pac<<1)-seed[i].pos-seed[i].length-1-chr_start]])!=(nst_nt4_table[(int)read->seq[x]])) break;
                else seed[i].length++;
            }
            else
            {
                if((seed[i].pos+seed[i].length < opt->chr->list[chr_order].start_site)
                    ||(seed[i].pos+seed[i].length >=opt->chr->list[chr_order].start_site+opt->chr->list[chr_order].length))
                    break;
                if(nst_nt4_table[(int)opt->chr->list[chr_order].seq[seed[i].pos+seed[i].length-chr_start]]!=nst_nt4_table[(int)read->seq[x]]) break;
                else seed[i].length++;
            }
        }
        while((seed[i].start >seed[start].start)&&(start<(*seed_num))) start++;
        for(j = start;j<i;j++)
        {
            if(seed[j].length == 0)continue;
            if((seed[i].start == seed[j].start)&&(seed[i].pos == seed[j].pos))
            {seed[i].start = 0;seed[i].pos = 0;seed[i].length = 0;seed[i].abs = 0;seed[i].score = 0;seed[i].lnum = 0;break;}
        }
        if((seed[i].start == 0)&&(seed[i].pos == 0)&&(seed[i].length == 0)) continue;
        if(i<init_num)
        {
            extend_seed_tail_forward_seed(opt,chr_order,strand,read->seq,seed,seed_num,seed[i].pos,seed[i].start);
            extend_seed_tail_backward_seed(opt,chr_order,strand,read->seq,read->length,seed,seed_num,seed[i].pos+seed[i].length,seed[i].start+seed[i].length);
        }
    }
    //qsort(seed,*seed_num,sizeof(struct seed_t),seed_cmp);
    qsort(seed,*seed_num,sizeof(struct seed_t),seed_pos_cmp);
    int k = 0;
    while((k<(*seed_num))&&(seed[k].start == 0)&&(seed[k].pos == 0)&&(seed[k].length == 0)) k++;
    start = k;

    for(i = k; i<*seed_num; i++)
    {
        if((seed[i].start == 0)&&(seed[i].pos == 0)&&(seed[i].length == 0)) continue;
        seed[i].abs = seed[i].pos-seed[i].start;
        RB=seed[i].pos+opt->change_length;//Insert
        //LB
        if(seed[i].pos>=opt->idx->bns->l_pac)
            LB = max(((opt->idx->bns->l_pac<<1)-opt->chr->list[seed[i].last[0]].start_site-opt->chr->list[seed[i].last[0]].length),(seed[i].pos>opt->area?seed[i].pos-opt->area:0));
        else
            LB = max(opt->chr->list[seed[i].last[0]].start_site,(seed[i].pos>opt->area?seed[i].pos-opt->area:0));

        //从头开始对每个seed赋分，即本seed长度+前置最优路径长度
        max_o = -1;
        max_s = -1000;
        seed[i].lnum = 0;
        seed[i].lorder = 0;
        seed[i].score = 0;
        while((LB >seed[start].pos)&&(start<(*seed_num))) start++;
        for(j = start;j<i;j++)
        //for(j = k;j<i;j++)
        {
            if(seed[j].length == 0)continue;
            //if((seed[i].start == seed[j].start)&&(seed[i].pos == seed[j].pos))
            //{seed[i].start = 0;seed[i].pos = 0;seed[i].length = 0;seed[i].score = 0;seed[i].lnum = 0;break;}

            if((seed[i].start > seed[j].start)
                &&((seed[i].start+seed[i].length)>(seed[j].start+seed[j].length))
                &&(RB>=(seed[j].pos+seed[i].start-seed[j].start))
                &&(LB<=(seed[j].pos+seed[i].start-seed[j].start))
             )
            {
                if(seed[i].start<=seed[j].start+seed[j].length)
                    max_s = seed[j].score + (seed[i].start+seed[i].length-seed[j].start-seed[j].length);
                else
                    max_s = seed[j].score + seed[i].length;

                if(seed[i].abs<seed[j].abs)
                    max_s -=(seed[j].abs-seed[i].abs)+2;
                else if((seed[i].abs>seed[j].abs)&&(seed[i].abs-seed[j].abs<opt->change_length))
                    max_s -=opt->gap;
                else if((seed[i].abs>seed[j].abs)&&(seed[i].abs-seed[j].abs>opt->change_length))
                    max_s -=opt->splice;

                if(seed[i].score < max_s)
                {
                    seed[i].score = max_s;
                    seed[i].lnum = 0;
                    seed[i].last[seed[i].lnum] = j;
                    seed[i].lnum++;
                }
                else if((seed[i].score == max_s)&&(seed[i].lnum<MAX_SEED_LAST))
                {
                    seed[i].last[seed[i].lnum] = j;
                    seed[i].lnum++;
                }
            }
        }
        if(seed[i].lnum == 0)
            seed[i].score = seed[i].length;
    }
    //回溯
    max_o = 0;
    max_s = -1000;
    int flag = 0;

    {
        int best[MAX_CAND_NUM];
        int best_num = 0;
        int mix_o = 0,mix_s =200;

        for(i = *seed_num-1; i>=0; i--)
        {
            if(seed[i].score>max_s) max_s = seed[i].score;
            if(seed[i].score>=score_T)
            {
                if(best_num<MAX_CAND_NUM)
                {
                    if(seed[i].score<mix_s){mix_s = seed[i].score;mix_o = best_num;}
                    best[best_num] = i;best_num++;
                }
                else
                {
                    if(seed[i].score<mix_s)continue;
                    else
                    {
                        best[mix_o] = i;
                        mix_s = seed[i].score;
                        for(j = 0;j<MAX_CAND_NUM;j++)
                        {
                            if(seed[best[j]].score<mix_s)
                            {mix_o = j;mix_s =seed[best[j]].score;}
                        }
                    }
                }
            }
        }
        read->out_flag = 0;
        for(i = 0;i<best_num;i++)
        {
            mix_o = best[0];
            for(j = 0;j<best_num;j++)
            {
                if((seed[best[j]].score>seed[mix_o].score)) mix_o = best[j];
            }
            seed_order_n = 0;

            //if(max_s-seed[mix_o].score>20){seed[mix_o].score = 0;continue;}

            max_o = mix_o;
            seed_order[seed_order_n] = max_o;
            seed_order_n++;
            while(seed_order_n>0)
            {
                max_o = seed_order[seed_order_n-1];
                flag = 1&seed[max_o].flag;
                seed[max_o].flag = 1;
                while(seed[max_o].lorder!=seed[max_o].lnum)
                {
                    seed[max_o].lorder++;
                    max_o = seed[max_o].last[seed[max_o].lorder-1];
                    flag = flag&seed[max_o].flag;
                    seed[max_o].flag = 1;

                    seed_order[seed_order_n] = max_o;
                    seed_order_n++;
                    if(seed_order_n>=30) break;
                }
                if(read->cand_num>=MAX_CAND_NUM) break;
                if(flag==0)
                generate_cigar(opt,read,seed,seed_order,seed_order_n,exon,snp);

                while((seed_order_n>0)&&(seed[seed_order[seed_order_n-1]].lorder==seed[seed_order[seed_order_n-1]].lnum))
                    {seed[seed_order[seed_order_n-1]].lorder=0;seed_order_n--;}
                flag = 1;
            }
            seed[mix_o].score = 0;
        }
    }
}
void find_pair1(struct m_opt *opt,struct read_inf_t *read1,struct read_inf_t *read2)
{
    int i = 0,j = 0;

    if((opt->output_mode==OUTPUT_BEST)&&(read1->out_flag==1)&&(read2->out_flag==1))
    {
        for(i = 0;i<read2->cand_num;i++)
        {
            read2->cand[i].pair = -1;
            read2->cand[i].pair_score = 0;
        }
        for(i = 0;i<read1->cand_num;i++)
        {
            read1->cand[i].pair = -1;
            read1->cand[i].pair_score = 0;

            if(read1->cand[i].out_put_flag==0)continue;

            for(j = 0;j<read2->cand_num;j++)
            {
                if(read2->cand[j].out_put_flag==0)continue;

                if(read1->cand[i].chr_order!=read2->cand[j].chr_order) continue;
                if(read1->cand[i].strand==read2->cand[j].strand) continue;

                if((read1->cand[i].pos>=read2->cand[j].pos-opt->area)&&(read1->cand[i].pos<=read2->cand[j].pos+opt->area))
                {
                    if(read2->cand[j].score > read1->cand[i].pair_score)
                    {
                        read1->cand[i].pair = j;
                        read1->cand[i].pair_score = read2->cand[j].score;
                    }
                    if(read1->cand[i].score > read2->cand[j].pair_score)
                    {
                        read2->cand[j].pair = i;
                        read2->cand[j].pair_score = read1->cand[i].score;
                    }
                }
            }
        }
    }
    else
    {
        for(i = 0;i<read2->cand_num;i++)
        {
            read2->cand[i].pair = -1;
            read2->cand[i].pair_score = 0;
        }
        for(i = 0;i<read1->cand_num;i++)
        {
            read1->cand[i].pair = -1;
            read1->cand[i].pair_score = 0;

            for(j = 0;j<read2->cand_num;j++)
            {
                if(read1->cand[i].chr_order!=read2->cand[j].chr_order) continue;
                if(read1->cand[i].strand==read2->cand[j].strand) continue;

                if((read1->cand[i].pos>=read2->cand[j].pos-opt->area)&&(read1->cand[i].pos<=read2->cand[j].pos+opt->area))
                {
                    if(read2->cand[j].score > read1->cand[i].pair_score)
                    {
                        read1->cand[i].pair = j;
                        read1->cand[i].pair_score = read2->cand[j].score;
                    }
                    if(read1->cand[i].score > read2->cand[j].pair_score)
                    {
                        read2->cand[j].pair = i;
                        read2->cand[j].pair_score = read1->cand[i].score;
                    }
                }
            }
        }
    }
}
void find_pair(struct m_opt *opt,struct read_inf_t *read1,struct read_inf_t *read2)
{
    int i = 0,j = 0;

    if((read1->out_flag==1)&&(read2->out_flag==1))
    {
        for(i = 0;i<read2->cand_num;i++)
        {
            read2->cand[i].pair = -1;
            read2->cand[i].pair_score = 0;
        }
        for(i = 0;i<read1->cand_num;i++)
        {
            read1->cand[i].pair = -1;
            read1->cand[i].pair_score = 0;

            for(j = 0;j<read2->cand_num;j++)
            {
                if(read1->cand[i].chr_order!=read2->cand[j].chr_order) continue;
                if(read1->cand[i].strand==read2->cand[j].strand) continue;

                if(((read1->cand[i].end>=read2->cand[j].pos-1000)&&(read1->cand[i].end<=read2->cand[j].pos))
                    ||((read1->cand[i].pos>=read2->cand[j].end)&&(read1->cand[i].pos<=read2->cand[j].end+1000)))
                //if(((read1->cand[i].pos>=read2->cand[j].pos-1000)&&(read1->cand[i].pos<=read2->cand[j].end+1000))
                    //||((read1->cand[i].end>=read2->cand[j].pos-1000)&&(read1->cand[i].end<=read2->cand[j].end+1000)))
                {
                    if(read2->cand[j].score > read1->cand[i].pair_score)
                    {
                        read1->cand[i].pair = j;
                        read1->cand[i].pair_score = read2->cand[j].score;
                    }
                    if(read1->cand[i].score > read2->cand[j].pair_score)
                    {
                        read2->cand[j].pair = i;
                        read2->cand[j].pair_score = read1->cand[i].score;
                    }
                }
            }
        }
    }
    else
    {
        for(i = 0;i<read2->cand_num;i++)
        {
            read2->cand[i].pair = -1;
            read2->cand[i].pair_score = 0;
        }
        for(i = 0;i<read1->cand_num;i++)
        {
            read1->cand[i].pair = -1;
            read1->cand[i].pair_score = 0;

            for(j = 0;j<read2->cand_num;j++)
            {
                if(read1->cand[i].chr_order!=read2->cand[j].chr_order) continue;
                if(read1->cand[i].strand==read2->cand[j].strand) continue;

                if((read1->cand[i].pos>=read2->cand[j].pos-opt->area)&&(read1->cand[i].pos<=read2->cand[j].pos+opt->area))
                {
                    if(read2->cand[j].score > read1->cand[i].pair_score)
                    {
                        read1->cand[i].pair = j;
                        read1->cand[i].pair_score = read2->cand[j].score;
                    }
                    if(read1->cand[i].score > read2->cand[j].pair_score)
                    {
                        read2->cand[j].pair = i;
                        read2->cand[j].pair_score = read1->cand[i].score;
                    }
                }
            }
        }
    }
}
void find_pair_cand(struct m_opt *opt,struct read_inf_t *read1,struct read_inf_t *read2,struct exon_array *exon,struct snp_list_t *snp)
{
    uint64_t pos;
    uint64_t LB,RB;

    int i = 0,j = 0,x = 0;//back_step = 0;;

    int c = 0;
    bwtint_t k, l, ok, ol;
    struct seed_a seed_w[MAX_READ_LENGTH];

    int strand;

    struct cigar_t cigar;
    cigar.c = 'M';
    cigar.l = read2->length;
    int flag = 0;

            //bwt_find_seed;
    k = 0; l = opt->idx->bwt->seq_len;
    for (j = read2->length-1;j>=0;j--)
    {
        if(nst_nt4_table[(int)read2->seq[j]]>=4) {break;}
        c = nst_nt4_table[(int)read2->seq[j]];
        bwt_2occ(opt->idx->bwt, k - 1, l, c, &ok, &ol);
        k = opt->idx->bwt->L2[c] + ok + 1;
        l = opt->idx->bwt->L2[c] + ol;
        seed_w[j].start = k;
        seed_w[j].end = l;
    }
    if(j == -1)
    {
        flag = 0;
        if (k <= l)
        {
            for (x = 0;x<=seed_w[j+1].end-seed_w[j+1].start;x++)
            {
                pos = bwt_sa(opt->idx->bwt,seed_w[j+1].start + x);
                strand = 0;
                if(pos>=opt->idx->bns->l_pac) {pos = (opt->idx->bns->l_pac<<1)-read2->length-pos;strand = 1;}
                for(i = 0;i<read1->cand_num;i++)
                {
                    LB = read1->cand[i].pos-opt->area;
                    RB = read1->cand[i].pos+opt->area;
                    if((pos>LB)&&(pos<RB))
                    {
                        for(k = 0;k<opt->chr->total;k++)
                        {
                            if((pos >= opt->chr->list[k].start_site)
                                &&(pos <opt->chr->list[k].start_site+opt->chr->list[k].length))
                                break;
                        }
                        read2->cand[read2->cand_num].chr_order = k;
                        sprintf(read2->cand[read2->cand_num].cigar,"100M");
                        sprintf(read2->cand[read2->cand_num].TAG,"100");
                        read2->cand[read2->cand_num].out_put_flag = 1;
                        read2->cand[read2->cand_num].strand = strand;
                        read2->cand[read2->cand_num].score = read2->length;
                        read2->cand[read2->cand_num].pos = pos;
                        read2->cand[read2->cand_num].dis = 0;

                        find_snp2(opt,read2,exon,snp,&cigar,1,pos);
                        read2->cand_num++;
                        if(read2->cand_num>=MAX_CAND_NUM)
                        {
                            flag = 1;
                            read2->out_flag = 0;
                            read2->cand_num = 0;
                            break;
                        }
                        break;
                    }
                }
                if(flag) break;
            }
        }
    }
    if(read2->cand_num==0) read2->out_flag = 0;
}
int find_s(struct snp_list_t *snp,unsigned int start)
{
    int left = 0, right = snp->total-1, middle;

    if (right == -1)
        return -1;
    while (left <= right)
    {
        middle = (left + right)/2;

        if(snp->snp[middle].start==start)
        {
            while((snp->snp[middle-1].start==start)&&(middle>0)) middle--;
            return middle;
        }
        else if (snp->snp[middle].start>start)
        {
            right = middle -1;
        }
        else
            left = middle + 1;
    }
    return left;
}
int find_exon(struct exon_array *exon,unsigned int start)
{
    int left = 0, right = exon->total-1, middle;

    if (right == -1)
        return -1;
    while (left <= right)
    {
        middle = (left + right)/2;

        if((exon->exon[middle].start<=start)&&(exon->exon[middle].end>=start))
        { return middle; }
        else if (exon->exon[middle].start > start)
        {
            right = middle -1;
        }
        else
            left = middle + 1;
    }
    return left;
}
void out_put_SAM_head(struct m_opt *opt,FILE *file)
{
    fprintf(file,"@PG\tID:DEEP\tVN:2018-10-10\n");
    int i = 0;
    for(i = 0;i<opt->chr->total;i++)
        fprintf(file,"@SG\tSN:%s\tLN:%u\n",opt->chr->list[i].name,opt->chr->list[i].length);
}
void *seed_align_core(void* arg)
{
    struct job_seed *job = (struct job_seed *)arg;

    struct snp_list_t *snp = (struct snp_list_t *)calloc(1,sizeof(struct snp_list_t));
    snp->snp = (struct snp_t *)calloc(MAX_SNP_NUM,sizeof(struct snp_t));
    snp->total = 0;

    struct exon_array *exon = (struct exon_array *)calloc(1,sizeof(struct exon_array));
    exon->exon = (struct exon_inf_t *)calloc(READ_BUF_LENGTH*3,sizeof(struct exon_inf_t));
    exon->total = 0;

	struct seed_t *seed = (struct seed_t *)calloc(SEED_BUF_LENGTH+1,sizeof(struct seed_t));
	int seed_num;
	struct read_inf_t *read = (struct read_inf_t *)calloc(READ_BUF_LENGTH+1,sizeof(struct read_inf_t));
	int read_num = 0;


	int i = 0;
	int break_flag = 0;

	while (1)
	{
        //read_file;
        pthread_mutex_lock(&ReadLock);
        read_num = read_file(read,opt,job->in_file_1,job->in_file_2);
		while (read_num == 0)
		{
            if(job->opt->file_flag >= (job->opt->input_file_1->total-1)){break_flag = 1;break;}
            else
            {
                job->opt->file_flag++;

                fclose(job->in_file_1);
                if(job->opt->pair == 1) fclose(job->in_file_2);

                job->in_file_1 = fopen(job->opt->input_file_1->file[job->opt->file_flag].name,"r");
                if(job->opt->pair == 1)
                {job->in_file_2 = fopen(job->opt->input_file_2->file[job->opt->file_flag].name,"r");}

                read_num = read_file(read,job->opt,job->in_file_1,job->in_file_2);
            }
		}
		pthread_mutex_unlock(&ReadLock);

		if(break_flag) break;

        //seed_align
        for (i = 0;i<read_num;i++)
        {
            if(job->opt->pair)
            {
                seed_num = 0;
                find_seed(job->opt,read+i,seed,&seed_num);
                if(read[i].out_flag==0)
                    find_cand(job->opt,read+i,seed,&seed_num,exon,snp,0.5*strlen(read->seq));
                if(read[i].out_flag==0)
                {
                    read[i].cand_num=0;
                    find_seed_r(job->opt,read+i,seed,&seed_num);
                    if(read[i].out_flag==0)
                    find_cand(job->opt,read+i,seed,&seed_num,exon,snp,0.35*strlen(read->seq));
                }

                //if(read[i].out_flag==0)//if(read[i].cand_num==0)
                //{
                    //find_seed_s(job->opt,read+i,seed,&seed_num);
                    //find_cand(job->opt,read+i,seed,&seed_num,exon,snp,0.35*strlen(read->seq));
                //}

                seed_num = 0;
                find_seed_r(job->opt,read+i+1,seed,&seed_num);
                if(read[i+1].out_flag==0)
                find_cand(job->opt,read+i+1,seed,&seed_num,exon,snp,0.5*strlen(read->seq));
                if(read[i+1].out_flag==0)
                {
                    read[i+1].cand_num=0;
                    find_seed(job->opt,read+i+1,seed,&seed_num);
                    if(read[i+1].out_flag==0)
                    find_cand(job->opt,read+i+1,seed,&seed_num,exon,snp,0.35*strlen(read->seq));
                }

                //if(read[i+1].out_flag==0)//if(read[i+1].cand_num==0)
                //{
                    //find_seed_s(job->opt,read+i+1,seed,&seed_num);
                    //find_cand(job->opt,read+i+1,seed,&seed_num,exon,snp,0.35*strlen(read->seq));
                //}
                if((read[i].out_flag==1)&&(read[i].cand_num ==0))
                {
                    if(read[i+1].cand_num!=0)find_pair_cand(job->opt,read+i+1,read+i,exon,snp);
                    else read[i].out_flag=0;
                }

                if((read[i+1].out_flag==1)&&(read[i+1].cand_num ==0))
                {
                    if(read[i].cand_num!=0) find_pair_cand(job->opt,read+i,read+i+1,exon,snp);
                    else read[i+1].out_flag = 0;
                }


                find_pair(job->opt,read+i,read+i+1);
                i++;
            }
            else
            {
                seed_num = 0;
                find_seed(job->opt,read+i,seed,&seed_num);
                find_cand(job->opt,read+i,seed,&seed_num,exon,snp,0.5*strlen(read->seq));
                if(read[i].out_flag==0)
                {
                    read[i].cand_num=0;
                    find_seed_r(job->opt,read+i,seed,&seed_num);
                    find_cand(job->opt,read+i,seed,&seed_num,exon,snp,0.35*strlen(read->seq));
                }
            }
        }

        pthread_mutex_lock(&OutputLock);
        out_put_read(job->opt,read,read_num,job->out_file_1,job->unmap_file_1);
        Update_Exon(job->opt,exon);
        pthread_mutex_unlock(&OutputLock);
	}
    Update_Snp(job->opt,snp);

    free(snp->snp);
    free(snp);

    free(exon->exon);
    free(exon);

	free(seed);
	free(read);

    return (void*)(0);
}
int seed_align(struct m_opt *opt)
{
    //open_file input output unmap
    FILE *in_file_1,*in_file_2;
    FILE *out_file_1;
    FILE *unmap_file_1;
    char out_file_name_1[MAX_NAME_LENGTH+50];
    char unmap_file_name_1[MAX_NAME_LENGTH+50];

    in_file_1 = fopen(opt->input_file_1->file[0].name,"r");
    if (in_file_1 == NULL)
        return 1;
    if(opt->pair == 1)
    {
        in_file_2 = fopen(opt->input_file_2->file[0].name,"r");
        if (in_file_2 == NULL)
            return 1;
    }

    sprintf(out_file_name_1,"%s/map_out.sam",opt->Output_path);
    out_file_1 = fopen(out_file_name_1,"w");
    if (out_file_1 == NULL)
        return 1;

    sprintf(unmap_file_name_1,"%s/unmap_out.sam",opt->Temp_path);
    unmap_file_1 = fopen(unmap_file_name_1,"w");
    if (unmap_file_1 == NULL)
        return 1;

    opt->snp_fn = 0;

    out_put_SAM_head(opt,out_file_1);

    struct job_seed *job = (struct job_seed *)calloc(1,sizeof(struct job_seed));
    job->opt = opt;
    job->in_file_1 = in_file_1;
    job->in_file_2 = in_file_2;
    job->out_file_1 = out_file_1;
    job->unmap_file_1 = unmap_file_1;

	pthread_t *pthreads = malloc(sizeof(pthread_t) * opt->thread_num);

	int i = 0;
    for (i = 0; i < opt->thread_num; i++) pthread_create(&pthreads[i], NULL, seed_align_core, job);
	for (i = 0; i < opt->thread_num; i++) pthread_join(pthreads[i], NULL);

    fclose(in_file_1);
    if(opt->pair == 1)fclose(in_file_2);
    fclose(out_file_1);
    fclose(unmap_file_1);

    free(pthreads);
    free(job);

    opt->input_file_1->total = 1;
    strcpy(opt->input_file_1->file[0].name,unmap_file_name_1);
    return 0;
}

int read_file_sam_sta(struct read_map_inf *read,FILE *file)
{
    int i = 0,j;
    char f_line[MAX_STRING_LENGTH];
    char chr[40];

    while((i<READ_BUF_LENGTH)&&fgets(f_line,MAX_STRING_LENGTH,file)!=NULL)
    {
        if(f_line[0]=='@') continue;
        sscanf(f_line,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%*d\t%s\t%s",
                read[i].name,&read[i].flag,chr,&read[i].start_site,&read[i].MapQ,read[i].cigar,
                read[i].pname,&read[i].psite,read[i].seq,read[i].qual);
        read[i].length = strlen(read[i].seq);
        read[i].porder = 0;
        for(j = 0;j<opt->chr->total;j++)
        {
            if(strcmp(opt->chr->list[j].name,chr)==0) {read[i].chr_order=j;break;}
        }
        i++;
    }

    return i;
}
void out_put_read2(struct read_map_inf *read,int read_num,FILE *out,FILE *un)
{
    int i = 0;

    for (i = 0;i<read_num;i++)
    {
        if(read[i].porder==1)
            fprintf(un,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%d\n",
                read[i].name,read[i].flag,opt->chr->list[read[i].chr_order].name,(read[i].start_site==0)?0:read[i].start_site-1+opt->chr->list[read[i].chr_order].start_site,read[i].MapQ,read[i].cigar,
                read[i].pname,read[i].psite,read[i].seq,read[i].qual,read[i].score);
        else
            fprintf(out,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t0\t%s\t%s\t%d\n",
                read[i].name,read[i].flag,opt->chr->list[read[i].chr_order].name,read[i].start_site,read[i].MapQ,read[i].cigar,
                read[i].pname,read[i].psite,read[i].seq,read[i].qual,read[i].score);
    }
}
void find_snp(struct m_opt *opt,struct read_map_inf *read,struct exon_array *exon,struct snp_list_t *snp)
{
    struct snp_t snp_r;
    unsigned int exon_start = 0,exon_end = 0;

    unsigned int ref_site;
    int read_site = 0;

    unsigned int i = 0,j = 0,k = 0;
    struct cigar_t cigar[50];
    int cigar_total = 0;
    char temp[50];

    struct cigar_t t_cigar[50];
    int t_total = -1;

    int length = 0;
    read->score = 0;

    if((read->start_site!=0)&&(read->cigar[0]!='*'))
    {
        j = 0;
        cigar_total = 0;
        length = strlen(read->cigar);
        for(i = 0;i<length;i++)
        {
            if((read->cigar[i]>='0')&&(read->cigar[i]<='9'))
            {
                temp[j] = read->cigar[i];
                j++;
            }
            else
            {
                temp[j] = '\0';
                cigar[cigar_total].c = read->cigar[i];
                cigar[cigar_total].l = atoi(temp);
                cigar_total++;
                j = 0;
            }
        }

        ref_site = exon_start = read->start_site+opt->chr->list[read->chr_order].start_site-1;
        read_site = 0;
        for(i = 0;i<cigar_total;i++)
        {
            switch(cigar[i].c)
            {
            case 'M':
                read->score+=cigar[i].l;
                k = 0;
                for(j = 0;j<cigar[i].l;j++)
                    {
                        if(nst_nt4_table[(int)read->seq[read_site+j]]!=nst_nt4_table[(int)opt->chr->list[read->chr_order].seq[ref_site+j-opt->chr->list[read->chr_order].start_site]])
                        {
                            if(k==0)
                            {
                                t_total++;
                                t_cigar[t_total].c='X';
                                t_cigar[t_total].l = 1;

                                snp_r.type = 'X';
                                snp_r.start = ref_site+j;
                                snp_r.seq = 0;
                                snp_r.seq = (snp_r.seq<<2)|nst_nt4_table[(int)read->seq[read_site+j]];
                                snp_r.length = 1;
                                snp_r.num = 1;
                                snp_r.dep = 0;
                                k = 1;
                            }
                            else
                            {snp_r.seq = (snp_r.seq<<2)|nst_nt4_table[(int)read->seq[read_site+j]];
                            snp_r.length++;
                            t_cigar[t_total].l++;}
                        }
                        else
                        {
                            if((t_total==-1)||(t_cigar[t_total].c!='M')){t_total++;t_cigar[t_total].c='M';t_cigar[t_total].l = 1;}
                            else t_cigar[t_total].l++;
                            if(k)
                            {
                                read->score-=(opt->miss+opt->match)*snp_r.length;
                                snp_r.end = ref_site+j-1;
                                insert_snp(snp,snp_r);
                                k = 0;
                            }
                        }
                    }

                    ref_site+=cigar[i].l;
                    read_site+=cigar[i].l;
                    break;
                case 'X':
                    t_total++;
                    t_cigar[t_total].c='X';
                    t_cigar[t_total].l=cigar[i].l;

                    read->score-=cigar[i].l*opt->miss;
                    snp_r.start = ref_site;
                    snp_r.end = ref_site+cigar[i].l-1;
                    snp_r.type = 'X';
                    snp_r.seq = 0;
                    snp_r.length = cigar[i].l;
                    snp_r.num = 1;
                    snp_r.dep = 0;

                    for(j = 0;j<cigar[i].l;j++)
                        snp_r.seq = (snp_r.seq<<2)|nst_nt4_table[(int)read->seq[read_site+j]];

                    insert_snp(snp,snp_r);

                    ref_site+=cigar[i].l;
                    read_site+=cigar[i].l;
                    break;
                case 'I':
                    t_total++;
                    t_cigar[t_total].c='I';
                    t_cigar[t_total].l=cigar[i].l;

                    read->score-=opt->miss;
                    snp_r.start = ref_site;
                    snp_r.end = ref_site;
                    snp_r.type = 'I';
                    snp_r.seq = 0;
                    snp_r.length = cigar[i].l;
                    snp_r.num = 1;
                    snp_r.dep = 0;

                    for(j = 0;j<cigar[i].l;j++)
                        snp_r.seq = (snp_r.seq<<2)|nst_nt4_table[(int)read->seq[read_site+j]];

                    insert_snp(snp,snp_r);

                    read_site+=cigar[i].l;
                    break;
                case 'D':
                    t_total++;
                    t_cigar[t_total].c='D';
                    t_cigar[t_total].l=cigar[i].l;

                    read->score-=opt->gap;
                    snp_r.start = ref_site;
                    snp_r.end = ref_site+cigar[i].l-1;
                    snp_r.type = 'D';
                    snp_r.seq = 0;
                    snp_r.length = cigar[i].l;
                    snp_r.num = 1;
                    snp_r.dep = 0;
                    insert_snp(snp,snp_r);

                    ref_site+=cigar[i].l;
                    break;
                case 'N':
                    t_total++;
                    t_cigar[t_total].c='N';
                    t_cigar[t_total].l=cigar[i].l;

                    read->score-=opt->gap;
                    snp_r.start = ref_site;
                    snp_r.end = ref_site+cigar[i].l-1;
                    snp_r.type = 'N';
                    snp_r.seq = 0;
                    snp_r.length = cigar[i].l;
                    snp_r.num = 1;
                    snp_r.dep = 0;
                    insert_snp(snp,snp_r);

                    exon_end = ref_site-1;
                    insert_exon(exon,snp,exon_start,exon_end);
                    ref_site+=cigar[i].l;
                    exon_start = ref_site;
                    break;
                case 'U':
                    t_total++;
                    t_cigar[t_total].c='U';
                    t_cigar[t_total].l=cigar[i].l;

                    ref_site+=cigar[i].l;
                    break;
                case 'S':
                    t_total++;
                    t_cigar[t_total].c='S';
                    t_cigar[t_total].l=cigar[i].l;

                    read->porder = 1;
                    read_site+=cigar[i].l;
                    break;
                }
        }
        t_total++;
            exon_end = ref_site-1;
            insert_exon(exon,snp,exon_start,exon_end);

            read->cigar[0] = '\0';
            for(i = 0;i<t_total;i++)
            {
                sprintf(temp,"%d%c",t_cigar[i].l,t_cigar[i].c);
                strcat(read->cigar,temp);
            }
        }
    else read[i].porder = 1;
}
void *stat_core(void* arg)
{
    struct job_seed *job = (struct job_seed *)arg;

    struct snp_list_t *snp = (struct snp_list_t *)calloc(1,sizeof(struct snp_list_t));
    snp->snp = (struct snp_t *)calloc(MAX_SNP_NUM,sizeof(struct snp_t));
    snp->total = 0;

    struct exon_array *exon = (struct exon_array *)calloc(1,sizeof(struct exon_array));
    exon->exon = (struct exon_inf_t *)calloc(READ_BUF_LENGTH*3,sizeof(struct exon_inf_t));
    exon->total = 0;

	struct read_map_inf *read = (struct read_map_inf *)calloc(READ_BUF_LENGTH,sizeof(struct read_map_inf));
	int read_num = 0;

	int i = 0;

	while (1)
	{
        //read_file;
        pthread_mutex_lock(&ReadLock);
        read_num = read_file_sam_sta(read,job->in_file_1);
		pthread_mutex_unlock(&ReadLock);

		if(read_num==0) break;

        //seed_align
        for (i = 0;i<read_num;i++)
        {
            find_snp(opt,read+i,exon,snp);
        }

        pthread_mutex_lock(&OutputLock);
        out_put_read2(read,read_num,job->out_file_1,job->unmap_file_1);
        Update_Exon(job->opt,exon);
        pthread_mutex_unlock(&OutputLock);
	}
    Update_Snp(job->opt,snp);

    free(snp->snp);
    free(snp);

    free(exon->exon);
    free(exon);

	free(read);

    return (void*)(0);
}
int stat_read(struct m_opt *opt)
{
    //open_file input output unmap
    FILE *in_file_1;
    FILE *out_file_1;
    FILE *unmap_file_1;
    char out_file_name_1[MAX_NAME_LENGTH+50];
    char unmap_file_name_1[MAX_NAME_LENGTH+50];

    in_file_1 = fopen(opt->input_file_1->file[0].name,"r");
    if (in_file_1 == NULL)
        return 1;

    sprintf(out_file_name_1,"%s/map_out.sam",opt->Output_path);
    out_file_1 = fopen(out_file_name_1,"w");
    if (out_file_1 == NULL)
        return 1;

    sprintf(unmap_file_name_1,"%s/unmap_out.sam",opt->Temp_path);
    unmap_file_1 = fopen(unmap_file_name_1,"w");
    if (unmap_file_1 == NULL)
        return 1;

    opt->snp_fn = 0;

    out_put_SAM_head(opt,out_file_1);

    struct job_seed *job = (struct job_seed *)calloc(1,sizeof(struct job_seed));
    job->opt = opt;
    job->in_file_1 = in_file_1;
    job->out_file_1 = out_file_1;
    job->unmap_file_1 = unmap_file_1;

	pthread_t *pthreads = malloc(sizeof(pthread_t) * opt->thread_num);

	int i = 0;
    for (i = 0; i < opt->thread_num; i++) pthread_create(&pthreads[i], NULL, stat_core, job);
	for (i = 0; i < opt->thread_num; i++) pthread_join(pthreads[i], NULL);

    fclose(in_file_1);
    fclose(out_file_1);
    fclose(unmap_file_1);

    free(pthreads);
    free(job);

    opt->input_file_1->total = 1;
    strcpy(opt->input_file_1->file[0].name,unmap_file_name_1);
    return 0;
}

