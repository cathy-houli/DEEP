#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <zlib.h>
//#include <io.h>
#include <ctype.h>//大小写转换 toupper

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

struct read_t
{
    char name[100];
    char cigar[1000];
    char chr[40];
    unsigned int pos;
    unsigned int site[1000];
    int flag;
    int ab;
    int full;
    char part_flag[1000];
    int part_score;
    int length;
};
struct cigar_t
{
    char c;
    int l;
};

void cigar2site(struct read_t *read)
{
    unsigned int ref_site = read->pos;
    int read_site = 0;

    unsigned int i = 0,j = 0;
    struct cigar_t cigar[50];
    int cigar_total = 0;
    char temp[20];

    if((read->pos!=0)&&(read->cigar[0]!='*'))
    {
        j = 0;
        cigar_total = 0;
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
                cigar[cigar_total].c = read->cigar[i];
                cigar[cigar_total].l = atoi(temp);
                cigar_total++;
                j = 0;
            }
        }
        for(i = 0;i<cigar_total;i++)
        {
            switch(cigar[i].c)
            {
            case 'M':case 'X':case 'H':
                for(j = 0;j<cigar[i].l;j++)
                {
                    read->site[j+read_site] = ref_site+j;
                }
                ref_site+=cigar[i].l;
                read_site+=cigar[i].l;
                break;
            case 'I':
                for(j = 0;j<cigar[i].l;j++)
                {
                    read->site[j+read_site] = ref_site-1;
                }
                read_site+=cigar[i].l;
                break;
            case 'D':case 'N':case 'U':
                ref_site+=cigar[i].l;
                break;
            case 'S':
                for(j = 0;j<cigar[i].l;j++)
                {
                    read->site[j+read_site] = 0;
                }
                read_site+=cigar[i].l;
                break;
            }
        }
    }

}
int read_sam(struct read_t *read,FILE *in)
{
    char f_line[10000];
    char seq[10000];
    int flag = 0;
    while((fgets(f_line,10000,in)!=NULL)&&(f_line[0]!='\n'))
    {
        while(f_line[0]=='@')
        {
            if((fgets(f_line,10000,in)==NULL))return 1;
        }
        if(f_line[0]!='@') {flag = 1;break;}
    }
    if(flag)
    {
        sscanf(f_line,"%s\t%d\t%s\t%u\t%*d\t%s\t%*s\t%*u\t%*d\t%s",read->name,&(read->flag),read->chr,&(read->pos),read->cigar,seq);
        read->name[strlen(read->name)-1] = '\0';
        read->full = 0;
        read->length = strlen(seq);

        if(((read->flag)&0x0040)==0x0040) read->ab=1;
        else read->ab=2;

        return 1;
    }
    return 0;
}
int read_rum(struct read_t *read,FILE *in)
{
    int i;
    char f_line[10000];
    char seq[1000];
    char *p;

    if((fgets(f_line,10000,in)!=NULL)&&(f_line[0]!='\n'))
    {
        sscanf(f_line,"%s\t%s\t%u\t%s\t%*s\t%*c\t%s",read->name,read->chr,&(read->pos),read->cigar,seq);
        p = strtok(f_line,"\t");
        i = 0;
        while(p!=NULL)
        {
            if(i==7) {strcpy(seq,p);break;}
            p=strtok(NULL,"\t");
            i++;
        }
        read->name[strlen(read->name)-1] = '\0';
        read->flag = 0;
        read->full = 0;
        read->part_score = 0;
        read->length = strlen(seq);//strlen(seq);
        if(seq[read->length-1] == '\n'){read->length--;}

        memset(read->part_flag,0,read->length);
        return 1;
    }
    return 0;
}
int read_hum(struct read_t *read,FILE *in)
{
    int i;
    char f_line[10000];
    char seq[1000];
    char *p;

    if((fgets(f_line,10000,in)!=NULL)&&(f_line[0]!='\n'))
    {
        sscanf(f_line,"%s\t%s\t%u\t%*u\t%s\t%*s\t%*c\t%s",read->name,read->chr,&(read->pos),read->cigar,seq);
        p = strtok(f_line,"\t");
        i = 0;
        while(p!=NULL)
        {
            if(i==7) {strcpy(seq,p);break;}
            p=strtok(NULL,"\t");
            i++;
        }
        read->name[strlen(read->name)-1] = '\0';
        read->flag = 0;
        read->full = 0;
        read->part_score = 0;
        read->length = strlen(seq);//strlen(seq);
        if(seq[read->length-1] == '\n'){read->length--;}

        memset(read->part_flag,0,read->length);
        return 1;
    }
    return 0;
}
int main(int argc, char *argv[])
{
    uint64_t A_base = 0,T_base = 0,RA_base = 0,RT_base = 0,U_base = 0;

    FILE *T,*A;
    A = fopen(argv[1],"r");
    T = fopen(argv[2],"r");

    char name[1000];

    FILE *Right,*A_Part,*T_Part,*T_Wrong,*A_Wrong,*Un;

    sprintf(name,"%s/Right",argv[3]); Right = fopen(name,"w");
    sprintf(name,"%s/A_Part",argv[3]); A_Part = fopen(name,"w");
    sprintf(name,"%s/A_Wrong",argv[3]); A_Wrong = fopen(name,"w");
    sprintf(name,"%s/T_Part",argv[3]); T_Part = fopen(name,"w");
    sprintf(name,"%s/T_Wrong",argv[3]); T_Wrong = fopen(name,"w");
    sprintf(name,"%s/Un",argv[3]); Un = fopen(name,"w");

    struct read_t answer1,answer2;
    struct read_t test;

    read_sam(&test,T);
    cigar2site(&test);
    T_base+=test.length;

    read_hum(&answer1,A);
    cigar2site(&answer1);
    A_base+=answer1.length;

    read_hum(&answer2,A);
    cigar2site(&answer2);
    A_base+=answer2.length;


    int flag = 0;
    int score = 0;
    char part_flag[1000];
    int i;

    int x;
    int num = 0;

    while(1)
    {
        if((strcmp(test.name,answer1.name)==0))//&&(test.cigar[0]!='*'))
        {
            if (test.cigar[0]=='*')
            num++;
            else if(test.ab==1)
            {
                answer1.flag = 1;
                if((strcmp(test.chr,answer1.chr)==0)&&(test.pos==answer1.pos)&&(strcmp(test.cigar,answer1.cigar)==0))
                    {RT_base+=test.length;answer1.full = 1;}
                else
                {
                    flag = 0;
                    score = 0;
                    memset(part_flag,0,answer1.length);
                    if(test.site[i]==0) part_flag[i] = 0;
                    else if(test.length==answer1.length)
                    {
                        for(i = 0;i<answer1.length;i++)
                        {
                            if(test.site[i]==answer1.site[i])
                            {
                                flag = 1;
                                RT_base++;
                                part_flag[i] = 1;
                                score++;
                            }
                            else part_flag[i] = 2;
                        }
                    }
                    if(score>answer1.part_score) {answer1.part_score = score;memcpy(answer1.part_flag,part_flag,sizeof(char)*answer1.length);}

                    if(flag) fprintf(T_Part,"%sa\t%s\t%u\t%s\n",test.name,test.chr,test.pos,test.cigar);
                    else fprintf(T_Wrong,"%sa\t%s\t%u\t%s\n",test.name,test.chr,test.pos,test.cigar);
                }
            }
            else if(test.ab ==2)
            {
                answer2.flag = 1;
                if((strcmp(test.chr,answer2.chr)==0)&&(test.pos==answer2.pos)&&(strcmp(test.cigar,answer2.cigar)==0))
                    {RT_base+=test.length;answer2.full = 1;}
                else
                {
                    flag = 0;
                    score = 0;
                    memset(part_flag,0,answer2.length);
                    if(test.length==answer2.length)
                    {
                        for(i = 0;i<answer2.length;i++)
                        {
                            if(test.site[i]==0) part_flag[i] = 0;
                            else if(test.site[i]==answer2.site[i])
                            {
                                flag = 1;
                                RT_base++;
                                part_flag[i] = 1;
                                score++;
                            }
                            else part_flag[i] = 2;
                        }
                    }
                    if(score>answer2.part_score) {answer2.part_score = score;memcpy(answer2.part_flag,part_flag,sizeof(char)*answer2.length);}

                    if(flag) fprintf(T_Part,"%sb\t%s\t%u\t%s\n",test.name,test.chr,test.pos,test.cigar);
                    else fprintf(T_Wrong,"%sb\t%s\t%u\t%s\n",test.name,test.chr,test.pos,test.cigar);
                }
            }
            if(read_sam(&test,T)==0) break;
            cigar2site(&test);
            T_base+=test.length;
        }
        else
        {
            if(answer1.flag == 0)
                {fprintf(Un,"%sa\t%s\t%u\t%s\n",answer1.name,answer1.chr,answer1.pos,answer1.cigar);U_base+=answer1.length;}
            else if(answer1.full == 1)
                {fprintf(Right,"%sa\t%s\t%u\t%s\n",answer1.name,answer1.chr,answer1.pos,answer1.cigar);RA_base+=answer1.length;}
            else
            {
                if(answer1.part_score!=0)
                    fprintf(A_Part,"%sa\t%s\t%u\t%s\n",answer1.name,answer1.chr,answer1.pos,answer1.cigar);
                else
                    fprintf(A_Wrong,"%sa\t%s\t%u\t%s\n",answer1.name,answer1.chr,answer1.pos,answer1.cigar);

                for(x = 0;x<answer1.length;x++)
                {
                    if(answer1.part_flag[x]==1) RA_base++;
                    else if (answer1.part_flag[x]==0) U_base++;
                }
            }

            if(answer2.flag == 0)
                {fprintf(Un,"%sb\t%s\t%u\t%s\n",answer2.name,answer2.chr,answer2.pos,answer2.cigar);U_base+=answer2.length;}
            else if(answer2.full == 1)
                {fprintf(Right,"%sb\t%s\t%u\t%s\n",answer2.name,answer2.chr,answer2.pos,answer2.cigar);RA_base+=answer2.length;}
            else
            {
                if(answer2.part_score!=0)
                    fprintf(A_Part,"%sb\t%s\t%u\t%s\n",answer2.name,answer2.chr,answer2.pos,answer2.cigar);
                else
                    fprintf(A_Wrong,"%sb\t%s\t%u\t%s\n",answer2.name,answer2.chr,answer2.pos,answer2.cigar);

                for(x = 0;x<answer2.length;x++)
                {
                    if(answer2.part_flag[x]==1) RA_base++;
                    else if (answer2.part_flag[x]==0) U_base++;
                }
            }

            if(read_hum(&answer1,A)==0)break;
            cigar2site(&answer1);
            A_base+=answer1.length;

            if(read_hum(&answer2,A)==0)break;
            cigar2site(&answer2);
            A_base+=answer2.length;
            continue;
        }
    }
        {
            if(answer1.flag == 0)
                {fprintf(Un,"%sa\t%s\t%u\t%s\n",answer1.name,answer1.chr,answer1.pos,answer1.cigar);U_base+=answer1.length;}
            else if(answer1.full == 1)
                {fprintf(Right,"%sa\t%s\t%u\t%s\n",answer1.name,answer1.chr,answer1.pos,answer1.cigar);RA_base+=answer1.length;}
            else
            {
                if(answer1.part_score!=0)
                    fprintf(A_Part,"%sa\t%s\t%u\t%s\n",answer1.name,answer1.chr,answer1.pos,answer1.cigar);
                else
                    fprintf(A_Wrong,"%sa\t%s\t%u\t%s\n",answer1.name,answer1.chr,answer1.pos,answer1.cigar);

                for(x = 0;x<answer1.length;x++)
                {
                    if(answer1.part_flag[x]==1) RA_base++;
                    else if (answer1.part_flag[x]==0) U_base++;
                }
            }

            if(answer2.flag == 0)
                {fprintf(Un,"%sb\t%s\t%u\t%s\n",answer2.name,answer2.chr,answer2.pos,answer2.cigar);U_base+=answer2.length;}
            else if(answer2.full == 1)
                {fprintf(Right,"%sb\t%s\t%u\t%s\n",answer2.name,answer2.chr,answer2.pos,answer2.cigar);RA_base+=answer2.length;}
            else
            {
                if(answer2.part_score!=0)
                    fprintf(A_Part,"%sb\t%s\t%u\t%s\n",answer2.name,answer2.chr,answer2.pos,answer2.cigar);
                else
                    fprintf(A_Wrong,"%sb\t%s\t%u\t%s\n",answer2.name,answer2.chr,answer2.pos,answer2.cigar);

                for(x = 0;x<answer2.length;x++)
                {
                    if(answer2.part_flag[x]==1) RA_base++;
                    else if (answer2.part_flag[x]==0) U_base++;
                }
            }
        }


    while(read_hum(&answer1,A)!=0)
    {
        A_base+=answer1.length;
        U_base+=answer1.length;
        {fprintf(Un,"%sa\t%s\t%u\t%s\n",answer1.name,answer1.chr,answer1.pos,answer1.cigar);U_base+=answer1.length;}

        if(read_hum(&answer2,A)==0)break;
        A_base+=answer2.length;
        U_base+=answer2.length;
        {fprintf(Un,"%sb\t%s\t%u\t%s\n",answer2.name,answer2.chr,answer2.pos,answer2.cigar);U_base+=answer2.length;}
    }
    fprintf(stdout,"Answer Total base: %lu\n",A_base);
    fprintf(stdout,"Test Total base: %lu\n",T_base);
    fprintf(stdout,"Answer Right base: %lu\n",RA_base);
    fprintf(stdout,"Test Right base: %lu\n",RT_base);
    fprintf(stdout,"Answer Wrong base: %lu\n",A_base-RA_base-U_base);
    //fprintf(stdout,"Test Wrong base: %lu\n",T_base-RT_base-TU_base);
    fprintf(stdout,"S base: %lu\n",U_base);

    fclose(T);
    fclose(A);

    fclose(A_Part);
    fclose(T_Part);
    fclose(T_Wrong);
    fclose(A_Wrong);
    fclose(Un);
    fclose(Right);

    return 0;
}

