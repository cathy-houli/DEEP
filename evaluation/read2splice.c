#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAX_STRING_LENGTH 400000

struct cigar_t
{
    int l;
    char c;
};
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
int main(int argc, char *argv[])
{
    FILE *file;
    FILE *write;

    file = fopen(argv[1],"r");
    write = fopen(argv[2],"w");

    char f_line[MAX_STRING_LENGTH];
    char cigar[MAX_STRING_LENGTH];
    char chr[40];
    unsigned int ref_pos,end_pos;
    int read_pos;

    struct cigar_t t_cigar[5000];
    int t_num = 0;
    int i;

    while(fgets(f_line,MAX_STRING_LENGTH,file)!=NULL)
    {
        if(f_line[0]=='@') continue;
        if(argv[3][0]=='r') sscanf(f_line,"%*s\t%*d\t%s\t%u\t%*d\t%s\t%*s",chr,&ref_pos,cigar);
	else if(argv[3][0]=='a') sscanf(f_line,"%*s\t%s\t%u\t%*u\t%s\t%*s",chr,&ref_pos,cigar);

        process_cigar(cigar,t_cigar,&t_num,&end_pos);

        for(i = 0;i<t_num;i++)
        {
            switch(t_cigar[i].c)
            {
                case 'M':
                case 'X': read_pos+=t_cigar[i].l; ref_pos+=t_cigar[i].l;break;
                case 'I':
                case 'S':
                case 'H': read_pos+=t_cigar[i].l;break;
                case 'D': ref_pos+=t_cigar[i].l;break;
                case 'N':
                    fprintf(write,"%s\t%u\t%u\n",chr,ref_pos,ref_pos+t_cigar[i].l-1);
                    ref_pos+=t_cigar[i].l;break;
                default:break;
            }
        }
    }

    fclose(file);
    fclose(write);

    return 0;
}
