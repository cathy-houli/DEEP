#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAX_NAME_LENGTH 400
#define READ_MAX_LENGTH 400
#define MAX_STRING_LENGTH 40000

#define MAX_SPLICE 1000000

struct exon_t
{
    char chr[40];
    unsigned int start;
    unsigned int end;
};
int main(int argc, char *argv[])
{
    FILE *file;
    FILE *write;

    file = fopen(argv[1],"r");
    write = fopen(argv[2],"w");

    char f_line[MAX_STRING_LENGTH];
    char name[400];
    char temp_name[400];
    char chr[40];
    char exonf[400];
    unsigned int start,end,length;

    struct exon_t exon[500];
    int exon_num = 0;

    int i;
    int total = 0;

    name[0] = '\0';
    length = 0;

    while(fgets(f_line,MAX_STRING_LENGTH,file)!=NULL)
    {
        sscanf(f_line,"%s\t%*s\t%s\t%u\t%u\t%*s\t%*c\t%*c\t%*s\t%s",chr,exonf,&start,&end,temp_name);
        if(strcmp(exonf,"exon")!=0) continue;

        if(strcmp(temp_name,name)!=0)
        {
            //if((length>100)&&(total<9998))
            {
		total++;
                for(i = 1;i<exon_num;i++)
                {
                    fprintf(write,"%s\t%u\t%u\n",chr,exon[i-1].end+1,exon[i].start-1);
                }
            }
            strcpy(name,temp_name);
            length = 0;
            exon_num = 0;
        }
        exon[exon_num].start = start;
        exon[exon_num].end = end;
        strcpy(exon[exon_num].chr,chr);
        exon_num++;

        length += end-start+1;
    }

    fclose(file);
    fclose(write);

    return 0;
}
