#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAX_NAME_LENGTH 400
#define READ_MAX_LENGTH 400
#define STRING_MAX_LENGTH 4000

#define MAX_SPLICE 2000000
struct splice_t
{
    char chr[50];
    unsigned int start;
    unsigned int end;
    int flag;
};
int main(int argc, char *argv[])
{
    FILE *file_1;
    FILE *file_2;

    file_1 = fopen(argv[1],"r");
    file_2 = fopen(argv[2],"r");

    //FILE *write_1;
    //FILE *write_2;
    //write_1 = fopen(argv[3],"w");
    //write_2 = fopen(argv[4],"w");

    int answer_total = 0;
    int answer_found = 0;
    int answer_unfound = 0;
    float found_per = 0;
    float unfound_per = 0;

    int test_right = 0;
    int test_wrong = 0;
    float right_per = 0;
    float wrong_per = 0;

    char f_line[STRING_MAX_LENGTH];

    struct splice_t *test = (struct splice_t *)calloc(MAX_SPLICE,sizeof(struct splice_t));
    struct splice_t *anwser = (struct splice_t *)calloc(1,sizeof(struct splice_t));

    int th = 5;

    int back_line = -1;
    int test_line = 0;

    int total_test = 0;

    while(fgets(f_line,STRING_MAX_LENGTH,file_2)!=NULL)
    {
        if(total_test<MAX_SPLICE)
            sscanf(f_line,"%s %u %u",test[total_test].chr,&(test[total_test].start),&(test[total_test].end));
        else return 0;
        test[total_test].flag = 0;
        total_test++;
    }
    while (fgets(f_line,STRING_MAX_LENGTH,file_1)!=NULL)//循环处理read,直到所有的文件都读完了
    {
        sscanf(f_line,"%s %u %u",anwser->chr,&(anwser->start),&(anwser->end));
        anwser->flag = 0;
        answer_total++;

        while((strcmp(anwser->chr,test[test_line].chr)>0)&&(test_line<total_test))
        {
            if(test_line<total_test) test_line++;
        }
        while((strcmp(anwser->chr,test[test_line].chr)==0)&&(anwser->start>test[test_line].start+th)&&(test_line<total_test))
        {
            if(test_line<total_test) test_line++;
        }
        back_line = test_line;
        while((strcmp(anwser->chr,test[test_line].chr)==0)&&(anwser->start>=test[test_line].start-th)&&(anwser->start<=test[test_line].start+th)&&(test_line<total_test))
        {
            if((anwser->end<test[test_line].end-th)&&(test_line<total_test))
            {
                if(test_line<total_test) test_line++;
                continue;
            }
            if((anwser->end>=test[test_line].end-th)&&(anwser->end<=test[test_line].end+th)&&(test_line<total_test))
            {
                anwser->flag = 1;
                test[test_line].flag = 1;
                if(test_line<total_test) test_line++;
            }
            else
            {
                if(test_line<total_test) test_line++;
            }
        }
        test_line = back_line;
        if(anwser->flag == 0 )
        {
            answer_unfound++;
            //fprintf(write_1,"%s",f_line);
        }
        else answer_found++;
    }
    int i = 0;
    for(i = 0;i< total_test;i++)
    {
        if(test[i].flag == 0 )
        {
            test_wrong++;
            //fprintf(write_2,"%s\t%u\t%u\n",test[i].chr,test[i].start,test[i].end);
        }
        else test_right++;
    }

    printf("all splice in gtf file:%d\n",answer_total);
    printf("found/unfound splice in gtf file:%d/%d\n",answer_found,answer_unfound);
    found_per = (float)answer_found/answer_total;
    unfound_per = (float)answer_unfound/answer_total;
    printf("percentage of found/unfound splice in gtf file:%.2f/%.2f\n\n",found_per,unfound_per);

    printf("all splice in input sam file:%d\n",total_test);
    printf("right/wrong splice in input sam file:%d/%d\n",test_right,test_wrong);
    right_per = (float)test_right/total_test;
    wrong_per = (float)test_wrong/total_test;
    printf("percentage of found/unfound splice in gtf file:%.2f/%.2f\n\n",right_per,wrong_per);

    free(anwser);
    free(test);

    fclose(file_1);
    fclose(file_2);

    //fclose(write_1);
    //fclose(write_2);

    return 0;
}
