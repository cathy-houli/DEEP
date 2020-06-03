#include "main.h"

int usage()
{
    fprintf(stdout, "\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "Program: use RNA-seq to find junctions, use align file to find more junctions\n");
    fprintf(stdout, "Version: %s\n", PACKAGE_VERSION);
    fprintf(stdout, "Contact: Hou.L <ye.wenzhu@gmail.com>\n\n");

    fprintf(stdout, "Usage:   DEEP <command>\n\n");

    fprintf(stdout, "command :\n");
    fprintf(stdout, "   index:         build reference hash index\n");
    fprintf(stdout, "   micro:         run micro mode\n");
    fprintf(stdout, "   complete:      run complete mode\n");
    fprintf(stdout, "   version:       show current version\n\n");

    fprintf(stdout, "index:  MAP index input_path output_path\n\n");

    fprintf(stdout, "   input_path:    reference files path, each need file's name must be *.fa\n");
    fprintf(stdout, "   output_path:   need enough space\n\n");


    fprintf(stdout, "complete option :\n");
    fprintf(stdout, "<basic> :\n");
    fprintf(stdout, "   -B <STRING>:   reference BWT index path\n");
    fprintf(stdout, "   -H <STRING>:   hash files path, each need file's name must be *.hash with *.ann files\n");
    fprintf(stdout, "   -f         :   fa formate\n");
    fprintf(stdout, "   -q         :   fq formate\n");
    fprintf(stdout, "   -1 <STRING>:   input files, need *.fa,*.fq files,split with ','\n");
    fprintf(stdout, "   -2 <STRING>:   input pair-end files, need *.fa,*.fq files,split with ','\n");
    fprintf(stdout, "   -O <STRING>:   output path,need enough space\n\n");

    fprintf(stdout, "<alignment parameter> :\n");
    fprintf(stdout, "   -t <INT>:      output read score threshold(90)\n");
    fprintf(stdout, "   -r <INT>:      search area(500000)\n");
    fprintf(stdout, "   -a :           print all read sites while coverage score and alignment score over set threld will be used(-b)\n");
    fprintf(stdout, "   -b :           print best read sites while coverage score and alignment score over set threld will be used(-b)\n");
    fprintf(stdout, "   -p <INT>:      thread(1)\n");
    fprintf(stdout, "   -m <INT>:      match score(1)\n");
    fprintf(stdout, "   -s <INT>:      miss score(-1)\n");
    fprintf(stdout, "   -g <INT>:      gap score(-3)\n");
    fprintf(stdout, "   -e <INT>:      min exon length(12)\n");
    fprintf(stdout, "   -d             deep mode(0)\n\n");

    return 1;
}
char base2char[5] = {'A','C','G','T'};
unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
int check_inputfile(struct file_list *input_list,char *input_file_name)
{
    char *p;
    int i = 0;
    p=strtok(input_file_name,",");
    while(p)
    {
        if((access(p,0))== -1)
        {
            fprintf(stdout,"[DEEP]error!!! cannot find file %s...\n",p);
            return 1;
        }

        if(i<MAX_FILE_NUM)
        {
            strncpy(input_list->file[i].name,p,MAX_NAME_LENGTH);
            i++;
        }
        else
        {
            fprintf(stdout,"[DEEP]error!!! too many input files...\n");
            return 1;
        }
        p = strtok(NULL,",");
    }
    input_list->total = i;
    return 0;
}
int check_outpath(struct m_opt *opt)
{
    char system_order[MAX_STRING_LENGTH];
    if((access(opt->Output_path,0))== -1)
    {
        fprintf(stdout,"[DEEP]cannot access %s, create output files folder %s...\n",opt->Output_path,opt->Output_path);
        sprintf(system_order,"mkdir %s",opt->Output_path);
        system(system_order);

        if((access(opt->Output_path,0))== -1)
        {
            fprintf(stdout,"[DEEP]error!!! cannot create output files folder...\n");
            return 1;
        }
    }
    //create temp folds
    sprintf(opt->Temp_path,"%s/temp",opt->Output_path);
    if((access(opt->Temp_path,0))== -1)
    {
        fprintf(stdout,"[DEEP]cannot access %s, create output files folder %s...\n",opt->Temp_path,opt->Temp_path);
        sprintf(system_order,"mkdir %s",opt->Temp_path);
        system(system_order);
        if((access(opt->Temp_path,0))== -1)
        {
            fprintf(stdout,"[DEEP]error!!! cannot create output files folder...\n");
            return 1;
        }
    }
    return 0;
}
void opt_init(struct m_opt **opt)
{
    (*opt)->idx = (bwaidx_t*)calloc(1, sizeof(bwaidx_t));

    (*opt)->input_file_1 = (struct file_list *)calloc(1,sizeof(struct file_list));
	(*opt)->input_file_1->file = (struct file_name *)calloc(MAX_FILE_NUM,sizeof(struct file_name));

	(*opt)->input_file_2 = (struct file_list *)calloc(1,sizeof(struct file_list));
	(*opt)->input_file_2->file = (struct file_name *)calloc(MAX_FILE_NUM,sizeof(struct file_name));
	(*opt)->pair = 0;
	(*opt)->file_flag = 0;

	(*opt)->chr = (struct chr_list *)calloc(1,sizeof(struct chr_list));
	(*opt)->chr->list = (struct chr_t *)calloc(MAX_FILE_NUM,sizeof(struct chr_t));

	(*opt)->exon = NULL;
	(*opt)->snp = NULL;
	(*opt)->snp_fn = 0;
	(*opt)->snp_num = 0;
	(*opt)->find_fn = 0;

	(*opt)->hash_front = 10;
	(*opt)->hash_back = 10;

	(*opt)->e_hash = (struct exon_hash **)calloc(pow(4,(*opt)->hash_front),sizeof(struct exon_hash *));
	(*opt)->e_num= (unsigned int *)calloc(pow(4,(*opt)->hash_front),sizeof(unsigned int));

	unsigned int i=0;
    for (i=0; i<pow(4,(*opt)->hash_front); i++)
    {
        (*opt)->e_hash[i] = NULL;
        (*opt)->e_num[i] = 0;
    }

	(*opt)->read_length = 0;
	(*opt)->change_length = 6;

	(*opt)->output_mode = OUTPUT_ALL;
	(*opt)->output_mode = FQ_FILE;
	(*opt)->area = 300000;
	(*opt)->min_exon = 12;
	(*opt)->match = 1;
	(*opt)->miss = 1;
	(*opt)->gap = 3;
	(*opt)->splice = 5;

	(*opt)->score_t= 92;

	(*opt)->deep_mode = 0;
	(*opt)->hard = 0;

	(*opt)->thread_num = 1;

	(*opt)->step_flag = SEED_STEP;

	(*opt)->total_read=0;
	(*opt)->unmapped_read=0;
}
void opt_free(struct m_opt **opt)
{
    free((*opt)->idx);

    free((*opt)->input_file_1->file);
    free((*opt)->input_file_1);
    free((*opt)->input_file_2->file);
    free((*opt)->input_file_2);

    int i = 0;
    for(i = 0;i<(*opt)->chr->total;i++)
        free((*opt)->chr->list[i].seq);
    free((*opt)->chr->list);
	free((*opt)->chr);
}
int check_hash_index(struct m_opt *opt)
{
    int i = 0;
    char name[MAX_STRING_LENGTH];
    for (i = 0;i<opt->chr->total;i++)
    {
        sprintf(name,"%s/%s.hash",opt->Hash_path,opt->chr->list[i].name);
        if((access(name,0))== -1)
        {
            fprintf(stdout,"[DEEP]error!!! cannot find %s...\n",name);
            return 1;
        }
    }
    return 0;
}
int main(int argc, char *argv[])
{
    time_t present_time;
    struct tm *present_tm;

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"[DEEP]begin micro program...%s\n",asctime(present_tm));

    if(strcmp(argv[1],"index")==0) return index_main(argc-1, argv+1);

    opt = (struct m_opt *)calloc(1,sizeof(struct m_opt));
    opt_init(&opt);

    char input_file_name_1[MAX_STRING_LENGTH],input_file_name_2[MAX_STRING_LENGTH];
    input_file_name_1[0]='\0';
    input_file_name_2[0]='\0';

    int i = 0;
    for(i = 2;i<argc;i++)
    {
        if(argv[i][0]=='-')
        {
        switch (argv[i][1]) {
            case 'B': strncpy(opt->BWTpath,argv[++i],MAX_NAME_LENGTH);break;
            case 'H': strncpy(opt->Hash_path,argv[++i],MAX_NAME_LENGTH);break;
            case 'O': strncpy(opt->Output_path,argv[++i],MAX_NAME_LENGTH);break;
            case '1': strncpy(input_file_name_1,argv[++i],MAX_NAME_LENGTH);break;
            case '2': strncpy(input_file_name_2,argv[++i],MAX_NAME_LENGTH);break;
            case 'r': opt->area = atoi(argv[++i]); break;
            case 't': opt->score_t = atoi(argv[++i]); break;
            case 's': opt->miss = atoi(argv[++i]); break;
            case 'm': opt->match = atoi(argv[++i]); break;
            case 'g': opt->gap = atoi(argv[++i]); break;
            case 'p': opt->thread_num = atoi(argv[++i]); break;
            case 'a': opt->output_mode = OUTPUT_ALL; break;
            case 'b': opt->output_mode = OUTPUT_BEST; break;
            case 'f': opt->input_mode = FA_FILE; break;
            case 'q': opt->input_mode = FQ_FILE; break;
            case 'e': opt->min_exon = atoi(argv[++i]); break;
            case 'd': opt->deep_mode = 1; break;
            default: return usage();
            }
        }
		else return usage();
    }

    if(input_file_name_1[0]=='\0') return usage();
    else {if(check_inputfile(opt->input_file_1,input_file_name_1)) return 1;}
    if(input_file_name_2[0]!='\0')
    {
        opt->pair = 1;
        if(check_inputfile(opt->input_file_2,input_file_name_2)) return 1;
    }
    else opt->pair = 0;
    if((opt->input_file_2->total!=0)&&(opt->input_file_1->total!=opt->input_file_2->total))
    {
        fprintf(stdout,"[DEEP]wrong input files...\n");
        return usage();
    }
    if(check_outpath(opt)) return 1;
    if(check_hash_index(opt)) return 1;

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"[DEEP]load index...%s\n",asctime(present_tm));
    if(load_index(opt)) return 1;

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"[DEEP]begin align...%s\n",asctime(present_tm));

    opt->dep =(char *)calloc(opt->idx->bns->l_pac,sizeof(char));

    if(strcmp(argv[1],"micro")==0)
    {
        if(stat_read(opt)) return 1 ;
    }
    else if(strcmp(argv[1],"complete")==0)
    {
        if(seed_align(opt)) return 1;
    }

    if (opt->idx->bwt) bwt_destroy(opt->idx->bwt);

    heap = (struct heap_array *)calloc(1,sizeof(struct heap_array));
    heap->heap = (unsigned int *)calloc(MAP_READ_BUF_LENGTH,sizeof(unsigned int));
    if(heap->heap == NULL)
    {
        fprintf(stdout,"[DEEP core]need more space...\n");
        return 1;
    }
    heap->heap_num = (unsigned int *)calloc(MAP_READ_BUF_LENGTH,sizeof(unsigned int));
    if(heap->heap_num == NULL)
    {
        fprintf(stdout,"[DEEP core]need more space...\n");
        return 1;
    }
    heap->order = (unsigned int *)calloc(MAP_READ_BUF_LENGTH,sizeof(unsigned int));
    if(heap->order == NULL)
    {
        fprintf(stdout,"[DEEP core]need more space...\n");
        return 1;
    }
    heap->heap_length = 0;

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"[DEEP]begin build hash...%s\n",asctime(present_tm));

    if(build_hash(opt)) return 1;

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"[DEEP]begin find site...%s\n",asctime(present_tm));

    if(find_site(opt)) return 1;

    if (opt->idx->bns) bns_destroy(opt->idx->bns);
    if (opt->idx->pac) free(opt->idx->pac);

    for (i=0; i<pow(4,opt->hash_front); i++)
    {
        if(opt->e_num[i]!=0) free(opt->e_hash[i]);
    }
    free(opt->e_num);
    free(opt->e_hash);

    if(opt->deep_mode==1)
    {
        if(opt->unmapped_read>(opt->total_read/20)) opt->hard=1;
        time(&present_time);
        present_tm = localtime(&present_time);
        fprintf(stdout,"[DEEP-D]begin map align...%s\n",asctime(present_tm));

        read_buf = (struct read_map_inf *)calloc(MAP_READ_BUF_LENGTH,sizeof(struct read_map_inf));
        if(map_align(opt)) return 1;
        free(read_buf);
    }

    free(heap->heap);
    free(heap->heap_num);
    free(heap->order);
    free(heap);

    opt_free(&opt);
    free(opt);

    time(&present_time);
    present_tm = localtime(&present_time);
    fprintf(stdout,"[DEEP]exit micro program...%s\n",asctime(present_tm));
    return 0;
}
