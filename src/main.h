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
#include <ctype.h>

#include "kseq.h"

#define PACKAGE_VERSION "2.0"

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })


/**bwa index struct**/
typedef uint64_t bwtint_t;
typedef unsigned char ubyte_t;
#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL   (1LL<<OCC_INTV_SHIFT)
#define OCC_INTV_MASK  (OCC_INTERVAL - 1)

#define bwt_bwt(b, k) ((b)->bwt[((k)>>7<<4) + sizeof(bwtint_t) + (((k)&0x7f)>>4)])
#define bwt_occ_intv(b, k) ((b)->bwt + ((k)>>7<<4))

/* retrieve a character from the $-removed BWT string. Note that
 * bwt_t::bwt is not exactly the BWT string and therefore this macro is
 * called bwt_B0 instead of bwt_B */
#define bwt_B0(b, k) (bwt_bwt(b, k)>>((~(k)&0xf)<<1)&3)

typedef struct {
	bwtint_t primary; // S^{-1}(0), or the primary index of BWT
	bwtint_t L2[5]; // C(), cumulative count
	bwtint_t seq_len; // sequence length
	bwtint_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT
	uint32_t cnt_table[256];
	int sa_intv;
	bwtint_t n_sa;
	bwtint_t *sa;
} bwt_t;

typedef struct {
	int64_t offset;
	int32_t len;
	int32_t n_ambs;
	uint32_t gi;
	char *name, *anno;
} bntann1_t;

typedef struct {
	int64_t offset;
	int32_t len;
	char amb;
} bntamb1_t;

typedef struct {
	int64_t l_pac;
	int32_t n_seqs;
	uint32_t seed;
	bntann1_t *anns; // n_seqs elements
	int32_t n_holes;
	bntamb1_t *ambs; // n_holes elements
	FILE *fp_pac;
} bntseq_t;

typedef struct {
	bwt_t    *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
} bwaidx_t;


#define OUTPUT_ALL        0x01
#define OUTPUT_BEST       0x02

#define FA_FILE 0x01
#define FQ_FILE 0x02

#define SEED_STEP 0x01
#define SNP_STEP 0x02
#define MAP_STEP 0x03

#define KMER_FRONT 8

#define MAX_NAME_LENGTH 400
#define MAX_FILE_NUM 2000
#define MAX_STRING_LENGTH 10000
#define MAX_READ_LENGTH 400

#define MAX_EXON_NUM 1000000
#define MAX_SNP_NUM 2000000

#define MAP_READ_BUF_LENGTH (256*1024)

struct file_name
{
    char name[MAX_NAME_LENGTH];
};
struct file_list
{
    struct file_name *file;
    int total;
};

struct chr_t
{
    char name[MAX_NAME_LENGTH];
    unsigned int length;
    unsigned int start_site;
    char *seq;
    unsigned int **c_hash;
    unsigned int *c_num;
    int thread_num;
};
struct chr_list
{
    struct chr_t *list;
    int total;
};
struct exon_hash
{
    uint32_t back;
    unsigned int start;
    unsigned int end;
    uint8_t c[2];
    uint8_t l[2];
};
struct hash_temp
{
    uint32_t front;
    uint32_t back;
    unsigned int start;
    unsigned int end;
    uint8_t c[2];
    uint8_t l[2];
};
struct heap_array
{
    unsigned int *heap;
    unsigned int *sort_buf;
    unsigned int *heap_num;
    unsigned int *order;
    unsigned int heap_length;
};
struct heap_array *heap;
struct m_opt{
	char BWTpath[MAX_NAME_LENGTH];
	bwaidx_t *idx;
	char Hash_path[MAX_NAME_LENGTH];

	char Output_path[MAX_NAME_LENGTH];
	char Temp_path[MAX_NAME_LENGTH];
	struct file_list *input_file_1,*input_file_2;
	char pair;
	int file_flag;
	struct chr_list *chr;

	struct exon_array *exon;
	struct snp_list_t *snp;
	int snp_fn;
	int snp_num;
	char *dep;
	uint64_t total;

	int find_fn;
	int result_fn;

	int hash_front,hash_back;

	struct exon_hash **e_hash;
	unsigned int *e_num;

	int read_length,change_length;

	int area;
	int min_exon;
	int score_t;
	int match,miss,gap,splice;

	int thread_num;
	int output_mode;
	int input_mode;

	int deep_mode;
	int hard;

	int step_flag;

	int total_read;
	int unmapped_read;
};
struct m_opt *opt;

struct read_map_inf
{
    char name[MAX_NAME_LENGTH];
    unsigned int ref_site;
    unsigned int start_site;
    unsigned int end_site;
    int front;
    int back;
    char seq[MAX_READ_LENGTH];
    char qual[MAX_READ_LENGTH];
    char cigar[MAX_READ_LENGTH];
    unsigned int *site;
    unsigned int *node;
    int sstart,send;

    int contig_order;
    int contig_start;

    char pname[MAX_NAME_LENGTH];
    unsigned int psite;
    int porder;

    char output_flag;

    int length;
    int flag;
    int MapQ;
    char read_order;
    int chr_order;
    int score;
    int dis;
    char strand;
    char TAG[MAX_READ_LENGTH];
};
struct read_map_inf *read_buf;
#define MAP_READ_BUF_LENGTH (256*1024)
struct candiate_t
{
    unsigned int pos;
    unsigned int end;
    int chr_order;
    char cigar[MAX_READ_LENGTH];
    int score;
    int dis;
    char strand;
    char TAG[MAX_READ_LENGTH];
    char out_put_flag;

    int pair;
    int pair_score;
};
#define MAX_SEED_LAST 10
struct seed_t
{
    int start;
    int length;
    uint64_t pos;
    char flag;
    int score;
    uint64_t abs;

    int last[MAX_SEED_LAST];
    int lnum;
    int lorder;
};
#define MAX_CAND_NUM 20
struct read_inf_t
{
    char name[MAX_NAME_LENGTH];
    int length;
    char seq[MAX_READ_LENGTH];
    char qual[MAX_READ_LENGTH];
    char rseq[MAX_READ_LENGTH];
    char rqual[MAX_READ_LENGTH];
    struct candiate_t cand[MAX_CAND_NUM];
    int cand_num;
    int MapQ;
    char out_flag;
    char pair_flag;
};

struct exon_inf_t
{
    unsigned int start;
    unsigned int end;
    uint16_t *dep;
};
struct exon_array
{
    struct exon_inf_t *exon;
    uint64_t total;
};

struct snp_t
{
    unsigned int start;
    unsigned int end;
    int length;
    char type;
    uint16_t seq;
    int num;
    int dep;
};
//struct snp_t *snp_buf;
struct snp_seq_t
{
    unsigned int start;
    unsigned int end;
    int snp_num;
    struct snp_t snp[4];
};
//struct snp_seq_t *snpb_buf;
struct snp_branch_t
{
    struct snp_seq_t *snp_seq;
    uint64_t total;
};
struct snp_list_t
{
    struct snp_t *snp;
    uint64_t total;
};

#define READ_BUF_LENGTH 10000//20000
#define SEED_BUF_LENGTH 40000
#define SEED_CAND_NUM 100


#define MAX_CIGAR_BUF 100
struct seed_a
{
    bwtint_t start,end;
};
struct cigar_t
{
    int l;
    char c;
};
#define AREA 200000
#define SCORE_BALANCE 1


extern char base2char[5];
extern unsigned char nst_nt4_table[256];
extern int usage();
//index.c
extern void SwapHeap(unsigned int a,unsigned int b,unsigned int *heap,unsigned int *heap_num,unsigned int *order);
extern int index_main(int argc, char *argv[]);

//bwa_index.c
extern void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol);
extern bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k);
extern int load_index(struct m_opt *opt);
extern void bwt_destroy(bwt_t *bwt);
extern void bns_destroy(bntseq_t *bns);

#define MAX_K 35
#define CountTrailingZeroes(x, ans) {ans = __builtin_ctzll(x);}
enum CigarFormat
{
    COMPACT_CIGAR_STRING = 0,
    EXPANDED_CIGAR_STRING = 1,
    COMPACT_CIGAR_BINARY = 2,
};
extern int landau2site(struct m_opt *opt,int length_r,int length_t,uint64_t r_start,char *cigarBuf,unsigned int *siteBuf,int mode,int scoreMode);
extern int computeEditDistanceWithCigar(const char* text, int textLen,const char* pattern, int patternLen,int k,char *cigarBuf, int cigarBufLen, int useM,enum CigarFormat format);


//tool.c
extern int generate_MapQ(struct read_inf_t *read);
extern int generate_flag_single(struct m_opt *opt,struct read_inf_t *read,int cand_order);
extern int generate_flag_paired(struct m_opt *opt,struct read_inf_t *read,struct read_inf_t *read_p,int cand_order,int read_order);
extern void write_MapQ(struct read_map_inf *read,int num);
extern void write_TAG(struct read_map_inf *read);
extern void process_cigar(char *cigarBuf,struct cigar_t *cigar,int *cigar_total,unsigned int *end_pos);


//stat.c
extern void insert_seed(struct seed_t *seed,int *seed_num,int max,struct seed_t seed_r);
extern int seed_cmp(const void *a,const void *b);
extern int seed_cmp_r(const void *a,const void *b);
extern int find_s(struct snp_list_t *snp,unsigned int start);
extern int find_exon(struct exon_array *exon,unsigned int start);
extern void check_read_splice(char *read_seq,char *read_rseq,struct cigar_t *cigar,int *cigar_num,uint64_t start_pos,int chr_order,char strand);
extern int middle_align_seed(struct m_opt *opt,int chr,char strand,uint64_t start,uint64_t end,char *seq,int length,struct cigar_t *cigar,int *cigar_num,int max);
extern int tail_align_seed(struct m_opt *opt,int chr,char strand,uint64_t pos,char *seq,int length,struct cigar_t *cigar,int *cigar_num,int max,uint64_t *start_pos,int mode);
extern void find_cand(struct m_opt *opt,struct read_inf_t *read,struct seed_t *seed,int *seed_num,struct exon_array *exon,struct snp_list_t *snp,float score_T);
extern int stat_read(struct m_opt *opt);
extern int seed_align(struct m_opt *opt);

//variation
extern int hash_e_find(struct exon_hash **hash_front,unsigned int *hash_num,uint32_t f_base,uint32_t b_base);
extern int build_hash(struct m_opt *opt);
extern int find_site(struct m_opt *opt);

//map_align.c
extern int map_align(struct m_opt *opt);

