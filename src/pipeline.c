

#include <pthread.h>

#include "htslib/thread_pool.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "pipeline.h"
#include "taxonomy.h"
#include "sequence.h"
#include "kmer_hash.h"
#include "pylink.h"

static void *pipe_in2tax(void *arg);
static void *process_taxonomy(void *arg);
static void *pipe_tax2seq(void *arg);
static void *process_sequence(void *arg);



typedef struct {

    hts_tpool *p;
    hts_tpool_process *q1;
    hts_tpool_process *q2;
    PARAM *param;
    uint32_t n_seq;
    char **id_list;
    TAX **blast_table;
    TAX *check_table;
    pthread_mutex_t *lock;

} pipe_opt;



typedef struct {

    pipe_opt *o;
    SEQDATA *seqdata;
    uint32_t id;

} pipe_job;


// called by run_pipeline() below; calls newSEQDATA in sequence.c, also calls process_taxonomy() in pipeline.c 
static void *pipe_in2tax(void *arg) {
//    printf("inside pipe_in2tax pipeline.c\n");

    pipe_opt *o = (pipe_opt *)arg;

    for (uint32_t i = 0; i < o->n_seq; i++) {

        pipe_job *job = malloc(sizeof(pipe_job));
        assert(job != 0);
        job->o = o;
        job->seqdata = newSEQDATA(o->param->n_kmers);
        job->id = i;

        if (hts_tpool_dispatch(o->p, o->q1, process_taxonomy, job) != 0) {
            free(job);
            pthread_exit((void *)1);
        }
    }

    pthread_exit(0);

}


// called by pipe_in2tax (this file: pipline.c) cals update_name() from sequencingdata.c
static void *process_taxonomy(void *arg) {
//     printf("inside process_taxonomy pipeline.c\n");

    pipe_job *job = (pipe_job *)arg;

    update_name(job->seqdata, job->o->id_list[job->id]); //calling update_name from sequencingdata.c
    free(job->o->id_list[job->id]);

    if(job->o->blast_table[job->id] != 0) {

      update_blast(job->seqdata, true); //calling update_blast from sequencingdata.c

      for(uint32_t i = 0; i < numHits(job->o->blast_table[job->id]); i++) {  //calling numhits in taxonomy.c

          uint32_t curr_TaxID = getTaxID(job->o->blast_table[job->id], i); //calling getTaxID in taxonomy.c

          // first check to see if current ID was already identified
          pthread_mutex_lock(job->o->lock);
          if(checkTAX(job->o->check_table, curr_TaxID)) { //calls checkTAX in taxonomy.c
              pthread_mutex_unlock(job->o->lock);
              update_tax(job->seqdata, true);
              break;
          }
          pthread_mutex_unlock(job->o->lock);

          // if current ID is unidentified then parse taxdump
          uint32_t check = curr_TaxID;
          if(check_delnodes(job->o->param->delnodes, check)) continue;
          check_merged(job->o->param->merged, &check);
          check_nodes(job->o->param->nodes, &check, check, job->o->param->rank);
          if(check_names(job->o->param->names, check, job->o->param->classification)) {
              pthread_mutex_lock(job->o->lock);
              insertTAX(job->o->check_table, curr_TaxID);
              pthread_mutex_unlock(job->o->lock);
              update_tax(job->seqdata, true);
          }
      }

      freeTAX(job->o->blast_table[job->id]);

    }

    return (void *)job;

}


// called by run_pipeline() below in this file pipeline.c
static void *pipe_tax2seq(void *arg) {
//    printf("inside pipe_tax2seq pipeline.c\n");

    pipe_opt *o = (pipe_opt *)arg;
    hts_tpool_result *result;

    while ((result = hts_tpool_next_result_wait(o->q1))) {

        pipe_job *job = (pipe_job *)hts_tpool_result_data(result);
        hts_tpool_delete_result(result, 0);

        if (hts_tpool_dispatch(job->o->p, job->o->q2, process_sequence, job) != 0) pthread_exit((void *)1);

        if (hts_tpool_process_empty(o->q1)) {
            free(o->id_list);
            free(o->blast_table);
            freeTAX(o->check_table);
            break;
        }
    }

    pthread_exit(0);

}




static void *process_sequence(void *arg) {
//    printf("inside process_sequence pipeline.c\n");
    pipe_job *job = (pipe_job *)arg;



    //  ================ GC CONTENT ================ //

    faidx_t *fa_idx = fai_load(job->o->param->assembly);
    assert(fa_idx != 0);
    uint32_t seq_len;
    char *seq = fai_fetch(fa_idx, get_name(job->seqdata), &seq_len);

    update_length(job->seqdata, seq_len);

    uint64_t gc_count = 0;
    SEQCODE *seqcode = encode(seq, seq_len, &gc_count);

    update_gc(job->seqdata, gc_count);

    free(seq);
    fai_destroy(fa_idx);



    // ================ KMER DISTRIBUTION ================ //

    HASH *map = newHASH(seq_len, seed0(seqcode), seed1(seqcode));
    for(uint32_t k = 0; k < job->o->param->n_kmers; k++) {
        uint32_t kmer_len = job->o->param->kmer_list[k];
        for(int i = 0; i <= (signed)(seq_len - kmer_len); i++) insertHASH(map, get_kmer(seqcode, kmer_len, i));
        update_kdist(job->seqdata, k, distHASH(map, job->o->param->max_kmer_freq));
        clearHASH(map);
    }

    freeHASH(map);
    freeSEQCODE(seqcode);



    // ================ READ COVERAGE ================ //

    samFile *sam = sam_open(job->o->param->alignment, "r");
    assert(sam != 0);
    hts_idx_t *hts_idx = sam_index_load(sam, job->o->param->alignment);
    assert(hts_idx != 0);
    bam_hdr_t *bam_hdr = sam_hdr_read(sam);
    assert(bam_hdr != 0);
    assert(seq_len == bam_hdr->target_len[job->id]);
    hts_itr_t *itr = sam_itr_querys(hts_idx, bam_hdr, get_name(job->seqdata));
    assert(itr != 0);
    bam1_t *itr_data = bam_init1();
    assert(itr_data != 0);

    uint64_t bp_count = 0;
    while (sam_itr_next(sam, itr, itr_data) >= 0) {
        bam1_core_t *c = &(itr_data->core);
        if((c->pos + c->l_qseq) < seq_len) bp_count += c->l_qseq;
        else bp_count += (seq_len - c->pos);
    }

    update_cov(job->seqdata, bp_count);

    bam_destroy1(itr_data);
    hts_itr_destroy(itr);
    bam_hdr_destroy(bam_hdr);
    hts_idx_destroy(hts_idx);
    sam_close(sam);

    export_seqdata(job->o->param, job->seqdata);

    free(job);

}




int run_pipeline(PARAM *param) {
//    printf("inside run_pipeline in pipeline.c \n");

    // Load quantity and names of sequences
    faidx_t *fa_idx = fai_load(param->assembly);
    assert(fa_idx != 0);
    uint32_t n_seq = faidx_nseq(fa_idx);
    char **id_list = malloc(n_seq * sizeof(char *));
    for(uint32_t i = 0; i < n_seq; i++){
        const char *seq_name = faidx_iseq(fa_idx, i);
        id_list[i] = malloc(strlen(seq_name) + 1);
        strcpy(id_list[i], seq_name);
    }
    fai_destroy(fa_idx);

    // calls functions in taxonomy.c
    // Load blast hits
    TAX **blast_table = parse_blast(param->blast, id_list, n_seq);
    TAX *check_table = newTAX(16);
    pthread_mutex_t lock;
    pthread_mutex_init(&lock, 0);

    // Initialize thread pool
    hts_tpool *p = hts_tpool_init(param->n_threads);
    hts_tpool_process *q1 = hts_tpool_process_init(p, n_seq*2, 0);    // taxonomy queue
    hts_tpool_process *q2 = hts_tpool_process_init(p, n_seq, 1);      // sequence queue
    pipe_opt o = {p, q1, q2, param, n_seq, id_list, blast_table, check_table, &lock};


    // calls functions pipe_in2tax and pipe_tax2seq within pipeline.c (this file)
    // Launch our data source and sink threads
    pthread_t tidIto1, tid1to2;
    pthread_create(&tidIto1, 0, pipe_in2tax, &o);
    pthread_create(&tid1to2, 0, pipe_tax2seq, &o);

    // Join pipe threads
    void *retv;
    int ret = 0;
    pthread_join(tidIto1, &retv); ret |= (retv != 0);
    pthread_join(tid1to2, &retv); ret |= (retv != 0);

    // Clean up thread pool allocations
    pthread_mutex_destroy(&lock);
    hts_tpool_process_destroy(q1);
    hts_tpool_process_flush(q2);
    hts_tpool_process_destroy(q2);
    hts_tpool_destroy(p);

    return ret;

}
