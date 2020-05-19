struct hapnode{
	int allele;
	void *nextlocus; /* make this a void pointer so it can point to either another hapnode or a lastnode */
	struct hapnode *nextallele;
};
struct lastnode{
	char *hapstring;
	int allele;
	int status;  /* stores e.g. whether the result for this node has been printed out */
	float count[5]; 
    // 7/12 assume we only need GS or EM, not both at once
    // No!  currently the tree scan is called for distributing unamb cts 
    // in any case everything else is free
	// GS: count[0] - stores GS unambig count -- Note this is not used directly (as of 7/2/12)
	// GS: count[1] - stores GS calculated count, calculated one at a time so this is also old count  
        // for true Gibbs sampler this should be the only count used in the calculation
        // -but for stability the moving average count is used
	// GS: count[2] - stores moving average count // is this safe? Should be, note below "count[2] ... utility storage after inference is complete."
	// GS: count[3] - 
	// GS: count[4] - utility for diagnostics TEMP use for old count for GS calc; old count should just be diagnostic
	// Use of counters for EM calculation; used initially in GS calculation
	/* not integers bcs we are fractionally distributing counts*/ /* reduce to two? just distribute unambig to get rid of [0]; do we really need [3] at all?
	 	count[0] = abscount:  distributed count of unambiguous haplotypes.
		count[1] - stores count from last round of calculation; stores final count. 
		count[2] - accumulates new count during this round; utility storage after inference is complete.
		count[3] = product of oldcount for this haplotype and oldcount 
			for complementary haplotype relative to current genotype [THIS APPEARS TO BE NOT NEEDED, CHECK];
			[No, I think this is needed unless we want to calc cross hap product again...]
			used as utility storage in output but count[2] better for this;
    // To transition from GS calc to output: assume all code after inference uses old (EM) count[] assignments
    // 1) move moving avg ct to count[1] 
    */
    float upperCI;
    float lowerCI;
	int hapnumber; // this is reset to give rank: descending order of frequency. Should it be an unchanging ID? Since ranked hap pointer will save the rank info.
	struct lastnode *nextallele;
	// 12/11 setting up for gibbs sampler calc:
	// vbl to record how many genotypes contain this hap
	int n_hap_gtypes; // can save a space by combining these, after malloc'ing this_haps_gtypes below reusing a counter, but dangerous 
	int hap_gtcount;  // this increments as we add gtypes carrying the haps, above carries final count
	struct gtype_allhaps **this_haps_gtypes; // will need to know the number of gtypes containing this hap, and malloc this vbl
	struct gtype_allhaps **this_haps_unseen_gtypes; // will need to know the number of gtypes containing this hap, and malloc this vbl
	// diagnostic and output:  pointer will be malloc'd to array with hap printout 
    float case_count;
    float ctl_count;
    float coalesce_p;
	char *hapgraphic;
};
struct grokSubBlock_data
{
	// data to gsb
	int *is_monosome;
	int exist_monosomes;
	// data back  
	int subblock_met_criteria;  
};
/* program will search first down each locus, then across to the next, and create 
a new series of nodes when it finds no match.  Not necessary if the genotypes are 
perfectly sorted in allele order-- but is that possible anyway? 

The last node must have a count (counts!).  For now add a count[3] vbl to all nodes, 
this is wasteful but simple.  Maybe one of the pointers is void, can point to counts at end?

But note also I could just use the simple struct haplotype, and search through them
to find a match or absence of a match.  But this way I need either a fixed hapset or 
another system of pointers.  Also this way I would have to follow each haplotype
that had a match at the first (or the first n) loci.*/
struct gs_stats{
    float gsample_freq;
    int rtn_0;
    float rand_p;
    float max_p_freq;
    float upper_limit;
    float lower_limit;
    float p_mvg_avg_freq;
    float int_mvg_avg_freq;
    // simple way to compare the new and old integration
    struct gs_stats *old_int_stats; // has to be pointer, otherwise infinite regress (how does compiler know?)
    int ncalls;
};

struct gtype_hap_pair{
	struct lastnode *hap_node[2];
	// ... but possibly faster to have specific pointers to needed "registers"
	/*float *oldhapct[2]
	float *newhapct[2]*/
	// these aren't suited for GS algorithm:
	float crossproduct; 
	float realcount;
	float countsum;
	float countsquaresum;
};
struct gtype_allhaps{ // all the haps pairs for one gtype
	int gtno; 
	struct gtype_hap_pair *gtype_happair;  // malloc this later....
	int max_n_pairs;  // DROP? keep re missings issues
	int n_happairs;  // FOR NO MISSINGS, JUST ASSUME 2^(nhet - 1) happairs
	int times_gt_appears; // how often the genotype appears in the resampling KEEP
	int hetcat;  // hetcat >= 2 means n hets >= 2 or n missing >0  --DROP? --Try to redo calc to treat unabiguous the same as others
};
// doesn't look like all_gtypes_haps is actually used. No, wrong I think, search
struct all_gtypes_haps{
	int status;   // originally, 0 if gtype_allhaps not filled yet; now 12/06 0 if neither malloc'd or filled
	struct gtype_allhaps *gtype_haps; // malloc this 
	struct missings_gt  **gts_with_missings_data;
	int hapdata_n_gts; // to record if n_gts changes (i.e. unseen gts are added)
	int thisgt;  // isn't this provided by calc params?
	// following for sum of hap freqs from boot calc
	int n_boot_haps;
	int *boot_haplist;
	float *hapfreq; // hap freqs for all genotypes
	// or could use the hap_nodes to sum all genotypes
};
//  if we are considering genotypes with missings:  we will malloc an array of pointers to struct missings_gt.;
// there will be a pointer for each subgenotype, with or without missings.
// For genotypes with missings, we will malloc the structure that it points to. Wastes space for the pointers, 
// and also we waste space for the all_gtypes_haps for genotypes that contain missings. Probably no matter.
struct deltas
{
    float GSsumdelta_steps;
    float GS_avg_sumdelta;
    float GS_avg_fract_sumdelta;
    float ratio; // what ratio?
    float sum_int_x_freq;
};
struct missings_gt
{
	int n_matching_gts;
	int *matching_gts;
};
struct output_spec{
	int inferred_from_loci[2];
	int inferred_loci[2];  /* after testing make this a pointer for batch output */
	int flank_inference;
}; 
struct locus_set{
	char locusname[LOCUS_NAMELENGTH];
	int position;
	int used; //added 3/22/05 by AS, flag to check if all snps in info
			  //file are actually in genotype file
};
/* building a struct to hold all the input genotype data! */
struct gtype_data{
	int nloci;
};
struct ind_gtype{ 
	char indiv_name[INDIV_NAMELENGTH];
	char (*gtype)[GT_STRINGLENGTH]; /* malloc'd to n_loci*GT_STRINGLENGTH */
	int gtypenum;
    int disease_status;
	struct ind_gtype *nextind;
} ;
/* struct to hold genotypes for an individual at all loci, points to another copy so we store all 
indiv. gtypes in a chain */
struct allele_namect{
	int count;
	char name[ALLELE_NAME_LENGTH];
};
struct hap_calc_args{ // set in order of the longer list, the one passed to node_scan.
                      // list passed to hap_node_scan is a subset (I think, except gt/gtset params)
                      // which of these have values passed as pointers?
    struct hapnode *node_ptr_start;
    int startlocus;
    int (**gtype_array)[2];
    int *hapvector;
    int LOCUS_ZERO;
    struct gtset_params *set_params;
    struct gtype_params *gt_params; // currently 1/15 contains set_params, to be fixed
    struct hap_calc_params *calc_params;
    int *tasks;
    int gtno;
};
/*  hap_calc_args issues:
 1) when to wrap and unwrap?
 */
struct gtset_params{/* way to pass all params used in hap calcs, replacing global vars */
	int n_loci;
	int n_gts; // for mc method dimensioning for real n_gts + n supp gts; use this var to distinguish?
	int *gt_count; 
	int *n_msng;
	int *hetct;
    int *allele_ct;
    int *hapct;
	struct allele_namect **allelelist;
	int *is_monosome; // "monosome" = male X; (or Y, but why would we ever look at Y?) //only for seen gts
	int n_indivs_seen; // need this because "indiv_ct" (next) is kludged (remove if fixed) TEMP
	int indiv_ct;
	int indiv_ct_nomissing;
	int indiv_ct_output; // count of indivs for output, including non infer gtypes (with too many missings)
	int chrom_ct;
	//int chrom_ct_sqrd;
	int chrom_ct_output;
	int monosome_ct;
	float obs_pred_factor;
	int *ordered_gt; // only for seen gts
	int startlocus; /* startlocus, endlocus (= startlocus + n_loci - 1) are start and end of 
		subgenotype, for inference of subhaplotypes from  subgenotypes. */	
	struct gtype_allhaps gthaps;
};
// 12/31/14: working to simplify calculation up and down trees (and other accesses)
// This struct can contain genotype number, locus number, etc.
// But the haplotype tree itself really determines where we are; the variable hapvector carries this knowledge
struct place_vars{ // variables to keep track of where we are in the calculation
    int gtype_nmbr;
};
struct gtype_params{ // passed down the tree.  (instance of) gtset_params carries all the data, but this has the gtype number and points to gtset_params
	struct gtset_params *setparams; /* break this off inside called function ? */
	int gtype_nmbr;
	/* double gtype_count;  ??? not an int?  not redundant */
};
/* THIS IS BACKWARDS, gtset_params SHOULD CONTAIN gtype_params, NOT VICE VERSA. But probably works... 
More confusing is the fact that gtype_array is independent, not included in gtset_params*/
struct subhap_ent{ /* make this simple by considering one gtype at a time... No! consider all! */
	int *subhap;
	double subhap_freqprod;
	float subhap_count;
};
struct subhap_ent_calc{
	struct subhap_ent *subhapdat; 
	int firstloc;
	int lastloc;
	int subhapmax_n;
	int n_subhaps_seen;
};
struct tree_connections{
	struct lastnode crossnode; // could have a forwardnode, backnode ("crossnode" is a backnode), here as these are passed by a calling task to a called task here. 
};
struct NRfunc_passer
{
	struct hap_calc_params *calcparams;
	struct gtset_params *setparams;
	struct gtype_params *gt_params;
	struct hapnode *firstnode_ptr;
	int (**gtype_array)[2];
	int *hapvector;
};
struct hap_calc_params{
	struct lastnode **temp_haplist;
	struct lastnode **ranked_haplist;
	int tot_n_haps;
	int hap_counter; // replacing tot_n_haps when just used as counter
	int n_haps_considered; //future use when we don't consider all possible haps
	int this_hap;
	int use_gtypehapstruct;
	float *hapfreqs;
	double sum_delta;
	float *gtfreq;
	float this_gtfreq;
	float crossfreq; /* to replace global variable */
	float binom_p;
	double dip_mutual_info;
	double dip_max_info; 
	double mutual_info;
	double max_info; 
	double this_dip_rowmarginal;
	float obs_pred_factor; 
	float max_unambig_hapcount;
	struct hapnode *firstnode_ptr;  /* needed at end of tree for cross hap call */
	// variables for MC calc
	float pred_gtype_freq;
	float logP;
	int calc_subhapent;
	struct subhap_ent_calc *subent_calc;
	float util_ct;
	void *util_ptr; // to attach whatever structure, for special purposes 
	char **hapoutstrings;
    int calc_w;
    double w_sum;
};
struct output_params{
	int baseoutput; // set for any output, so output function is called (not used this way 3-06)
	int startoutput; /* not used yet */
	int hapoutput;
	int renumber_haps_by_freq;
	int indivoutput;
	int sasoutput;
	int ld_output;
	int block_output; // basically, block output shouldn't happen in same run as other output, but need a safety
	int bestscore_output;
	int distant_locus_diseq_calc;
	int distant_loc1;
	int distant_loc2;
	int distant_real1;
	int distant_real2;
	int calc_LD;
	int sum_boot_LD;
	int calc_boot_LD;
	int print_boot_LD;
	int free_LD_matrices;
	int print_reg_LD;
	int batch;
};
struct all_hap_inferences{
	struct all_hapvet_results *all_vetresults; /* so random hce results are available */
	/*unsigned dim; /* not used; we just give inferred and set the 2 dimensions tot_n_loci */
	int **inferred; 
	struct hap_inference_set ***set; /*  malloc this ptr to dim*dim (2D array of pointers!), dim - tot_n_loci; points to nothing (NULL?) unless inferred */
};
struct hap_vet_ptrs{ /* for whatever has to be passed to vet_hapcalc */
	struct all_hap_inferences *all_inf_ptr;
};
struct hap_inference{
	float pf;
	float hce_pct;
	float rand_pct;
	int subhap_start; /* this and _end will be redundant but helpful for orientation/checking */
	int subhap_end;
};
struct hap_inference_set{ /* */
	unsigned dim[2]; /* or are dimensions necessarily equal? No, extent of subhaps containing this subhap can be constrained by start/end of sequence */
	int minlocus;
	int subhap_start;
	int subhap_end;	
	/* subhap_end is offset for rows, minlocus is offset for columns (add offset to index to get (true) locusnumber).*/
	int **inferred;
	struct hap_inference ***inference; /* for each set malloc this ptr to dims, points to nothing (NULL?) unless inferred */
	/*add vbls for best calls here...*/
	float min_inference;
};
/* redundancy now between hap_inference and hce structures, fix? */
struct hce_result{
	int could_calc;
	float hce;
	float boot_hce_pct;
	float equil_hce1; /* equil_hce1 considers only entropy of call of ambiguous gts; */
	float equil_hce2; /* equil_hce2 considers entropy of call of all gts. */
	float pf;
};
struct hap_diagnostics{
	int could_calc; 
	struct hce_result hce_data; 
	struct hce_result *rand_hce_data; /* need malloc for this dimension, fix */
	struct hap_inference temp_best; /* FIX this in when function is done */
	float analhap_pct;
	float hap1freq;	
};
struct all_hapvet_results{
	int dim;
	int **inferred;
	struct hap_diagnostics ***diagnostics;
};
struct gthap_pair{ /* not used, currently  */
	double *hap1; /*...points to oldcount? offset to get newcount? */
	double *hap2;
	double count_product; /* maybe could eliminate this and use count[3] for hap1 (&2?) */
} /* *gthaps[MAX_N_GTS]*/ ;

