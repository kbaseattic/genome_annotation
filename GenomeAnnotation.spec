/*
  API Access to the Genome Annotation Service.

  Provides support for gene calling, functional annotation, re-annotation. Use to extract annotation in
formation about an existing genome, or to create new annotations.

 */
module GenomeAnnotation
{
    typedef int bool;
    typedef string md5;
    typedef list<md5> md5s;
    typedef string genome_id;
    typedef string feature_id;
    typedef string contig_id;
    typedef string feature_type;

    /* A region of DNA is maintained as a tuple of four components:

		the contig
		the beginning position (from 1)
		the strand
		the length

	   We often speak of "a region".  By "location", we mean a sequence
	   of regions from the same genome (perhaps from distinct contigs).

	   Strand is either '+' or '-'.
        */
    typedef tuple<contig_id, int begin, string strand,int length> region_of_dna;

    /*
	a "location" refers to a sequence of regions
    */
    typedef list<region_of_dna> location;

    typedef string analysis_event_id;
    typedef structure {
	analysis_event_id id;
	string tool_name;
	float execution_time;
	list<string> parameters;
	string hostname;
    } analysis_event;

    typedef tuple<string comment, string annotator, int annotation_time, analysis_event_id> annotation;

    typedef structure {
	bool truncated_begin;
	bool truncated_end;
	/* Is this a real feature? */
	float existence_confidence;

	bool frameshifted;
	bool selenoprotein;
	bool pyrrolysylprotein;

	bool overlap_allowed;

	float hit_count;
	float weighted_hit_count;
    } feature_quality_measure;

    /* A feature object represents a feature on the genome. It contains 
       the location on the contig with a type, the translation if it
       represents a protein, associated aliases, etc. It also contains
       information gathered during the annotation process that is involved
       in stages that perform overlap removal, quality testing, etc.
    */
    typedef structure {
	feature_id id;
	location location;
	feature_type type;
	string function;
	string protein_translation;
	list<string> aliases;
	list<annotation> annotations;
	feature_quality_measure quality;
	analysis_event_id feature_creation_event;
    } feature;

    /* Data for DNA contig */
    typedef structure {
	contig_id id;
	string dna;
	int genetic_code;
	string cell_compartment;
	string replicon_type;
	/* circular / linear */
	string replicon_geometry;
	bool complete;
    } contig;

    typedef structure {
	genome_id genome;
	float closeness_measure;
    } close_genome;

    typedef structure
    {
	float frameshift_error_rate;
	float sequence_error_rate;
    } genome_quality_measure;

    /* All of the information about particular genome */
    typedef structure {
	genome_id id;
	string scientific_name;
	string domain;
	int genetic_code;
	string source;
	string source_id;

	genome_quality_measure quality;
	
	list<contig> contigs;
	list<feature> features;

	list<close_genome> close_genomes;

	list <analysis_event> analysis_events;
    } genomeTO;

    typedef string subsystem;
    typedef string variant;
    typedef tuple<subsystem,variant> variant_of_subsystem;
    typedef list<variant_of_subsystem> variant_subsystem_pairs;
    typedef string fid;
    typedef string role;
    typedef string function;
    typedef tuple<fid,role> fid_role_pair;
    typedef list<fid_role_pair> fid_role_pairs;
    typedef tuple<fid,function> fid_function_pair;
    typedef list<fid_function_pair> fid_function_pairs;

    /* Metabolic reconstruction
       represents the set of subsystems that we infer are present in this genome
    */
    typedef structure {
	variant_subsystem_pairs subsystems;
	fid_role_pairs bindings;
	fid_function_pairs assignments;
    } reconstructionTO;

    typedef tuple<fid,md5,location,function> fid_data_tuple;
    typedef list<fid_data_tuple> fid_data_tuples;

    funcdef genomeTO_to_reconstructionTO (genomeTO) returns (reconstructionTO);
    funcdef genomeTO_to_feature_data (genomeTO) returns (fid_data_tuples);
    funcdef reconstructionTO_to_roles (reconstructionTO) returns (list<role>);
    funcdef reconstructionTO_to_subsystems(reconstructionTO) returns (variant_subsystem_pairs);

    /*
     * Given a genome object populated with contig data, perform gene calling
     * and functional annotation and return the annotated genome.
     *
     *  NOTE: Many of these "transformations" modify the input hash and
     *        copy the pointer.  Be warned.
     */
    funcdef annotate_genome(genomeTO) returns (genomeTO);
    funcdef call_selenoproteins(genomeTO) returns (genomeTO);
    funcdef call_pyrrolysoproteins(genomeTO) returns (genomeTO);

    /* [ validate.enum("5S", "SSU", "LSU", "ALL") ] */
    typedef string rna_type;
    
    /*
     * Given a genome typed object, find instances of ribosomal RNAs in
     * the genome.
     *
     * The types parameter is used to select the types of RNAs to
     * call. It is a list of strings where each value is one of
     *
     *    "5S"
     *    "SSU"
     *    "LSU"
     *
     * or "ALL" to choose all available rRNA types.
     */
    funcdef call_features_rRNA_SEED(genomeTO genome_in, list<rna_type> types) returns (genomeTO genome_out);

    /*
     * Given a genome typed object, find instances of tRNAs in
     * the genome.
     */
    funcdef call_features_tRNA_trnascan(genomeTO genome_in) returns (genomeTO genome_out);

    /*
     * Given a genome typed object, find instances of all RNAs we currently
     * have support for detecting.
     */
    funcdef call_RNAs(genomeTO genome_in) returns (genomeTO genome_out);

    typedef structure
    {
	int min_training_len;
    } glimmer3_parameters;
    
    funcdef call_features_CDS_glimmer3(genomeTO, glimmer3_parameters params) returns (genomeTO);
    
    funcdef call_features_CDS_prodigal(genomeTO) returns (genomeTO);
    funcdef call_features_CDS_SEED_projection(genomeTO) returns (genomeTO);
    funcdef call_features_CDS_FragGeneScan(genomeTO) returns (genomeTO);

    typedef structure
    {
	float min_identity;
	int min_length;
    } repeat_region_SEED_parameters;
    funcdef call_features_repeat_region_SEED(genomeTO genome_in, repeat_region_SEED_parameters params) returns (genomeTO genome_out);

    funcdef call_features_prophage_phispy(genomeTO genome_in) returns (genomeTO genome_out);

    funcdef call_features_scan_for_matches(genomeTO genome_in, string pattern, string feature_type) returns (genomeTO genome_out);
    
    typedef structure
    {
	int kmer_size;
	string dataset_name;
	int return_scores_for_all_proteins;
	int score_threshold;
	int hit_threshold;
	int sequential_hit_threshold;
	int detailed;
	int min_hits;
	int min_size;
	int max_gap;
    } kmer_v1_parameters;

    funcdef annotate_proteins_kmer_v1(genomeTO, kmer_v1_parameters params) returns (genomeTO);

    typedef structure {
	int min_hits;
	int max_gap;
    } kmer_v2_parameters;
    
    funcdef annotate_proteins_kmer_v2(genomeTO genome_in, kmer_v2_parameters params) returns (genomeTO genome_out);

    funcdef call_features_ProtoCDS_kmer_v1(genomeTO, kmer_v1_parameters params) returns (genomeTO);		/* RAST-style kmers */
    funcdef call_features_ProtoCDS_kmer_v2(genomeTO genome_in, kmer_v2_parameters params) returns (genomeTO genome_out);		/* Ross's new kmers */

    funcdef annotate_proteins(genomeTO) returns (genomeTO);

    /* Determine close genomes. */
    funcdef estimate_crude_phylogenetic_position_kmer(genomeTO) returns (string position_estimate);

    funcdef find_close_neighbors(genomeTO) returns (genomeTO);

    /*
     * Interface to Strep repeats and "boxes" tools
     */
    funcdef call_features_strep_suis_repeat(genomeTO) returns (genomeTO);
    funcdef call_features_strep_pneumo_repeat(genomeTO) returns (genomeTO);
    funcdef call_features_crispr(genomeTO) returns (genomeTO);

    /*
     * Export genome typed object to one of the supported output formats:
     * genbank, embl, or gff.
     * If feature_types is a non-empty list, limit the output to the given
     * feature types.
     */
    funcdef export_genome(genomeTO genome_in, string format, list<string> feature_types) returns (string exported_data);
};
