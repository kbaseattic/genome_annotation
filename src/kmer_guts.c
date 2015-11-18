/*

kmer_guts.c can be compiled into either a 5-mer or an 8-mer version.  I have 
labeled the critical changes with 

     ### CHANGE THIS FOR DIFFERING Ks ###

Note: to run use

   build_data_directory KmerData 8

Then use

   kmer_guts -w -D KmerData < input.contigs

to build an memory map -- it can take 10-15 minutes.

Then use

   kmer_guts -D KmerData < input.contigs

to run.  Eventually, the memory-mapped hash table, which is large,
will become resident, and these should go fairly quickly.

   kmer_guts -D KmerData -a  < amino_acid_sequences.fasta > hits

can be used to call translated PEGs.

where KmersData is a directory that must contain

       final.kmers     [a file of [kmer,avg-offset-from-end,function-index,otu-index]]
       function.index  [a file of [index,function] pairs]
       otu.index       [a file of [index,otu] pairs]
----------------

Conceptually, the data associated with each signature Kmer is

        1. the protein kmer
        2. the average offset from the end of the protein
        3. a set of numeric values that include

                function index
                OTU index

The program takes as input a  file that should be thought of as a sequence
of fasta files in which a single line

          >FLUSH

occurs after each infut file.  That is,

          >contig1
	  acgtacgt
	  >FLUSH
	  >contig2
	  acgtacgt

is an input file containing two input fasta files (each containing a single, short
contig).  The '>FLUSH' causes all files to be flushed.

In effect, the program processes a sequence of requests.  The output for each request
is a piece of stadout terminated by a line

          //

The lines in the output for a request are of four forms.  The first type of line begins
the processing of a contig and looks like

         processing contig[length]

Then, there will be a message like

         TRANSLATION contig length-contig frame offset-of-start  

before beginning output for each of the six frames.

Thus, there will be six such lines printed for each input contig.  Following
each of these "TRANSLATION..." lines, there will be a set of lines like

         CALL Start-in-protein-seq end-in-protein-seq number-hits function-index function

Finally, after all six frames have been processed (for DNA - just one sequence for aa input),
you get a line of the form

         OTU-COUNTS cnt1-otu1 cnt2-otu2 ...

This code uses a table indicating which K-mers are signatures. I call this the 
"kmer_bits" table.  You can run this code for K ranging from 5 to 8.

#################
The line

    #define K 8

defines the kmer size.  
The code is intended to work with only 5-mers or 8-mers.  That is:

YOU MUST SET K TO EITHER 5 OR 8 (sorry).
############################################

COMMAND LINE ARGUMENTS:

    -a         means amino acid input sequence (defaults to DNA)

    -d Level   sets debugging level (1 shows hits; after that it gets intense

    -m MinHits minimum number of hits to get CALLed

    -M MinWeighted  Kmers are now weighted; this is the min of the sum of the weights,

    -O         use order contraint

    -g MaxGap  sets maximum allowed gap between HITS

    -D Data    sets the Data directory where the memory map lives

    -s HashSize make sure that the value is the same when you save the memory map
                and when you use it to search

    -w          write the memory map (means Data must contain final.kmers and the indexes

    -l port	Run in server mode, listening on the given port. If port = 0, pick a port

    -L pfile	When running in server mode, write the port number into the given file
*/


#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

/* parameters to main -- accessed globally */
int debug = 0;
int aa    = 0;
int hits_only = 0;
long long size_hash =  1400303159; /* 1400303159  tot_lookups=13474100 retry=2981020 for 5.contigs 4.684 sec */
                                   /* 2147483648  tot_lookups=13474100 retry=1736650  */
			           /* 1073741824  tot_lookups=13474100 retry=4728020  */
int write_mem_map = 0;
char *data_dir;

#define K 8
#define MAX_SEQ_LEN 500000000

#if K == 5 
const CORE = 20L*20L*20L*20L;
#endif
#if K == 8
const unsigned long long CORE = 20L*20L*20L*20L*20L*20L*20L;
#endif

#define MAX_ENCODED CORE*20L 

static int tot_lookups = 0;
static int retry  = 0;


const  char genetic_code[64] = {
      'K','N','K','N','T','T','T','T','R','S','R','S','I','I','M','I',
      'Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L',
      'E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V',
      '*','Y','*','Y','S','S','S','S','*','C','W','C','L','F','L','F'
    };

const  char prot_alpha[20] = {
  'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y' 
    };


typedef struct sig_kmer {
  unsigned long long  which_kmer;
  int  otu_index;
  unsigned short  avg_from_end;
  int  function_index;
  float function_wt;
} sig_kmer_t;

#define VERSION 1
typedef struct kmer_memory_image {
  unsigned long long num_sigs;
  unsigned long long entry_size;
  long long  version;
} kmer_memory_image_t;

typedef struct kmer_handle {
  sig_kmer_t *kmer_table;
  unsigned long long num_sigs;
  char **function_array;   /* indexed by fI */
  char **otu_array;        /* OTU indexes point at a representation of multiple OTUs */
} kmer_handle_t;

/* the following stuff was added to condense sets of hits to specific calls.
   The basic strategy is to use a set of global variables to retain state and flush
   calls (bad, bad, bad...).
*/

typedef struct hit {
  unsigned int   oI;
  unsigned int from0_in_prot;      /* offset from start of protein sequence */
  unsigned short avg_off_from_end;  /* average offset from the end */
  unsigned int fI;
  float function_wt;
} hit_t;

#define MAX_HITS_PER_SEQ 40000
static hit_t hits[MAX_HITS_PER_SEQ]; 
static int   num_hits = 0;

#define OI_BUFSZ 5
struct otu_count {
    int oI;
    int count;
} oI_counts[OI_BUFSZ];
static int num_oI = 0;

static int   current_fI;
static char  current_id[300];
static int   current_length_contig;
static char  current_strand;
static short current_prot_off;
static int   order_constraint = 0;
static int   min_hits = 5;
static int   min_weighted_hits = 0;
static int   max_gap  = 200;

void run_accept_loop(kmer_handle_t *kmersH, in_port_t port, char *port_file, pid_t parent);
void run_from_filehandle(kmer_handle_t *kmersH, FILE *fh_in, FILE *fh_out);

/* =========================== end of reduction global variables ================= */

unsigned char to_amino_acid_off(char c) {
  switch (c)
    {
      case 'A':
	return 0;

      case 'C':
	return 1;

      case 'D':
	return 2;

      case 'E':
	return 3;

      case 'F':
	return 4;

      case 'G':
	return 5;

      case 'H':
	return 6;

      case 'I':
	return 7;

      case 'K':
	return 8;

      case 'L':
	return 9;

      case 'M':
	return 10;

      case 'N':
	return 11;

      case 'P':
	return 12;

      case 'Q':
	return 13;

      case 'R':
	return 14;

      case 'S':
	return 15;

      case 'T':
	return 16;

      case 'V':
	return 17;

      case 'W':
	return 18;

      case 'Y':
	return 19;

      default:
	return 20;
    }
}

char compl(c)
char c;
{
    switch (c)
    {
      case 'a':
	return 't';
      case 'A':
	return 'T';

      case 'c':
	return 'g';
      case 'C':
	return 'G';

      case 'g':
	return 'c';
      case 'G':
	return 'C';

      case 't':
      case 'u':
	return 'a';
      case 'T':
      case 'U':
	return 'A';

      case 'm':
	return 'k';
      case 'M':
	return 'K';

      case 'r':
	return 'y';
      case 'R':
	return 'Y';

      case 'w':
	return 'w';
      case 'W':
	return 'W';

      case 's':
	return 'S';
      case 'S':
	return 'S';

      case 'y':
	return 'r';
      case 'Y':
	return 'R';

      case 'k':
	return 'm';
      case 'K':
	return 'M';

      case 'b':
	return 'v';
      case 'B':
	return 'V';

      case 'd':
	return 'h';
      case 'D':
	return 'H';

      case 'h':
	return 'd';
      case 'H':
	return 'D';

      case 'v':
	return 'b';
      case 'V':
	return 'B';

      case 'n':
	return 'n';
      case 'N':
	return 'N';

      default:
	return c;
    }
}


void rev_comp(char *data,char *cdata) {

  int n = strlen(data);
  char *p  = data + (n-1);
  char *pc = cdata;
  while (n--) {
    *(pc++) = compl(*(p--));
  }
  *pc = 0;
}

unsigned long long encoded_kmer(unsigned char *p) {
  unsigned long long encodedK = *p;
  int i;
  for (i=1; (i <= K-1); i++) {
    encodedK = (encodedK * 20) + *(p+i);
  }

  if (encodedK > MAX_ENCODED) {
    fprintf(stderr,"bad encoding - input must have included invalid characters\n");
    for (i=0; (i < K); i++) {
      fprintf(stderr,"%d ",*(p+i));
    }
    fprintf(stderr,"\n");
    exit(2);
  }
  return encodedK;
}

unsigned long long encoded_aa_kmer(char *p)
{
  unsigned char aa_off[K];
  int j;
  for (j=0; (j < K); j++) {
    int prot_c = *(p+j);
    aa_off[j] = to_amino_acid_off(prot_c);
  }
  return encoded_kmer(aa_off);
}

void decoded_kmer(unsigned long long encodedK,char *decoded) {
  
  int i;
  *(decoded+K) = '\0';
  unsigned long long x = encodedK;

  for (i=K-1; (i >= 0); i--) {
    *(decoded+i) = prot_alpha[x % 20];
    x = x / 20;
  }
}


int dna_char(c)
char c;
{
    switch (c)
    {
      case 'a':
      case 'A':
	return 0;

      case 'c':
      case 'C':
	return 1;

      case 'g':
      case 'G':
	return 2;

      case 't':
      case 'u':
      case 'T':
      case 'U':
	return 3;

      default:
	return 4;;
    }
}

void translate(char *seq,int off,char *pseq, unsigned char *pIseq) {

  int i;
  int max = strlen(seq) - 3;
  char *p = pseq;
  unsigned char *pI = pIseq;
  for (i=off; (i <= max); ) {
    int c1 = dna_char(seq[i++]);
    int c2 = dna_char(seq[i++]);
    int c3 = dna_char(seq[i++]);
    if ((c1 < 4) && (c2 < 4) && (c3 < 4)) {
      int I = (c1 * 16) + (c2 * 4) + c3;
      char prot_c = genetic_code[I];
      *(p++) = prot_c;
      *(pI++) = to_amino_acid_off(prot_c);
    }
    else {
      *(p++)  = 'x';
      *(pI++) = 20;
    }
  }
  *p = 0;
  *pI = 21;
  if (debug >= 3) {
    fprintf(stderr,"len-seq=%d max=%d p=%ld\n",(int) strlen(seq),max,p-pseq);
  }
}

#define MAX_FUNC_OI_INDEX 1000000
#define MAX_FUNC_OI_VALS  100000000

char **load_indexed_ar(char *filename,int *sz) {
  char **index_ar = malloc(MAX_FUNC_OI_INDEX * sizeof(char *));
  char *vals      = malloc(MAX_FUNC_OI_VALS);
  char *p         = vals;
  FILE *ifp      = fopen(filename,"r");
  if (ifp == NULL) { 
    fprintf(stderr,"could not open %s\n",filename);
    exit(1);
  }

  *sz = 0;
  int j;
  while ((fscanf(ifp,"%d\t",&j) == 1) && fgets(p,1000,ifp)) {
    if (*sz != j) {
      fprintf(stderr,"Your index must be dense and in order (see line %ld, should be %d)\n",p-vals,*sz);
      exit(1);
    }
    /* fprintf(stderr,"%d is %s\n",*sz,p); */
    index_ar[*sz] = p;               /* the fgets leaves the \n at the end of each line */
    p += strlen(index_ar[*sz]) -1;
    *(p++) = '\0';
    if ((*sz >= MAX_FUNC_OI_INDEX) || ((p-vals) > (MAX_FUNC_OI_VALS - 1000))) {
      fprintf(stderr,"Your function or oI index arrays are too small; bump MAX_FUNC_OI_INDEX and MAX_FUNC_OI_VALS\n");
      exit(1);
    }

    *sz += 1;
  }
  return index_ar;
}

char **load_functions(char *file) {
  int sz;
  return load_indexed_ar(file,&sz);
}

char **load_otus(char *file) {
  int sz;
  return load_indexed_ar(file,&sz);
}

long long find_empty_hash_entry(sig_kmer_t sig_kmers[],unsigned long long encodedK) {
    long long hash_entry = encodedK % size_hash;
    while (sig_kmers[hash_entry].which_kmer <= MAX_ENCODED)
      hash_entry = (hash_entry+1)%size_hash;
    return hash_entry;
}

long long lookup_hash_entry(sig_kmer_t sig_kmers[],unsigned long long encodedK) {
    long long  hash_entry = encodedK % size_hash;
    // printf("%lld\n", size_hash);
    if (debug >= 2)
      tot_lookups++;
    while ((sig_kmers[hash_entry].which_kmer <= MAX_ENCODED) && (sig_kmers[hash_entry].which_kmer != encodedK)) {
      if (debug >= 2)
	retry++;
      hash_entry++;
      if (hash_entry == size_hash)
	hash_entry = 0;
    }
    if (sig_kmers[hash_entry].which_kmer > MAX_ENCODED) {
      return -1;
    }
    else {
      return hash_entry;
    }
}

kmer_memory_image_t *load_raw_kmers(char *file,unsigned long long num_entries, unsigned long long *alloc_sz) {
  /*
   * Allocate enough memory to hold the kmer_memory_image_t header plus the hash table itself.
   */
  *alloc_sz = sizeof(kmer_memory_image_t) + (sizeof(sig_kmer_t) * num_entries);

  kmer_memory_image_t *image = malloc(*alloc_sz);

  /*
   * Initialize our table pointer to the first byte after the header.
   */
  image->num_sigs = num_entries;
  image->entry_size = sizeof(sig_kmer_t);
  image->version = (long long) VERSION;

  sig_kmer_t *sig_kmers = (sig_kmer_t *) (image + 1);

  FILE *ifp      = fopen(file,"r");
  if (ifp == NULL) { 
    fprintf(stderr,"could not open %s",file);
    exit(1);
  }

  long long i;
  for (i=0; (i < size_hash); i++)
    sig_kmers[i].which_kmer = MAX_ENCODED + 1;

  char kmer_string[K+1];
  int end_off;
  int fI;
  float f_wt;
  int oI;
  long long loaded = 0;
  while (fscanf(ifp,"%s\t%d\t%d\t%f\t%d",
		kmer_string,&end_off,&fI,&f_wt,&oI) >= 4) {
    unsigned long long encodedK = encoded_aa_kmer(kmer_string);
    long long hash_entry = find_empty_hash_entry(sig_kmers,encodedK);
    loaded++;
    if (loaded >= (size_hash / 2)) {
      fprintf(stderr,"Your Kmer hash is half-full; use -s (and -w) to bump it\n");
      exit(1);
    }
    sig_kmers[hash_entry].which_kmer     = encodedK;
    sig_kmers[hash_entry].avg_from_end   = end_off;
    sig_kmers[hash_entry].function_index = fI;
    sig_kmers[hash_entry].otu_index      = oI;
    sig_kmers[hash_entry].function_wt    = f_wt;
  }
  if (debug >= 2)
    fprintf(stderr,"loaded %lld kmers\n",loaded);

  return image;
}

kmer_handle_t *init_kmers(char *dataD) {
  kmer_handle_t *handle = malloc(sizeof(kmer_handle_t));

  kmer_memory_image_t *image;

  char file[300];
  strcpy(file,dataD);
  strcat(file,"/function.index");
  handle->function_array = load_functions(file);

  strcpy(file,dataD);
  strcat(file,"/otu.index");
  handle->otu_array      = load_otus(file);

  char fileM[300];
  strcpy(fileM,dataD);
  strcat(fileM,"/kmer.table.mem_map");

  if (write_mem_map) {
    unsigned long long sz, table_size;
    strcpy(file,dataD);
    strcat(file,"/final.kmers");
    
    unsigned long long image_size;

    image = load_raw_kmers(file, size_hash, &image_size);

    handle->kmer_table = (sig_kmer_t *) (image + 1);
    handle->num_sigs   = image->num_sigs;

    FILE *fp = fopen(fileM,"w");
    if (fp == NULL) { 
      fprintf(stderr,"could not open %s for writing: %s ",fileM, strerror(errno));
      exit(1);
    }
    fwrite(image, image_size, 1, fp);
    fclose(fp);

    strcpy(fileM,dataD);
    strcat(fileM,"/size_hash.and.table_size");
    fp = fopen(fileM,"w");
    fprintf(fp,"%lld\t%lld\n",sz,table_size);
    fclose(fp);
  }
  else {
    int fd;
    if ((fd = open(fileM, O_RDONLY)) == -1) {
      perror("open");
      exit(1);
    }

    /*
     * Set up for creating memory image from file. Start by determining file size
     * on disk with a stat() call.
     */
    struct stat sbuf;
    if (stat(fileM, &sbuf) == -1) {
      fprintf(stderr, "stat %s failed: %s\n", fileM, strerror(errno));
      exit(1);
    }
    unsigned long long file_size = sbuf.st_size;

    /* 
     * Memory map.
     */
    int flags = MAP_SHARED;
    #ifdef MAP_POPULATE
    flags |= MAP_POPULATE;
    #endif
    
    image = (kmer_memory_image_t *) mmap((caddr_t)0, file_size, PROT_READ, flags, fd, 0);

    if (image == (kmer_memory_image_t *)(-1)) {
      fprintf(stderr, "mmap of kmer_table %s failed: %s\n", fileM, strerror(errno));
      exit(1);
    }

    /* 
     * Our image is mapped. Validate against the current version of this code.
     */
    if (image->version != (long long) VERSION) {
      fprintf(stderr, "Version mismatch for file %s: file has %lld code has %lld\n", 
	      fileM, image->version, (long long) VERSION);
      exit(1);
    }

    if (image->entry_size != (unsigned long long) sizeof(sig_kmer_t)) {
      fprintf(stderr, "Version mismatch for file %s: file has entry size %lld code has %lld\n",
	      fileM, image->entry_size, (unsigned long long) sizeof(sig_kmer_t));
      exit(1);
    }

    size_hash = image->num_sigs;
    handle->num_sigs = size_hash;
    handle->kmer_table = (sig_kmer_t *) (image + 1);

    /* Validate overall file size vs the entry size and number of entries */
    if (file_size != ((sizeof(sig_kmer_t) * image->num_sigs) + sizeof(kmer_memory_image_t))) {
      fprintf(stderr, "Version mismatch for file %s: file size does not match\n", fileM);
      exit(1);
    }

    fprintf(stderr, "Set size_hash=%lld from file size %lld\n", size_hash, file_size);

  }
  return handle;
}

void advance_past_ambig(unsigned char **p,unsigned char *bound) {

  if (K == 5) {
    while (((*p) < bound) &&
	   ((*(*p) == 20)     || 
            (*((*p)+1) == 20) || 
            (*((*p)+2) == 20) || 
            (*((*p)+3) == 20) || 
	    (*((*p)+4) == 20) )) {
      (*p)++;
    }
  }
  else {   /*  ##### ASSUMING K == 8 #### */
    int bad = 1;
    while ((*p < bound) && (bad == 1)) {
      bad = 0;
      if      (*((*p)+7) == 20) {
	bad = 1;
	(*p) += 8;
      }
      else if (*((*p)+6) == 20) {
	bad = 1;
	(*p) += 7;
      }
      else if (*((*p)+5) == 20) {
	bad = 1;
	(*p) += 6;
      }
      else if (*((*p)+4) == 20) {
	bad = 1;
	(*p) += 5;
      }
      else if (*((*p)+3) == 20) {
	bad = 1;
	(*p) += 4;
      }
      else if (*((*p)+2) == 20) {
	bad = 1;
	(*p) += 3;
      }
      else if (*((*p)+1) == 20) {
	bad = 1;
	(*p) += 2;
      }
      else if (*((*p)+0) == 20) {
	bad = 1;
	(*p) += 1;
      }
    } 
  }
}

void display_hits(FILE *fh) {
  fprintf(fh, "hits: ");
  int i;
  for (i=0; (i < num_hits); i++) {
    fprintf(fh, "%d/%f/%d ", hits[i].from0_in_prot,hits[i].function_wt,hits[i].fI);
  }
  fprintf(fh, "\n");
}


void process_set_of_hits(kmer_handle_t *kmersH, FILE *fh) {
  int fI_count = 0;
  float weighted_hits = 0;
  int last_hit=0;
  int i=0;
  while (i < num_hits) {
    if (hits[i].fI == current_fI) {
      last_hit = i;
      fI_count++;
      weighted_hits += hits[i].function_wt;
    }
    i++;
  }
  if ((fI_count >= min_hits) && (weighted_hits >= min_weighted_hits)) {
      if (!hits_only)
	  fprintf(fh, "CALL\t%d\t%d\t%d\t%d\t%s\t%f\n",hits[0].from0_in_prot,
		  hits[last_hit].from0_in_prot+(K-1),
		  fI_count,
		  current_fI,
		  kmersH->function_array[current_fI],
		  weighted_hits);

      if (debug > 1) {
	  fprintf(fh, "after-call: ");
	  display_hits(fh);
    }
    /* once we have decided to call a region, we take the kmers for fI and
       add them to the counts maintained to assign an OTU to the sequence */
    for (i=0; (i <= last_hit); i++) {
      if (hits[i].fI == current_fI) {
	int j;
	for (j=0; (j < num_oI) && (oI_counts[j].oI != hits[i].oI); j++) {}
	if (j == num_oI) {
	  if (num_oI == OI_BUFSZ) {
	    j--;   /* we overwrite the last entry */
	  }
	  else
	    num_oI++;
	  oI_counts[j].oI    = hits[i].oI;
	  oI_counts[j].count = 1;
	}
	else {
	  oI_counts[j].count++;
	}
	/* now we bubble the count back, allowing it to establish */
	while ((j > 0) && (oI_counts[j-1].count <= oI_counts[j].count)) {
	  int oI_tmp    = oI_counts[j-1].oI;
	  int count_tmp = oI_counts[j-1].count;
	  oI_counts[j-1].oI    = oI_counts[j].oI;
	  oI_counts[j-1].count = oI_counts[j].count;
	  oI_counts[j].oI    = oI_tmp;
	  oI_counts[j].count = count_tmp;
	  j--;
	}
      }
    }
  }

  if ((hits[num_hits-2].fI != current_fI) && (hits[num_hits-2].fI == hits[num_hits-1].fI)) {
      /*
    fprintf(stderr, "Copying two entries cur=%d %d: %d, %d: %d\n",
	    current_fI,
	    num_hits - 2, hits[num_hits-2].fI, 
	    num_hits - 1, hits[num_hits-1].fI);
      */
    current_fI = hits[num_hits-1].fI;
    /* now copy the last two entries to the start of the hits array.  Sorry this is so clumsy */
    hits[0].oI               = hits[num_hits-2].oI;
    hits[0].from0_in_prot    = hits[num_hits-2].from0_in_prot;
    hits[0].avg_off_from_end = hits[num_hits-2].avg_off_from_end;
    hits[0].fI               = hits[num_hits-2].fI;
    hits[0].function_wt      = hits[num_hits-2].function_wt;
      
    hits[1].oI               = hits[num_hits-1].oI;
    hits[1].from0_in_prot    = hits[num_hits-1].from0_in_prot;
    hits[1].avg_off_from_end = hits[num_hits-1].avg_off_from_end;
    hits[1].fI               = hits[num_hits-1].fI;
    hits[1].function_wt      = hits[num_hits-1].function_wt;
    num_hits                 = 2;
  }
  else {
    num_hits = 0;
  }
}

void gather_hits(int ln_DNA, char strand,int prot_off,char *pseq,
		 unsigned char *pIseq, kmer_handle_t *kmersH, FILE *fh) {
  
  if (debug >= 3) {
      fprintf(fh, "translated: %c\t%d\t%s\n",strand,prot_off,pseq);
  }

  unsigned char *p = pIseq;
 /* pseq and pIseq are the same length */

  unsigned char *bound = pIseq + strlen(pseq) - K;
  advance_past_ambig(&p,bound);
  unsigned long long encodedK=0;
  if (p < bound) {
    encodedK = encoded_kmer(p);
  }
  while (p < bound) {
    long long  where = lookup_hash_entry(kmersH->kmer_table,encodedK);
    // printf("%lu %lld\n", p - pIseq, where);
    if (where >= 0) {
      sig_kmer_t *kmers_hash_entry = &(kmersH->kmer_table[where]);
      int avg_off_end = kmers_hash_entry->avg_from_end;
      int fI        = kmers_hash_entry->function_index;
      int oI          = kmers_hash_entry->otu_index;
      float f_wt      = kmers_hash_entry->function_wt;
      if (debug >= 1) {
	  if (hits_only)
	      fprintf(fh, "%ld\t%s\n",encodedK, current_id);
	  else
	      fprintf(fh, "HIT\t%ld\t%lld\t%d\t%d\t%0.3f\t%d\n",p-pIseq,encodedK,avg_off_end,fI,f_wt,oI);
      }

      if ((num_hits > 0) && (hits[num_hits-1].from0_in_prot + max_gap) < (p-pIseq)) {
	if (num_hits >= min_hits) {
	    // fprintf(stderr, "pset from %d cur=%d\n",  __LINE__, current_fI);
	    process_set_of_hits(kmersH, fh);
	}
	else {
	  num_hits = 0;
	}
      }

      if (num_hits == 0) {
	current_fI = fI;   /* if this is the first, set the current_fI */
      }

      if ((! order_constraint) || (num_hits == 0) ||
	  ((fI == hits[num_hits-1].fI) &&
	   (abs(((p-pIseq) - hits[num_hits-1].from0_in_prot) - 
                (hits[num_hits-1].avg_off_from_end - avg_off_end)
		) <= 20))) {
          /* we have a new hit, so we add it to the global set of hits */
	hits[num_hits].oI = oI;
	hits[num_hits].fI = fI;
	hits[num_hits].from0_in_prot = p-pIseq;
	hits[num_hits].avg_off_from_end = avg_off_end;
	hits[num_hits].function_wt = f_wt;
        if (num_hits < MAX_HITS_PER_SEQ - 2) 
	  num_hits++;
	if (debug > 1) {
	    fprintf(fh, "after-hit: ");
	    display_hits(fh);
	}
	if ((num_hits > 1) && (current_fI != fI) &&           /* if we have a pair of new fIs, it is time to */
	    (hits[num_hits-2].fI == hits[num_hits-1].fI)) {   /* process one set and initialize the next */
	    // fprintf(stderr, "pset from %d cur=%d\n",  __LINE__, current_fI);
	    process_set_of_hits(kmersH, fh);
	}
      }
    }
    p++;
    if (p < bound) {
      if (*(p+K-1) < 20) {
	encodedK = ((encodedK % CORE) * 20L) + *(p+K-1);
      }
      else {
	p += K;
	advance_past_ambig(&p,bound);
	if (p < bound) {
	  encodedK = encoded_kmer(p);
	}
      }
    }    
  }
  if (num_hits >= min_hits) {
      // fprintf(stderr, "pset from %d cur=%d\n",  __LINE__, current_fI);
      process_set_of_hits(kmersH, fh);
  }
  num_hits = 0;
}

void tabulate_otu_data_for_contig(FILE *fh) {
  int i;
  if (!hits_only)
  {
      fprintf(fh, "OTU-COUNTS\t%s[%d]",current_id,current_length_contig);
      for (i=0; (i < num_oI); i++) {
	  fprintf(fh, "\t%d-%d",oI_counts[i].count,oI_counts[i].oI);
      }
      fprintf(fh, "\n");
  }
  num_oI = 0;
}

void process_aa_seq(char *id,char *pseq,size_t ln,kmer_handle_t *kmersH, FILE *fh) {
  static unsigned char *pIseq = 0;
  if (pIseq == 0)
  {
    pIseq = malloc(MAX_SEQ_LEN / 3);
  }
  //static unsigned char pIseq[MAX_SEQ_LEN / 3];

  strcpy(current_id,id);
  if (!hits_only)
      fprintf(fh, "PROTEIN-ID\t%s\t%d\n",id,ln);

  current_length_contig = ln;
  current_strand        = '+';
  current_prot_off      = 0;
  int i;
  for (i=0; (i < ln); i++)
    pIseq[i] = to_amino_acid_off(*(pseq+i));
  gather_hits(ln,'+',0,pseq,pIseq,kmersH,fh);  
  tabulate_otu_data_for_contig(fh);
}

void process_seq(char *id,char *data,kmer_handle_t *kmersH, FILE *fh) {

  static char *cdata = 0;
  static char *pseq = 0;
  static unsigned char *pIseq = 0;

  if (cdata == 0)
  {
    cdata = malloc(MAX_SEQ_LEN);
    pseq = malloc(MAX_SEQ_LEN / 3);
    pIseq = malloc(MAX_SEQ_LEN / 3);
  }
     
//  static char cdata[MAX_SEQ_LEN];
//  static char pseq[MAX_SEQ_LEN / 3];
//  static unsigned char pIseq[MAX_SEQ_LEN / 3];

  strcpy(current_id,id);
  int ln = strlen(data);
  current_length_contig = ln;
  fprintf(fh, "processing %s[%d]\n",id,ln);
  int i;
  for (i=0; (i < 3); i++) {
    
    translate(data,i,pseq,pIseq);
    current_strand   = '+';
    current_prot_off = i;
    if (!hits_only)
	fprintf(fh, "TRANSLATION\t%s\t%d\t%c\t%d\n",current_id,
		current_length_contig,
		current_strand,
		current_prot_off);
    gather_hits(ln,'+',i,pseq,pIseq,kmersH, fh);
  }
  rev_comp(data,cdata);
  for (i=0; (i < 3); i++) {
    translate(cdata,i,pseq,pIseq);

    current_strand   = '-';
    current_prot_off = i;
    if (!hits_only)
	fprintf(fh, "TRANSLATION\t%s\t%d\t%c\t%d\n",current_id,
		current_length_contig,
		current_strand,
		current_prot_off);
    gather_hits(ln,'-',i,pseq,pIseq,kmersH, fh);
  }
  tabulate_otu_data_for_contig(fh);
}

int main(int argc,char *argv[]) {
  int c;
  char *past;
  char file[300];
  int is_server = 0;
  in_port_t port;
  char port_file[1024];
  pid_t parent = -1;

  port_file[0] = 0;
  file[0] = 0;

  while ((c = getopt (argc, argv, "ad:s:wD:m:g:OM:l:L:P:H")) != -1) {
    switch (c) {
    case 'a':
      aa = 1;
      break;
    case 'H':
      hits_only = 1;
      break;
    case 'd':
      debug = strtol(optarg,&past,0);
      break;
    case 'l':
	port = atoi(optarg);
	is_server = 1;
	break;

    case 'L':
	strncpy(port_file, optarg, sizeof(port_file) -1 );
	port_file[sizeof(port_file) - 1] = 0;
	break;

    case 'P':
	parent = atoi(optarg);
	break;
	
    case 'm':
      min_hits = strtol(optarg,&past,0);
      break;
    case 'M':
      min_weighted_hits = strtol(optarg,&past,0);
      break;
    case 'O':
      order_constraint = 1;
      break;
    case 'g':
      max_gap = strtol(optarg,&past,0);
      break;
    case 'D':
      strcpy(file,optarg);
      break;
    case 's':
      size_hash = strtol(optarg,&past,0);
      break;
    case 'w':
      write_mem_map = 1;
      break;
    default:
      fprintf(stderr,"arguments: [-a] [-d level] [-s hash-size] [-w] [-m min_hits] -D DataDir \n");
      abort ();
    }
  }

  kmer_handle_t *kmersH = init_kmers(file);

  if (is_server)
  {
      run_accept_loop(kmersH, port, port_file, parent);
  }
  else
  {
      run_from_filehandle(kmersH, stdin, stdout);
  }
  return 0;
}

void run_from_filehandle(kmer_handle_t *kmersH, FILE *fh_in, FILE *fh_out)
{
  static char *data = 0;
  if (data == 0)
  {
    data = malloc(MAX_SEQ_LEN);
  }
  
  // static char data[MAX_SEQ_LEN];
  int got_gt = 0;
  int i;
  char *p;
  char id[2000];

  while (((!got_gt && (fscanf(fh_in,">%s",id) == 1)) ||
	  (got_gt && (fscanf(fh_in,"%s",id) == 1)))) {
    while (getc(fh_in) != '\n')
      ;
    if ((id[0] == 'F') &&
	(id[1] == 'L') &&
	(id[2] == 'U') &&
	(id[3] == 'S') &&
	(id[4] == 'H')) {   /* ugly compare, since \n is included in id */

/*      int rc = fflush(fh_in);
      if (rc != 0) {
	fprintf(stderr,"fflush did not seem to work\n");
      }
*/
	fprintf(fh_out, "//\n");
	got_gt = 0;
    }
    else {

      for (p=data; ((i = getc(fh_in)) != -1) && (i != '>');) {
	if ((i != ' ') && (i != '\n'))
	    *(p++) = toupper(i);
      }

      if (i == '>')
	got_gt = 1;
      else
	got_gt = 0;

      *p=0;
      size_t len = p - data;
      if ((len) > MAX_SEQ_LEN) {
	fprintf(stderr,"The contig size exceeds %d; bump MAX_SEQ_LEN\n",MAX_SEQ_LEN);
	exit(1);
      }

      /* fprintf(stderr,"%d bytes read\nends with %s\n",(p-data),p-50); */

      if (! aa)
	  process_seq(id,data,kmersH, fh_out);
      else
	  process_aa_seq(id,data, len, kmersH, fh_out);
      fflush(fh_out);
    }
  }

  if (debug >= 2)
      fprintf(fh_out, "tot_lookups=%d retry=%d\n",tot_lookups,retry);
}

void run_accept_loop(kmer_handle_t *kmersH, in_port_t port, char *port_file, pid_t parent)
{
    int listenfd = 0, connfd = 0;
    struct sockaddr_in serv_addr;

    signal(SIGPIPE, SIG_IGN);

    listenfd = socket(AF_INET, SOCK_STREAM, 0);
    memset(&serv_addr, 0, sizeof(serv_addr));

    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
    serv_addr.sin_port = htons(port);

    if (bind(listenfd, (struct sockaddr*)&serv_addr, sizeof(serv_addr)) < 0)
    {
	perror("bind failed");
	exit(1);
    }
 
    struct sockaddr_in my_addr;
    socklen_t my_len = sizeof(my_addr);
    if (getsockname(listenfd, (struct sockaddr *) &my_addr, &my_len) < 0)
    {
	perror("getsockname failed");
	exit(1);
    }
    in_port_t my_port = ntohs(my_addr.sin_port);
    printf("Listening on %d\n", my_port);
    if (port_file[0])
    {
	FILE *fp = fopen(port_file, "w");
	if (!fp)
	{
	    fprintf(stderr, "error opening %s for writing: %s\n", port_file, strerror(errno));
	    exit(1);
	}
	fprintf(fp, "%d\n", my_port);
	fclose(fp);
    }

    listen(listenfd, 10);

    int save_aa = aa;
    int save_hits_only = hits_only;
    int save_debug = debug;
    int save_min_hits = min_hits;
    int save_min_weighted_hits = min_weighted_hits;
    int save_order_constraint = order_constraint;
    int save_max_gap = max_gap;
    
    while(1)
    {
	/*
	 * If we have a parent set, and parent doesn't exist, exit.
	 */
	if (parent > 0)
	{
	    int rc = kill(parent, 0);
	    if (rc < 0)
	    {
		fprintf(stderr, "Parent process %d does not exist any more, exiting\n", parent);
		exit(0);
	    }
	}
	
      aa = save_aa;
      hits_only = save_hits_only;
      debug = save_debug;
      min_hits = save_min_hits;
      min_weighted_hits = save_min_weighted_hits;
      order_constraint = save_order_constraint;
      max_gap = save_max_gap;
      
      connfd = accept(listenfd, (struct sockaddr*)NULL, NULL);
      struct sockaddr_in peer;
      socklen_t peer_len = sizeof(peer);
      memset(&peer, 0, sizeof(peer));
      
      getpeername(connfd, (struct sockaddr *) &peer, &peer_len);

      char *who = inet_ntoa(peer.sin_addr);
      // fprintf(stderr, "connection from %s\n", who);
	
      FILE *fh_in = fdopen(connfd, "r");
      FILE *fh_out = fdopen(connfd, "w");

      /*
       * If the first line starts with a '-', then it's
       * a set of options to be set for this document.
       */

      int c = fgetc(fh_in);
      if (c == '-')
      {
	char linebuf[1024];
	linebuf[0] = c;
	  
	if (fgets(linebuf + 1, sizeof(linebuf) - 2, fh_in) == 0)
	{
	  fprintf(stderr, "Error reading options line from %s\n", who);
	  fclose(fh_in);
	  fclose(fh_out);
	  continue;
	}

	char *s = strchr(linebuf, '\n');
	if (s)
	  *s = 0;
	  
	/* parse into an argv */
	const int max_args = 20;
	char *argv[max_args + 2];
	int n = 0;
	// fprintf(stderr, "Parsing option line '%s'\n", linebuf);
	argv[n++] = "nothing";

	s = strtok(linebuf, " \t");
	while (s)
	{
	    // fprintf(stderr, "argv[%d] = '%s'\n", n, s);
	  argv[n++] = s;
	  if (n >= max_args)
	  {
	    fprintf(stderr, "too many args in connection from %s\n", who);
	    fprintf(fh_out, "ERR too many args\n");
	    fflush(fh_out);
	    fclose(fh_out);
	    fclose(fh_in);
	    s = 0;
	    break;
	  }
	  s = strtok(0, " \t");
	}
	if (n >= max_args)
	{
	  continue;
	}
	argv[n] = 0;

	char *past;
	int arg_error = 0;

	optind = 1;
	while ((c = getopt(n, argv, "ad:m:M:Og:")) != -1)
	{
	  switch (c) {
	  case 'a':
	    aa = 1;
	    break;
	  case 'd':
	    debug = strtol(optarg,&past,0);
	    break;
	  case 'm':
	    min_hits = strtol(optarg,&past,0);
	    break;
	  case 'M':
	    min_weighted_hits = strtol(optarg,&past,0);
	    break;
	  case 'O':
	    order_constraint = 1;
	    break;
	  case 'g':
	    max_gap = strtol(optarg,&past,0);
	    break;
	  default:
	    fprintf(fh_out, "ERR invalid argument %c\n", c);
	    fflush(fh_out);
	    fclose(fh_out);
	    fclose(fh_in);
	    arg_error = 1;
	    break;
	  }
	}
	if (arg_error)
	  continue;
	if (!hits_only)
	    fprintf(fh_out, "OK aa=%d debug=%d min_hits=%d min_weighted_hits=%d order_constraint=%d max_gap=%d\n",
		    aa, debug, min_hits, min_weighted_hits, order_constraint, max_gap);
      }
      else
      {
	ungetc(c, fh_in);
      }

      
      run_from_filehandle(kmersH, fh_in, fh_out);

      fflush(fh_out);
      fclose(fh_in);
      fclose(fh_out);

    }
}
