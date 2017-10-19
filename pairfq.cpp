#include <fstream>
#include <iostream>
#include <deque>
#include <string>
#include <stdexcept>
#include <cstring>

using namespace std;


struct fastq_entry {
    /// read name
    string qname;
    /// read sequence
    string seq;
    /// read quality scores
    string qual;
};

/**
 * Natural string comparison, account for numbers within a string.
 * This is different from lexigraphical comparison (e.g. std::string::compare).
 * Taken from https://github.com/samtools/samtools/blob/develop/bam_sort.c#L13
 */
int strnum_cmp(const char *_a, const char *_b) {
    const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
    const unsigned char *pa = a, *pb = b;
    while (*pa && *pb) {
        if (isdigit(*pa) && isdigit(*pb)) {
            while (*pa == '0') ++pa;
            while (*pb == '0') ++pb;
            while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
            if (isdigit(*pa) && isdigit(*pb)) {
                int i = 0;
                while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
                return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
            } else if (isdigit(*pa)) return 1;
            else if (isdigit(*pb)) return -1;
            else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
        } else {
            if (*pa != *pb) return (int)*pa - (int)*pb;
            ++pa; ++pb;
        }
    }
    return *pa? 1 : *pb? -1 : 0;
}

/**
 * Compare two read names by lexical or natural order.
 */
int compare_qnames(const string& a, const string& b, bool natural) {
    size_t a_end = a.length();
    size_t b_end = b.length();
    if (a.substr(a_end - 2) == "/1") {
        a_end -= 2;
    }
    if (b.substr(b_end - 2) == "/2") {
        b_end -= 2;
    }
    
    if (natural) {
        return strnum_cmp(a.substr(0, a_end).c_str(), b.substr(0, b_end).c_str());
    }
    return a.substr(0, a_end).compare(b.substr(0, b_end));
}

/** 
 * Read one fastq entry from file.
 *
 * Assume that each sequence and quality score entry `r1.fq` and `r2.fq` is single-line.
 */
bool read_fastq_entry(istream& f, fastq_entry& x) {
    getline(f, x.qname);
    if (x.qname.empty()) return false;
    
    getline(f, x.seq);

    string marker;
    getline(f, marker);
    if (marker != "+") {
        throw runtime_error("fastq entry is malformed");
    }

    getline(f, x.qual);
    
    return true;
}

/** 
 * Write one fastq entry to file.
 */
void write_fastq_entry(ostream& f, const fastq_entry& x) {
    f << x.qname << endl;
    f << x.seq << endl;
    f << '+' << endl;
    f << x.qual << endl;
}

/** 
 * Flush all fastq entries from deque and writea them to file.
 */
void flush_fastq_entries(ostream& f, deque<fastq_entry>& q) {
    while(!q.empty()) {
        write_fastq_entry(f, q.back());
        q.pop_back();
    }
}

/**
 * Read the next fastq entry from file and push onto front of the deque.
 */
inline bool read_next_fastq(istream& f, deque<fastq_entry>& q) {
    q.push_front(fastq_entry());
    if (read_fastq_entry(f, q.front())) {
        return true;
    }
    // read failed: remove empty entry
    q.pop_front();
    return false;
}

/**
 * Pair first-read and second-second fastq files together, separating unpaired reads.
 *
 * Assume that each sequence and quality score entry `r1.fq` and `r2.fq` is single-line.
 *
 * Assume that `r1.fq` and `r2.fq` are both sorted according `sort_type`
 * (`lexical` or `natural`).
 * If `sort_type == natural`, then the reads in the source BAM file should 
 * be been sorted by `sambamba sort -N` or `samtools sort -n`).
 * If `sort_type == lexical`, then the reads in should be been sorted
 * by `sambamba sort -n` or `java -jar Picard.jar RevertSam
 * SORT_ORDER=queryname`.
 *
 * The output file `out.fq` will contain interleaved read pairs with 
 * unpaired reads removed, and this file should be suitable for alignment 
 * with `bwa mem -p`.
 *
 *
 * Unpaired reads will be written to `unpaired.fq` if it has not been set to
 * `null`.
 */
int main(int argc, char* argv[]) {

    char* r1_fname;
    char* r2_fname;
    char* out_fname;
    char* unpaired_fname = NULL;
    bool natural_sort = false;

    int nargs = argc - 1;
    if (nargs < 5) {
        cerr << "usage: pairfq <r1.fq> <r2.fq> <out.fq> <sort_type> <unpaired.fq>" << endl;
        return 1;
    } else {
        r1_fname = argv[1];
        r2_fname = argv[2];
        out_fname = argv[3];
        if (strcmp(argv[4], "natural") == 0) {
            natural_sort = true;
        }
        if (strcmp(argv[5], "null") != 0) {
            unpaired_fname = argv[5];
        }
    }
    
    // open input and output files
    
    ifstream r1f(r1_fname);
    ifstream r2f(r2_fname);
    ofstream outf(out_fname);
    ofstream unpairedf;
    bool write_unpaired = false;
    if (unpaired_fname != NULL) {
        unpairedf.open(unpaired_fname);
        write_unpaired = true;
    }
    
    // process reads
    
    deque<fastq_entry> r1q, r2q;

    while (true) {
        if (r1q.empty()) {
            // get the next fastq entry from list 1
            if (!read_next_fastq(r1f, r1q)) break;
        }

        if (r2q.empty()) {
            // get the next fastq entry from list 2
          if (!read_next_fastq(r2f, r2q)) break;
        }
        
        int cmp = compare_qnames(r1q.front().qname, r2q.front().qname, natural_sort);

        if (cmp < 0) {
            // current element of list 1 is smaller:
            // get the next element of list 2 for comparison
            if (!read_next_fastq(r1f, r1q)) break; 
        } else if (cmp > 0) {
            // current element of list 2 is smaller:
            // get the next element of list 1 for comparison
            if (!read_next_fastq(r2f, r2q)) break; 
        } else {
            // elements of list 1 and list 2 match: write matched pair to file
            write_fastq_entry(outf, r1q.front());
            r1q.pop_front();
            write_fastq_entry(outf, r2q.front());
            r2q.pop_front();
            // deques now possibly contain unpaired reads: write them out
            if (write_unpaired) {
                flush_fastq_entries(unpairedf, r1q);
                flush_fastq_entries(unpairedf, r2q);
            } else {
                r1q.clear();
                r2q.clear();
            }
        }
    }
    
    // flush the remaining reads
    if (write_unpaired) {
        flush_fastq_entries(unpairedf, r1q);
        flush_fastq_entries(unpairedf, r2q);
    } else {
        r1q.clear();
        r2q.clear();
    }

    // clean up

    r1f.close();
    r2f.close();
    outf.close();
    if (write_unpaired) {
        unpairedf.close();
    }

    return 0;
}
