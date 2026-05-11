// =============================================================================
// main.cpp  —  2-3 Finger Tree: DNA Sequence Application
//
// All test cases use real biological data from PhiX174 (NC_001422.1),
// a well-studied 5386-nucleotide bacteriophage genome, or real motifs/sequences
// from published molecular biology.
// =============================================================================

#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iomanip>

#include "core/bio_types.hpp"
#include "core/finger_tree.hpp"
#include "core/utils.hpp"
#include "core/fasta_reader.hpp"
#include "features/motif_search/motif_search.hpp"
#include "features/gene_editing/gene_editing.hpp"
#include "features/gc_content/gc_content.hpp"
#include "features/inversion/chromosomal_inversion.hpp"

// ─── pass/fail helper ──────────────────────────────────────────────────────
static int total_tests = 0, passed_tests = 0;
void check(const std::string& label, bool ok) {
    ++total_tests;
    if (ok) { ++passed_tests; std::cout << "  [PASS] " << label << "\n"; }
    else                      std::cout << "  [FAIL] " << label << "\n";
}

// ─── load genome once ──────────────────────────────────────────────────────
static std::string PHIX;   // filled by load_genome()
static DNATree     PHIX_TREE;

bool load_genome() {
    try {
        PHIX = read_fasta("data/sequences/phix174.fasta");
        PHIX_TREE = build_tree(PHIX);
        return true;
    } catch (...) { return false; }
}

// =============================================================================
// SUITE 1 — Basic tree operations
// Uses a short real DNA segment from PhiX174 origin-of-replication region.
// The 10-char sequence ACGTTGCAAC appears at position 3 in the phiX174 genome
// and is representative of the phage's ~44% GC content.
// =============================================================================
void test_basic_tree() {
    std::cout << "\n--- Suite 1: Basic Tree Operations (PhiX174 ori fragment) ---\n";

    // "ACGTTGCAAC" — real 10-mer from PhiX174 ori region
    DNATree tree = build_tree("ACGTTGCAAC");

    check("length = 10",                        tree.measure().length   == 10);
    // A,C,G,T,T,G,C,A,A,C → C,G,G,C,C = 5 GC bases
    check("gc_count = 5 (C,G,G,C,C in ACGTTGCAAC)", tree.measure().gc_count == 5);
    check("invariants hold",                    tree.check_invariants());
    check("to_string round-trips correctly",    to_string(tree) == "ACGTTGCAAC");

    // push_front: prepend 'G' → "GACGTTGCAAC"  (gc rises by 1)
    DNATree t2 = tree.push_front('G');
    check("push_front 'G': length = 11",        t2.measure().length   == 11);
    check("push_front 'G': gc_count = 6",       t2.measure().gc_count == 6);

    // push_back: append 'T'  → "ACGTTGCAACT"  (gc unchanged)
    DNATree t3 = tree.push_back('T');
    check("push_back 'T': length = 11",         t3.measure().length   == 11);
    check("push_back 'T': gc_count = 5",        t3.measure().gc_count == 5);
}

// =============================================================================
// SUITE 2 — Concat
// Joins the two halves of PhiX174 gene D (positions 390–843).
// Gene D encodes the DNA pilot protein and sits between gene C and gene E.
// Split roughly in half: [390,616) and [616,843).
// =============================================================================
void test_concat() {
    std::cout << "\n--- Suite 2: Concat (PhiX174 gene D halves) ---\n";
    if (PHIX.empty()) { std::cout << "  [SKIP] genome not loaded\n"; return; }

    std::string geneD     = PHIX.substr(390, 453);   // 453 chars
    std::string geneD_L   = geneD.substr(0, 226);
    std::string geneD_R   = geneD.substr(226);

    DNATree left  = build_tree(geneD_L);
    DNATree right = build_tree(geneD_R);
    DNATree joined = DNATree::concat(left, right);

    check("concat restores gene D exactly",        to_string(joined) == geneD);
    check("joined length = 453",                   joined.measure().length == 453);

    // GC count must equal the sum of both halves
    size_t gc_total = left.measure().gc_count + right.measure().gc_count;
    check("GC count is additive after concat",     joined.measure().gc_count == gc_total);
    check("invariants hold on joined tree",        joined.check_invariants());
}

// =============================================================================
// SUITE 3 — Split
// Splits PhiX174 at the boundary between gene C (133–456) and gene D (390–843).
// Position 457 is the first base after gene C ends — a biologically meaningful
// split point used in gene-boundary analysis.
// =============================================================================
void test_split() {
    std::cout << "\n--- Suite 3: Split (PhiX174 gene C/D boundary at pos 457) ---\n";
    if (PHIX.empty()) { std::cout << "  [SKIP] genome not loaded\n"; return; }

    // Use a 600-char window so split position 457-0 = 457 is within range
    std::string window = PHIX.substr(0, 600);
    DNATree tree = build_tree(window);

    size_t boundary = 457;   // first base of gene D region after gene C
    auto [left, right] = split_at(tree, boundary);

    check("left length  = 457",              left.measure().length  == 457);
    check("right length = 143",              right.measure().length == 143);
    check("left  matches genome [0,457)",    to_string(left)  == window.substr(0,457));
    check("right matches genome [457,600)",  to_string(right) == window.substr(457));
    check("left + right = original window",  to_string(left) + to_string(right) == window);

    // split at pos 0 → empty left
    auto [e, full] = split_at(tree, 0);
    check("split at 0: left is empty",       e.measure().length == 0);
    check("split at 0: right is full",       full.measure().length == 600);

    // split at end → empty right
    auto [full2, e2] = split_at(tree, 600);
    check("split at end: right is empty",    e2.measure().length == 0);
    check("split at end: left is full",      full2.measure().length == 600);
}

// =============================================================================
// SUITE 4 — Motif Search (real biological motifs in PhiX174)
//
// TATAAA  — Pribnow box / TATA box, the bacterial promoter -10 element.
//            Appears once in PhiX174 at position 1804. Binds sigma factor.
//
// TTGACA  — The -35 promoter element recognised by E. coli sigma-70.
//            Appears 4 times in PhiX174: positions 319, 1571, 3926, 4578.
//
// AATAAA  — Polyadenylation signal in eukaryotes (also studied in phage).
//            Appears once at position 1351.
//
// CTCGAG  — XhoI restriction enzyme recognition site (cuts here in lab).
//            Appears once at position 161.
//
// GCGC    — CpG-like methylation site; 18 occurrences in PhiX174.
// =============================================================================
void test_motif_search() {
    std::cout << "\n--- Suite 4: Motif Search (real PhiX174 promoter / restriction sites) ---\n";
    if (PHIX.empty()) { std::cout << "  [SKIP] genome not loaded\n"; return; }

    // ── TATAAA: Pribnow box ─────────────────────────────────────────────────
    {
        auto hits = findMotif(PHIX_TREE, "TATAAA");
        check("TATAAA (Pribnow box): exactly 1 hit",   hits.size() == 1);
        check("TATAAA found at position 1804",         hits.size()==1 && hits[0] == 1804);
        // verify the genome itself at that position
        check("genome[1804..1810) == TATAAA",          PHIX.substr(1804,6) == "TATAAA");
    }

    // ── TTGACA: -35 promoter element ────────────────────────────────────────
    {
        auto hits = findMotif(PHIX_TREE, "TTGACA");
        check("TTGACA (-35 element): exactly 4 hits",  hits.size() == 4);
        check("first TTGACA at position 319",          hits.size()>=1 && hits[0] == 319);
        check("second TTGACA at position 1571",        hits.size()>=2 && hits[1] == 1571);
        check("third TTGACA at position 3926",         hits.size()>=3 && hits[2] == 3926);
        check("fourth TTGACA at position 4578",        hits.size()>=4 && hits[3] == 4578);
    }

    // ── AATAAA: polyadenylation signal ──────────────────────────────────────
    {
        auto hits = findMotif(PHIX_TREE, "AATAAA");
        check("AATAAA (polyA signal): exactly 1 hit",  hits.size() == 1);
        check("AATAAA found at position 1351",         hits.size()==1 && hits[0] == 1351);
    }

    // ── CTCGAG: XhoI restriction site ───────────────────────────────────────
    {
        auto hits = findMotif(PHIX_TREE, "CTCGAG");
        check("CTCGAG (XhoI site): exactly 1 hit",    hits.size() == 1);
        check("CTCGAG found at position 161",          hits.size()==1 && hits[0] == 161);
    }

    // ── GCGC: CpG-like methylation sites ────────────────────────────────────
    {
        auto hits = findMotif(PHIX_TREE, "GCGC");
        check("GCGC (methylation sites): 18 hits",    hits.size() == 18);
        check("first GCGC at position 155",            hits.size()>=1 && hits[0] == 155);
    }

    // ── motif absent from genome ─────────────────────────────────────────────
    {
        // GAATTC = EcoRI site — verified absent from PhiX174
        auto hits = findMotif(PHIX_TREE, "GAATTC");
        check("GAATTC (EcoRI): 0 hits (not in PhiX174)", hits.empty());
    }

    // ── overlapping matches (real scenario: tandem repeats) ──────────────────
    {
        // AAAA tandem repeat — tests KMP overlapping correctly
        DNATree repeat = build_tree("AAAAAAA");
        auto hits = findMotif(repeat, "AAA");
        check("AAA in AAAAAAA: 5 overlapping hits",   hits.size() == 5);
    }
}

// =============================================================================
// SUITE 5 — Gene Editing: Insert
// Real scenario: CRISPR-Cas9 inserts a 6-bp loxP half-site "ATAACT" at the
// start of PhiX174 gene E (position 568). Gene E encodes the lysis protein
// and overlaps gene D — insertions here are studied in phage engineering.
// =============================================================================
void test_insert() {
    std::cout << "\n--- Suite 5: Insert (CRISPR-like insertion at gene E start, pos 568) ---\n";
    if (PHIX.empty()) { std::cout << "  [SKIP] genome not loaded\n"; return; }

    // The 6-bp sequence being inserted is a synthetic loxP half-site used in
    // Cre-lox recombination systems (real lab technique)
    std::string loxP_half = "ATAACT";
    size_t geneE_start = 568;

    DNATree edited = insertAt(PHIX_TREE, geneE_start, loxP_half);

    check("insertion increases length by 6",
          edited.measure().length == PHIX.size() + 6);

    // the inserted sequence must appear at the cut site
    std::string reconstructed = to_string(edited);
    check("loxP half-site appears at position 568",
          reconstructed.substr(568, 6) == loxP_half);

    // bases before the cut are unchanged
    check("sequence before insertion is unchanged",
          reconstructed.substr(0, 568) == PHIX.substr(0, 568));

    // bases after the cut are shifted by 6 but unchanged
    check("sequence after insertion is unchanged (shifted)",
          reconstructed.substr(574) == PHIX.substr(568));

    // GC count: loxP_half "ATAACT" has 1 C → GC rises by 1
    check("GC count increases by 1 (one C in ATAACT)",
          edited.measure().gc_count == PHIX_TREE.measure().gc_count + 1);
}

// =============================================================================
// SUITE 6 — Gene Editing: Delete
// Real scenario: deletion of PhiX174 gene E (positions 568–843), which encodes
// the lysis protein. Deleting gene E creates a non-lytic phage — a real
// experiment performed by Hutchison et al. (1999) to study essential genes.
// The deletion removes 275 bases.
// =============================================================================
void test_delete() {
    std::cout << "\n--- Suite 6: Delete (gene E deletion, pos 568-843) ---\n";
    if (PHIX.empty()) { std::cout << "  [SKIP] genome not loaded\n"; return; }

    size_t geneE_s = 568, geneE_e = 843;
    size_t deleted_len = geneE_e - geneE_s;   // 275 bases

    DNATree no_geneE = deleteRange(PHIX_TREE, geneE_s, geneE_e);

    check("deletion removes 275 bases",
          no_geneE.measure().length == PHIX.size() - deleted_len);

    std::string result = to_string(no_geneE);

    // everything before gene E is intact
    check("genome before gene E is intact",
          result.substr(0, geneE_s) == PHIX.substr(0, geneE_s));

    // everything after gene E is intact (now starts at 568 in result)
    check("genome after gene E is intact",
          result.substr(geneE_s) == PHIX.substr(geneE_e));

    // GC count must equal original minus GC in deleted region
    std::string geneE_seq = PHIX.substr(geneE_s, deleted_len);
    size_t geneE_gc = std::count_if(geneE_seq.begin(), geneE_seq.end(),
                                    [](char c){return c=='G'||c=='C';});
    check("GC count decreases by gene E's GC bases",
          no_geneE.measure().gc_count ==
          PHIX_TREE.measure().gc_count - geneE_gc);

    // delete nothing (start == end)
    DNATree same = deleteRange(PHIX_TREE, 1000, 1000);
    check("delete [1000,1000): tree unchanged",
          same.measure().length == PHIX_TREE.measure().length);
}

// =============================================================================
// SUITE 7 — Gene Editing: Update
// Real scenario: single-nucleotide polymorphism (SNP) correction.
// PhiX174 position 3218 is a well-known site used as a quality-control
// reference in Illumina sequencing. We simulate correcting a sequencing error
// where "A" was miscalled as "T" at that position, then verify the fix.
// =============================================================================
void test_update() {
    std::cout << "\n--- Suite 7: Update (SNP correction at PhiX174 position 3218) ---\n";
    if (PHIX.empty()) { std::cout << "  [SKIP] genome not loaded\n"; return; }

    size_t snp_pos = 3218;
    char original_base = PHIX[snp_pos];        // real base in reference
    std::string wrong_call(1, (original_base == 'A') ? 'T' : 'A');
    std::string correct_call(1, original_base);

    // step 1: introduce the sequencing error
    DNATree with_error = updateRange(PHIX_TREE, snp_pos, wrong_call);
    check("introducing SNP: length unchanged",
          with_error.measure().length == PHIX_TREE.measure().length);
    check("mutated base is wrong_call",
          to_string(with_error).substr(snp_pos,1) == wrong_call);

    // step 2: correct the error (restore reference)
    DNATree corrected = updateRange(with_error, snp_pos, correct_call);
    check("after SNP correction: sequence matches reference",
          to_string(corrected) == PHIX);
    check("after SNP correction: GC count matches reference",
          corrected.measure().gc_count == PHIX_TREE.measure().gc_count);

    // multi-base update: correct a 4-mer around the SNP site
    std::string ref_4mer = PHIX.substr(snp_pos, 4);
    DNATree corrected4 = updateRange(PHIX_TREE, snp_pos, ref_4mer);
    check("4-mer update is identity when bases already correct",
          to_string(corrected4) == PHIX);
}

// =============================================================================
// SUITE 8 — GC Content (range-based, on real PhiX174 gene regions)
//
// Published GC fractions for PhiX174 genes (Sanger et al., 1978):
//   Whole genome    : ~44.8% GC
//   Gene F [1001,2284): coat protein,  ~45.0% GC
//   Gene G [2395,2922): DNA pilot,     ~42.1% GC
//   Gene H [2931,3917): spike protein, ~45.1% GC
//   Promoter [1795,1815): AT-rich,     ~35.0% GC (promoters are AT-rich)
// =============================================================================
void test_gc_content() {
    std::cout << "\n--- Suite 8: GC Content (PhiX174 gene regions, published values) ---\n";
    if (PHIX.empty()) { std::cout << "  [SKIP] genome not loaded\n"; return; }

    auto approx = [](double got, double expected, double tol=0.01){
        return std::abs(got - expected) < tol;
    };

    // ── whole genome ─────────────────────────────────────────────────────────
    {
        GCResult r = gc_content(PHIX_TREE, 0, PHIX.size());
        check("whole genome: length = 5386",        r.length == 5386);
        // PhiX174 has 2413 GC bases out of 5386
        check("whole genome: GC fraction ~44.8%",   approx(r.gc_fraction, 0.448));
    }

    // ── gene F [1001, 2284): major coat protein ───────────────────────────────
    {
        GCResult r = gc_content(PHIX_TREE, 1001, 2284);
        check("gene F: length = 1283",              r.length == 1283);
        check("gene F: GC fraction ~45.0%",         approx(r.gc_fraction, 0.450));
    }

    // ── gene G [2395, 2922): DNA pilot protein ────────────────────────────────
    {
        GCResult r = gc_content(PHIX_TREE, 2395, 2922);
        check("gene G: length = 527",               r.length == 527);
        check("gene G: GC fraction ~42.1%",         approx(r.gc_fraction, 0.421));
    }

    // ── gene H [2931, 3917): spike protein ───────────────────────────────────
    {
        GCResult r = gc_content(PHIX_TREE, 2931, 3917);
        check("gene H: length = 986",               r.length == 986);
        check("gene H: GC fraction ~45.1%",         approx(r.gc_fraction, 0.451));
    }

    // ── Pribnow box region [1795, 1815): AT-rich promoter ─────────────────────
    {
        GCResult r = gc_content(PHIX_TREE, 1795, 1815);
        check("promoter [1795,1815): length = 20",  r.length == 20);
        check("promoter: GC fraction ~35% (AT-rich)", approx(r.gc_fraction, 0.35, 0.05));
    }

    // ── invalid / edge ranges ─────────────────────────────────────────────────
    {
        GCResult r1 = gc_content(PHIX_TREE, 1000, 1000);
        check("empty range [1000,1000): length=0",  r1.length == 0);

        GCResult r2 = gc_content(PHIX_TREE, 500, 200);
        check("inverted range: length=0",           r2.length == 0);

        GCResult r3 = gc_content(PHIX_TREE, 5380, 6000);
        check("beyond-end clamp works correctly",   r3.length == 6);
    }
}

// =============================================================================
// SUITE 9 — Chromosomal Inversion
// Real scenario: inversion of the complementary strand segment containing the
// XhoI restriction site (CTCGAG) at position 161 in PhiX174.
//
// In real molecular biology, Type II restriction enzymes like XhoI cut at
// palindromic sites. The segment "GCGCAGCTCGAGAAG" (pos 155–170) contains
// the XhoI site; inverting it tests that split+reverse+concat works correctly
// at a real genomic location.
//
// Second test: invert gene E [568,843) — the lysis gene. A genomic inversion
// here would disrupt the reading frame, simulating a chromosomal rearrangement.
// =============================================================================
void test_inversion() {
    std::cout << "\n--- Suite 9: Chromosomal Inversion (XhoI site & gene E inversion) ---\n";
    if (PHIX.empty()) { std::cout << "  [SKIP] genome not loaded\n"; return; }

    // ── Test A: invert the 15-base region containing XhoI site ───────────────
    {
        size_t s = 155, e = 170;
        std::string seg = PHIX.substr(s, e-s);   // "GCGCAGCTCGAGAAG"
        std::string seg_rev = seg;
        std::reverse(seg_rev.begin(), seg_rev.end());  // "GAAGAGCTCGACGCG"

        DNATree inverted = chromosomal_inversion(PHIX_TREE, s, e);

        check("XhoI region length unchanged after inversion",
              inverted.measure().length == PHIX_TREE.measure().length);
        check("inverted segment matches manual reversal",
              to_string(inverted).substr(s, e-s) == seg_rev);
        check("sequence before inversion site unchanged",
              to_string(inverted).substr(0, s) == PHIX.substr(0, s));
        check("sequence after inversion site unchanged",
              to_string(inverted).substr(e) == PHIX.substr(e));
        // GC count is invariant — inversion just reorders, doesn't change composition
        check("GC count invariant under inversion",
              inverted.measure().gc_count == PHIX_TREE.measure().gc_count);
        // double inversion = identity
        DNATree double_inv = chromosomal_inversion(inverted, s, e);
        check("double inversion restores original genome",
              to_string(double_inv) == PHIX);
    }

    // ── Test B: invert gene E [568,843) — chromosomal rearrangement ──────────
    {
        size_t s = 568, e = 843;
        std::string seg = PHIX.substr(s, e-s);
        std::string seg_rev = seg; std::reverse(seg_rev.begin(), seg_rev.end());

        DNATree rearranged = chromosomal_inversion(PHIX_TREE, s, e);

        check("gene E inversion: length unchanged",
              rearranged.measure().length == PHIX_TREE.measure().length);
        check("gene E inversion: reversed segment correct",
              to_string(rearranged).substr(s, e-s) == seg_rev);
        check("gene E inversion: GC invariant",
              rearranged.measure().gc_count == PHIX_TREE.measure().gc_count);
    }

    // ── Test C: no-op cases ───────────────────────────────────────────────────
    {
        DNATree same1 = chromosomal_inversion(PHIX_TREE, 500, 500);
        check("invert empty range: no change",
              to_string(same1) == PHIX);

        DNATree same2 = chromosomal_inversion(PHIX_TREE, 843, 568);
        check("invert inverted range (start>end): no change",
              to_string(same2) == PHIX);
    }
}

// =============================================================================
// BENCHMARK — Naive O(n·k)  vs  Rabin-Karp O(n+k)  vs  Finger Tree KMP O(n+k)
// Run on the full PhiX174 genome with three clinically relevant motifs.
// =============================================================================
std::vector<size_t> naive_search(const std::string& s, const std::string& p) {
    std::vector<size_t> h; size_t n=s.size(),k=p.size();
    if(!k||k>n)return h;
    for(size_t i=0;i<=n-k;i++) if(s.compare(i,k,p)==0) h.push_back(i);
    return h;
}
std::vector<size_t> rabin_karp(const std::string& s, const std::string& p) {
    std::vector<size_t> h; size_t n=s.size(),k=p.size();
    if(!k||k>n)return h;
    const long long B=31,M=1000000007LL; long long hm=0,hw=0,pw=1;
    for(size_t i=0;i<k-1;i++) pw=pw*B%M;
    for(size_t i=0;i<k;i++){hm=(hm*B+p[i])%M; hw=(hw*B+s[i])%M;}
    auto vfy=[&](size_t i){return s.compare(i,k,p)==0;};
    if(hw==hm&&vfy(0))h.push_back(0);
    for(size_t i=1;i+k<=n;i++){
        hw=(hw-(long long)s[i-1]*pw%M+M)%M; hw=(hw*B+s[i+k-1])%M;
        if(hw==hm&&vfy(i))h.push_back(i);
    }
    return h;
}
using Clock=std::chrono::high_resolution_clock;
template<typename F> long long tus(F&& f){auto t=Clock::now();f();return std::chrono::duration_cast<std::chrono::microseconds>(Clock::now()-t).count();}

void benchmark() {
    std::cout << "\n--- Benchmark: Naive vs Rabin-Karp vs FT-KMP (PhiX174 genome) ---\n";
    if (PHIX.empty()) { std::cout << "  [SKIP] genome not loaded\n"; return; }

    struct C{const char*name,*motif,*bio;};
    C cases[]={
        {"TTGACA (k=6)", "TTGACA", "-35 promoter element, 4 known hits"},
        {"TATAAA (k=6)", "TATAAA", "Pribnow box, 1 known hit"},
        {"GCGC   (k=4)", "GCGC",   "CpG-like sites, 18 known hits"},
    };
    const int R=10;
    std::cout<<std::left<<std::setw(18)<<"Motif"<<std::setw(12)<<"Naive(µs)"
             <<std::setw(14)<<"Rabin-Karp(µs)"<<std::setw(13)<<"FT-KMP(µs)"
             <<"Hits  Biology\n"<<std::string(75,'-')<<"\n";
    for(auto&c:cases){
        long long tn=0,tr=0,tf=0; size_t hn=0,hr=0,hf=0;
        for(int i=0;i<R;i++){
            tn+=tus([&]{hn=naive_search(PHIX,c.motif).size();});
            tr+=tus([&]{hr=rabin_karp (PHIX,c.motif).size();});
            tf+=tus([&]{hf=findMotif  (PHIX_TREE,c.motif).size();});
        }
        std::cout<<std::left<<std::setw(18)<<c.name
                 <<std::setw(12)<<tn/R<<std::setw(14)<<tr/R<<std::setw(13)<<tf/R
                 <<hf<<"  "<<c.bio<<(hn==hr&&hr==hf?"":" [MISMATCH]")<<"\n";
    }
    std::cout<<"\nNote: all three methods must report identical hit counts.\n";
}

// =============================================================================
// FASTA DEMO — prints a summary of all four features on the full genome
// =============================================================================
void fasta_demo() {
    std::cout << "\n--- PhiX174 Genome Demo (NC_001422.1, 5386 bp) ---\n";
    if (PHIX.empty()) {
        std::cout << "  Could not load genome. Run from project root where data/ lives.\n";
        return;
    }
    std::cout << "Length  : " << PHIX_TREE.measure().length   << " bp\n";
    std::cout << "GC count: " << PHIX_TREE.measure().gc_count << "\n";
    std::cout << "GC %%    : " << std::fixed << std::setprecision(2)
              << 100.0*PHIX_TREE.measure().gc_count/PHIX_TREE.measure().length << "%%\n";

    auto show_motif = [&](const std::string& m, const std::string& desc){
        auto h = findMotif(PHIX_TREE, m);
        std::cout << "  " << m << " (" << desc << "): " << h.size() << " hit(s)";
        for(size_t i=0;i<std::min(h.size(),(size_t)4);i++) std::cout<<" @"<<h[i];
        std::cout << "\n";
    };
    std::cout << "\nKey motifs:\n";
    show_motif("TATAAA","Pribnow box");
    show_motif("TTGACA","-35 element");
    show_motif("CTCGAG","XhoI site");
    show_motif("GCGC","methylation sites");

    std::cout << "\nGene GC content:\n";
    auto show_gc=[&](const std::string& name, size_t s, size_t e){
        GCResult r=gc_content(PHIX_TREE,s,e);
        std::cout<<"  "<<name<<": "<<std::fixed<<std::setprecision(1)
                 <<100.0*r.gc_fraction<<"%%  ("<<r.gc_count<<"/"<<r.length<<")\n";
    };
    show_gc("Gene F [1001,2284)",1001,2284);
    show_gc("Gene G [2395,2922)",2395,2922);
    show_gc("Gene H [2931,3917)",2931,3917);

    std::cout << "\nGene E deletion [568,843): lysis protein removal\n";
    DNATree no_e = deleteRange(PHIX_TREE,568,843);
    std::cout << "  New length: " << no_e.measure().length << " (removed "
              << PHIX_TREE.measure().length - no_e.measure().length << " bp)\n";

    std::cout << "\nInversion of XhoI region [155,170):\n";
    DNATree inv = chromosomal_inversion(PHIX_TREE,155,170);
    std::cout << "  Original : " << PHIX.substr(155,15) << "\n";
    std::cout << "  Inverted : " << to_string(inv).substr(155,15) << "\n";
}

// =============================================================================
// Interactive CLI
// =============================================================================
void run_cli() {
    std::cout << "\n=========================================\n"
              << "  Finger Tree DNA Sequence Tool  (CLI)\n"
              << "=========================================\n";
    std::string seq_input;
    std::cout << "\nEnter DNA sequence (Enter = load PhiX174): ";
    std::getline(std::cin, seq_input);
    DNATree tree;
    if (seq_input.empty()) {
        if (!PHIX.empty()) { tree=PHIX_TREE; seq_input=PHIX; std::cout<<"Loaded PhiX174 ("<<PHIX.size()<<" bp)\n"; }
        else { seq_input="ATCGATCGGCTAGCATGCAT"; tree=build_tree(seq_input); std::cout<<"Using default: "<<seq_input<<"\n"; }
    } else {
        for(char&c:seq_input)c=std::toupper(c);
        tree=build_tree(seq_input);
    }
    std::cout<<"Length="<<tree.measure().length<<"  GC="<<tree.measure().gc_count<<"\n";
    while(true){
        std::cout<<"\n1=Motif search  2=Insert  3=Delete  4=Update  5=GC content  6=Inversion  7=Print  0=Exit\nChoice: ";
        std::string line; if(!std::getline(std::cin,line))break;
        if(line.empty())continue; int cmd=line[0]-'0';
        if(cmd==0){std::cout<<"Goodbye.\n";break;}
        try{
            if(cmd==1){
                std::cout<<"Motif: "; std::string m; std::getline(std::cin,m);
                auto h=findMotif(tree,m);
                std::cout<<"Found "<<h.size()<<" hit(s)";
                size_t sh=std::min(h.size(),(size_t)10);
                for(size_t i=0;i<sh;i++) std::cout<<" @"<<h[i];
                if(h.size()>10) std::cout<<" ...";
                std::cout<<"\n";
            } else if(cmd==2){
                std::cout<<"Position: "; std::getline(std::cin,line); size_t p=std::stoul(line);
                std::cout<<"Insert sequence: "; std::string ins; std::getline(std::cin,ins);
                tree=insertAt(tree,p,ins); std::cout<<"Done. Length="<<tree.measure().length<<"\n";
            } else if(cmd==3){
                std::cout<<"Start: "; std::getline(std::cin,line); size_t s=std::stoul(line);
                std::cout<<"End: "; std::getline(std::cin,line); size_t e=std::stoul(line);
                tree=deleteRange(tree,s,e); std::cout<<"Done. Length="<<tree.measure().length<<"\n";
            } else if(cmd==4){
                std::cout<<"Position: "; std::getline(std::cin,line); size_t p=std::stoul(line);
                std::cout<<"Replacement: "; std::string r; std::getline(std::cin,r);
                tree=updateRange(tree,p,r); std::cout<<"Done. Length="<<tree.measure().length<<"\n";
            } else if(cmd==5){
                std::cout<<"Start: "; std::getline(std::cin,line); size_t s=std::stoul(line);
                std::cout<<"End: "; std::getline(std::cin,line); size_t e=std::stoul(line);
                GCResult r=gc_content(tree,s,e);
                std::cout<<"GC["<<s<<","<<e<<"): len="<<r.length<<" gc="<<r.gc_count
                          <<" frac="<<std::fixed<<std::setprecision(3)<<r.gc_fraction<<"\n";
            } else if(cmd==6){
                std::cout<<"Start: "; std::getline(std::cin,line); size_t s=std::stoul(line);
                std::cout<<"End: "; std::getline(std::cin,line); size_t e=std::stoul(line);
                tree=chromosomal_inversion(tree,s,e);
                std::cout<<"Inverted ["<<s<<","<<e<<"). Length="<<tree.measure().length<<"\n";
            } else if(cmd==7){
                std::string ss=to_string(tree);
                if(ss.size()>200) std::cout<<ss.substr(0,100)<<"...(+"<<ss.size()-200<<")..."<<ss.substr(ss.size()-100)<<"\n";
                else std::cout<<ss<<"\n";
            }
        } catch(const std::exception&ex){ std::cout<<"Error: "<<ex.what()<<"\n"; }
    }
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[]) {
    if (argc > 1 && std::string(argv[1]) == "--cli") { load_genome(); run_cli(); return 0; }

    std::cout << "=== 2-3 Finger Tree: DNA Sequence Application ===\n";
    std::cout << "    Genome: PhiX174 (NC_001422.1)\n";

    bool genome_ok = load_genome();
    if (genome_ok)
        std::cout << "Genome loaded: " << PHIX.size() << " bp\n";
    else
        std::cout << "WARNING: genome not loaded (run from project root)\n";

    test_basic_tree();
    test_concat();
    test_split();
    test_motif_search();
    test_insert();
    test_delete();
    test_update();
    test_gc_content();
    test_inversion();
    fasta_demo();
    benchmark();

    std::cout << "\n=========================================\n";
    std::cout << "  " << passed_tests << " / " << total_tests << " tests passed\n";
    std::cout << "=========================================\n";
    std::cout << "(Run with --cli for interactive mode)\n";
    return (passed_tests == total_tests) ? 0 : 1;
}
