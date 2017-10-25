# Blocked Pattern Matching with Priority Search Trees

## Compile
> make

## Usage

### Unit tests
> ./UnitTest

### (Exact) Blocked Pattern Matching
We first create the index data structure for our database:
> ./CreateIndex  sample/sample.fasta

The file *sample/patterns.txt* contains three patterns. The first pattern contains the mass of the string ASV, the second pattern additionally the mass of AN, and the last pattern additionally the mass of GLP.
The algorithm reads the patterns from stdin and outputs the matching substrings.

> ./BPM sample/sample.fasta.db < sample/patterns.txt

### Modification-Tolerant Blocked Pattern Matching
We create the index:
> ./CreateIndex_Mod sample/sample.fasta.db cfg/modifications.cfg
The file *cfg/modifications.cfg* specifies two modifications, namely +16 Da at M (all sites) and +1 Da at A (n-terminal, only allowed at the first character).

*sample/patterns-mod.txt* contains two pattern. The first pattern has the masses of A(+1)VS,AN, and PGL. The second pattern has the masses of M(+16)FS,SV,M(+16)V
The algorithm outputs all substrings T with a modified string T' matching one of the patterns.

> ./BPM_Mod sample/sample.fasta.db cfg/modifications.cfg < sample/patterns-mod.txt

### Mutation-Tolerant Blocked Pattern Matching
We create the index for our database:
> ./CreateIndex_Mut  sample/sample.fasta

*sample/patterns-mut.txt* contains a single pattern with the masses of ASV,A, and GLP.
The algorithm outputs all substrings T with a mutated string T' matching the pattern.

> ./BPM_Mut sample/sample.fasta.db < sample/patterns-mut.txt

1. AVSANPGL - deletion of N
2. VSANPGL  - substitution of N by A
3. AVSANPG  - substitution of L by N
