# Extensions developed for the bluetongue virus project

The primalscheme software was modified to support design of an amplicon scheme for bluetongue virus (BTV) genome sequencing, as part of a collaborative project between UTS and NSW DPI under the AusGEM3 program.
The design builds upon concepts that were originally introduced by Monahan et al for [massively parallel diagnosis and genome sequencing of SARS-CoV-2](https://www.protocols.io/view/a-protocol-for-massively-parallel-diagnosis-and-ge-betrjem6), DOI 10.17504/protocols.io.betrjem6 .
This document describes the extensions to the primalscheme software to support the design, with text that could be included in brief or in full in manuscript methods descriptions.

At a high level, the following features were introduced to primalscheme for application to BTV and similar diverse viruses:

* accept multiple fasta input alignments, one per viral genome segment (could also represent different virus genomes)
* design oligos at ends of each linear genome segment to pair with endswitch oligos for full length genomes. requires bed file input
* degenerate primer designs to target diverse virus populations
* exclusion amplification design for single pool reverse transcription - inhibits amplification of short overlaps
* short UMIs for reverse transcription
* sample barcodes for early pooling & deep multiplexing
* output to IDT oPools spreadsheet format

## Design approach

The amplicon design approach leverages IDT's oPools product, whereby a single pool of thousands of defined oligonucleotide sequences can be synthesized. The general structure of oligos is as follows.

Reverse transcription primers:

```
[Nextera 7 tail]-[4nt UMI]-[6nt sample barcode]-[target sequence]

or

[TruSeq 5 tail]-[4nt UMI]-[6nt sample barcode]-[target sequence]
```

Second strand synthesis primers:

```
Pool 1: [TruSeq 5 tail]-[target sequence]
Pool 2: [Nextera 7 tail]-[target sequence]
```

As with the original primalscheme approach, overlapping amplicons are designed to tile the entire genome. Unlike the original approach, the primer tail used for each successive amplicon target alternates between the Nextera 7 and the TruSeq 5 stub. During 2nd strand synthesis this creates products with four topologies: 

1. `[TS5]---[TS5']`
2. `[N7]---[N7']`
3. `[TS5]---[N7']`
4. `[N7]---[TS5']`

where ' denotes the reverse complement sequence. Of those four, topologies 1 and 2 are expected to amplify poorly because the 3' and 5' ends are reverse complement and expected to occlude the priming sites used for subsequent enrichment PCR. This type of exclusion occurs during standard Nextera library prep as well, preventing the amplification of library fragments with the N7 on both ends, or N5 on both ends.

In principle it may be possible to carry out the 2nd strand synthesis in a single reaction rather than two separate reactions, although there is risk of primer dimers forming that contain Nextera 7 + TruSeq 5 adapters that would be active on the flowcell.

Enrichment primer structure:

```
[TruSeq 5]
[Nextera 7]
```

These oligos are used to enrich the product of 2nd strand synthesis and can include sample barcode sequences.


### degenerate target sequences

If a viral genome alignment is provided as input, each candidate target primer is evaluated for similarity across the entire set of aligned sequences. Where mismatches are found between the target primer and one of the aligned sequences, a greedy approach is used to introduce degeneracy to the oligo to maximize the number of sequence identities. Up to 5 target oligo bases are iteratively replaced with N to improve the oligo target sequence match to the aligned sequences. In practice on the bluetongue virus, we observe that this approach frequently results in 3rd codon position bases represented as N, as expected.

## Possible design advantages

There are several possible advantages to the above design: 

1. the use of the alternating primer tails enables overlapping amplicon targets to be reverse-transcribed in a single RT reaction
2. By introducing sample barcodes via the RT primers, it may be possible to anneal the RT primers per-sample, then pool samples for reverse transcription, enabling a large number of samples to be processed in a single reverse transcription reaction, saving on RT costs
3. The UMIs introduced during RT might identify PCR duplicates, helping to distinguish between true within-host variants and error induced by RNA/DNA damage and polymerase error occuring during sample prep
4. The use of oPools enables low cost synthesis of large amplicon sets and supports the use of up to 10 degenerate sites (IUPAC N or K during synthesis). Currently up to 384 oligos can be synthesized at 50nM scale.
5. The approach can be scaled to extremely large sample batches via combinatorial barcoding.


Some potential concerns:

* It is not possible to normalize individual oligos in an oPools oligo set. It is possible to add the same oligo multiple times to a single oPools set, so it may be possible to boost yields of poorly amplifying target regions in this manner, but requires re-ordering the entire oligo set.
* The algorithm used to design degenerate oligos is not phylogeny-aware. Variants will of course be linked to each other via common evolutionary history. It should be possible to design a small number of separate target sequences, which capture fixed differences in different parts of the phylogeny, and in so doing, introduce a smaller number of degeneracies to the design. This would have the advantage of increasing the molarity of oligo that closely matches target sequences. A "manual" approach to achieve this would be to split the alignment into separate alignments by major phylogenetic groups prior to processing with primalscheme. An automated procedure would be preferable, however.


## Example usage

The following is an example invocation to design an amplicon scheme for bluetongue virus:

```
source venv/bin/activate
cd bluetongue
primalscheme multiplex --prefix btv --unaligned-bed bluetongue/all.bed --output-path btv_scheme bluetongue/*.fa
```
