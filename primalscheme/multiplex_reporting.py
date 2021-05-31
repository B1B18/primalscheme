import logging
import os
import pickle

import pandas as pd

from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation, SeqFeature
from reportlab.lib import colors

from .SMARTplex import SMARTplex
from .multiplex import MultiplexScheme
from . import settings

logger = logging.getLogger('Primal Log')


class MultiplexReporter(MultiplexScheme):
    """Reporting methods to extend MultiplexScheme"""

    def write_bed(self, path='./'):
        logger.info('Writing BED')
        filepath = os.path.join(path, '{}.scheme.bed'.format(self.prefix))
        with open(filepath, 'w') as bedhandle:
            for r in self.segment_regions:
                print(*map(str,
                           [self.primary_reference(r.segment_id).id, r.top_pair.left.start,
                            r.top_pair.left.end, r.top_pair.left.name, r.pool]
                           ), sep='\t', file=bedhandle)
                print(*map(str,
                           [self.primary_reference(r.segment_id).id, r.top_pair.right.end,
                            r.top_pair.right.start, r.top_pair.right.name,
                            r.pool]),
                      sep='\t', file=bedhandle)

    def write_tsv(self, path='./'):
        logger.info('Writing TSV')
        filepath = os.path.join(path, '{}.tsv'.format(self.prefix))
        with open(filepath, 'w') as tsvhandle:
            print(*['name', 'pool', 'seq', 'length', '%gc', 'tm (use 65)'],
                  sep='\t', file=tsvhandle)
            for r in self.segment_regions:
                left = r.top_pair.left
                right = r.top_pair.right
                print(*map(str,
                           [left.name, r.pool, left.seq, left.length,
                            '%.2f' % left.gc, '%.2f' % left.tm]),
                      sep='\t', file=tsvhandle)
                print(*map(str,
                           [right.name, r.pool, right.seq, right.length,
                            '%.2f' % right.gc, '%.2f' % right.tm]),
                      sep='\t', file=tsvhandle)
                if r.alternates:
                    for alt in r.alternates:
                        print(*map(str,
                                   [alt.name, r.pool, alt.seq, alt.length,
                                    '%.2f' % alt.gc, '%.2f' % alt.tm]),
                              sep='\t', file=tsvhandle)

    def write_oPools(self, path='./'):
        logger.info('Writing oPools')
        opool_oligos = []
        for bc in self.good_barcodes:
            filepath = os.path.join(path, '{}.{}.RToligos.oPools.tsv'.format(
                    self.prefix, bc))
            with open(filepath, 'w') as tsvhandle:
                print(*['name', 'seq'], sep='\t', file=tsvhandle)
                for r in self.segment_regions:
                    right = r.top_pair.right
                    r_seq = r.right_tail + 'NNNN' + bc + right.seq
                    print(*map(str, [right.name, r_seq]),
                          sep='\t', file=tsvhandle)
                    opool_oligos.append(['RT_pool_'+bc,r_seq])
                    # ensure both endswitch region primers go into RT
                    if r.is_endswitch_region:
                        left = r.top_pair.left
                        l_seq = r.left_tail + 'NNNN' + bc + left.seq
                        print(*map(str, [left.name, l_seq]),
                              sep='\t', file=tsvhandle)
                        opool_oligos.append(['RT_pool_'+bc,l_seq])


        filepath = os.path.join(path, '{}.2ndstrand_1.oPools.tsv'.format(
                self.prefix))
        with open(filepath, 'w') as tsvhandle:
            print(*['name', 'seq'], sep='\t', file=tsvhandle)
            for r in self.segment_regions:
                if r.region_num % 2 == 0: continue
                if r.is_endswitch_region: continue
                left = r.top_pair.left
                l_seq = r.left_tail + left.seq
                print(*map(str, [left.name, l_seq]),
                      sep='\t', file=tsvhandle)
                opool_oligos.append(['2ndstrand_1_pool',l_seq])

        filepath = os.path.join(path, '{}.2ndstrand_2.oPools.tsv'.format(
                self.prefix))
        with open(filepath, 'w') as tsvhandle:
            print(*['name', 'seq'], sep='\t', file=tsvhandle)
            for r in self.segment_regions:
                if r.region_num % 2 == 1: continue
                if r.is_endswitch_region: continue
                left = r.top_pair.left
                l_seq = r.left_tail + left.seq
                print(*map(str, [left.name, l_seq]),
                      sep='\t', file=tsvhandle)
                opool_oligos.append(['2ndstrand_2_pool',l_seq])

        filepath = os.path.join(path, '{}.opools_for_idt.xlsx'.format(
                self.prefix))
        opool_df = pd.DataFrame(opool_oligos, columns=['Pool name','Sequence'])
        opool_df.to_excel(filepath, sheet_name='Sheet1', index=False)

        filepath = os.path.join(path, '{}.tsoligo.oPools.tsv'.format(
                self.prefix))
        with open(filepath, 'w') as tsvhandle:
            print(*['name', 'seq'], sep='\t', file=tsvhandle)
            print('ts_oligo_P5\t/5Me-isodC//iisodG//iMe-isodC/'+
                    settings.TAIL_P5+'rGrGrG', file=tsvhandle)

        filepath = os.path.join(path, '{}.enrichment.oPools.tsv'.format(
                self.prefix))
        # columns for IDT bulk entry are name, sequence, scale, purification
        with open(filepath, 'w') as tsvhandle:
            # write the N7 enrichment barcode oligos
            for bc in settings.ENRICHMENT_BARCODES_N7:
                oligo_name = self.prefix+"_N7_"+bc
                n7_oligo = settings.ENRICHMENT_F7_STUB + \
                    settings.ENRICHMENT_BARCODES_N7[bc] + \
                    settings.ENRICHMENT_N7_STUB
                print(oligo_name+'\t'+n7_oligo+'\t4nmU\tSTD', file=tsvhandle)
            # and the P5 enrichment barcode oligos
            for bc in settings.ENRICHMENT_BARCODES_P5:
                oligo_name = self.prefix+"_P5_"+bc
                p5_oligo = settings.ENRICHMENT_F5_STUB + \
                    settings.ENRICHMENT_BARCODES_P5[bc] + \
                    settings.ENRICHMENT_P5_STUB
                print(oligo_name+'\t'+p5_oligo+'\t4nmU\tSTD', file=tsvhandle)

    def write_SMARTplex(self, path='./'):
        logger.info('Writing SMARTplex')
        filepath = os.path.join(path, '{}_SMARTplex.tsv'.format(self.prefix))
        with open(filepath, 'w') as tsvhandle:
            print(*['name', 'fullseq', 'tm (use 52)', 'seq', 'subseq',
                    'lensubseq', 'lenmatch'], sep='\t', file=tsvhandle)
            for r in self.segment_regions:
                right = r.top_pair.right
                name = '{}_{}_SMARTplex'.format(self.prefix, r.region_num)
                RTprimer, thermo, subseq, lensubseq, lenmatch = SMARTplex(
                    right)
                print(*map(str,
                           [name, RTprimer, '%.2f' % thermo, right.seq, subseq,
                            lensubseq, lenmatch]), sep='\t', file=tsvhandle)

    def write_pickle(self, path='./'):
        logger.info('Writing pickles')
        filepath = os.path.join(path, '{}.pickle'.format(self.prefix))
        with open(filepath, 'wb') as pickleobj:
            pickle.dump(self.segment_regions, pickleobj)

    def write_refs(self, path='./'):
        logger.info('Writing references')
        filepath = os.path.join(path, '{}.reference.fasta'.format(self.prefix))
        with open(filepath, 'w'):
            for refset in self.references:
                SeqIO.write(refset, filepath, 'fasta')

    def apply_to_window(self, sequence, window_size, function, step=None):
        """
        Modified from
        https://github.com/biopython/biopython/blob/master/Tests/test_GenomeDiagram.py
        Apply function to windows of the given sequence.
        Returns a list of (position, value) tuples for fragments of the passed
        sequence of length window_size (stepped by step), calculated by the
        passed function.  Returned positions are the midpoint of each window.
        Arguments:
        - sequence - Bio.Seq.Seq object.
        - window_size - an integer describing the length of sequence
          to consider.
        - step - an integer describing the step to take between windows
          (default = window_size//2).
        - function - Method or function that accepts a Bio.Seq.Seq object
          as its sole argument and returns a single value.
        apply_to_window(sequence, window_size, function) ->
            [(int, float),(int, float),...]
        """
        seqlen = len(sequence)  # Total length of sequence to be used
        if step is None:  # Use half window-width or 1 if larger
            step = max(window_size // 2, 1)
        else:  # Use specified step, or 1 if greater
            step = max(step, 1)

        results = []  # Holds (position, value) results

        # Perform the passed function on as many windows as possible, short of
        # overrunning the sequence
        pos = 0
        while pos < seqlen - window_size + 1:
            # Obtain sequence fragment
            start = pos
            middle = (pos + window_size + pos) // 2
            end = pos + window_size
            fragment = sequence[start:end]
            # Apply function to the sequence fragment
            value = function(fragment)
            results.append((middle, value))  # Add results to list
            # Advance to next fragment
            pos += step

        # Use the last available window on the sequence, even if it means
        # re-covering old ground
        if pos != seqlen - window_size:
            # Obtain sequence fragment
            pos = seqlen - window_size
            start = pos
            middle = (pos + window_size + pos) // 2
            end = pos + window_size
            fragment = sequence[start:end]
            # Apply function to sequence fragment
            value = function(fragment)
            results.append((middle, value))  # Add results to list

        return results      # Return the list of (position, value) results

    def calc_gc_skew(self, sequence):
        """Return the (G-C)/(G+C) GC skew in a passed sequence.
        Arguments:
            - sequence   - a Bio.Seq.Seq object.
        calc_gc_skew(sequence)
        """
        g = sequence.count('G') + sequence.count('g')
        c = sequence.count('C') + sequence.count('c')
        if g + c == 0:
            return 0.0  # TODO - return NaN or None here?
        else:
            return (g - c) / float(g + c)

    def write_schemadelica_plot(self, path='./'):
        logger.info('Writing plot')

        gd_diagram = GenomeDiagram.Diagram("Primer Scheme", track_size=0.15)
        primer_feature_set = GenomeDiagram.FeatureSet()

        # make the gc track
        window = 50
        gc_set = GenomeDiagram.GraphSet('GC skew')
        graphdata1 = self.apply_to_window(
            self.primary_reference.seq, window, self.calc_gc_skew)
        gc_set.new_graph(graphdata1, 'GC Skew', style='line',
                         color=colors.violet, altcolor=colors.purple)
        gc_track = GenomeDiagram.Track('GC Skew', height=1.5, greytrack=0,
                                       scale_largetick_interval=1e3)
        gc_track.add_set(gc_set)

        # make the primer track
        for r in self.regions:
            region = str(r.region_num)
            strand = 1 if r.region_num % 2 else -1

            fwd_feature = SeqFeature(
                FeatureLocation(r.top_pair.left.start,
                                r.top_pair.left.end, strand=strand))
            rev_feature = SeqFeature(
                FeatureLocation(r.top_pair.right.end,
                                r.top_pair.right.start, strand=strand))
            region_feature = SeqFeature(
                FeatureLocation(r.top_pair.left.start,
                                r.top_pair.right.start, strand=strand))

            primer_color = colors.red
            region_color = colors.palevioletred

            primer_feature_set.add_feature(
                region_feature, color=region_color, name=region, label=True,
                label_position="middle",
                label_angle=0 if strand == 1 else -180)
            primer_feature_set.add_feature(fwd_feature, color=primer_color,
                                           name=region)
            primer_feature_set.add_feature(rev_feature, color=primer_color,
                                           name=region)

        primer_track = GenomeDiagram.Track(name="Annotated Features", height=1)
        primer_track.add_set(primer_feature_set)

        gd_diagram.add_track(primer_track, 2)
        gd_diagram.add_track(gc_track, 1)

        rows = max(2, int(round(len(self.primary_reference) / 10000.0)))
        gd_diagram.draw(format='linear', pagesize=(300 * rows, 200 * rows),
                        fragments=rows, start=0,
                        end=len(self.primary_reference))

        pdf_filepath = os.path.join(path, '{}.pdf'.format(self.prefix))
        svg_filepath = os.path.join(path, '{}.svg'.format(self.prefix))
        gd_diagram.write(pdf_filepath, 'PDF', dpi=300)
        gd_diagram.write(svg_filepath, 'SVG', dpi=300)
