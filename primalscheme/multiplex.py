import logging
import primer3

from . import settings
from .exceptions import NoSuitablePrimers
from .models import CandidatePrimer, CandidatePrimerPair, Region

logger = logging.getLogger('Primal Log')


class MultiplexScheme(object):
    """A complete multiplex primer scheme."""

    def __init__(self, references, amplicon_length, min_overlap, max_gap,
                 max_alts, max_candidates, step_size, max_variation,
                 prefix='PRIMAL_SCHEME', flanks=None):
        self.references = references
        self.amplicon_length = amplicon_length
        self.min_overlap = min_overlap
        self.max_gap = max_gap
        self.max_alts = max_alts
        self.max_candidates = max_candidates
        self.step_size = step_size
        self.max_variation = max_variation
        self.prefix = prefix
        self.segment_regions = []
        self.flanks = flanks

        self.run()

    def primary_reference(self, segment_id):
        return self.references[segment_id][0]

    def run(self):
        for s in range(len(self.references)):
            self.segment_regions.extend(self.run_segment(s))


        regions = self.segment_regions
        logger.info("Checking for hairpins in "+str(len(regions))+" regions")
        # check for hairpins with tails
        for region in regions:
            new_candidates = []
            best_candidate = region.top_pair
            best_candidate_tm_ssq = 10000000


            for c in region.candidate_pairs:
                l_seq = region.left_tail + c.left.seq
                r_seq = region.right_tail + c.right.seq
                if len(l_seq) > 60: l_seq = l_seq[-60:]
                if len(r_seq) > 60: r_seq = r_seq[-60:]
                left_hp = primer3.bindings.calcHairpin(l_seq)
                right_hp = primer3.bindings.calcHairpin(r_seq)
                hairpin_ok = False
                if left_hp.tm < settings.MAX_HAIRPIN_TM and right_hp.tm < settings.MAX_HAIRPIN_TM:
                    hairpin_ok = True
                if left_hp.tm**2 + right_hp.tm**2 < best_candidate_tm_ssq:
                    best_candidate = c
                    best_candidate_tm_ssq = left_hp.tm**2 + right_hp.tm**2
                # now check for self-dimers
                left_homo = primer3.bindings.calcHomodimer(l_seq)
                right_homo = primer3.bindings.calcHomodimer(r_seq)
                if left_homo.tm < settings.MAX_HOMODIMER_TM and right_homo.tm < settings.MAX_HOMODIMER_TM and hairpin_ok:
                    new_candidates.append(c)
                if left_homo.tm > settings.MAX_HOMODIMER_TM:
                    print(left_homo)
                if right_homo.tm > settings.MAX_HOMODIMER_TM:
                    print(right_homo)
            if len(new_candidates)==0:
                logger.info("Warning, no primers found without hairpin for region "+str(region.region_num))
                left_hp = primer3.bindings.calcHairpin(l_seq)
                right_hp = primer3.bindings.calcHairpin(r_seq)
                logger.info("Best candidate left T_m: " + str(left_hp.tm) + " right T_m: " + str(right_hp.tm))
                new_candidates.append(best_candidate)
            region.candidate_pairs = new_candidates

        logger.info("Checking for heterodimers")
        filt_regions_0 = self._remove_heterodimers(regions,0)
        filt_regions_1 = self._remove_heterodimers(filt_regions_0,1)


        # select some barcodes
        self.good_barcodes = []
        need_barcodes = 8
        bc_index = 0
        while len(self.good_barcodes) < need_barcodes:
            cur_bc = settings.INLINE_BARCODES[bc_index]
            bc_index+=1
            worse = 0
            for region in filt_regions_1:
                right_seq = region.right_tail + 'NNN' + cur_bc + region.top_pair.right.seq
                right_seq_no_bc = region.right_tail + region.top_pair.right.seq
                if len(right_seq) > 60: right_seq = right_seq[-60:]
                if len(right_seq_no_bc) > 60: right_seq_no_bc = right_seq_no_bc[-60:]
                right_hp = primer3.bindings.calcHairpin(right_seq)
                right_hp_no_bc = primer3.bindings.calcHairpin(right_seq_no_bc)
                if right_hp.tm > right_hp_no_bc.tm and right_hp.tm > settings.MAX_HAIRPIN_TM:
                    worse+=1
            if worse < 10:
                self.good_barcodes.append(cur_bc)
            else:
                logger.info("barcode "+cur_bc+" made "+str(worse)+" hairpins worse")

        logger.info(self.good_barcodes)

        # diversify the oligos to better match the entire target set
        for region in regions:
            for c in region.candidate_pairs:
                c.left.diversify()
                c.right.diversify()

    def run_segment(self, segment_id):
        regions = []
        region_num = 0
        is_last_region = False
        is_endswitch_region = False

        while True:
            region_num += 1

            # Get the previous region in each pool
            prev_pair = (regions[-1].candidate_pairs[0]
                         if len(regions) >= 1 else None)
            prev_pair_same_pool = (regions[-2].candidate_pairs[0]
                                   if len(regions) > 2 else None)

            # If there are two regions or more
            if prev_pair_same_pool:
                # Gap opened between -1 and -2 regions
                if prev_pair.left.start > prev_pair_same_pool.right.start:
                    # If gap, left primer cannot overlap with -1 region
                    left_primer_left_limit = prev_pair.left.end + 1
                else:
                    # Left primer cannot overlap -2 region
                    left_primer_left_limit = (prev_pair_same_pool.right.start
                                              + 1)
            # If there is more than one region
            elif prev_pair:
                # Left primer cannot overlap with -1 region or you don't move
                left_primer_left_limit = prev_pair.left.end + 1
            else:
                # Region one only limit is 0
                left_primer_left_limit = 0

            # Right start limit maintains the minimum_overlap
            left_primer_right_limit = (prev_pair.right.end
                                       - self.min_overlap - 1
                                       if prev_pair else self.max_gap)

            # Last region if less than one amplicon length remaining
            if prev_pair:
                if ((len(self.primary_reference(segment_id)) - prev_pair.right.end)
                        < self.amplicon_length):
                    is_last_region = True
                    logger.debug(
                        'Region {}: is last region'.format(region_num))

            if is_endswitch_region:
                is_last_region = False
                logger.info('Segment {} region {} designed for template switching'.format(segment_id,region_num))
                left_primer_left_limit = regions[0].candidate_pairs[0].left.end

            # Log limits
            logger.debug('Region {}: forward primer limits {}:{}'.format(
                region_num, left_primer_left_limit, left_primer_right_limit))

            # Find primers or handle no suitable error
            try:
                region = self._find_primers(region_num, left_primer_left_limit,
                                            left_primer_right_limit,
                                            is_last_region, is_endswitch_region,
                                            segment_id)
                regions.append(region)
            except NoSuitablePrimers:
                logger.debug('Region {}: no suitable primer error'.format(
                    region_num))
                break

            # Handle the end
            if is_endswitch_region:
                logger.debug('Region {}: ending normally'.format(region_num))
                break

            if is_last_region:
                is_endswitch_region = True
                is_last_region = False

            # Report scores and alignments
            for i in range(0, len(self.references[segment_id])):
                # Don't display alignment to reference
                logger.debug(regions[-1].candidate_pairs[0].left.alignments[i]
                             .formatted_alignment)
            logger.debug('Identities for sorted left candidates: ' + ','.join(
                ['%.2f' % each.left.percent_identity
                 for each in regions[-1].candidate_pairs]))
            logger.debug('Left start for sorted candidates: ' + ','.join(
                ['%i' % each.left.start
                 for each in regions[-1].candidate_pairs]))
            logger.debug('Left end for sorted candidates: '
                         + ','.join(['%i' % each.left.end
                                     for each in regions[-1].candidate_pairs]))
            logger.debug('Left length for sorted candidates: ' + ','.join(
                ['%i' % each.left.length
                 for each in regions[-1].candidate_pairs]))

            for i in range(0, len(self.references[segment_id])):
                logger.debug(regions[-1].candidate_pairs[0].right.alignments[i]
                             .formatted_alignment)
            logger.debug('Identities for sorted right candidates: ' + ','.join(
                ['%.2f' % each.right.percent_identity
                 for each in regions[-1].candidate_pairs]))
            logger.debug('Right start for sorted candidates: ' + ','.join(
                ['%i' % each.right.start
                 for each in regions[-1].candidate_pairs]))
            logger.debug('Right end for sorted candidates: ' + ','.join(
                ['%i' % each.right.end
                 for each in regions[-1].candidate_pairs]))
            logger.debug('Right length for sorted candidates: ' + ','.join(
                ['%i' % each.right.length
                 for each in regions[-1].candidate_pairs]))
            logger.debug('Totals for sorted pairs: ' + ','.join(
                ['%.2f' % each.mean_percent_identity
                 for each in regions[-1].candidate_pairs]))

            if len(regions) > 1:
                # Remember, results now include this one, so -2
                # is the other pool
                trimmed_overlap = (regions[-2].candidate_pairs[0].right.end
                                   - regions[-1].candidate_pairs[0].left.end
                                   - 1)
                logger.info(
                    "Segment %i Region %i: highest scoring product %i:%i, length %i, "
                    "trimmed overlap %i" % ( segment_id,
                        region_num, regions[-1].candidate_pairs[0].left.start,
                        regions[-1].candidate_pairs[0].right.start,
                        regions[-1].candidate_pairs[0].product_length,
                        trimmed_overlap))
            else:
                logger.info(
                    "Segment %i Region %i: highest scoring product %i:%i, length %i" % (
                        segment_id,
                        region_num, regions[-1].candidate_pairs[0].left.start,
                        regions[-1].candidate_pairs[0].right.start,
                        regions[-1].candidate_pairs[0].product_length))

        # Return regions
        return regions


    def _remove_heterodimers(self,regions,pool):
        """
        Removes heterodimers from a pool. Valid pool IDs are 0 or 1.

        Returns a filtered list of Region objects.
        """
        # check for heterodimers...
        h_network = self._calc_heterodimers(regions,pool)
        while len(h_network)>0:
            # find the most troublesome primer and remove it
            max_heterodimers = 0
            max_pair = -1
            for h in h_network:
                if len(h_network[h]) > max_heterodimers:
                    max_heterodimers = len(h_network[h])
                    max_pair = h
            # take the next pair in list
            if(len(regions[max_pair].candidate_pairs)>0):
                regions[max_pair].candidate_pairs.pop(0)
                logger.info("Swapped in alt primer pair for "+str(max_pair)+" with "+str(max_heterodimers)+" heterodimer interactions")
            else:
                logger.info("Warning: no alt primer pair for "+str(max_pair)+" with "+str(max_heterodimers)+" heterodimer interactions")

            h_network = {}
            h_network = self._calc_heterodimers(regions,pool)

        filtered_regions = []
        for r in regions:
            if len(r.candidate_pairs)>0:
                filtered_regions.append(r)
        return filtered_regions

    def _calc_heterodimers(self,regions,pool):
        """
        Compute heterodimers among top candidate primer pairs.

        Return a hash of region numbers with a list of heterodimer interactions
        that exceed the Tm threshold.
        """
        # first build a heterodimer Tm network
        het_network = {}
        for i in range(len(regions)):
            if len(regions[i].candidate_pairs)==0: continue
            # there are two pools - only find heterodimers within a single pool
            if regions[i].region_num % 2 != pool: continue

            # the left oligos get used together with the enrichment oligos
            left_het_n7 = primer3.bindings.calcHeterodimer(
                regions[i].top_pair.left.seq,
                settings.ENRICHMENT_N7)
            left_het_p5 = primer3.bindings.calcHeterodimer(
                regions[i].top_pair.left.seq,
                settings.ENRICHMENT_P5)
            if left_het_n7.tm > settings.MAX_HETERODIMER_TM or \
                left_het_p5.tm > settings.MAX_HETERODIMER_TM:
                lh_tm = max(left_het_n7.tm,left_het_p5.tm)
                het_network[i].append([lh_tm,-1])

            for j in range(i+1,len(regions)):
                if len(regions[j].candidate_pairs)==0: continue
                if regions[j].region_num % 2 != pool: continue
                # the left oligos get used together in a reaction, and the
                # right oligos get used in a separate reaction. So no need
                # to test left against right
                left_het = primer3.bindings.calcHeterodimer(
                    regions[i].top_pair.left.seq,
                    regions[j].top_pair.left.seq)

                right_het = primer3.bindings.calcHeterodimer(
                    regions[i].top_pair.right.seq,
                    regions[j].top_pair.right.seq)

                # get the max T_m
                max_tm = left_het.tm
                if max_tm < right_het.tm: max_tm = right_het.tm
                if max_tm > settings.MAX_HETERODIMER_TM:
                    if not i in het_network:
                        het_network[i] = []
                    if not j in het_network:
                        het_network[j] = []
                    het_network[i].append([max_tm,j])
                    het_network[j].append([max_tm,i])

        return het_network

    def _find_primers(self, region_num, left_primer_left_limit,
                      left_primer_right_limit, is_last_region,
                      is_endswitch_region, segment_id):
        """
        Find primers for a given region.

        Return a list of Region objects containing candidate
        primer pairs sorted by mean percent identity against all references.
        """

        size_mult = 1  # target size multiplier - used for endswitch
        # Calculate where to slice the reference
        if region_num == 1:
            chunk_start = 0
            chunk_end = (int((1 + self.max_variation / 2)
                             * self.amplicon_length))
        elif is_last_region:
            # Last time work backwards
            chunk_start = (int(len(self.primary_reference(segment_id))
                               - ((1 + self.max_variation / 2)
                               * self.amplicon_length)))
            chunk_end = len(self.primary_reference(segment_id))
        elif is_endswitch_region:
            size_mult = 2  # these are in fact two separate amplicons
            seq_len = len(self.primary_reference(segment_id).seq)
            chunk_start = (int(left_primer_right_limit
                               - (self.max_variation/10 * self.amplicon_length)
                               - settings.global_args['PRIMER_MAX_SIZE']))
            chunk_end = int(left_primer_left_limit
                               + (self.max_variation/10 * self.amplicon_length)
                               + settings.global_args['PRIMER_MAX_SIZE'])
            chunk_start = min(chunk_start, seq_len - self.amplicon_length)
            chunk_end = max(chunk_end, self.amplicon_length)
        else:
            # right limit
            # - min overlap
            # - diff max min product length
            # - max primer length
            chunk_start = (int(left_primer_right_limit
                               - (self.max_variation/10 * self.amplicon_length)
                               - settings.global_args['PRIMER_MAX_SIZE']))
            chunk_end = (int(chunk_start
                             + ((1 + self.max_variation/2)
                                 * self.amplicon_length)))
        initial_chunk_start = chunk_start
        initial_chunk_end = chunk_end

        # Primer3 setup
        p3_global_args = settings.global_args
        p3_seq_args = settings.seq_args
        flank_size = 0
        p3_global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [
            [int(size_mult * self.amplicon_length * (1 - self.max_variation / 10)),
             int(size_mult * self.amplicon_length * (1 + self.max_variation / 2))]]
        p3_global_args['PRIMER_NUM_RETURN'] = self.max_candidates
        logger.info("target amplicon size range:")
        logger.info(p3_global_args['PRIMER_PRODUCT_SIZE_RANGE'][0])

        # Run primer3 until primers are found
        hit_left_limit = False
        while True:
            # Slice primary reference
            if is_endswitch_region: print("chunk_start " + str(chunk_start))
            if is_endswitch_region: print("chunk_end " + str(chunk_end))
            seq = str(self.primary_reference(segment_id).seq[chunk_start:chunk_end])

            if is_endswitch_region:
                seq = str(self.primary_reference(segment_id).seq[chunk_start:])
                if self.flanks is not None:
                    seq += "N" * self.flanks[segment_id][1]
                    seq += "N" * self.flanks[segment_id][0]
                    logger.info("added flank N")
                seq += str(self.primary_reference(segment_id).seq[0:chunk_end])
                if chunk_start < chunk_end:
                    logger.info("WARNING!! overlap in endswitch!!")

            p3_seq_args['SEQUENCE_TEMPLATE'] = seq
            p3_seq_args['SEQUENCE_INCLUDED_REGION'] = [0, len(seq) - 1]
            logger.debug("Region %i: reference chunk %i:%i, length %i" % (
                region_num, chunk_start, chunk_end, len(seq)))
            # only run the design if the target sequence is big enough, otherwise need to step out first
            if len(seq) >= p3_global_args['PRIMER_PRODUCT_SIZE_RANGE'][0][0]:
                primer3_output = primer3.bindings.designPrimers(p3_seq_args,
                                                            p3_global_args)

                candidate_pairs = []

                for cand_num in range(self.max_candidates):
                    lenkey = 'PRIMER_LEFT_%i' % (cand_num)
                    left_name = '%s_%i_%i_%s' % (self.prefix, segment_id,
                                                    region_num, 'LEFT')
                    right_name = '%s_%i_%i_%s' % (self.prefix, segment_id,
                                                    region_num, 'RIGHT')

                    if lenkey not in primer3_output:
                        break

                    left_seq = str(
                        primer3_output['PRIMER_LEFT_%i_SEQUENCE' % (cand_num)])
                    right_seq = str(
                        primer3_output['PRIMER_RIGHT_%i_SEQUENCE' % (cand_num)])

                    left_start = int(
                        primer3_output['PRIMER_LEFT_%i' % (cand_num)][0]
                        + chunk_start)
                    right_start = int(
                        primer3_output['PRIMER_RIGHT_%i' % (cand_num)][0]
                        + chunk_start + 1)

                    left = CandidatePrimer('LEFT', left_name, left_seq, left_start)
                    right = CandidatePrimer('RIGHT', right_name, right_seq,
                                            right_start)

                    candidate_pairs.append(CandidatePrimerPair(left, right))

                set_left = set(pair.left.seq for pair in candidate_pairs)
                set_right = set(pair.right.seq for pair in candidate_pairs)

                logger.info(
                    "Region %i: current position returned %i left and %i "
                    "right unique" % (region_num, len(set_left), len(set_right)))

                if len(set_left) > 2 and len(set_right) > 2:
                    return Region(region_num, chunk_start, candidate_pairs, segment_id,
                                  self.references, self.prefix, self.max_alts,
                                  is_endswitch_region)
            else:
                logger.info("target amplicon size range does not include sequence length:")

            # Move right if first region or to open gap
            if region_num == 1 or hit_left_limit:
                logger.debug("Region %i: stepping right, position %i"
                             % (region_num, chunk_start))
                chunk_start += self.step_size
                chunk_end += self.step_size
                # Hit end of reference
                if chunk_end > len(self.primary_reference(segment_id)):
                    logger.debug("Region %i: hit right limit %i"
                                 % (region_num, len(self.primary_reference(segment_id))))
                    raise NoSuitablePrimers("No suitable primers in region")
            elif is_endswitch_region:
                # increase distance to end of target seq
                logger.debug(
                    "Region %i: stepping inward, position %i, limit %s"
                    % (region_num, chunk_start, left_primer_left_limit))
                chunk_start -= self.step_size
                chunk_end += self.step_size
            else:
                # Move left for all other regions
                logger.debug(
                    "Region %i: stepping left, position %i, limit %s"
                    % (region_num, chunk_start, left_primer_left_limit))
                chunk_start -= self.step_size
                chunk_end -= self.step_size
                if chunk_start <= left_primer_left_limit:
                    # Switch direction to open gap
                    logger.debug("Region %i: hit left limit" % (region_num))
                    chunk_start = initial_chunk_start
                    chunk_end = initial_chunk_end
                    hit_left_limit = True
