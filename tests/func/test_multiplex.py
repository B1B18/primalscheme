import pytest

from primalscheme.align import align_primer
from primalscheme.cli import process_fasta
from primalscheme.multiplex import MultiplexScheme


@pytest.fixture(scope="session")
def default_chikv_scheme(chikv_input):
    return get_scheme(chikv_input)


@pytest.fixture(scope="session")
def ebola_scheme(ebola_input):
    return get_scheme(ebola_input)


def get_scheme(fasta, **kwargs):
    references = process_fasta(fasta)
    scheme = MultiplexScheme(references, **kwargs)
    scheme.design_scheme()
    return scheme


def no_collisions(regions):
    pool_ids = (1, 2)
    for pool_id in pool_ids:
        pool = [r for r in regions if r.pool == pool_id]
        last_coord = -1
        for r in pool:
            if r.left.start <= last_coord:
                return False
            last_coord = r.right.start
    return True


def test_chikv_scheme_has_no_gaps(default_chikv_scheme):
    regions = default_chikv_scheme.regions
    inserts = [(r.left.end + 1, r.right.end) for r in regions]
    covered_coords = set([x for insert in inserts for x in range(*insert)])
    all_coords = set(range(regions[0].left.end + 1, regions[-1].right.end))

    assert all_coords.issubset(covered_coords)


def test_no_collisions_in_any_test_scheme(all_stored_inputs):
    scheme = get_scheme(all_stored_inputs)
    assert no_collisions(scheme.regions)


@pytest.mark.parametrize("amplicon_size", range(200, 801, 200))
def test_scheme_varying_amplicon_sizes(amplicon_size, chikv_input):
    variation = int(0.1 * amplicon_size / 2)
    amplicon_size_min = amplicon_size - variation
    amplicon_size_max = amplicon_size + variation
    scheme = get_scheme(
        chikv_input,
        amplicon_size_min=amplicon_size_min,
        amplicon_size_max=amplicon_size_max,
    )

    regions = scheme.regions
    for region in regions[:4]:
        product_size = region.product_size
        assert product_size <= amplicon_size_max
        assert product_size >= amplicon_size_min


def test_large_target_overlap_does_not_result_in_collision(chikv_input):
    scheme = get_scheme(chikv_input, target_overlap=300)

    assert no_collisions(scheme.regions)


def test_left_candidates_are_correctly_sorted(default_chikv_scheme):
    region = default_chikv_scheme.regions[0]
    penalty = region.left.combined_penalty
    for candidate in region.left_candidates:
        if candidate.combined_penalty < penalty:
            pytest.fail("Candidates are not sorted by combined penalty")
        penalty = candidate.combined_penalty


def test_right_candidates_are_correctly_sorted(default_chikv_scheme):
    region = default_chikv_scheme.regions[0]
    penalty = region.right.combined_penalty
    for candidate in region.right_candidates:
        if candidate.combined_penalty < penalty:
            pytest.fail("Candidates are not sorted by combined penalty")
        penalty = candidate.combined_penalty


def test_left_primer_seq_matches_some_ref_slice(default_chikv_scheme):
    left = default_chikv_scheme.regions[0].left
    ref_slices = left.reference_msa

    assert any(ref_slice.seq == left.seq for ref_slice in ref_slices)


def test_right_primer_seq_matches_some_ref_slice(default_chikv_scheme):
    right = default_chikv_scheme.regions[0].right
    ref_slices = right.reference_msa

    assert any(ref_slice.seq == right.seq for ref_slice in ref_slices)


def test_scheme_with_single_reference(chikv_input):
    references = process_fasta(chikv_input)[:1]
    scheme = MultiplexScheme(references)
    scheme.design_scheme()

    assert len(scheme.regions) > 0


def test_inferred_mismatch_alignments(ebola_scheme):
    scheme = ebola_scheme
    primers = scheme.primers
    for primer in primers:
        mismatch_counts = [align_primer(primer, ref)[0] for ref in scheme.references]
        print(primer.reference_msa)
        assert primer.mismatch_counts == mismatch_counts
