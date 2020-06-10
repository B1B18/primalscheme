import primer3

from primalscheme import config
from primalscheme.primer import Primer, Direction


def test_base_penalty_matches_primer3_penalty(random_reference_slice):
    """Check internal penalty calculation matches primer3 output"""

    p3_global = {
        "PRIMER_NUM_RETURN": 10,
        "PRIMER_PRODUCT_SIZE_RANGE": [
            [config.AMPLICON_SIZE_MIN, config.AMPLICON_SIZE_MIN]
        ],
        "PRIMER_OPT_SIZE": config.PRIMER_OPT_SIZE,
        "PRIMER_MIN_SIZE": config.PRIMER_SIZE_MIN,
        "PRIMER_MAX_SIZE": config.PRIMER_SIZE_MAX,
        "PRIMER_OPT_TM": config.PRIMER_OPT_TM,
        "PRIMER_MIN_TM": config.PRIMER_MIN_TM,
        "PRIMER_MAX_TM": config.PRIMER_MAX_TM,
        "PRIMER_MIN_GC": config.PRIMER_MIN_GC,
        "PRIMER_MAX_GC": config.PRIMER_MAX_GC,
        "PRIMER_MAX_POLY_X": config.PRIMER_MAX_HOMO,
        "PRIMER_SALT_MONOVALENT": config.MV_CONC,
    }

    p3_seq = {
        "SEQUENCE_TEMPLATE": random_reference_slice,
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [-1, -1, -1, -1],
        "SEQUENCE_INCLUDED_REGION": [-1, -1],
    }

    p3_output = primer3.bindings.designPrimers(p3_seq, p3_global)
    text_dir = ("LEFT", "RIGHT")

    for d in range(2):
        num_returned = p3_output[f"PRIMER_{text_dir[d]}_NUM_RETURNED"]
        for i in range(num_returned):
            seq = str(p3_output[f"PRIMER_{text_dir[d]}_{i}_SEQUENCE"])
            p3_penalty = float(p3_output[f"PRIMER_{text_dir[d]}_{i}_PENALTY"])
            primer = Primer(
                seq, 0, Direction.LEFT if text_dir[d] == "LEFT" else Direction.RIGHT, 1
            )

            assert primer.base_penalty == p3_penalty
