from collections import namedtuple
from functools import cache
from riot_na import create_riot_aa, Scheme, Organism, RiotNumberingAA
from Bio.PDB.Polypeptide import aa1

Position = namedtuple("Position", ["number", "insertion"])
NumberedResidue = namedtuple("NumberedResidue", ["position", "amino_acid"])
NumberingOutput = list[NumberedResidue]


SET_AMINO_ACIDS = set(aa1)
ASCII_SHIFT = 64

SCHEME_SHORT_TO_LONG = {
    "m": "martin",
    "c": "chothia",
    "k": "kabat",
    "imgt": "imgt",
    "kabat": "kabat",
    "chothia": "chothia",
    "martin": "martin",
    "i": "imgt",
}


@cache
def get_riot_aa(allowed_species: tuple[str, ...]) -> RiotNumberingAA:
    return create_riot_aa(
        allowed_species=[Organism(species) for species in allowed_species]
    )


def map_position_to_tuple(pos: str) -> tuple[int, str]:
    if "." in pos:
        position, insertion = pos.split(".")
        insertion_letter = chr(int(insertion) + ASCII_SHIFT)
        return Position(number=int(position), insertion=insertion_letter)
    return Position(number=int(pos), insertion=" ")


def validate_sequence(sequence: str):
    """
    Check whether a sequence is a protein sequence or if someone has submitted something nasty.
    """
    assert len(sequence) < 10000, "Sequence too long."
    assert not (
        set(sequence.upper()) - SET_AMINO_ACIDS
    ), "Unknown amino acid letter found in sequence: %s" % ", ".join(
        list((set(sequence.upper()) - SET_AMINO_ACIDS))
    )


def get_raw_output(output: NumberingOutput) -> NumberingOutput:
    raw_output = [output[0]]
    for residue in output[1:]:
        prev_residue = raw_output[-1]
        if prev_residue.position.number >= residue.position.number:
            residue = NumberedResidue(
                Position(number=prev_residue.position.number + 1, insertion=" "),
                amino_acid=residue.amino_acid,
            )
        raw_output.append(residue)
    return raw_output


def number_single_sequence(
    sequence: str,
    chain: str,
    scheme: str = "imgt",
    allowed_species: list[str] = ["human", "mouse"],
):
    validate_sequence(sequence)

    try:
        if scheme != "raw":
            scheme = SCHEME_SHORT_TO_LONG[scheme.lower()]
    except KeyError:
        raise NotImplementedError(f"Unimplemented numbering scheme: {scheme}")

    assert (
        len(sequence) > 70
    ), f"Sequence too short to be an Ig domain. Please give whole sequence:\n{sequence}"

    allow = [chain]
    if chain == "L":
        allow.append("K")

    # Use imgt scheme for numbering sanity checks
    riot_aa = get_riot_aa(tuple(allowed_species))
    airr = riot_aa.run_on_sequence(header="", query_sequence=sequence)

    assert (
        airr.locus and airr.locus[-1].upper() in allow
    ), f"Sequence provided as an {chain} chain is not recognised as an {chain} chain."

    output = [
        NumberedResidue(position=map_position_to_tuple(pos), amino_acid=res)
        for pos, res in airr.scheme_residue_mapping.items()
    ]
    numbers = [residue.position.number for residue in output]

    # Check for missing residues assuming imgt numbering
    assert (max(numbers) > 120) and (
        min(numbers) < 8
    ), f"Sequence missing too many residues to model correctly. Please give whole sequence:\n{sequence}"
    # Renumber once sanity checks done

    if scheme == "raw":
        output = get_raw_output(output)
    elif scheme != "imgt":
        airr = riot_aa.run_on_sequence(
            header="", query_sequence=sequence, scheme=Scheme(scheme)
        )
        output = [
            (map_position_to_tuple(pos), res)
            for pos, res in airr.scheme_residue_mapping.items()
        ]

    return output


def number_sequences(seqs, scheme="imgt", allowed_species=["human", "mouse"]):
    return {
        chain: number_single_sequence(
            seqs[chain], chain, scheme=scheme, allowed_species=allowed_species
        )
        for chain in seqs
    }
