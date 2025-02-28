from collections import namedtuple
from functools import cache
from typing import Literal
from riot_na import (
    AirrRearrangementEntryAA,
    create_riot_aa,
    Scheme,
    Organism,
    RiotNumberingAA,
)
from Bio.PDB.Polypeptide import aa1
from riot_na.api.utils import int_to_str_insertion

Position = namedtuple("Position", ["number", "insertion"])
NumberedResidue = namedtuple("NumberedResidue", ["position", "amino_acid"])
NumberingOutput = list[NumberedResidue]

SET_AMINO_ACIDS = set(aa1)

SCHEME_SHORT_TO_LONG = {
    "m": "martin",
    "c": "chothia",
    "k": "kabat",
    "i": "imgt",
    "imgt": "imgt",
    "kabat": "kabat",
    "chothia": "chothia",
    "martin": "martin",
}


@cache
def get_riot_aa(allowed_species: tuple[str, ...]) -> RiotNumberingAA:
    return create_riot_aa(
        allowed_species=[Organism(species) for species in allowed_species]
    )


def map_position_to_tuple(pos: str) -> tuple[int, str]:
    if "." in pos:
        position, insertion = pos.split(".")
        insertion_letter = int_to_str_insertion(int(insertion))
        return Position(number=int(position), insertion=insertion_letter)
    return Position(number=int(pos), insertion=" ")


def validate_sequence(sequence: str):
    """
    Check whether a sequence is a protein sequence or if someone has submitted something nasty.
    """
    assert (
        len(sequence) > 70
    ), f"Sequence too short to be an Ig domain. Please give whole sequence:\n{sequence}"
    assert len(sequence) < 10000, "Sequence too long."
    assert not (
        set(sequence.upper()) - SET_AMINO_ACIDS
    ), "Unknown amino acid letter found in sequence: %s" % ", ".join(
        list((set(sequence.upper()) - SET_AMINO_ACIDS))
    )


def airr_to_numbering_output(airr: AirrRearrangementEntryAA) -> NumberingOutput:
    assert airr.scheme_residue_mapping
    return [
        NumberedResidue(position=map_position_to_tuple(pos), amino_acid=res)
        for pos, res in airr.scheme_residue_mapping.items()
    ]


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


def get_continuous_output(output: NumberingOutput) -> NumberingOutput:
    continuous_output = []
    for i, residue in enumerate(output, start=1):
        residue = NumberedResidue(
                Position(number=i, insertion=" "),
                amino_acid=residue.amino_acid,
            )
        continuous_output.append(residue)
    return continuous_output


def validate_numbering(
    sequence: str, chain: str, airr: AirrRearrangementEntryAA
) -> None:

    allow = [chain]
    if chain == "L":
        allow.append("K")
    assert (
        airr.locus and airr.locus[-1].upper() in allow
    ), f"Sequence provided as an {chain} chain is not recognised as an {chain} chain."

    assert airr.scheme_residue_mapping, f"Could not number sequence: {sequence}"

    imgt_positions = list(airr.scheme_residue_mapping.keys())
    first_imgt_position = int(float(imgt_positions[0]))
    last_imgt_position = int(float(imgt_positions[-1]))

    # Check for missing residues assuming imgt numbering
    assert (
        first_imgt_position < 8 and last_imgt_position > 120
    ), f"Sequence missing too many residues to model correctly. Please give whole sequence:\n{sequence}"


def number_single_sequence(
    sequence: str,
    chain: str,
    scheme: str = "imgt",
    allowed_species: list[str] = ["human", "mouse"],
):
    validate_sequence(sequence)

    try:
        if scheme not in ["raw", "continuous"]:
            scheme = SCHEME_SHORT_TO_LONG[scheme.lower()]
    except KeyError:
        raise NotImplementedError(f"Unimplemented numbering scheme: {scheme}")

    # Use imgt scheme for numbering sanity checks
    riot_aa = get_riot_aa(tuple(allowed_species))
    airr = riot_aa.run_on_sequence(header="", query_sequence=sequence)

    validate_numbering(sequence=sequence, chain=chain, airr=airr)

    output = airr_to_numbering_output(airr)

    match scheme:
        case "imgt":
            return output
        case "raw":
            return get_raw_output(output)
        case "continuous":
            return get_continuous_output(output)
        case _:
            airr = riot_aa.run_on_sequence(
                header="", query_sequence=sequence, scheme=Scheme(scheme)
            )
            return airr_to_numbering_output(airr)


def number_sequences(
    seqs: dict[str, str], scheme="imgt", allowed_species=["human", "mouse"]
) -> dict[str, NumberingOutput]:
    return {
        chain: number_single_sequence(
            seqs[chain], chain, scheme=scheme, allowed_species=allowed_species
        )
        for chain in seqs
    }


def heavy_light_airr_to_numbering_output(
    seqs: dict[str, str],
    airr_dict: dict[str, AirrRearrangementEntryAA],
    scheme: str = "imgt",
) -> dict[str, NumberingOutput]:
    # support only imgt numbering for now
    for airr in airr_dict.values():
        assert (
            airr.numbering_scheme == "imgt"
        ), "Only IMGT scheme numbering is supported when passing airr objects"
    validate_sequence(seqs["H"])
    validate_sequence(seqs["L"])
    validate_numbering(sequence=seqs["H"], chain="H", airr=airr_dict["H"])
    validate_numbering(sequence=seqs["L"], chain="L", airr=airr_dict["L"])
    numbering_outputs: dict[str, NumberingOutput] = {
        chain: airr_to_numbering_output(airr) for chain, airr in airr_dict.items()
    }
    match scheme:
        case "raw":
            return {
            chain: get_raw_output(output) for chain, output in numbering_outputs.items()
            }
        case "continuous":
            return {
            chain: get_continuous_output(output) for chain, output in numbering_outputs.items()
            }
        case "imgt":
            return numbering_outputs
        case _:
            print(f"Requested numbering scheme '{scheme}' but passed AIRR objects contain {airr_dict['H'].numbering_scheme} numbering scheme, which will be used instead.")
            return numbering_outputs
