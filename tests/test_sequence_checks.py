from dis import disco
import pytest
from riot_na import AirrRearrangementEntryAA, create_riot_aa
from ImmuneBuilder.sequence_checks import (
    ChainType,
    heavy_light_airr_to_numbering_output,
    number_sequences,
    number_single_sequence,
)
from data import IMGT_OUTPUT, RAW_OUTPUT, MARTIN_OUTPUT


def test_number_single_sequence_imgt():

    seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYNMNWVRQAPGKGLEWVSYISSSSSTIYYADSVKGRFTISRDNAKNSLSLQMNSLRDEDTAVYYCARAYYAAAAAAAYGMDVWGQGTTVTVSS"

    output = number_single_sequence(seq, "H", scheme="imgt")

    assert output == IMGT_OUTPUT


def test_number_single_sequence_martin():

    seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYNMNWVRQAPGKGLEWVSYISSSSSTIYYADSVKGRFTISRDNAKNSLSLQMNSLRDEDTAVYYCARAYYAAAAAAAYGMDVWGQGTTVTVSS"

    output = number_single_sequence(seq, "H", scheme="martin")

    assert output == MARTIN_OUTPUT


def test_number_single_sequence_raw():

    seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYNMNWVRQAPGKGLEWVSYISSSSSTIYYADSVKGRFTISRDNAKNSLSLQMNSLRDEDTAVYYCARAYYAAAAAAAYGMDVWGQGTTVTVSS"

    output = number_single_sequence(seq, "H", scheme="raw")

    assert output == RAW_OUTPUT


def test_number_single_sequence_exceptions():

    seq = "LVQPGGSLRLSCAASGFTFSSYNMNWVRQAPGKGLEWVSYISSSSSTIYYADSVKGRFTISRDNAKNSLSLQMNSLRDEDTAVYYCARAYYAAAAAAAYGMDVWGQGTTVTVSS"

    with pytest.raises(
        AssertionError,
        match="Sequence missing too many residues to model correctly. Please give whole sequence",
    ):
        number_single_sequence(seq, "H")

    seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYNMNWVRQAPGKGLEWVSYISSSSSTIYYADSVKGRFTISRDNAKNSLSLQMNSLRDEDTAVYYCARAYYAAAAAAAYGMDVWG"

    with pytest.raises(
        AssertionError,
        match="Sequence missing too many residues to model correctly. Please give whole sequence",
    ):
        number_single_sequence(seq, "H")

    seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYNMNWVRQAPGKGLEWVSYISSSSSTIYYADSVKGRFT"

    with pytest.raises(
        AssertionError,
        match="Sequence too short to be an Ig domain. Please give whole sequence",
    ):
        number_single_sequence(seq, "H")

    seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYNMNWVRQAPGKGLEWVSYISSSSSTIYYADSVKGRFTISRDNAKNSLSLQMNSLRDEDTAVYYCARAYYAAAAAAAYGMDVWGQGTTVTVSS"

    with pytest.raises(
        AssertionError,
        match="Sequence provided as an L chain is not recognised as an L chain.",
    ):
        number_single_sequence(seq, "L")

    seq = "EVQLVESGGGLVQPGGSLRLSCAA!GFTFSSYNMNWVRQAPGKGLEWVSYISSSSSTIYYADSVKGRFTISRDNAKNSLSLQMNSLRDEDTAVYYCARAYYAAAAAAAYGMDVWGQGTTVTVSS"

    with pytest.raises(
        AssertionError, match="Unknown amino acid letter found in sequence: !"
    ):
        number_single_sequence(seq, "H")


def test_airr_to_numbering_output():
    sequence_dict: dict[ChainType, str] = {
        "H": "EVQLVESGGGLVQPGGSLRLSCAASGFSLTIYGAHWVRQAPGKGLEWVSVIWAGGSTNYNSALMSRFTISKDNSKNTVYLQMNSLRAEDTAVYYCARDGSSPYYYSMEYWGQGTTVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKRVEPKSC",
        "L": "EIVLTQSPATLSLSPGERATLSCSATSSVSYMHWFQQKPGQAPRLLIYSTSNLASGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQRSSYPFTFGPGTKLDIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC",
    }
    riot_aa = create_riot_aa()
    airr_dict: dict[ChainType, AirrRearrangementEntryAA] = {chain_type: riot_aa.run_on_sequence(header="", query_sequence=seq) for chain_type, seq in sequence_dict.items()}
    assert number_sequences(sequence_dict) == heavy_light_airr_to_numbering_output(sequence_dict, airr_dict)
