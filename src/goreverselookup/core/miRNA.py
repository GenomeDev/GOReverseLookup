from typing import Dict


class miRNA:
    def __init__(
        self,
        id: str,
        sequence: str = None,
        mRNA_overlaps: Dict[str, float] = None,
        scores: Dict[str, float] = None,
    ) -> None:
        """
        Initializes an instance of miRNA class.

        Args:
        - id: a string that uniquely identifies this miRNA.
        - sequence: an optional string that represents the sequence of this miRNA.
        - mRNA_overlaps: an optional dictionary that represents the overlaps of this miRNA with mRNA sequences.
        - scores: an optional dictionary that represents the scores of this miRNA.

        Returns: None
        """
        self.id = id
        self.sequence = sequence
        self.mRNA_overlaps = {} if mRNA_overlaps is None else mRNA_overlaps.copy()
        self.scores = {} if scores is None else scores.copy()

    @classmethod
    def from_dict(cls, d: dict) -> "miRNA":
        """
        Creates a new instance of miRNA class based on the values in the input dictionary.

        Args:
        - d: a dictionary that represents the values of the miRNA.

        Returns:
        - A new instance of miRNA class based on the values in the input dictionary.
        """
        return cls(d["id"], d.get("sequence"), d.get("mRNA_overlaps"), d.get("scores"))
