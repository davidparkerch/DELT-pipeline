"""RDKit helpers shared by the shortlist, the model and the explorer, so they all fingerprint alike."""

from rdkit import RDLogger
from rdkit.Chem import Mol, MolFromSmiles, rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

FP_RADIUS = 2
FP_BITS = 2048
"""The Morgan fingerprint KinDEL uses (radius 2, 2048 bits); RDKit.js computes the same one in the browser."""


def morgan_generator() -> rdFingerprintGenerator.FingeprintGenerator64:
    """A Morgan fingerprint generator with the pipeline's settings."""
    return rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)


def parse_smiles(text: str) -> Mol:
    """Parse a SMILES string.

    Raises:
        ValueError: If RDKit cannot parse it.
    """
    molecule = MolFromSmiles(text)
    if molecule is None:
        msg = f"cannot parse SMILES {text!r}"
        raise ValueError(msg)
    return molecule
