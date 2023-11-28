"""
A convenience wrapper for the ProteinMPNN package.

# TODO:
* make_tied_positions_dict support
* homooligomer

"""
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


def read_fasta_sequences(fasta_file: str) -> list[str]:
    """Retrieve sequences from a fasta file."""
    sequences = []
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            sequences.append(line.strip().replace("/", ""))
    return sequences


class ProteinMPNNWorkflow:
    def __init__(self,
                 *,
                 pdb_in: str,
                 pdb_chains: Optional[str] = None,
                 results_dir: str,
                 num_seq_per_target: int = 10,
                 sampling_temp: float = 0.5,
                 score_only: bool = False,
                 seed: int = 13,
                 batch_size: int = 32,
                 fix_positions: Optional[str] = None,
                 inverse_fix_positions: bool = False,
                 verbose: bool = True,
                 path_to_model_weights: str = None,
                 use_antibody_model: bool = False,
                 ):

        pdb_in = Path(pdb_in)
        if not pdb_in.exists():
            raise FileNotFoundError(f"pdb_in file {pdb_in} does not exist.")

        results_dir = Path(results_dir)
        results_dir.mkdir(exist_ok=True)

        self.pdb_in = pdb_in
        self.fasta_in = None
        self.pdb_chains = pdb_chains
        self.results_dir = results_dir
        self.verbose = verbose

        # sampling parameters
        self.num_seq_per_target = num_seq_per_target
        self.sampling_temp = sampling_temp
        self.score_only = score_only
        self.seed = seed
        if batch_size > num_seq_per_target:
            batch_size = num_seq_per_target
        self.batch_size = batch_size

        # model
        self.path_to_model_weights = Path(path_to_model_weights) if path_to_model_weights else None
        self.use_antibody_model = use_antibody_model

        # proteinmpnn uses these default folders
        self.sequence_folder = self.results_dir / "seqs"
        self.score_folder = self.results_dir / "score_only"

        # intermediate files
        self.parsed_pdbs_jsonl = None
        self.assigned_pdbs_jsonl = None
        self.fixed_pdbs_jsonl = None
        self.fix_positions = fix_positions
        self.inverse_fix_positions = inverse_fix_positions
        self.calls = []

    def make_inputs(self):
        pass

    def score(self, *, sequences: Optional[list[str]] = None, fasta_in: Optional[str] = None):
        """
        Use ProtienMPNN to score the sequence against the prepared pdb.
        """
        # use temporary file to score sequences
        with tempfile.TemporaryDirectory() as tmpdirname:
            # fasta file
            if sequences is None and fasta_in is None:
                raise ValueError("No sequences provided and no fasta file provided.")

            elif sequences is None and fasta_in is not None:
                print("Using provided fasta file.")
                self.fasta_in = fasta_in
                sequences = read_fasta_sequences(fasta_in)

            elif sequences is not None and fasta_in is None:
                # directory for temporary files
                tmpdir = Path(tmpdirname)
                tmpdir.mkdir(exist_ok=True)
                seq_file = tmpdir / "seqs.fasta"

                print(f"write sequences to: {tmpdir}")
                self.fasta_in = seq_file
                for idx, seq in enumerate(sequences):
                    with open(tmpdir / f"seq_{idx}.fasta", "w") as f:
                        f.write(f">seq_{str(idx).zfill(4)}\n{seq}\n")
            else:
                raise ValueError("Provide either sequences or fasta_in.")

            # run scoring
            self._call_proteinmpnn(score=True)

        sequence_df = self._parse_scores()
        sequence_df["seq_id"] = sequence_df["sequence"].apply(str)
        sequence_df = sequence_df.set_index("seq_id")
        sequence_df = sequence_df.loc[sequences].reset_index(drop=True)
        return sequence_df

    def _parse_scores(self) -> pd.DataFrame:
        """Parse the npz files from ProteinMPNN."""
        npz_files = list(Path(self.score_folder).glob("*.npz"))
        results_data = {"sequence": [],
                        "score_mean": [],
                        "score_std": [],
                        "global_score_mean": [],
                        "global_score_std": [],
                        "sample_size": []}
        # numpy data has the keys: ['score', 'global_score', 'S', 'seq_str']
        for npz_file in npz_files:
            # this is the PDB npz file from the input
            if "fasta" not in npz_file.name:
                continue

            # sampled sequences
            tmp_data = np.load(npz_file)
            scores = tmp_data["score"]
            global_scores = tmp_data["global_score"]
            sequence = tmp_data["seq_str"]
            results_data["sequence"].append(sequence)
            results_data["score_mean"].append(scores.mean())
            results_data["score_std"].append(scores.std())
            results_data["global_score_mean"].append(global_scores.mean())
            results_data["global_score_std"].append(global_scores.std())
            results_data["sample_size"].append(len(scores))
        return pd.DataFrame(results_data)

    @property
    def chain_string(self):
        """Return the chain string for the proteinmpnn command."""
        return f'"{self.pdb_chains}"'

    def _run_cli(self, command: str, tool: str):
        """Run the command line interface."""
        command_joined = " ".join(command)

        if self.verbose:
            print(command_joined)

        self.calls.append((tool, command_joined))
        log = subprocess.check_output(command_joined, shell=True)
        print(log)

    def _call_proteinmpnn(self, score: bool = False):
        """Call ProteinMPNN"""
        command = ["proteinmpnn",
                   "--out_folder", str(self.results_dir),
                   "--num_seq_per_target", str(self.num_seq_per_target),
                   "--sampling_temp", str(self.sampling_temp),
                   "--score_only", str(int(self.score_only)),
                   "--seed", str(self.seed),
                   "--batch_size", str(self.batch_size)]

        # chains are needed for score and sample, if provided
        if self.pdb_chains:
            command.extend(["--pdb_path_chains", self.chain_string])

        if self.use_antibody_model:
            command.extend(["--use_antibody_model"])
            command.extend(["--model_name", "abmpnn"])

        if self.path_to_model_weights:
            command.extend(["--path_to_model_weights", str(self.path_to_model_weights)])

        if score:
            # requires a fasta but no parsing
            command.extend(["--pdb_path", str(self.pdb_in)])
            command.extend(["--path_to_fasta", str(self.fasta_in)])
        else:
            command.extend(["--jsonl_path", str(self.parsed_pdbs_jsonl)])

            if self.pdb_chains:
                # parsed
                command.extend(["--chain_id_jsonl", str(self.assigned_pdbs_jsonl)])

                if self.fix_positions:
                    command.extend(["--fixed_positions_jsonl", str(self.fixed_pdbs_jsonl)])

        self._run_cli(command, "proteinmpnn")

    def _call_proteinmpnn_parse(self):
        """Call ProteinMPNN"""
        self.parsed_pdbs_jsonl = self.results_dir / "parsed_pdbs.jsonl"
        command = ["proteinmpnn-parse",
                   "--input_path", str(self.pdb_in),
                   "--output_path", str(self.parsed_pdbs_jsonl)]
        self._run_cli(command, "proteinmpnn-parse")

    def _call_proteinmpnn_assign_chains(self):
        """Call fix positions / chains"""
        self.assigned_pdbs_jsonl = self.results_dir / "assigned_pdbs.jsonl"
        command = ["proteinmpnn-assign",
                   "--input_path", str(self.parsed_pdbs_jsonl),
                   "--output_path", str(self.assigned_pdbs_jsonl),
                   "--chain_list", self.chain_string,
                   ]
        self._run_cli(command, "proteinmpnn-assign")

    def __call_proteinmpnn_fix_positions(self):
        """Call fix positions / chains"""
        self.fixed_pdbs_jsonl = self.results_dir / "fixed_pdbs.jsonl"
        command = ["proteinmpnn-fix-positions",
                   "--input_path", str(self.parsed_pdbs_jsonl),
                   "--output_path", str(self.fixed_pdbs_jsonl),
                   "--chain_list", self.chain_string,
                   "--position_list", f'"{self.fix_positions}"',
                   ]

        if self.inverse_fix_positions:
            command.append("--specify_non_fixed")

        self._run_cli(command, "proteinmpnn-fix-positions")

    def __call__(self):
        # parse_multiple_chains
        if self.verbose:
            print("Parsing chains...")
        self._call_proteinmpnn_parse()

        # assign_fixed_chains
        if self.verbose:
            print("Assigning fixed chains...")
        if self.fix_positions or self.pdb_chains:
            self._call_proteinmpnn_assign_chains()

        if self.verbose:
            print("Fixing positions...")
        if self.fix_positions:
            self.__call_proteinmpnn_fix_positions()

        # protein_mpnn_run
        if self.verbose:
            print("Running ProteinMPNN...")
        self._call_proteinmpnn()

    def get_sequences(self, pdb: str):
        return read_fasta_sequences(str(self.sequence_folder / f"{pdb}.fa"))


def get_pssm_example():
    pssm_ar = np.load("../examples/inputs/PSSM_inputs/3HTN.npz")
    print(list(pssm_ar.keys()))
    A_coef = pd.DataFrame(pssm_ar["A_coef"])
    A_bias = pd.DataFrame(pssm_ar["A_bias"])
    A_odds = pd.DataFrame(pssm_ar["A_odds"])


def create_parser():
    import argparse
    parser = argparse.ArgumentParser(description="Run ProteinMPNN on a PDB file.")

    # input parameters
    group_input = parser.add_argument_group('General: Control in- and output')
    group_input.add_argument("--pdb_in", type=str, required=True,
                             help="Input PDB directory.")
    group_input.add_argument("--results_out", type=str, required=True,
                             help="Results directory.")

    # model group
    group_model = parser.add_argument_group('Model Parameters')
    group_model.add_argument("--path_to_model_weights", type=str, default=None,
                             help="Path to model weights (default: None).")
    group_model.add_argument("--use_antibody_model", action="store_true", default=False,
                             help="Use antibody model (default: False).")
    group_model.add_argument("--model_name", type=str, default=None,
                             help="Model name (default: None).")

    # constrain parameters
    constrain_group = parser.add_argument_group('Sampling constrain parameters')
    constrain_group.add_argument("--fix_positions", type=str, default=None,
                                 help="Positions to fix (default: None). Example: '1 2 3 4 5'")
    constrain_group.add_argument("--pdb_chains", type=str, default=None,
                                 help="PDB chains to use (default: None). Example: 'A B'")
    constrain_group.add_argument("--inverse_fix_positions", action="store_true", default=False,
                                 help="Inverse the fixed positions (default: False).")

    # sampling parameters for scoring and sequence sampling
    group_sampling = parser.add_argument_group('Sampling parameters')
    group_sampling.add_argument("--num_seq_per_target", type=int, default=10,
                                help="Number of sequences to sample per target (default: 10).")
    group_sampling.add_argument("--sampling_temp", type=float, default=0.25,
                                help="Sampling temperature (default: 0.25).")
    group_sampling.add_argument("--score_only", action="store_true", default=False,
                                help="Only score the sequences (default: False)")

    # parameters influencing the run environment
    group_run_parameters = parser.add_argument_group('Run Parameters')
    group_run_parameters.add_argument("--seed", type=int, default=13,
                                      help="Random seed (default: 13).")
    group_run_parameters.add_argument("--batch_size", type=int, default=32,
                                      help="Batch size (default: 32).")
    group_run_parameters.add_argument("--verbose", action="store_true", default=False,
                                      help="Verbose output (default: False).")

    return parser.parse_args()


def main():
    args = create_parser()
    print("Running ProteinMPNN...")
    pmpnn_model = ProteinMPNNWorkflow(pdb_in=args.pdb_in,
                                      results_dir=args.results_out,
                                      num_seq_per_target=args.num_seq_per_target,
                                      sampling_temp=args.sampling_temp,
                                      score_only=args.score_only,
                                      seed=args.seed,
                                      batch_size=args.batch_size,
                                      pdb_chains=args.pdb_chains,
                                      fix_positions=args.fix_positions,
                                      inverse_fix_positions=args.inverse_fix_positions,
                                      verbose=args.verbose,
                                      path_to_model_weights=args.path_to_model_weights,
                                      use_antibody_model=args.use_antibody_model,
                                      )
    pmpnn_model()
    print(f"Results Dir: {pmpnn_model.results_dir}")


if __name__ == '__main__':
    main()
