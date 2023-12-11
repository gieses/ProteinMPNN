from proteinmpnn import proteinmpnn_wf

if __name__ == '__main__':
    pmpnn = proteinmpnn_wf.ProteinMPNNWorkflow(pdb_in="fixture/monomer/",
                                               results_dir="fixture/monomer/result",
                                               num_seq_per_target=10)
    pmpnn()
