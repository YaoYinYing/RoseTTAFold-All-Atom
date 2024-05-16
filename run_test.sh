#! bash
set -e

export HYDRA_FULL_ERROR=1

# echo testing 7u7w_protein
# rf2aa_inference --config-name=base \
#     job_name=7u7w_protein \
#     +protein_inputs.A.fasta_file=examples/protein/7u7w_A.fasta \
#     output_path=examples/output/7u7w_protein/



echo testing 7u7w_protein_nucleic
rf2aa_inference --config-name=base \
    job_name=7u7w_protein_nucleic \
    +protein_inputs.A.fasta_file=examples/protein/7u7w_A.fasta \
    +na_inputs.B.fasta=examples/nucleic_acid/7u7w_B.fasta \
    +na_inputs.B.input_type=dna \
    +na_inputs.C.fasta=examples/nucleic_acid/7u7w_C.fasta \
    +na_inputs.C.input_type=dna \
    output_path=examples/output/7u7w_protein_nucleic/


echo testing 3fap_protein_sm
rf2aa_inference --config-name=base \
    job_name=3fap_protein_sm \
    +protein_inputs.A.fasta_file=examples/protein/3fap_A.fasta \
    +protein_inputs.B.fasta_file=examples/protein/3fap_B.fasta \
    +sm_inputs.C.input=examples/small_molecule/ARD_ideal.sdf \
    +sm_inputs.C.input_type=sdf \
    output_path=examples/output/3fap_protein_sm/


echo testing 3fap_protein_nucleic_sm
rf2aa_inference --config-name=base \
    job_name=7u7w_protein_nucleic_sm \
    +protein_inputs.A.fasta_file=examples/protein/7u7w_A.fasta \
    +na_inputs.B.fasta=examples/nucleic_acid/7u7w_B.fasta \
    +na_inputs.B.input_type=dna \
    +na_inputs.C.fasta=examples/nucleic_acid/7u7w_C.fasta \
    +na_inputs.C.input_type=dna \
    +sm_inputs.D.input=examples/small_molecule/XG4.sdf \
    +sm_inputs.D.input_type=sdf \
    output_path=examples/output/7u7w_protein_nucleic_sm/


echo testing covalently_modified_7s69_A
rf2aa_inference --config-name=base \
    job_name=covalently_modified_7s69_A \
    +protein_inputs.A.fasta_file=examples/protein/7s69_A.fasta \
    +sm_inputs.B.input=examples/small_molecule/7s69_glycan.sdf \
    +sm_inputs.B.input_type=sdf \
    '+covale_inputs="A,74,ND2:B,1:CW,null"' \
    output_path=examples/output/covalently_modified_7s69_A/

