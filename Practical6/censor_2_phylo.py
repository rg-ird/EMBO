import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

print("######## You need EMBOSS, MAFFT, FASTTREE and GENEWISE installed")

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    if isinstance(seq, Seq):
        return str(seq.reverse_complement())
    else:
        return str(Seq(seq).reverse_complement())

def extract_and_align_sequences(input_genome_fasta, ref_domain_fasta, input_map_file, output_folder, final_fasta_base, min_similarity=0.0, min_coverage=0.0):
    os.environ['WISECONFIGDIR'] = '/shared/home/rguyot/wisecfg'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    sequences = SeqIO.to_dict(SeqIO.parse(input_genome_fasta, "fasta"))
    domains = SeqIO.to_dict(SeqIO.parse(ref_domain_fasta, "fasta"))
    aggregated_records = []

    with open(input_map_file, 'r') as file:
        for line in file:
            elements = line.strip().split()
            if len(elements) < 8:
                print(f"Skipping malformed line: {line}")
                continue

            try:
                similarity_raw = float(elements[7])  # colonne 8 (entre 0 et 1)
                coverage = float(elements[-1])       # dernière colonne
            except (ValueError, IndexError):
                print(f"Skipping line with invalid similarity or coverage: {line}")
                continue

            similarity = similarity_raw * 100.0  # conversion en pourcentage

            if similarity < min_similarity or coverage < min_coverage:
                print(f"Skipping due to similarity ({similarity:.2f}%) < {min_similarity}% or coverage ({coverage}) < {min_coverage}: {line}")
                continue

            seq_id = elements[0]
            start = int(elements[1]) - 200
            end = int(elements[2]) + 200
            strand = elements[6]
            domain_id = elements[3]

            print(f"Processing {seq_id} from {start} to {end}.")

            start = max(start, 0)
            end = min(end, len(sequences[seq_id].seq))

            extracted_seq = sequences[seq_id].seq[start:end]
            if strand.lower() == 'c':
                extracted_seq = reverse_complement(extracted_seq)

            if not isinstance(extracted_seq, Seq):
                extracted_seq = Seq(extracted_seq)

            extracted_record = SeqRecord(extracted_seq,
                                         id=f"{seq_id}-{elements[1]}_{elements[2]}",
                                         description="")

            temp_nuc_file = f"{output_folder}/{seq_id}-{elements[1]}_{elements[2]}.fasta"
            domain_file = f"{output_folder}/{domain_id}.fasta"
            result_file = f"{output_folder}/{seq_id}-{elements[1]}_{elements[2]}.domain.fasta"

            SeqIO.write(extracted_record, temp_nuc_file, "fasta")

            if domain_id not in domains:
                print(f"Domain ID {domain_id} not found in {ref_domain_fasta}")
                continue
            SeqIO.write(domains[domain_id], domain_file, "fasta")

            print(f"Running genewise: {domain_file} {temp_nuc_file} -pep > {result_file}")
            cmd = f"genewise {domain_file} {temp_nuc_file} -pep > {result_file}"
            process = subprocess.run(cmd, shell=True, capture_output=True)
            if process.returncode != 0:
                print(f"Error running genewise: {process.stderr.decode()}")
                continue

            if not os.path.exists(result_file) or os.path.getsize(result_file) == 0:
                print(f"Genewise did not produce an output file or the file is empty: {result_file}")
                continue

            for record in SeqIO.parse(result_file, 'fasta'):
                if len(record.seq) >= 200:
                    aggregated_records.append(record)

    final_fasta_path = os.path.join(output_folder, f"{final_fasta_base}.genewise.fas")
    if aggregated_records:
        SeqIO.write(aggregated_records, final_fasta_path, 'fasta')
        print(f"All sequences saved to {final_fasta_path}.")
    else:
        print("No sequences were aggregated. The final fasta file will be empty.")

    return final_fasta_path

def cleanup_intermediate_files(output_folder, final_fasta_base):
    """Remove all files not starting with final_fasta_base in output_folder."""
    for f in os.listdir(output_folder):
        file_path = os.path.join(output_folder, f)
        if os.path.isfile(file_path) and not f.startswith(final_fasta_base):
            os.remove(file_path)
            print(f"Deleted intermediate file: {file_path}")

def concatenate_domains(final_fasta_path, domain_fasta):
    """Concatenate the final fasta with the domain fasta."""
    with open(final_fasta_path, 'a') as outfile:
        with open(domain_fasta, 'r') as infile:
            domain_content = infile.read()
            if domain_content.strip():
                outfile.write(domain_content)
                print(f"Appended domain sequences to {final_fasta_path}:")
                print(domain_content)
            else:
                print(f"No content found in {domain_fasta} to append.")
    print(f"Domains concatenated to {final_fasta_path}")

def clean_fasta_file(fasta_path):
    """Clean the FASTA file by removing '//' from sequence ends and '.pep' from IDs."""
    cleaned_records = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequence = str(record.seq).replace("//", "")
        record.id = record.id.replace(".pep", "")
        record.seq = Seq(sequence)
        cleaned_records.append(record)

    cleaned_fasta_path = fasta_path.replace('.fas', '.clean.fasta')
    SeqIO.write(cleaned_records, cleaned_fasta_path, "fasta")
    print(f"Cleaned FASTA file saved to {cleaned_fasta_path}")
    return cleaned_fasta_path

def run_mafft(input_fasta, num_threads=12):
    output_alignment = input_fasta.replace('.clean.fasta', '.clean.aln')
    cmd = f"mafft --auto --thread {num_threads} {input_fasta} > {output_alignment}"
    process = subprocess.run(cmd, shell=True, capture_output=True)
    if process.returncode != 0:
        print(f"Error running MAFFT: {process.stderr.decode()}")
        return None
    print(f"Alignment completed and saved to {output_alignment}")
    return output_alignment

def run_fasttree(input_alignment):
    output_tree = input_alignment.replace('.clean.aln', '.tree')
    cmd = f"FastTree -out {output_tree} {input_alignment}"
    process = subprocess.run(cmd, shell=True, capture_output=True)
    if process.returncode != 0:
        print(f"Error running FastTree: {process.stderr.decode()}")
        return None
    print(f"Phylogenetic tree completed and saved to {output_tree}")
    return output_tree

if __name__ == "__main__":
    if len(sys.argv) < 6 or len(sys.argv) > 8:
        print("Usage:\n  python script.py \\\n    input_genome_fasta.fasta \\\n    ref_domain_fasta.fasta \\\n    input_map_file.txt \\\n    output_folder \\\n    final_fasta_base \\\n    [min_similarity] [min_coverage]\n")
        print("Note: min_similarity should be given in %, even if values in map are between 0 and 1.")
    else:
        min_similarity = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0
        min_coverage = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0

        final_fasta_path = extract_and_align_sequences(
            sys.argv[1], sys.argv[2], sys.argv[3],
            sys.argv[4], sys.argv[5],
            min_similarity, min_coverage
        )

        if os.path.exists(final_fasta_path) and os.path.getsize(final_fasta_path) > 0:
            concatenate_domains(final_fasta_path, sys.argv[2])
            cleaned_fasta_path = clean_fasta_file(final_fasta_path)
            aligned_fasta = run_mafft(cleaned_fasta_path)
            if aligned_fasta:
                phylogenetic_tree = run_fasttree(aligned_fasta)
                if phylogenetic_tree:
                    print(f"Phylogenetic analysis completed and saved to {phylogenetic_tree}")
                else:
                    print(f"Phylogenetic analysis failed.")
            else:
                print(f"Alignment failed.")
        else:
            print(f"The final fasta file {final_fasta_path} is empty or was not created.")

        cleanup_intermediate_files(sys.argv[4], sys.argv[5])

