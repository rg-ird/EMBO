import sys
from Bio import SeqIO
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def identify_and_remove_outliers(fasta_file, output_file):
    #
    classification_sizes = defaultdict(list)
    sequences_by_classification = defaultdict(list)

    for record in SeqIO.parse(fasta_file, "fasta"):
        classification = record.id.split("#")[1] if "#" in record.id else "Unknown"
        classification_sizes[classification].append(len(record.seq))
        sequences_by_classification[classification].append(record)

    #
    outliers = defaultdict(list)
    for classification, sizes in classification_sizes.items():
        sizes_np = np.array(sizes)
        q1 = np.percentile(sizes_np, 25)
        q3 = np.percentile(sizes_np, 75)
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        for i, size in enumerate(sizes):
            if size < lower_bound or size > upper_bound:
                outliers[classification].append(sequences_by_classification[classification][i].id)

    #
    filtered_sequences = []
    filtered_classification_sizes = defaultdict(list)
    for classification, records in sequences_by_classification.items():
        for record in records:
            if record.id not in outliers[classification]:
                filtered_sequences.append(record)
                filtered_classification_sizes[classification].append(len(record.seq))

    #
    SeqIO.write(filtered_sequences, output_file, "fasta")

    #
    plt.figure(figsize=(12, 8))
    plt.boxplot(classification_sizes.values(), labels=classification_sizes.keys())
    plt.xlabel('Classification')
    plt.ylabel('Sequences length')
    plt.title('Distribution of Sequences length before filtering of outliers')
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig('LTR_Sequence_Sizes_Boxplot_Before_Filtering.pdf')
    plt.close()

    #
    plt.figure(figsize=(12, 8))
    plt.boxplot(filtered_classification_sizes.values(), labels=filtered_classification_sizes.keys())
    plt.xlabel('Classification')
    plt.ylabel('Sequences length')
    plt.title('Distribution of Sequences length after filtering of outliers')
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig('LTR_Sequence_Sizes_Boxplot_After_Filtering.pdf')
    plt.close()

    print(f"cleaning done --> '{output_file}'.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python count_ltr_classifications.py <path_to_fasta_file> <output_fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_file = sys.argv[2]
    identify_and_remove_outliers(fasta_file, output_file)
