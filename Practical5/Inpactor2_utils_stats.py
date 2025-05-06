import sys
from Bio import SeqIO
from collections import defaultdict, Counter
import matplotlib.pyplot as plt

#
if len(sys.argv) < 2:
    print("Usage: python count_ltr_classifications.py <path_to_fasta_file>")
    sys.exit(1)

fasta_file = sys.argv[1]

#
classification_sizes = defaultdict(list)
for record in SeqIO.parse(fasta_file, "fasta"):
    classification = record.id.split("#")[1] if "#" in record.id else "Unknown"
    classification_sizes[classification].append(len(record.seq))

#
classification_counts = {cls: len(sizes) for cls, sizes in classification_sizes.items()}
plt.figure(figsize=(10, 8))
plt.bar(classification_counts.keys(), classification_counts.values(), color='blue')
plt.xlabel('Classification')
plt.ylabel('Sequence number')
plt.title('Sequence number according to lineage classification')
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig('LTR_Classification_Counts.pdf')
plt.close()

# 
plt.figure(figsize=(12, 8))
plt.boxplot(classification_sizes.values(), labels=classification_sizes.keys())
plt.xlabel('Classification')
plt.ylabel('Sequence length')
plt.title('Distribution of  Sequence length according to lineage classification')
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig('LTR_Sequence_Sizes_Boxplot.pdf')
plt.close()
