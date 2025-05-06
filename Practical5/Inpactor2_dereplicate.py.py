import os
import subprocess
import sys
import shutil

def run_cmd(cmd):
    print("[commande]", ' '.join(cmd))
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print("❌ Erreur :")
        print(result.stderr.decode())
        sys.exit(1)

if len(sys.argv) < 2:
    print("Usage : python dereplicate_easycluster.py fichier.fasta [identity=0.8]")
    sys.exit(1)

# === Paramètres ===
input_fasta = sys.argv[1]
identity = float(sys.argv[2]) if len(sys.argv) > 2 else 0.8
basename = os.path.splitext(os.path.basename(input_fasta))[0]

output_dir = "mmseqs_output"
tmp_dir = "mmseqs_tmp"

rep_fasta = f"{basename}.clust"
tsv_output = f"{output_dir}_cluster.tsv"

# === Nettoyage préalable ===
shutil.rmtree(output_dir, ignore_errors=True)
shutil.rmtree(tmp_dir, ignore_errors=True)

# === Exécution de MMseqs2 easy-cluster ===
run_cmd([
    "mmseqs", "easy-cluster",
    input_fasta,
    output_dir,
    tmp_dir,
    "--min-seq-id", str(identity),
    "--cov-mode", "1",
    "-c", "0.8"
])

# === Renommer la sortie FASTA ===
rep_fasta_source = f"{output_dir}_rep_seq.fasta"
if os.path.exists(rep_fasta_source):
    shutil.copyfile(rep_fasta_source, rep_fasta)
    print(f"✔ Séquences non redondantes : {rep_fasta}")
else:
    print("❌ Fichier FASTA de sortie introuvable.")
    sys.exit(1)

# === Afficher l'existence du fichier TSV ===
if os.path.exists(tsv_output):
    print(f"✔ Fichier des clusters : {tsv_output}")
else:
    print("⚠️ Fichier des clusters introuvable.")

# === Nettoyage des fichiers MMseqs inutiles ===
for f in os.listdir("."):
    if f.startswith("mmseqs_output") and not f.endswith("_cluster.tsv"):
        os.remove(f)
shutil.rmtree(tmp_dir, ignore_errors=True)

print("🧹 Nettoyage terminé : fichiers temporaires supprimés.")

