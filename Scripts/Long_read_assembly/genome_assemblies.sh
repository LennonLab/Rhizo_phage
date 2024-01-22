#!/bin/bash
#SBATCH --job-name=genome_assembly
#SBATCH -p general
#SBATCH -o genome_assembly%j.txt
#SBATCH -e genome_assembly%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jgmcmull@iu.edu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time=23:59:00
#SBATCH --mem=100gb
#SBATCH -A r00324

#Command line arguments
initials="" #i
target_bp="" #t
coverage_subset="" #c
genome_size="" #g
reads_file=""#r
threads="" #n

# Function to display script usage
usage() {
    echo "Usage: $0 -i <initials> -t <target_bp> -c <coverage_subset> -g <genome_size> -r <reads_file> -n <threads>"
    echo "  -i <initials>    : Specify the initials of the user for the working directory name, e.g., jgm (mandatory)"
    echo "  -t <target_bp> : Specify the target bp value for filtlong, e.g., 1000000000 (mandatory)"
    echo "  -c <coverage_subset>    : Specify the coverage of the subsetted reads, e.g., 50 (mandatory)"
    echo "  -g <genome_size>    : Specify the genome size in mb, e.g., $genome_size (mandatory)"
    echo "  -r <reads_file>    : Specify the path to the directory of the raw reads (mandatory)"
    echo "  -n <threads>    : Specify the number of threads in the script allocated, e.g., 24 (mandatory)"
    echo "  -h           : Display this help message"
    exit 1
}

# Parse command-line options
while getopts ":i:t:c:g:r:n:h" opt; do
    case $opt in
        i)
            initials="$OPTARG"
            ;;
        t)
            target_bp="$OPTARG"
            ;;
        c)
            coverage_subset="$OPTARG"
            ;;
        g)
            genome_size="$OPTARG"
            ;;
        r)
            reads_file="$OPTARG"
            ;;
        n)
            threads="$OPTARG"
            ;;
        h)
            usage
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done

# Check if mandatory argument is provided
if [ -z "$initials" ] || [ -z "$target_bp" ] || [ -z "$coverage_subset" ] || [ -z "$genome_size" ] || [ -z "$reads_file" ] || [ -z "$threads" ]; then
    echo "Error: -n <name>, -t <target_bp>, -c <coverage_subset>, -g <genome_size>, -r <reads_file>, and -n <threads> are mandatory arguments. Please use -h if you need help."
    usage
fi

mkdir $initials

#Save the arguments into a log file
touch "$initials"/arguments_log.txt
echo "The input arguments were:
		-i = $initials
		-t = $target_bp
		-c = $coverage_subset
		-g = $genome_size
		-r = $reads_file
		-n = $threads" >> "$initials"/arguments_log.txt

mkdir "$initials"/clean

#filtlong will obtain the reads that are at least 1 kb, it will keep top 90% of the best reads (based on bp, with the total estimated genome size from NCBI), and over a scanning 250bp window, reads kept must have at least an average quality score of 10 (90% accuracy)
filtlong --min_length 1000 --keep_percent 90 --mean_q_weight 10 --target_bases $target_bp "$reads_file" | gzip > "$initials"/clean/clean_reads.fastq.gz

#Subsample to get <coverage_subset> x coverage of genomes 3 times
rasusa --input "$initials"/clean/clean_reads.fastq.gz --coverage $coverage_subset --genome-size $genome_size --seed 123 --output "$initials"/clean/clean_"$coverage_subset"x_1.fastq.gz
rasusa --input "$initials"/clean/clean_reads.fastq.gz --coverage $coverage_subset --genome-size $genome_size --seed 456 --output "$initials"/clean/clean_"$coverage_subset"x_2.fastq.gz
rasusa --input "$initials"/clean/clean_reads.fastq.gz --coverage $coverage_subset --genome-size $genome_size --seed 789 --output "$initials"/clean/clean_"$coverage_subset"x_3.fastq.gz
rasusa --input "$initials"/clean/clean_reads.fastq.gz --coverage $coverage_subset --genome-size $genome_size --seed 258 --output "$initials"/clean/clean_"$coverage_subset"x_4.fastq.gz
rasusa --input "$initials"/clean/clean_reads.fastq.gz --coverage $coverage_subset --genome-size $genome_size --seed 5 --output "$initials"/clean/clean_"$coverage_subset"x_5.fastq.gz
rasusa --input "$initials"/clean/clean_reads.fastq.gz --coverage $coverage_subset --genome-size $genome_size --seed 6 --output "$initials"/clean/clean_"$coverage_subset"x_6.fastq.gz
rasusa --input "$initials"/clean/clean_reads.fastq.gz --coverage $coverage_subset --genome-size $genome_size --seed 7 --output "$initials"/clean/clean_"$coverage_subset"x_7.fastq.gz
rasusa --input "$initials"/clean/clean_reads.fastq.gz --coverage $coverage_subset --genome-size $genome_size --seed 8 --output "$initials"/clean/clean_"$coverage_subset"x_8.fastq.gz
rasusa --input "$initials"/clean/clean_reads.fastq.gz --coverage $coverage_subset --genome-size $genome_size --seed 9 --output "$initials"/clean/clean_"$coverage_subset"x_9.fastq.gz
rasusa --input "$initials"/clean/clean_reads.fastq.gz --coverage $coverage_subset --genome-size $genome_size --seed 10 --output "$initials"/clean/clean_"$coverage_subset"x_10.fastq.gz
rasusa --input "$initials"/clean/clean_reads.fastq.gz --coverage $coverage_subset --genome-size $genome_size --seed 11 --output "$initials"/clean/clean_"$coverage_subset"x_11.fastq.gz
rasusa --input "$initials"/clean/clean_reads.fastq.gz --coverage $coverage_subset --genome-size $genome_size --seed 12 --output "$initials"/clean/clean_"$coverage_subset"x_12.fastq.gz

mkdir "$initials"/assemblies

source ~/miniconda3/bin/activate flye_29

#Assemble each subset with Flye, raven, and miniams+minipolish
flye --nano-hq "$initials"/clean/clean_"$coverage_subset"x_1.fastq.gz --threads $threads --out-dir assembly_01 && cp assembly_01/assembly.fasta "$initials"/assemblies/assembly_01.fasta && cp assembly_01/assembly_graph.gfa "$initials"/assemblies/assembly_01.gfa && rm -r assembly_01
bash miniasm_and_minipolish.sh "$initials"/clean/clean_"$coverage_subset"x_2.fastq.gz $threads > "$initials"/assemblies/assembly_02.gfa && any2fasta "$initials"/assemblies/assembly_02.gfa > "$initials"/assemblies/assembly_02.fasta
raven --threads $threads --disable-checkpoints --graphical-fragment-assembly "$initials"/assemblies/assembly_03.gfa "$initials"/clean/clean_"$coverage_subset"x_3.fastq.gz > "$initials"/assemblies/assembly_03.fasta

flye --nano-hq "$initials"/clean/clean_"$coverage_subset"x_4.fastq.gz --threads $threads --out-dir assembly_04 && cp assembly_04/assembly.fasta "$initials"/assemblies/assembly_04.fasta && cp assembly_04/assembly_graph.gfa "$initials"/assemblies/assembly_04.gfa && rm -r assembly_04
bash miniasm_and_minipolish.sh "$initials"/clean/clean_"$coverage_subset"x_5.fastq.gz $threads > "$initials"/assemblies/assembly_05.gfa && any2fasta "$initials"/assemblies/assembly_05.gfa > "$initials"/assemblies/assembly_05.fasta
raven --threads $threads --disable-checkpoints --graphical-fragment-assembly "$initials"/assemblies/assembly_06.gfa "$initials"/clean/clean_"$coverage_subset"x_6.fastq.gz > "$initials"/assemblies/assembly_06.fasta

flye --nano-hq "$initials"/clean/clean_"$coverage_subset"x_7.fastq.gz --threads $threads --out-dir assembly_07 && cp assembly_07/assembly.fasta "$initials"/assemblies/assembly_07.fasta && cp assembly_07/assembly_graph.gfa "$initials"/assemblies/assembly_07.gfa && rm -r assembly_07
bash miniasm_and_minipolish.sh "$initials"/clean/clean_"$coverage_subset"x_8.fastq.gz $threads > "$initials"/assemblies/assembly_08.gfa && any2fasta "$initials"/assemblies/assembly_08.gfa > "$initials"/assemblies/assembly_08.fasta
raven --threads $threads --disable-checkpoints --graphical-fragment-assembly "$initials"/assemblies/assembly_09.gfa "$initials"/clean/clean_"$coverage_subset"x_9.fastq.gz > "$initials"/assemblies/assembly_09.fasta

flye --nano-hq "$initials"/clean/clean_"$coverage_subset"x_10.fastq.gz --threads $threads --out-dir assembly_10 && cp assembly_10/assembly.fasta "$initials"/assemblies/assembly_10.fasta && cp assembly_10/assembly_graph.gfa "$initials"/assemblies/assembly_10.gfa && rm -r assembly_10
bash miniasm_and_minipolish.sh "$initials"/clean/clean_"$coverage_subset"x_11.fastq.gz $threads > "$initials"/assemblies/assembly_11.gfa && any2fasta "$initials"/assemblies/assembly_11.gfa > "$initials"/assemblies/assembly_11.fasta
raven --threads $threads --disable-checkpoints --graphical-fragment-assembly "$initials"/assemblies/assembly_12.gfa "$initials"/clean/clean_"$coverage_subset"x_12.fastq.gz > "$initials"/assemblies/assembly_12.fasta

source ~/miniconda3/bin/deactivate

#Cluster 12 assemblies with trycycler
trycycler cluster --assemblies "$initials"/assemblies/*.fasta --threads $threads --reads "$initials"/clean/clean_reads.fastq.gz --out_dir "$initials"/trycycler

