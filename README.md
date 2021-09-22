# Needleman-Wunsch algorithm
## Options:
+ `-g, --gap`       Penalty for the gap   
+ `-a`              Alphabet
+ \* `-i`           Paths to input sequences
+ `-o`              Path to output file
+ `-optimization`   Enable optimization


## Examples:
+ `-i ./seq1.fasta ./seq2.fasta -a BLOSUM62 -g -1 -optimization false`

+ `-i ./seq1.fasta ./seq2.fasta -a DNAFull -g -5 -o ./out.txt`

+ `-i ./seq1.fasta ./seq2.fasta -a Default`

+ `-i ./seq1.fasta ./seq2.fasta -a Default -optimization true`
