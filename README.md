# Reversals_identification

## About

### Identification of reversals (inversions) comparing multiple genomes

Input: Tab-delimited text file (.txt) This file should contain the gene order data (from gene 1 to gene 'n') for all the genomes

### Requirements

+ Java
+ python 2.7

## Run

### Create orthologous gene cluster table

Run the script runCreateOrthologousTable.py from the following link (https://github.com/tipputa/Circular-genome-visualizer) to create the gene cluster table.

### Remove the non-conserved gene clusters

Run the removing_nonconserved.py script to obtain only the almost conserved gene clusters.

python removing_nonconserved.py

Note: Name your input file as "FILE.xlsx" and place in the same directory as the script. Output file will be generated in the same directory named "removed_not_conserved.xlsx"

### Generate gene order Data file

Run the script gene_order(almost_conserved).py to obtain gene order data file.
python gene_order(almost_conserved).py

Note: Name your input file as "Genomic_Position.xlsx" and place in the same directory as the script. Output file will be generated in the same directory named "order_new_rotated.xlsx"
Genomic_Position is the file generated by removing_nonconserved.py script.
In Genomic_Position.xlsx file replace the first column (product_name) with 1 to n (up till the total gene clusters in the file)

### Identification of reversals

Run the Reversal_program.java to obtain the reversals information.

Note: Gene order file in tab-delimited text file named "gene_order.txt" should be given as an input the program. The output will be displayed on the console.
