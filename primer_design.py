# primer_design.py

from Bio import SeqIO
import primer3
import csv

def design_primers(input_fasta, output_csv):
    # Primer design settings
    primer_settings = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 22,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 58.0,
        'PRIMER_MAX_TM': 62.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[100, 130]],
    }

    # Open output CSV for writing
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Name", "Forward_Primer", "Reverse_Primer", "F_GC_Content", "R_GC_Content", "F_Tm", "R_Tm", "Product_Size"])

        # Read each sequence from the FASTA
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq_id = record.id
            sequence = str(record.seq)

            # Print to see the sequence and ID being processed
            print(f"Processing {seq_id}...")

            # Design primers
            primers = primer3.bindings.designPrimers(
                {
                    'SEQUENCE_ID': seq_id,
                    'SEQUENCE_TEMPLATE': sequence,
                },
                primer_settings
            )

            # Debug: print the primers to check
            print(f"Primers for {seq_id}: {primers}")

            # Try to extract designed primer info
            try:
                f_primer = primers['PRIMER_LEFT_0_SEQUENCE']
                r_primer = primers['PRIMER_RIGHT_0_SEQUENCE']
                product_size = primers['PRIMER_PAIR_0_PRODUCT_SIZE']

                # Calculate GC content and Tm for forward primer
                f_gc_content = (f_primer.count('G') + f_primer.count('C')) / len(f_primer) * 100
                f_tm = primer3.calcTm(f_primer)

                # Calculate GC content and Tm for reverse primer
                r_gc_content = (r_primer.count('G') + r_primer.count('C')) / len(r_primer) * 100
                r_tm = primer3.calcTm(r_primer)

                # Write to CSV
                writer.writerow([seq_id, f_primer, r_primer, f"{f_gc_content:.2f}", f"{r_gc_content:.2f}", f"{f_tm:.2f}", f"{r_tm:.2f}", product_size])

            except KeyError:
                print(f"⚠️ No suitable primers found for {seq_id}")

if __name__ == "__main__":
    input_fasta = "exampleGenomic.txt"    # <-- Change to your FASTA filename
    output_csv = "exampleprimers.csv"     # <-- Desired output CSV filename

    print(f"Starting primer design using {input_fasta}...")

    design_primers(input_fasta, output_csv)

    print("✅ Primer design finished! Results saved to", output_csv)
