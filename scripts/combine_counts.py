import pandas as pd

count_files = snakemake.input
output_file = snakemake.output[0]

combined = None

for file in count_files:
    sample = str(file).split("/")[-1].replace("_gene_counts.txt", "")

    df = pd.read_csv(file, sep="\t", comment="#")

    df = df[["Geneid", df.columns[-1]]]
    df.columns = ["Geneid", sample]

    if combined is None:
        combined = df
    else:
        combined = combined.merge(
            df,
            on="Geneid",
            how="inner",
            validate="one_to_one"
        )

combined.to_csv(output_file, sep="\t", index=False)

# Combine gene counts from all RNA-seq samples
