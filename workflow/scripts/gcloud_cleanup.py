from google.cloud import storage
import sys
from os.path import join
from pprint import pprint


def cleanup(bucket_name: str, org: str, ref: str):
    dirs_to_kill = [
        join("results", org, ref, x)
        for x in [
            "00_fastqFiltered",
            "01_mappedReads/preMerge",
            "01_mappedReads/postMerge",
        ]
    ]
    genome = join("data", "genome", ref, "*")
    fastq = join("data", "fastq", org, "*")
    dirs_to_kill.extend([genome, fastq])

    cls = storage.Client()
    bucket = cls.get_bucket(bucket_name)

    blob_list = []
    for d in dirs_to_kill:
        blobs = bucket.list_blobs(prefix=d)
        for blob in blobs:
            print(f"Would delete file: {blob}")
            blob_list.append(blob)
    
    outputname = join(bucket, "results", org, ref, "cleanup_confirmation.txt")

    bucket.delete_blobs(blob_list)

    with open(outputname, "w+") as writer:
        writer.write("cleanup script ran successfully")

def main():
    organism = snakemake.input[0].split("/")[2]
    ref = snakemake.input[0].split("/")[3]
    bucket = snakemake.input[0].split("/")[0]
    cleanup(bucket_name=bucket, org=organism, ref=ref)


if __name__ == "__main__":
    main()
