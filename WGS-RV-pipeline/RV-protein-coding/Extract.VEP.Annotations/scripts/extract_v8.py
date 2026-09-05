#!/usr/bin/env python3

"""
Extract transcript-level VEP annotations from GEL AggV3 functional annotation VCFs.

Design:
    - Input VCF is already restricted to a genomic subshard.
    - One output row corresponds to one VEP transcript consequence.
    - Designed for Nextflow execution on CloudOS.
    - Supports local paths and S3 paths through bcftools.

Author:
    Shiyu Zhang

Version:
    v1.0.0
"""

import argparse
import csv
import subprocess
import sys
import time
from pathlib import Path


##############################################
# VEP CSQ fields to retain
##############################################

INFO_FIELDS = [
    "AC",
    "AN",
    "AF",
]


CSQ_FIELDS_KEEP = [

    "Allele",
    "Consequence",
    "IMPACT",

    "SYMBOL",
    "Gene",

    "Feature_type",
    "Feature",
    "BIOTYPE",

    "EXON",
    "INTRON",

    "HGVSc",
    "HGVSp",

    "CANONICAL",

    "MANE",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",

    "CADD_PHRED",
    "CADD_RAW",

    "REVEL",

    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",

    "LoF",
    "LoF_filter",
    "LoF_flags",
    "LoF_info",

    "am_class",
    "am_pathogenicity",

    "MechPredict_pDN",
    "MechPredict_pGOF",
    "MechPredict_pLOF",
    "MechPredict_prediction",

    "GERP",
    "PhyloP",
]


##############################################
# Arguments
##############################################

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Extract transcript-level VEP annotations "
            "from GEL AggV3 functional annotation VCF."
        )
    )

    parser.add_argument(
        "--vcf",
        required=True,
        help="Input VCF.gz file (local or S3 path)"
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV file"
    )

    return parser.parse_args()



##############################################
# Read CSQ header
##############################################

def get_csq_fields(vcf):

    print(
        "[INFO] Reading VCF header...",
        flush=True
    )


    cmd = [
        "bcftools",
        "view",
        "-h",
        vcf
    ]


    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        text=True
    )


    csq_fields = None


    for line in p.stdout:

        if line.startswith("##INFO=<ID=CSQ"):

            if "Format: " not in line:
                continue


            fmt = line.split(
                "Format: ",
                1
            )[1]


            fmt = fmt.replace(
                '">',
                ""
            )

            fmt = fmt.strip()


            csq_fields = fmt.split("|")

            break


    p.stdout.close()


    if csq_fields is None:

        raise RuntimeError(
            "Cannot find CSQ annotation header"
        )


    print(
        f"[INFO] Detected CSQ fields: {len(csq_fields)}",
        flush=True
    )


    return csq_fields



##############################################
# Main extraction
##############################################

def extract(args):

    vcf = args.vcf
    out = args.out


    csq_fields = get_csq_fields(vcf)


    csq_index = {
        field: idx
        for idx, field in enumerate(csq_fields)
    }


    missing = [
        x
        for x in CSQ_FIELDS_KEEP
        if x not in csq_index
    ]


    if missing:

        print(
            "[WARNING] Missing CSQ fields:",
            ",".join(missing),
            flush=True
        )


    print(
        "[INFO] Starting extraction...",
        flush=True
    )


    cmd = [

        "bcftools",
        "query",

        "-f",

        "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n",

        vcf

    ]


    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        text=True
    )


    output_path = Path(out)

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )


    header = [

        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER"

    ]


    header.extend(INFO_FIELDS)

    header.extend(CSQ_FIELDS_KEEP)



    n_variant = 0
    n_transcript = 0


    start_time = time.time()



    with open(
        output_path,
        "w",
        newline=""
    ) as fout:


        writer = csv.writer(
            fout,
            delimiter="\t",
            lineterminator="\n"
        )


        writer.writerow(header)



        for line in process.stdout:


            line = line.rstrip("\n")


            if not line:
                continue



            cols = line.split("\t")


            chrom, pos, vid, ref, alt, qual, filt, info = cols



            info_dict = {}


            for item in info.split(";"):

                if "=" in item:

                    key, value = item.split(
                        "=",
                        1
                    )

                    info_dict[key] = value



            base = [

                chrom,
                pos,
                vid,
                ref,
                alt,
                qual,
                filt

            ]



            for field in INFO_FIELDS:

                base.append(
                    info_dict.get(
                        field,
                        "NA"
                    )
                )



            csq = info_dict.get(
                "CSQ",
                ""
            )


            if csq == "":

                writer.writerow(
                    base +
                    ["NA"] * len(CSQ_FIELDS_KEEP)
                )


            else:

                transcripts = csq.split(",")


                for transcript in transcripts:


                    values = transcript.split("|")


                    row = []


                    for field in CSQ_FIELDS_KEEP:


                        idx = csq_index.get(
                            field
                        )


                        if idx is None:

                            row.append("NA")


                        elif idx >= len(values):

                            row.append("NA")


                        else:

                            value = values[idx]

                            if value == "":
                                value = "NA"

                            row.append(value)



                    writer.writerow(
                        base + row
                    )


                    n_transcript += 1



            n_variant += 1



            if n_variant % 10000 == 0:


                elapsed = (
                    time.time()
                    -
                    start_time
                )


                print(

                    f"[INFO] "
                    f"Variants: {n_variant:,} | "
                    f"Transcripts: {n_transcript:,} | "
                    f"Elapsed: {elapsed/60:.1f} min",

                    flush=True

                )



    process.stdout.close()


    return_code = process.wait()


    if return_code != 0:

        raise RuntimeError(
            "bcftools query failed"
        )



    print(
        "\n[INFO] Extraction finished",
        flush=True
    )


    print(
        f"[INFO] Variants: {n_variant:,}",
        flush=True
    )


    print(
        f"[INFO] Transcripts: {n_transcript:,}",
        flush=True
    )


    print(
        f"[INFO] Output: {out}",
        flush=True
    )



##############################################
# Entry point
##############################################

def main():

    args = parse_args()

    extract(args)



if __name__ == "__main__":

    main()
