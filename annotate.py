#!/usr/bin/env python3
"""
annotate.py
This script reads a joint VCF file (the one with all patients from one clinic)
and creates a simple CSV table for each patient.
The CSV has fake medical info + MSUD symptom flags.
Easy to use for machine learning (like Elastic Net).
Now also adds control patients if you use --control_no.
"""

import argparse
import csv
import os
import random
from collections import defaultdict

def main():
    # 1. Read command line: user gives the path to the VCF file
    p = argparse.ArgumentParser(description="Annotate joint VCF for MSUD symptoms")
    p.add_argument("--path_to_vcf", required=True, help="Path to your joint VCF file (e.g. klinik1.vcf)")
    p.add_argument("--control_no", type=int, default=0, help="Number of fake healthy control patients to add (default: 0)")
    args = p.parse_args()
    vcf = args.path_to_vcf

    # 2. Count how many 1/1 variants each patient has
    patients = []
    counts = defaultdict(int)
    with open(vcf) as f:
        for line in f:
            if line.startswith("##"):      # skip header lines that start with ##
                continue
            if line.startswith("#CHROM"):  # this line has the patient names
                patients = line.strip().split("\t")[9:]
                continue
            fields = line.strip().split("\t")
            for i, pat in enumerate(patients, 9):
                gt = fields[i].split(":")[0]   # get genotype (1/1 = has variant)
                if gt == "1/1":
                    counts[pat] += 1

    # 3. Create output CSV in the same folder
    out = os.path.join(os.path.dirname(vcf) or ".", 
                      os.path.basename(vcf).replace(".vcf", "_annotated.csv"))

    # MSUD-related symptoms we care about
    symptoms = ["poor_feeding", "vomiting", "lethargy", "maple_syrup_odor", 
                "irritability", "seizures", "developmental_delay"]

    random.seed(42)  # so every run gives same fake data (good for testing)

    # 4. Write CSV with one row per patient
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, ["patient_id", "age", "sex", "variant_count"] + symptoms + ["msud_risk"])
        w.writeheader()

        # Real patients from the VCF
        for pat in patients:
            c = counts[pat]                    # how many variants this patient carries
            row = {
                "patient_id": pat,
                "age": random.randint(0, 60),   # fake age in years
                "sex": random.choice(["M", "F"]), # fake sex
                "variant_count": c,
                "msud_risk": "high" if c > 10 else "medium" if c > 5 else "low",
            }
            # fake symptoms (more variants = higher chance of symptoms)
            for s in symptoms:
                row[s] = 1 if random.random() < min(1.0, c / 20) else 0
            w.writerow(row)

        # Add control patients (healthy, zero variants)
        for i in range(1, args.control_no + 1):
            row = {
                "patient_id": f"control{i}",
                "age": random.randint(0, 60),
                "sex": random.choice(["M", "F"]),
                "variant_count": 0,
                "msud_risk": "low",
            }
            for s in symptoms:
                row[s] = 0   # controls have no symptoms
            w.writerow(row)

    print(f"done → {out}")   # tells you where the CSV was saved

if __name__ == "__main__":
    main()