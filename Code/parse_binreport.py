from bs4 import BeautifulSoup
import os
import math

# TODO run this on a report to make sure it still works.

# This is a beautiful soup program that parses a standard binning report downloaded as an html file
HealthyOrDiseased = "Healthy"

# Function to load representative genomes from a JSON file
def load_representative_genomes(file_path):
    representative_genomes = {}
    with open(file_path, 'r') as f:
        for line in f:
            genome_id, representative = line.strip().split('\t')
            representative_genomes[genome_id] = representative
    return representative_genomes

# Main function to parse binning reports
def parse_reports(input_path, representative_genomes):
    samples_to_genomes = {}
    bad_coverage_genomes = {}
    
    # Iterate through those files
    for report in os.listdir(input_path):
        full_path = os.path.join(input_path, report)
        with open(full_path) as fp:
            #make it into soup
            soup = BeautifulSoup(fp, 'html.parser')
            all_paragraph_tags = soup.find_all('p')
            for text in all_paragraph_tags:
                if "FIGCorePr" not in text.text:
                    continue
                else:
                    sample_text = text.text
            sample_text_list = sample_text.split("/")
            sample = sample_text_list[-2].lstrip(".").strip("\n")
            samples_to_genomes[sample] = []
            bad_coverage_genomes[sample] = []
            tables = soup.find_all("table")
            for table in tables:
                cells = soup.find_all("td")
                rows = math.floor((len(cells))/15)
                for i in range(rows):
                    genome = {'genome_id': '', 'dnasize': 0, 'coverage': 0}
                    row = i * 15
                    genome_raw = cells[row+3].text.strip("\n")
                    if "\n" in genome_raw:
                        genome_raw = genome_raw.split("\n")
                        genome_raw = genome_raw[0]
                    genome['genome_id'] = genome_raw
                    genome['dnasize'] = int(cells[row+10].text)
                    genome['coverage'] = float(cells[row+12].text)
                    if genome['coverage'] == 50 or genome['coverage'] < 1:
                        bad_coverage_genomes[sample].append(genome)
                    else:
                        samples_to_genomes[sample].append(genome)
                    # Map genome_id to representative genomes
                    if genome['genome_id'] in representative_genomes:
                        # Logic to handle mapping
                        pass  # Implement mapping logic here

    sample_to_coverage = []
    with open(f"{HealthyOrDiseased}_binReport_scrapped.tbl", "w+") as output:
        output.write(f"Sample\tReference Genome\tDNA Size\tMean Coverage\n")
        for sample in samples_to_genomes:
            if len(samples_to_genomes[sample]) >= 1:
                numerator_sum = 0
                denominator_sum = 0
                for list in samples_to_genomes[sample]:
                    if len(list) == 3:
                        numerator_sum += (list['dnasize'] * list['coverage'])
                        denominator_sum += list['dnasize']
                        output.write(f"{sample}\t{list['genome_id']}\t{list['dnasize']}\t{list['coverage']}\n")
                sample_mean_coverage = round((numerator_sum / denominator_sum), 2)
                sample_to_coverage.append([sample, sample_mean_coverage])

    with open(f"{HealthyOrDiseased}_binReport_coverage.tbl", "w+") as cov_output:
        cov_output.write(f"Sample\tMean Coverage For Entire Sample\n")
        for sample in sample_to_coverage:
            cov_output.write(f"{sample[0]}\t{sample[1]}\n")

    with open(f"{HealthyOrDiseased}_bad_coverage_genomes.tbl", "w+") as bad_out:
        bad_out.write(f"Sample\tReference Genome\tDNA Size\tMean Coverage\n")
        for sample in bad_coverage_genomes:
            for list in bad_coverage_genomes[sample]:
                bad_out.write(f"{sample}\t{list['genome_id']}\t{list['dnasize']}\t{list['coverage']}\n")

# Example usage
if __name__ == "__main__":
    input_path = input("Enter the path to the bin reports: ")
    representative_genomes_file = input("Enter the path to the representative genomes JSON file: ")
    representative_genomes = load_representative_genomes(representative_genomes_file)
    parse_reports(input_path, representative_genomes)
