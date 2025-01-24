########################################################################
#...Prompt that generated this program:
"""
Please write a program named 'hammer_compare_multirole.py'
that meets the following requirements. The program must:

* Accept two mandatory arguments and one optional argument, as follows:

-H or --hammers: a TSV-filename argument (mandatory).
-G or --genome-names: a TSV-filename argument (mandatory).
-F or --role-fraction: an optional argument with a default value of 0.8.

* When run with the --help argument, the program must use the following
as its "description" line:
    "program for multirole hammer analysis. Reads DNA FASTA from STDIN. Writes TSV report to STDOUT."

The program should:

* Open the hammer-file and perform the following:
    
    - Read the header-line and extract the column-names
    
    - Foreach subsequent line, read the "hammer" and "fid"
    columns as 'hammer' and 'feature_id' variables,
    respectively, and built a dictionary mapping
    hammers to feature_ids.

    - If a "role" column exists, read it as the 'role' variable,
    else the role should default to "Unknown",
    and build a dictionary mapping feature_ids to roles.

    - A 'feature_id' has the format 'fig|x.y.peg.z',
    where 'x', 'y', and 'z' are integers, and 'fig|' and '.peg.'
    are literal substrings; the portion 'x.y' is the 'genome_id'
    for that 'feature_id'. Extract the 'genome_id' from the 'feature_id'
    using a regex, and build a dictionary mapping feature_id to its genome_id.

    - Determine the Kmer-length 'K' of the hammers.
    (NOTE: all the hammers in the file will have the same length.)

    - Return the hammer-length 'K', the list of roles,
    the hammer-to-feature_id dictionary,
    the feature_id-to-genome_id dictionary,
    and the feature_id-to-role dictionary

* Open the 'genome-names' file and perform the following: 

    - Read the header-line and extract the column-names

    - For subsequent lines, read the 'genome_id' and 'genome_name' columns
    into the genome_id-to-genome_names dictionary.

* Create an empty dictionary whose keys will be genome-IDs,
and values will be an integer count of the number of hammers
found for that genome.

* Create an empty "dictionary of dictionaries" whose outer-keys will be genome_ids,
inner-keys will be roles, and inner-values will be the integer count
for the number of hammers found for that role in that genome.

* Use BioPython to read the sequences of the genome from 'STDIN';
the sequence should then be converted to lower-case.

* For each sequence, extract all possible Kmers, and if a Kmer
is a hammer, increment the hammer-count for its associated 'genome_id'
in the genome-counts dictionary, and also the hammer-count for that role
in the genome-to-role-to-counts directory-of-directories;
then repeat this operation on the reverse-complement of that sequence,
since a gene can face in either direction.

* Once all the sequences have been processed,
foreach genome_id in the genome-to-roles-to-counts directory
of directories, if the number of roles for that genome
is greater than or equal to the total number of roles
times the minimum fraction of roles,
print to STDOUT a TSV file of the genome_ids found
and their associated genome_name and score,
sorted by decreasing genome-to-counts score.
If the total number of roles exceeds 1,
also add a column for the number of roles for the genome.
Please handle missing genome-names gracefully;
if a genome_id does not have an associated genome-name,
display the genome_name as 'Unknown sp.' in the output TSV file,
and send a warning to STDERR that the name of genome_id was not in
the genome-names file.
"""

########################################################################
#...Pseudocode for this program:
"""
BEGIN
    PROCEDURE ParseHammers(hammersFile);
        DECLARE HammerToFeature, FeatureToGenome, FeatureToRole DICTIONARIES;
        DECLARE Roles SET;
        DECLARE KmerLength INTEGER;
        OPEN hammersFile FOR READING;
        READ HeaderLine;
        WHILE NOT EndOfFile DO
            READ Line;
            LET Hammer := Line["hammer"];
            LET FeatureID := Line["fid"];
            LET Role := IF "role" IN Line THEN Line["role"] ELSE "Unknown";
            ADD Hammer -> FeatureID TO HammerToFeature;
            ADD FeatureID -> Role TO FeatureToRole;
            ADD Role TO Roles;
            IF FeatureID MATCHES "fig|x.y.peg.z" THEN
                LET GenomeID := EXTRACT GenomeID FROM FeatureID;
                ADD FeatureID -> GenomeID TO FeatureToGenome;
            ENDIF;
        ENDWHILE;
        CLOSE hammersFile;
        LET KmerLength := LENGTH OF ANY_KEY IN HammerToFeature;
        RETURN KmerLength, Roles, HammerToFeature, FeatureToGenome, FeatureToRole;
    END;

    PROCEDURE ParseGenomeNames(genomeNamesFile);
        DECLARE GenomeIDToName DICTIONARY;
        OPEN genomeNamesFile FOR READING;
        READ HeaderLine;
        WHILE NOT EndOfFile DO
            READ Line;
            LET GenomeID := Line["genome_id"];
            LET GenomeName := Line["genome_name"];
            ADD GenomeID -> GenomeName TO GenomeIDToName;
        ENDWHILE;
        CLOSE genomeNamesFile;
        RETURN GenomeIDToName;
    END;

    PROCEDURE Main();
        DECLARE Args RECORD;
        DECLARE KmerLength, Roles, HammerToFeature, FeatureToGenome, FeatureToRole;
        DECLARE GenomeIDToName DICTIONARY;
        DECLARE GenomeCounts, GenomeRoleCounts DICTIONARY;
        DECLARE TotalRoles, RoleFraction REAL;
        PARSE Args FROM CommandLine;
        TRY
            KmerLength, Roles, HammerToFeature, FeatureToGenome, FeatureToRole := ParseHammers(Args.Hammers);
            GenomeIDToName := ParseGenomeNames(Args.GenomeNames);
            LET GenomeCounts := EMPTY_DICTIONARY;
            LET GenomeRoleCounts := EMPTY_DICTIONARY;

            FOR EACH Sequence IN ReadSequencesFromSTDIN() DO
                LET SequenceLower := LOWERCASE(Sequence);
                LET RevCompSequence := ReverseComplement(SequenceLower);

                FOR I := 0 TO LENGTH(SequenceLower) - KmerLength DO
                    LET Kmer := SUBSTRING(SequenceLower, I, KmerLength);
                    IF Kmer IN HammerToFeature THEN
                        LET FeatureID := HammerToFeature[Kmer];
                        LET GenomeID := FeatureToGenome[FeatureID];
                        LET Role := FeatureToRole[FeatureID];
                        INCREMENT GenomeCounts[GenomeID];
                        INCREMENT GenomeRoleCounts[GenomeID][Role];
                    ENDIF;

                    LET KmerRev := SUBSTRING(RevCompSequence, I, KmerLength);
                    IF KmerRev IN HammerToFeature THEN
                        LET FeatureID := HammerToFeature[KmerRev];
                        LET GenomeID := FeatureToGenome[FeatureID];
                        LET Role := FeatureToRole[FeatureID];
                        INCREMENT GenomeCounts[GenomeID];
                        INCREMENT GenomeRoleCounts[GenomeID][Role];
                    ENDIF;
                ENDFOR;
            ENDFOR;

            LET TotalRoles := LENGTH(Roles);
            LET RoleFraction := Args.RoleFraction;
            DECLARE Output ARRAY;

            FOR EACH GenomeID IN GenomeRoleCounts DO
                LET NumRoles := LENGTH(GenomeRoleCounts[GenomeID]);
                IF NumRoles >= TotalRoles * RoleFraction THEN
                    LET GenomeName := IF GenomeID IN GenomeIDToName THEN GenomeIDToName[GenomeID] ELSE "Unknown sp.";
                    IF GenomeName = "Unknown sp." THEN
                        PRINT "Warning: GenomeID {GenomeID} not found in genome-names file." TO STDERR;
                    ENDIF;
                    LET Score := GenomeCounts[GenomeID];
                    ADD (GenomeID, GenomeName, Score, NumRoles) TO Output;
                ENDIF;
            ENDFOR;

            SORT Output BY Score DESCENDING;
            PRINT Headers ("genome_id", "genome_name", "score", "num_roles" IF TotalRoles > 1);
            FOR EACH Row IN Output DO
                PRINT Row TO STDOUT;
            ENDFOR;

        CATCH Exception DO
            PRINT "Error: {Exception.Message}" TO STDERR;
            EXIT 1;
        END;
    END;

    CALL Main();
END;
"""

########################################################################
#...Code generated by Grimoire:
import argparse
import sys
import re
import csv
from Bio import SeqIO
from Bio.Seq import reverse_complement

def parse_hammers(hammers_file):
    hammer_to_feature = {}
    feature_to_genome = {}
    feature_to_role = {}
    roles = set()

    with open(hammers_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        headers = reader.fieldnames
        
        if 'hammer' not in headers or 'fid' not in headers:
            raise ValueError("The hammers file must contain 'hammer' and 'fid' columns.")

        for row in reader:
            hammer = row['hammer']
            feature_id = row['fid']
            role = row.get('role', 'Unknown')

            hammer_to_feature[hammer] = feature_id
            feature_to_role[feature_id] = role
            roles.add(role)

            genome_match = re.match(r'^fig\|(\d+\.\d+)\.peg\.\d+$', feature_id)
            if genome_match:
                genome_id = genome_match.group(1)
                feature_to_genome[feature_id] = genome_id

    kmer_length = len(next(iter(hammer_to_feature)))
    return kmer_length, roles, hammer_to_feature, feature_to_genome, feature_to_role

def parse_genome_names(genome_names_file):
    genome_id_to_name = {}
    with open(genome_names_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        headers = reader.fieldnames

        if 'genome_id' not in headers or 'genome_name' not in headers:
            raise ValueError("The genome-names file must contain 'genome_id' and 'genome_name' columns.")

        for row in reader:
            genome_id_to_name[row['genome_id']] = row['genome_name']

    return genome_id_to_name

def main():
    parser = argparse.ArgumentParser(description="program for multirole hammer analysis. Reads DNA FASTA from STDIN. Writes TSV report to STDOUT.")
    parser.add_argument('-H', '--hammers', required=True, help="Path to the hammers TSV file.")
    parser.add_argument('-G', '--genome-names', required=True, help="Path to the genome names TSV file.")
    parser.add_argument('-F', '--role-fraction', type=float, default=0.8, help="Minimum fraction of roles for genome inclusion (default: 0.8).")

    args = parser.parse_args()

    try:
        kmer_length, roles, hammer_to_feature, feature_to_genome, feature_to_role = parse_hammers(args.hammers)
        genome_id_to_name = parse_genome_names(args.genome_names)

        genome_counts = {}
        genome_role_counts = {}

        for record in SeqIO.parse(sys.stdin, 'fasta'):
            sequence = str(record.seq).lower()
            rev_comp_sequence = reverse_complement(sequence)

            for i in range(len(sequence) - kmer_length + 1):
                kmer = sequence[i:i + kmer_length]
                if kmer in hammer_to_feature:
                    feature_id = hammer_to_feature[kmer]
                    genome_id = feature_to_genome[feature_id]
                    role = feature_to_role[feature_id]

                    genome_counts[genome_id] = genome_counts.get(genome_id, 0) + 1

                    if genome_id not in genome_role_counts:
                        genome_role_counts[genome_id] = {}
                    genome_role_counts[genome_id][role] = genome_role_counts[genome_id].get(role, 0) + 1

                kmer = rev_comp_sequence[i:i + kmer_length]
                if kmer in hammer_to_feature:
                    feature_id = hammer_to_feature[kmer]
                    genome_id = feature_to_genome[feature_id]
                    role = feature_to_role[feature_id]

                    genome_counts[genome_id] = genome_counts.get(genome_id, 0) + 1

                    if genome_id not in genome_role_counts:
                        genome_role_counts[genome_id] = {}
                    genome_role_counts[genome_id][role] = genome_role_counts[genome_id].get(role, 0) + 1

        total_roles = len(roles)
        output = []

        for genome_id, role_counts in genome_role_counts.items():
            num_roles = len(role_counts)
            if num_roles >= total_roles * args.role_fraction:
                genome_name = genome_id_to_name.get(genome_id, 'Unknown sp.')
                if genome_name == 'Unknown sp.':
                    print(f"Warning: genome_id {genome_id} not found in genome-names file.", file=sys.stderr)

                score = genome_counts[genome_id]
                output.append((genome_id, genome_name, score, num_roles))

        output.sort(key=lambda x: x[2], reverse=True)

        writer = csv.writer(sys.stdout, delimiter='\t')
        headers = ['genome_id', 'genome_name', 'score']
        if total_roles > 1:
            headers.append('num_roles')
        writer.writerow(headers)

        for row in output:
            writer.writerow(row)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
