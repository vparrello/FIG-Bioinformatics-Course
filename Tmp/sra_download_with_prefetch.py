########################################################################
#...Prompt that generated this program:
"""
Write a Python script named "sra_download.py" that automates
the download of paired-end SRA data and converts it directly
into FASTA format without using external tools like seqtk.
The script should be fully functional and ready to run
on a system where the SRA Toolkit is installed.

The script must support the following command-line arguments:
- "-i" or "--id-file": A required argument specifying a file
   that contains a list of SRA IDs, one per line.

- "-d" or "--download-directory": A required argument specifying
   the directory where the processed FASTA files should be stored.

- "-c" or "--sra-cache": An optional argument specifying a directory
   to be used as a cache location for SRA downloads.

- "-D" or "--Debug": An optional flag to enable debug output.

### Script Behavior:
1. Read the list of SRA IDs from the specified file.
2. Use the "prefetch" command to download the corresponding SRA files
   into the cache directory (or the default "~/.ncbi/public/sra/" location
   if no cache directory is provided).
3. Ensure that the downloaded ".sra" file exists in the expected location
   before proceeding.
4. Use "fastq-dump --fasta 0 --split-files" to convert the downloaded
   SRA file into a pair of FASTA files (for paired-end reads)
   or a single FASTA file (for single-end reads).
5. Store the resulting FASTA files in a subdirectory under the specified
   download directory, named after the corresponding SRA ID.
6. Perform error handling at every step:
   - If "prefetch" fails, print a warning and move to the next SRA ID.
   - If "fastq-dump" fails, print a warning and move to the next SRA ID.
   - If no FASTA files are created, print a warning.
7. Ensure all informational messages (`[INFO]`) and debugging messages
   (`[DEBUG]`) are printed to STDERR.
8. Debug messages should only be printed if the script is invoked with the
   `-D` or `--Debug` flag.
9. After processing all SRA IDs, print a summary that includes:
   - The total number of SRA entries requested.
   - The number of successful conversions.
   - The count of single-end and paired-end datasets.

### Implementation Requirements:
- Use "subprocess.run()" to execute system commands.
- Ensure all output directories exist before writing files.
- Capture and log STDERR output for debugging.
- Use readable function names and modular code organization.
- Print all informational messages and warnings to STDERR.
"""

########################################################################
#...Pseudocode for this program:
"""
FUNCTION log(message, debug)
    IF debug
        PRINT "[DEBUG] " + message TO STDERR
    ELSE
        PRINT "[INFO] " + message TO STDERR
    ENDIF
END FUNCTION


FUNCTION ensure_directory_exists(directory)
    IF directory DOES NOT EXIST
        CREATE directory
    ENDIF
END FUNCTION


FUNCTION run_command(command, debug)
    TRY
        EXECUTE command, CAPTURE STDOUT AND STDERR
        IF debug
            log("Command: " + command, TRUE)
            log("STDOUT: " + STDOUT, TRUE)
            log("STDERR: " + STDERR, TRUE)
        ENDIF
        IF execution FAILED
            log("Error executing: " + command + "\n" + STDERR, FALSE)
            RETURN FALSE
        ENDIF
        RETURN TRUE
    CATCH exception
        log("Exception running command: " + command + ", Error: " + exception, FALSE)
        RETURN FALSE
    END TRY
END FUNCTION


FUNCTION process_sra_id(sra_id, download_dir, sra_cache, debug)
    sra_output_dir = download_dir + "/" + sra_id
    ensure_directory_exists(sra_output_dir)

    prefetch_cmd = "prefetch " + sra_id
    IF sra_cache IS NOT EMPTY
        prefetch_cmd = prefetch_cmd + " --output-directory " + sra_cache
    ENDIF

    IF NOT run_command(prefetch_cmd, debug)
        log("Failed to prefetch " + sra_id, FALSE)
        RETURN FALSE
    ENDIF

    sra_file = (sra_cache IF sra_cache ELSE DEFAULT_SRA_PATH) + "/" + sra_id + ".sra"
    IF sra_file DOES NOT EXIST
        log("SRA file " + sra_file + " not found after prefetch", FALSE)
        RETURN FALSE
    ENDIF

    fastq_dump_cmd = "fastq-dump --fasta 0 --split-files " + sra_file + " --outdir " + sra_output_dir
    IF NOT run_command(fastq_dump_cmd, debug)
        log("Failed to convert " + sra_id + " to FASTA", FALSE)
        RETURN FALSE
    ENDIF

    fasta_files = LIST FILES in sra_output_dir WHERE file EXTENSION IS ".fasta"
    IF fasta_files IS EMPTY
        log("No FASTA files were created for " + sra_id, FALSE)
        RETURN FALSE
    ENDIF

    log("Successfully processed " + sra_id + ": " + LENGTH(fasta_files) + " FASTA files created.", FALSE)
    RETURN TRUE
END FUNCTION


FUNCTION main()
    PARSE COMMAND-LINE ARGUMENTS:
        - id_file (REQUIRED): File containing SRA IDs
        - download_directory (REQUIRED): Directory for FASTA storage
        - sra_cache (OPTIONAL): Cache directory for SRA downloads
        - Debug (OPTIONAL FLAG): Enable debug output

    TRY
        READ LINES FROM id_file INTO sra_ids, REMOVE WHITESPACE
    CATCH FILE NOT FOUND
        log("Error: SRA ID file " + id_file + " not found.", FALSE)
        EXIT(1)
    END TRY

    log("Total SRA IDs to process: " + LENGTH(sra_ids), Debug)

    successful_conversions = 0
    single_end_count = 0
    paired_end_count = 0

    FOR EACH sra_id IN sra_ids
        success = process_sra_id(sra_id, download_directory, sra_cache, Debug)
        IF success
            successful_conversions += 1
            fasta_files = LIST FILES in download_directory + "/" + sra_id WHERE file EXTENSION IS ".fasta"
            IF LENGTH(fasta_files) == 1
                single_end_count += 1
            ELSE
                paired_end_count += 1
            ENDIF
        ENDIF
    ENDFOR

    log("\nProcessing Summary:", FALSE)
    log("Total SRA IDs: " + LENGTH(sra_ids), FALSE)
    log("Successful Conversions: " + successful_conversions, FALSE)
    log("Single-end datasets: " + single_end_count, FALSE)
    log("Paired-end datasets: " + paired_end_count, FALSE)

    IF successful_conversions == 0
        log("No successful conversions. Check logs for details.", FALSE)
        EXIT(1)
    ENDIF
END FUNCTION


IF EXECUTED AS MAIN PROGRAM
    main()
ENDIF
"""

########################################################################
#...Code generated by Grimoire:
import os
import sys
import argparse
import subprocess


def log(message, debug=False):
    """Prints messages to STDERR, handles debug flag."""
    if debug:
        print(f"[DEBUG] {message}", file=sys.stderr)
    else:
        print(f"[INFO] {message}", file=sys.stderr)


def ensure_directory_exists(directory):
    """Ensure that a directory exists."""
    if not os.path.exists(directory):
        os.makedirs(directory)


def run_command(command, debug=False):
    """Executes a shell command and captures output."""
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if debug:
            log(f"Command: {command}", debug=True)
            log(f"STDOUT: {result.stdout}", debug=True)
            log(f"STDERR: {result.stderr}", debug=True)
        if result.returncode != 0:
            log(f"Error executing: {command}\n{result.stderr}", debug=False)
            return False
        return True
    except Exception as e:
        log(f"Exception running command: {command}, Error: {e}", debug=False)
        return False


def process_sra_id(sra_id, download_dir, sra_cache, debug):
    """Download and convert SRA file to FASTA."""
    sra_output_dir = os.path.join(download_dir, sra_id)
    ensure_directory_exists(sra_output_dir)

    # Prefetch command
    prefetch_cmd = f"prefetch {sra_id}"
    if sra_cache:
        prefetch_cmd += f" --output-directory {sra_cache}"

    if not run_command(prefetch_cmd, debug):
        log(f"Failed to prefetch {sra_id}", debug=False)
        return False

    # Determine downloaded SRA file path
    sra_file = os.path.join(sra_cache if sra_cache else os.path.expanduser("~/.ncbi/public/sra"), f"{sra_id}.sra")
    if not os.path.exists(sra_file):
        log(f"SRA file {sra_file} not found after prefetch", debug=False)
        return False

    # Convert SRA to FASTA
    fastq_dump_cmd = f"fastq-dump --fasta 0 --split-files {sra_file} --outdir {sra_output_dir}"
    if not run_command(fastq_dump_cmd, debug):
        log(f"Failed to convert {sra_id} to FASTA", debug=False)
        return False

    # Check for output files
    fasta_files = [f for f in os.listdir(sra_output_dir) if f.endswith(".fasta")]
    if not fasta_files:
        log(f"No FASTA files were created for {sra_id}", debug=False)
        return False

    log(f"Successfully processed {sra_id}: {len(fasta_files)} FASTA files created.", debug=False)
    return True


def main():
    parser = argparse.ArgumentParser(description="Download and convert SRA files to FASTA.")
    parser.add_argument("-i", "--id-file", required=True, help="File containing list of SRA IDs.")
    parser.add_argument("-d", "--download-directory", required=True, help="Directory to store FASTA files.")
    parser.add_argument("-c", "--sra-cache", help="Optional cache directory for SRA downloads.")
    parser.add_argument("-D", "--Debug", action="store_true", help="Enable debug output.")
    args = parser.parse_args()

    # Read SRA IDs
    try:
        with open(args.id_file, 'r') as file:
            sra_ids = [line.strip() for line in file if line.strip()]
    except FileNotFoundError:
        log(f"Error: SRA ID file {args.id_file} not found.", debug=False)
        sys.exit(1)

    log(f"Total SRA IDs to process: {len(sra_ids)}", debug=args.Debug)

    successful_conversions = 0
    single_end_count = 0
    paired_end_count = 0

    for sra_id in sra_ids:
        success = process_sra_id(sra_id, args.download_directory, args.sra_cache, args.Debug)
        if success:
            successful_conversions += 1
            fasta_files = [f for f in os.listdir(os.path.join(args.download_directory, sra_id)) if f.endswith(".fasta")]
            if len(fasta_files) == 1:
                single_end_count += 1
            else:
                paired_end_count += 1

    log("\nProcessing Summary:", debug=False)
    log(f"Total SRA IDs: {len(sra_ids)}", debug=False)
    log(f"Successful Conversions: {successful_conversions}", debug=False)
    log(f"Single-end datasets: {single_end_count}", debug=False)
    log(f"Paired-end datasets: {paired_end_count}", debug=False)

    if successful_conversions == 0:
        log("No successful conversions. Check logs for details.", debug=False)
        sys.exit(1)

if __name__ == "__main__":
    main()
