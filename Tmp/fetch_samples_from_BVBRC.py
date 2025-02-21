import sys
import subprocess

def fetch_samples(base_path, target_dir):
    """
    Reads SRA-IDs from STDIN and fetches data using p3-cp commands.
    """
    for line in sys.stdin:
        sra_id = line.strip()
        if not sra_id:
            continue
        
        print(f"Fetching data for SRA-ID: {sra_id}", file=sys.stderr)
        
        file1_path = f"ws:{base_path}/.{sra_id}/{sra_id}_1_ptrim.fq.gz"
        file2_path = f"ws:{base_path}/.{sra_id}/{sra_id}_2_ptrim.fq.gz"
        
        try:
            subprocess.run(["p3-cp", file1_path, target_dir], check=True)
        except subprocess.CalledProcessError:
            print(f"Warning: Failed to fetch {file1_path}", file=sys.stderr)
        
        try:
            subprocess.run(["p3-cp", file2_path, target_dir], check=True)
        except subprocess.CalledProcessError:
            print(f"Warning: Failed to fetch {file2_path}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fetch_samples.py <base_path> <target_dir>")
        sys.exit(1)
    
    base_path = sys.argv[1]
    target_dir = sys.argv[2]
    fetch_samples(base_path, target_dir)
