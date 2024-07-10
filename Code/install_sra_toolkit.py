import subprocess
import sys
import os
from pathlib import Path

def install_sra_toolkit_user(top_dir):
    """Install the SRA Toolkit and EDirect in user mode if they are not already installed."""
    bin_dir = Path(top_dir)
    if not bin_dir.is_dir():
        print(f"Directory '{top_dir}' does not exist", file=sys.stderr)
        exit(1)
    
    sra_toolkit_dir = bin_dir / "sratoolkit"
    fastq_dump_path = sra_toolkit_dir / "bin" / "fastq-dump"
    edirect_dir = bin_dir / "edirect"
    esearch_path = edirect_dir / "esearch"

    if fastq_dump_path.exists():
        print("SRA Toolkit is already installed.")
    else:
        print("Installing SRA Toolkit in user mode...")
        if sys.platform.startswith("linux"):
            sra_toolkit_url = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz"
        elif sys.platform == "darwin":
            sra_toolkit_url = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz"
        elif sys.platform == "win32":
            sra_toolkit_url = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-win64.zip"
        else:
            print("Unsupported OS. Please install SRA Toolkit manually.")
            sys.exit(1)

        sra_toolkit_archive = sra_toolkit_dir / "sratoolkit.tar.gz" if sys.platform != "win32" else sra_toolkit_dir / "sratoolkit.zip"
        bin_dir.mkdir(parents=True, exist_ok=True)
        sra_toolkit_dir.mkdir(parents=True, exist_ok=True)

        # Download the toolkit
        subprocess.run(["curl", "-L", "-o", str(sra_toolkit_archive), sra_toolkit_url], check=True)

        # Extract the toolkit
        if sys.platform != "win32":
            subprocess.run(["tar", "-xzf", str(sra_toolkit_archive), "-C", str(sra_toolkit_dir), "--strip-components=1"], check=True)
        else:
            subprocess.run(["powershell", "Expand-Archive", str(sra_toolkit_archive), str(sra_toolkit_dir)], check=True)

        print("SRA Toolkit installed successfully.")

    if esearch_path.exists():
        print("EDirect tools are already installed.")
    else:
        print("Installing EDirect tools in user mode...")
        subprocess.run(["sh", "-c", f"mkdir -p {bin_dir} && cd {bin_dir} && curl -O https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.zip && unzip -o edirect.zip && rm edirect.zip"], check=True)
        
        print("EDirect tools installed successfully.")

def main():
    top_dir = sys.argv[1]
    install_sra_toolkit_user(top_dir)

if __name__ == "__main__":
    main()
