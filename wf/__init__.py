"""
Trimming tasks for illumina paired-end and single ended data.
"""

from asyncore import read
import subprocess
from pathlib import Path
from enum import Enum
from threading import local
from latch import small_task, workflow
from latch.types import LatchFile, LatchDir
from typing import Optional, Tuple
import os


class phred(Enum):
    phred33 = "-phred33"
    phred64 = "-phred64"


@small_task
def trimmomatic_task(read1: LatchFile, read2: LatchFile, out_dir: LatchDir, adapter_seq: LatchFile, Slid_wnd: str = "4:30:10", min_len: str = "30", threads: Optional[str] = "4",  PE_status: bool = False, phred: phred = phred.phred64) -> LatchDir:

    # defining outputs

    local_dir = "/root/latch_trim/"
    # local_prefix = os.path.join(local_dir, "latch")

    # Defining output files
    # orphaned output for reads
    orphaned_r1 = os.path.join(local_dir, "untrimmed_", read1)
    orphaned_r2 = os.path.join(local_dir, "untrimmed_", read2)

    # surviving pairs
    r1_surviving = os.path.join(local_dir, "trimmed_", read1)
    r2_surviving = os.path.join(local_dir, "trimmed_", read2)

    # Define the trimmomatic function
    if PE_status == True:
        _trimmomatic_cmd = [
            "java",
            "-jar",
            "trimmomatic-0.39.jar",
            "PE",
            "-threads",
            str(threads),
            phred.value,
            "-basein",
            read1.local_path,
            read2.local_path,
            "-baseout",
            str(r1_surviving),
            str(orphaned_r1),
            str(r2_surviving),
            str(orphaned_r2),
            "ILLUMINACLIP:",
            adapter_seq.local_path,
            "SLIDINGWINDOW:",
            str(Slid_wnd),
            "MINLEN:",
            str(min_len), ]
        subprocess.run(_trimmomatic_cmd, check=True)
        return LatchDir(str(local_dir), out_dir.remote_path)


@ workflow
def trimmomatic(read1: LatchFile, read2: LatchFile, out_dir: LatchDir, adapter_seq: LatchFile, Slid_wnd: str = "4:30:10", min_len: str = "30", threads: Optional[str] = "4",  PE_status: bool = False, phred: phred = phred.phred64) -> LatchDir:
    """Trimming tasks for illumina paired-end and single ended data.

    Trimmomatic

    __metadata__:
        display_name: Trimmomatic_latch

        author:
            name: GeOdette

            email: steveodettegeorge@gmail.com

            github:
        repository:

        license:
            id: MIT

    Args:


        PE_status:
          Tell the program whether you want to process paired end reads

          __metadata__:
            display_name: Paired End reads

        phred:
          phred. i.e default is phred64

          __metadata__:
            display_name: phred 

        read1:
          Paired-end read 1 file to be trimmed

          __metadata__:
            display_name: Read 1

        read2:
          Paired-end read 2 file to be trimmed

          __metadata__:
            display_name: Read 2

        threads:
          Tell the program how many threads you want it to run with

          __metadata__:
            display_name: Threads

        adapter_seq:
          Tell the program which adapter and other illumina-specific sequences to cut from the read.

          __metadata__:
            display_name: ILLUMINACLIP e.g NexteraPE-PE.fa

        Slid_wnd:
          Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.

          __metadata__:
            display_name: SLIDINGWINDOW e.g 4:30:10

        min_len:
          Tell the program to discard any reads that are less than set minimum length after trimming

          __metadata__:
            display_name: MINLEN e.g 30

        out_dir:
          Output directory

          __metadata__:
            display_name: Output directory
    """
    return trimmomatic_task(read1=read1, read2=read2, out_dir=out_dir, adapter_seq=adapter_seq, Slid_wnd=Slid_wnd, min_len=min_len, threads=threads, PE_status=PE_status, phred=phred)
