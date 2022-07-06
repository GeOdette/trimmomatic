"""
Trimming tasks for illumina paired-end and single ended data.
"""


from latch.types.glob import file_glob
import glob
import subprocess
from pathlib import Path
from enum import Enum
from latch import small_task, workflow
from latch.types import LatchFile, LatchDir
import typing
from typing import Optional, Tuple, List
import os
from dataclasses import dataclass
from dataclasses_json import dataclass_json, LetterCase
from flytekit import Resources, map_task


@dataclass_json(letter_case=LetterCase.CAMEL)
@dataclass
class SingleEndReads:
    r1: LatchFile


@dataclass_json(letter_case=LetterCase.CAMEL)
@dataclass
class PairedEndReads:
    r1: LatchFile
    r2: LatchFile


@dataclass_json(letter_case=LetterCase.CAMEL)
@dataclass
class ReverseEndReadsR1:
    trimmed_r1: str
    untrimmed_r1: str


@dataclass_json(letter_case=LetterCase.CAMEL)
@dataclass
class ReverseEndReadsR2:
    trimmed_r2: str
    untrimmed_r2: str


@dataclass_json(letter_case=LetterCase.CAMEL)
@dataclass
class EndReadsMapTuple:
    paired_end_reads: PairedEndReads
    reverse_end_reads_r1: ReverseEndReadsR1
    reverse_end_reads_r2: ReverseEndReadsR2

        
class ReadType(Enum):
    Single_end_reads = "SE"
    paired_end_reads = "PE"


class phred(Enum):
    phred33 = "-phred33"
    phred64 = "-phred64"


@small_task
# The aim here is to create file pairs from a list of files in a directory ending with either _r1 or _r2, for any case, _r1 is supposed to pair to the next _r2
def LatchFilePairs(input_dir: LatchDir) -> typing.List[EndReadsMapTuple]:
    pairs = []  # This should create a list to contain the forwad read file pairs
    for file in input_dir.local_path:
        r1 = glob.glob("*_r1")
        r2 = glob.glob("*_r2")
        trimmed_r1 = [read1.replace("_r1", "trimmed_r1") for read1 in r1]
        untrimmed_r1 = [read1.replace("_r1", "utrimmed_r1") for read1 in r1]
        trimmed_r2 = [read2.replace("_r2", "trimmed_r2") for read2 in r2]
        untrimmed_r2 = [read2.replace("_r2", "utrimmed_r1") for read2 in r2]

        paired_forwad_files = PairedEndReads(r1.local_path, r2.local_path)
        paired_reverse_read1 = ReverseEndReadsR1(trimmed_r1, untrimmed_r1)
        paired_reverse_read2 = ReverseEndReadsR2(trimmed_r2, untrimmed_r2)
        pairs.append(EndReadsMapTuple(paired_forwad_files, paired_reverse_read1, paired_reverse_read2))
        
    return pairs


@small_task
def trimmomatic_task(end_read_tuple: EndReadsMapTuple, out_dir: LatchDir, adapter_seq: LatchFile, Slid_wnd: str = "4:30:10", min_len: str = "30", threads: Optional[str] = "4", read_type: ReadType = ReadType.paired_end_reads, phred: phred = phred.phred64) -> Tuple[LatchFilePairs, LatchDir]:

    # local_prefix = os.path.join(local_dir, "latch")
    forwad = end_read_tuple.paired_end_reads
    reverse1 = end_read_tuple.reverse_end_reads_r1
    reverse2 = end_read_tuple.reverse_end_reads_r2
    # Define the trimmomatic function
    if read_type.value == "PE":
        _trimmomatic_cmd = [
            "java",
            "-jar",
            "trimmomatic-0.39.jar",
            read_type.value,
            "-threads",
            str(threads),
            phred.value,
            "-basein",
            *forwad, # /dir/file_r1 /dir/file_r2
#             forwad.r1,
#             forwad.r2,
            "-baseout",
            *reverse1,
            * reverse2,
            "ILLUMINACLIP:",
            adapter_seq.local_path,
            "SLIDINGWINDOW:",
            str(Slid_wnd),
            "MINLEN:",
            str(min_len), ]
        subprocess.run(_trimmomatic_cmd, check=True)
        return LatchDir(str(local_dir), out_dir.remote_path)

    
@dynamic
def my_own_loop(end_read_tuples: typing.List[EndReadsMapTuple], out_dir: LatchDir, adapter_seq: LatchFile, Slid_wnd: str = "4:30:10", min_len: str = "30", threads: Optional[str] = "4", read_type: ReadType = ReadType.paired_end_reads, phred: phred = phred.phred64) -> typing.List[Tuple[LatchFilePairs, LatchDir]]:
    res = []
    for end_read_tuple in end_read_tuples:
        o = trimmomatic_task(end_read_tuple=end_read_tuple, out_dir=out_dir, adapter_seq=adapter_seq,
                                Slid_wnd=Slid_wnd, min_len=min_len, threads=threads, read_type=read_type, phred=phred)
        res.append(o)
    return res


@ workflow
def trimmomatic(input_dir: LatchDir, out_dir: LatchDir, adapter_seq: LatchFile, Slid_wnd: str = "4:30:10", min_len: str = "30", threads: Optional[str] = "4", read_type: ReadType = ReadType.paired_end_reads, phred: phred = phred.phred64) -> LatchDir:
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
    pairs = LatchFilePairs(input_dir=input_dir)
    # map_task(partial(trimmomatic_task, out_dir=out_dir, adapter_seq=adapter_seq,
    #                            Slid_wnd=Slid_wnd, min_len=min_len, threads=threads, read_type=read_type, phred=phred))(end_read_tuple=pairs)
    
    trimming = my_own_loop(end_read_tuples=pairs, out_dir=out_dir, adapter_seq=adapter_seq,
                                Slid_wnd=Slid_wnd, min_len=min_len, threads=threads, read_type=read_type, phred=phred)
