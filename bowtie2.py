"""
Wrapper script for bowtie2-align and bowtie2-build with version control
VERSION by Scotty for improved bowtie2 version checking
"""
import subprocess
from itertools import izip_longest


def build(version, fasta):
    """
    Construct bowtie2 indices (.bt2 files)
    :param version: Enforces bowtie2 version number (e.g., '2.2.1')
    :param fasta: Path to FASTA containing reference sequences
    :return:
    """
    base_str = 'bowtie2-build'
    version_str = '%s-%s' % (base_str, version)
    build_cmd = version_str
    try:
        subprocess.check_output([build_cmd, '-h'])
    except OSError:
        build_cmd = base_str
        try:
            subprocess.check_output([build_cmd, '-h'])
        except:
            raise RuntimeError("""Neither '%s' nor '%s' found;
check if at least one of these is installed and in $PATH\n""" % (version_str, base_str))

    stdout = subprocess.check_output([build_cmd, '--version'])
    local_version = stdout.split('\n')[0].split()[-1]
    assert version == local_version, 'bowtie2-build version incompatibility %s != %s' % (version, local_version)
    subprocess.check_call([build_cmd, '-q', '-f', fasta, fasta])


def align_paired(version, refpath, fastq1, fastq2, nthreads, flags=('--quiet', '--no-unal', '--local')):
    """
    Call bowtie2-align on paired read data
    :param version: Enforces bowtie2 version number
    :param refpath: Path to bowtie2 index files (.bt2)
    :param fastq1: Files with #1 mates
    :param fastq2: Files with #2 mates
    :param flags: A tuple of bowtie2 flags such as --local
    :return:
    """
    # check that we have access to bowtie2
    base_str = 'bowtie2'
    version_str = '%s-%s' % (base_str, version)
    bt2_cmd = version_str
    try:
        subprocess.check_output([bt2_cmd, '-h'])
    except OSError:
        bt2_cmd = base_str
        try:
            subprocess.check_output([bt2_cmd, '-h'])
        except:
            raise RuntimeError("""Neither '%s' nor '%s' found;
check if at least one of these is installed and in $PATH\n""" % (version_str, base_str))

    # check that version is the expected version
    stdout = subprocess.check_output([bt2_cmd, '--version'])
    local_version = stdout.split('\n')[0].split()[-1]
    assert version == local_version, 'bowtie2 version incompatibility %s != %s' % (version, local_version)

    # stream output from bowtie2, two at a time
    bowtie_args = [bt2_cmd, '-x', refpath, '-1', fastq1, '-2', fastq2, '-p', str(nthreads)] + list(flags)
    p2 = subprocess.Popen(bowtie_args, stdout=subprocess.PIPE)
    for line, line2 in izip_longest(p2.stdout, p2.stdout, fillvalue=''):
        yield line, line2

    # exception handling
    return_code = p2.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, bowtie_args)


def main():
    iter = align_paired('2.2.8', 'gb-ref.fa', 'test1.fastq', 'test2.fastq', 1)
    for line in iter:
        print line

if __name__ == '__main__':
    main()
