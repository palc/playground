# playground
Tod Stuber's miscellaneous and experimental scripts

## OVERVIEW
**alignmentreport-latex.sh** Alignment and run LaTeX summary report as pdf.  Working directory must contain reference with .fasta extension and single or paired FASTQ reads.  Kraken can be called with -k option.

**assemble_verify.sh** Select assemble program, verify reads against local BLAST NT database and run Kraken.  Run with -h option for help.

**assemblereport-latex.sh** Assemble and run LaTeX summary report as pdf

**blast-contigs.sh** Give fasta file as argument.  Get a quick summary of read identifications against NT database.

**Bruc_MLST.sh** Get Brucella species MLST type.

**highqualitymaker.sh** Make a HQ VCF for GATK best practices.

**krakenreadfinder.sh** Pull reads based on Kraken identifications.  Pull reads by either identification name or report line number.  See script header for details.

**ksnp_validate.sh** For a specified node get SNP positions and annotation.
 
**ksnp3_run.sh** Run kSNP3 on working directory containing draft or finished assembled genomes.  Run with -h option for help.

**makeReference.sh** Make a reference from a reference guided assembly mimicking idvirus script.

**ncbi_download_assembled_genomes_from_xml.sh** Quickly select WGS from NCBI in XML format and download.

**oligo_identifier_test.sh** Count oligos in raw FASTQ files.  Used for primer/probe design.

**packagefastqs.sh** Package FASTQs into directories and echo one-liner to call on all directories.

**quickalignment.sh** Working directory to have paired FASTQs and reference with .fasta extension.

**readfindingbyends.sh** Gather reads from alignment/assembly and gather additional reads by pair read information.

**sixteen_s.sh** 16s and rpoB analysis and report.  Run without arguments for help.

**sra_download.sh** Download SRA using DNA Nexus urls.