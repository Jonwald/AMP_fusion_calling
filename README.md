# AMP fusion calling Pipeline

## Description

nextflow pipeline to perform gene fusion calling on anhcored multiplex PCR based RNA-seq libraries. Currently configured for use with the fusionplex solid tumor panel, but can be used as a generalised fusion calling pipeline for targetted RNA-seq libraries with minor modifications.

### stages:
  - Trim UMI sequence and space append UMIs to read headers (custom script)
  - Adapter / quality trim reads (bbduk.sh)
  - Align reads to reference genome (STAR)
  - Deduplicate Aligned reads using coordinates and UMI sequences (UMItools)
  - Re-align deduplicated reads (STAR-fusion)

## setup / requirements

Entire pipleine currently runs within a conda environemnt, future plans to containerise all steps and add a docker profile.

to set up the conda environemnt:
	conda create --name <env-name> --file env_list

Requires reference genome which is compatable with STAR-fusion v1.10.0, recommend using pre-build CTAT genome index listed in this table (https://github.com/STAR-Fusion/STAR-Fusion/wiki/STAR-Fusion-release-and-CTAT-Genome-Lib-Compatibility-Matrix)

## Usage

activate conda environment:
conda activate <env-name>

ensure nextflow is installed (https://www.nextflow.io/)
edit nextflow.config file and AMP_fusion_calling.nf to reflect available resources

run with nextflow run AMP_fusion_calling.nf

## Notes

## to do
- add test data / expected results
- remove hard paths from code
- update to Nextflow DSL2 syntax
- add docker profile
- modify for automatic resource detection and scaling
