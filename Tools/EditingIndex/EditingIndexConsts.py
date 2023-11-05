__author__ = 'Hillel'
"""
This module contains the consts for the editing index.
"""
# region Builtin Imports
import os
# endregion

# region Internal Imports
# endregion

# region Mismatches & Coverage Keys and Headers
#: Total coverage at positions with reference base %s (excluding SNPs, low quality and Ns)
COVERAGE_AT_N_POSITIONS_FORMAT = "TotalCoverageAt%sPositions"

#: The number of positions with reference base %s, covered (excluding SNPs, low quality and Ns)
NUM_OF_N_SITES_COVERED_FORMAT = "NumOf%sPositionsCovered"

NUM_OF_MM_FORMAT = "NumOf%sMismatches"
SNPS_NUM_OF_MM_FORMAT = "NumOf%sMismatchesAtSNPs"

NUM_OF_MM_SITES_FORMAT = "NumOf%sMismatchesSites"
SNPS_NUM_OF_MM_SITES_FORMAT = "NumOf%sMismatchesSitesAtSNPs"

NUM_OF_CANONICAL_FORMAT = "NumOf%s"

EDITING_INDEX_FORMAT = "%sEditingIndex"

INDEXED_EDITED_FORMAT = "IndexedMismatchesOf%s"
STRANDED_INDEXED_EDITED_FORMAT = "IndexedMismatchesOf%sOn%sStrand"
REGIONED_INDEXED_EDITED_FORMAT = "IndexedMismatchesOf%sFrom%sRegions"
SUMMERY_G_REGIONED_INDEXED_EDITED_FORMAT = "IndexedMismatchesFrom%sRegions"
SUMMERY_G_STRANDED_INDEXED_EDITED_FORMAT = "IndexedMismatchesOn%sStrand"

INDEXED_CANONICAL_FORMAT = "IndexedCanonicalOf%s"
STRANDED_INDEXED_CANONICAL_FORMAT = "IndexedCanonicalOf%sOn%sStrand"
REGIONED_INDEXED_CANONICAL_FORMAT = "IndexedCanonicalOf%sFrom%sRegions"
SUMMERY_G_STRANDED_INDEXED_CANONICAL_FORMAT = "IndexedCanonicalOf%sStrand"
SUMMERY_G_REGIONED_INDEXED_CANONICAL_FORMAT = "IndexedCanonicalOf%sRegions"

NUM_OF_INDEXED_MM_SITES_FORMAT = "NumOfIndexedMismatchesSitesOf%s"
STRANDED_NUM_OF_INDEXED_MM_SITES_FORMAT = "NumOfIndexedMismatchesSitesOf%sOn%sStrand"
REGIONED_NUM_OF_INDEXED_MM_SITES_FORMAT = "NumOfIndexedMismatchesSitesOf%sFrom%sRegions"
SUMMERY_G_REGIONED_NUM_OF_INDEXED_MM_SITES_FORMAT = "NumOfIndexedMismatchesSitesFrom%sRegions"
SUMMERY_G_STRANDED_NUM_OF_INDEXED_MM_SITES_FORMAT = "NumOfIndexedMismatchesSitesOn%sStrand"

NUM_OF_INDEXED_OVERALL_SITES_FORMAT = "NumOfIndexedOverallSitesOf%s"
STRANDED_NUM_OF_INDEXED_OVERALL_SITES_FORMAT = "NumOfIndexedOverallSitesOf%sOn%sStrand"
REGIONED_NUM_OF_INDEXED_OVERALL_SITES_FORMAT = "NumOfIndexedOverallSitesOf%sFrom%sRegions"
SUMMERY_G_STRANDED_NUM_OF_INDEXED_OVERALL_SITES_FORMAT = "NumOfIndexedOverallSitesOn%sStrand"
SUMMERY_G_REGIONED_NUM_OF_INDEXED_OVERALL_SITES_FORMAT = "NumOfIndexedOverallSitesFrom%sRegions"

NUM_OF_REGIONS = "NumOfRegionsWithCoverageFor%s"
STRANDED_NUM_OF_REGIONS = "NumOfRegionsWithCoverageFor%sOn%sStrand"
REGIONED_NUM_OF_REGIONS = "NumOfRegionsWithCoverageFor%sFrom%sRegions"
SUMMERY_G_STRANDED_NUM_OF_REGIONS = "NumOfRegionsWithCoverageOn%sStrand"
SUMMERY_G_REGIONED_NUM_OF_REGIONS = "NumOfRegionsWithCoverageFrom%sRegions"

STRAND_DECIDING_METHOD = "StrandDecidingMethod"

SUMMERY_RANK_FORMAT = "RankFor%sEditingIndex"

PER_REGIONS_OUTPUT_COUNTS_FORMAT = "%(count_type)s%(operation)s"

STRAND_OF_ORIGIN = "StrandOfOrigin"
# endregion

# region Paths Formats
HG38 = "hg38"
HG19 = "hg19"
MM10 = "mm10"
MM9 = "mm9"

RESOURCES_INI_GENOME_OPT = "Genome"
RESOURCES_INI_REGIONS_OPT = "RERegions"
RESOURCES_INI_CDS_REGIONS_OPT = "CDSRegions"
RESOURCES_INI_CDS_REGIONS_REFSEQ_ALL_OPT = "CDSRegionsRefSeqAll"
RESOURCES_INI_SNPS_OPT = "SNPs"
RESOURCES_INI_REFSEQ_OPT = "RefSeq"
RESOURCES_INI_REFSEQ_ALL_OPT = "RefSeqAll"
RESOURCES_INI_GENES_EXPRESSION_OPT = "GenesExpression"

RESOURCES_INI_SAMTOOLS_OPT = "SAMToolsPath"
RESOURCES_INI_BEDTOOLS_OPT = "BEDToolsPath"
RESOURCES_INI_BAM_UTILS_OPT = "BAMUtilsPath"
RESOURCES_INI_JAVA_OPT = "JavaHome"
RESOURCES_INI_JAVA_UTILS_OPT = "EIJavaUtils"

AVAILABLE_BUILTIN_GENOMES = [HG38, HG19, MM10, MM9]
NOT_BUILTIN_GENOME = "UserProvided"
AVAILABLE_GENOMES = AVAILABLE_BUILTIN_GENOMES[:] + [NOT_BUILTIN_GENOME]

RESOURCES_INI_FORMAT = os.path.join("%(script_dir)s", "Configs", "ResourcesPaths.ini")
DEFAULTS_CONFIG_PATH_FORMAT = os.path.join("%(script_dir)s", "Configs", "DefaultsConfig.ini")
FULL_CONFIG_PATH_FORMAT = os.path.join("%(script_dir)s", "Configs", "FullConfig.ini")
STRANDED_FULL_CONFIG_PATH_FORMAT = os.path.join("%(script_dir)s", "Configs", "FullConfigStranded.ini")
# endregion


# region EIConfig Consts
REGIONS_NAME_OPTION = "regions_name"
GENOME_FASTA_OPTION = "genome_fasta"
REGIONS_BED_OPTION = "regions_coordinates_bed"
BAM_FILE_SUFFIX_OPTION = "bam_file_suffix"
IN_BAM_FILE_OPTION = "input_bam_file"
ALIGNER_OUT_FORMAT = "aligner_output_format"
SAMTOOLS_PATH_OPTION = "samtools_path"
BEDTOOLS_PATH_OPTION = "bedtools_path"
BAM_UTILS_PATH_OPTION = "bam_utils_path"
JAVA_HOME_OPTION = "java_home_path"
JAVA_UTILS_HOME_OPTION = "ei_java_utils_jar"
REGIONS_PILEUP_COUNT_OPTION = "regions_pileup_with_count"
REGIONS_PILEUP_COUNT_STRAND_1_OPTION = "regions_pileup_with_count_strand_1"
REGIONS_PILEUP_COUNT_STRAND_2_OPTION = "regions_pileup_with_count_strand_2"
GENOME_INDEX_PATH_OPTION = "genome_index_path"
IS_PE_OPTION = "is_paired_end"
OVERALL_MAX_PROCESSES_NUM_OPTION = "overall_max_processes_num"
Q_THRESHOLD_OPTION = "PileupToCount_quality_threshold"
# endregion

# region EditingIndexJavaUtils Consts
ORIG_REG_SEP = "#"
# endregion
