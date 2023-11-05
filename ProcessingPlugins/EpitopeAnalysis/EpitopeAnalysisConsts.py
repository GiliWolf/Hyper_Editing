__author__ = 'Hillel'
"""This is the module containing the consts needed for the package"""

from Commons.data_structs import  AffinitiesEnum


OFFSET_HEADER = "Pos"
P_ID_HEADER = "ID"
WINDOW_CONTENT_HEADER = "Peptide"
NM_AFIINITY_HEADER = "nM"
RANK_HEADER = "Rank"
CORE_HEADER = "Core"
RANKS_AVG_HEADER = "H_Avg_Ranks"
BINDERS_NUM_HEADER = "N_binders"
NET_MHC_PAN_RANKS_AVG_HEADER = "Ave"
NET_MHC_PAN_BINDERS_NUM_HEADER = "NB"
CHR_I = 0
POSITION_I = 1
CHANGE_P_I = 4
STRONG_BINDING_PEPTIDES_THRESHOLD = 0.5
WEAK_BINDING_PEPTIDES_THRESHOLD = 2.0
STROG_BIND_VAL = AffinitiesEnum.STRONG
INTER_BIND_VAL = AffinitiesEnum.INTERMEDIATE
WEAK_BIND_VAL = AffinitiesEnum.WEAK
ORIG_PEP_TAG = "Orig"
NEO_PEP_TAG = "Neo"
NM_AFFINITY_FORMAT = "Affinity_%s_%s"
RANK_FORMAT = "Rank_%s_%s"
CORE_FORMAT = "Core_%s_%s"
OFFSET_FORMAT = "Offset"
BIND_STRENGTH_FORMAT = "BindStrength_%s_%s"
P_LENGTH = "PeptideLength"
FILTER_NET_MHC_WINDOWA_ERROR = "Peptides Metadata CSV File Is Not Compliant With NetMHC Results!"
LOG_FILE = "neo_epitopes_analysis.log"

OUTPUT_FILENAME = "%(prediction_tool)s_Analysis_Epitope%(summary_type)s.csv"
