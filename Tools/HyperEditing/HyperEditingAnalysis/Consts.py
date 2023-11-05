from Commons.general_functions import format_to_re
from Tools.HyperEditing.HyperEditingConsts import HE_DETECT_DIR_PER_SAMPLE_BEDS_FORMAT

SAMPLE_NAME_REGEX = format_to_re(HE_DETECT_DIR_PER_SAMPLE_BEDS_FORMAT)
POSITION = "Position"
MOTIF_OUTPUT_FORMAT = "%(sample_name)s.DeaminationMotif.tab"
