__author__ = 'Hillel'
"""
This module contains the constants needed for any configuration of the session
"""
# =====================imports=====================#
import os

# =====================constants===================#

# region Condition Constants
COND_STR_CMP_SEP = "="
NEXT_STEPS_SEPARATOR = ";"
STEP_CONDITION_SEPARATOR = ","
IN_CONDITION_SEPARATOR = " "


def nested_condition(cond):
    import ConfigFunctions
    return ConfigFunctions.condition(cond)

CONDITIONS = {
    "Exists": os.path.exists,
    "Bool": lambda x: (x == "True" or x == "False") and eval(x),  # for sanity and security reasons.
    "StrCmp": lambda x: x.split(COND_STR_CMP_SEP)[0] == x.split(COND_STR_CMP_SEP)[1],
    "Condition": nested_condition
}
# endregion


# region Sections and Options
MAX_STEP_PROCESSES_OPTION = 'max_processes_num'
OVERALL_MAX_PS_OPTION = 'overall_max_processes_num'
NEXT_STEP_OPTION = "next_step"
ERROR_STEP_OPTION = "error_step"
STEP_TYPE_OPTION = "Type"
STEP_NAME_OPTION = "name"
PROGRAM_NAME_OPTION = "program"
PROGRAM_PARAMS_OPTION = "parameters"
CONSTRAIN_COND_OPTION = "constraint"
WALL_TIME_OPTION = "wall_time"

ENABLE_SECTION_OPTION = "enable"
RECOVER_SECTION = "Recover"

# The default name for the first step and error step
FIRST_STEP = "Step_1"
ERROR_STEP = "Step_-1"
STEP_SEMI_CONDITION = None
# endregion
