__author__ = 'Hillel'
# =====================imports=====================#

# =====================constants===================#

# region Step Types
STEP_MODULE_TYPE_PIPELINE = "pipeline"
STEP_MODULE_TYPE_WHITEBOX = "whitebox"
STEP_MODULE_TYPE_BLACKBOX = "blackbox"
#: These are the possible values of step type, change according to implemented steps.
STEP_MODULE_TYPES_ENUM = (STEP_MODULE_TYPE_WHITEBOX, STEP_MODULE_TYPE_BLACKBOX, STEP_MODULE_TYPE_PIPELINE)
# endregion

# region Logging
# region Error Msgs
GOT_INIT_PARAMS_ERRS = "Step %s Of Type %s Should Not Have Init Parameters!"
GOT_TARGET_FUNC_ERRS = "Step %s Of Type %s Should Not Have Target Function!"
INVALID_FUNC_ERR = "Step %s's Function %s Of Module %s Doesn't Exist!!"
INVALID_MODULE_ERR = "Step %s's Module Doesn't Exist!! (error msg: %s)"
INVALID_STEP_TYPE_ERR = "Step Type Of %s Is Invalid!"
MISSING_RUN_PARAMS_ERRS = "Step %s Of Type %s Must Have Run Parameters!"
NO_MODULE_ERR = "No Module Was Specified For Step %s"
NO_TARGET_WHITE_BOX_ERR = "No Target Function Was Specified For Step %s"
SESSION_MULTI_ARGS_ERR = "Step %s (Type %s) Have Malformed Multiple Args (Args: %s)! Message: %s"
SESSION_CONF_FORMAT_ERR = "Could Not Load Conf (Failed On Step: %s) It Has Malformed Formatting! %s"
SESSION_CONF_STRUCT_ERR = "Could Not Load Conf (%s) It Is Malformed! (Many times this is due to looped steps)"
SESSION_STEP_FAILURE_ERR = "Step %s (Type %s) Have Failed To Run Properly!"
STEP_RETURNED_VALUE_ERR = "Step %s of type %s, Has Invalid Return Value(s)! All Must Be Instances Of ReturnedValue"
# endregion

# region Warning Msgs
ERROR_STEP_WRN = "Session Going To Error Step: %(step)s"
SESSION_STEP_FAILURE_REVERT_WRN = "Step %s (Type %s) Reverting!"
# endregion

# region Info & Debug Msgs
END_STEP_MSG = "Session Finished Step: %s"
NEXT_STEP_MSG = "Session Going From Step %s To Step %s"
RET_VAL_LOG_MSG = "Step %(name)s of Type %(type)s returned:\r\n %(ret_val)s"
RUNNING_STEP_MSG = "Session Running Step: %(step)s"
SESSION_END_MSG = "Session Finished!"
SESSION_START_MSG = "Started Session With:" + "\r\n\t".join(["Input Dir: %(input_dir)s",
                                                             "Outputer Dir: %(output_dir)s",
                                                             "Logging Dir: %(logging_dir)s",
                                                             "Session Config: %(config)s",
                                                             "Session Config Args: %(c_args)s",])
SKIPPED_STEP_MSG = "Session Skipped On Step %s"
STARTED_RUNNING_MSG = "Started Running"
START_STEP_MSG = "Session Started Step: %(step)s - %(name)s"
# endregion

SESSION_LOG_FILE_FORMAT = "%s_SessionLog.log"
SESSION_STEP_LOG_FILE_FORMAT = "%s_SessionLog_Step%s.log"

# endregion

# region Session Option Names and Special Chars

#: A flag option. If set the session step will receive passed on values from previous steps.
SESSION_ACCEPT_PASSED_VALUES_OPTION = "AcceptPassedValues"

#: The option that holds the type of the step
SESSION_STEP_TYPE_OPTION = "step_type"

#: This is the separator between a name of the parameter and the value to send to it.
SESSION_PARAMS_SEP = "<<"

#: This is value for paramters of a cmd line that do not have a name (e.g. <script> some_param)
SESSION_BLACKBOX_EMPTY_PARAM_OPTION = "no_option"

#: The name of the param to use for sending args for the a piplne configuration itself
PIPELINE_PARAM_CONFIG_ARGS = "config_args"

#: options for StepRunners
LOG_RETURNED_VALUE_OPTION = "log_returned_value"
PASS_RETURNED_VALUE_OPTION = "pass_returned_value"


#: parameter to init a class with
SESSION_INIT_PARAM_OPTIONS = "init_params" #: TODO - combine later (not functional now)
#: parameter to run a function with
SESSION_RUN_PARAM_OPTIONS = "run_params"


#: The full reference of the module to import (module name for white boxes, run path for black boxes)
SESSION_STEP_MODULE_OPTION = "module"
#: The name of the  function to use (for white boxes)
SESSION_STEP_TARGET_FUNC_OPTION = "target_function"

#: If an arg ends with this suffix, it will not be sent if it's value is '' (empty str)
SESSION_OVERRIDE_ONLY_PARAM_FORMAT = '_OVERRIDE_ONLY_FORMAT'
# endregion


STEP_SUCCESS_RT_NAME = "Success"

STEP_RT_VAL_SKIPPED = "SKIPPED"

STEP_PROCESS_NAME_FORMAT = "%(timestamp)s_%(type)s_%(name)s_%(desc)s"

STEP_PRINT_ONLY_FORMAT = "%(step)s (%(name)s) of Type %(type)s - '%(run_str)s'"