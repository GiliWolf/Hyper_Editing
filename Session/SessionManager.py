import Session.ConfigConsts
import Session.ConfigFunctions

__author__ = 'Hillel'
# =====================imports=====================#
# region Builtin Python Imports
import abc
import argparse
import importlib
import logging
import multiprocessing
import os
import re
import subprocess
from ConfigParser import InterpolationError
from collections import namedtuple
from datetime import date
from pprint import pformat
# endregion

# region Internal Imports
if __name__ == '__main__':
    import sys
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import Config
import PipelineManger
from Commons.consts import ERROR_VAL
from Commons.general_functions import remove_files, init_logging_dict, convert_args_to_dict
from ConfigConsts import FIRST_STEP, ENABLE_SECTION_OPTION, STEP_NAME_OPTION
from SessionConsts import *
# endregion

# =====================constants===================#
#: These are all the option that contain name of a target (e.g. class or function name) update accordingly
STEP_MODULE_NAME_PARAMS_OPTIONS = [SESSION_STEP_TARGET_FUNC_OPTION, SESSION_STEP_MODULE_OPTION]

#: These are all the option that contain multiple params under them (using SESSION_PARAMS_SEP) update accordingly
STEP_MULTI_PARAMS_OPTIONS = [SESSION_INIT_PARAM_OPTIONS, SESSION_RUN_PARAM_OPTIONS]
#: This regex is used to parse the multi params
STEP_MULTI_PARAMS_RE = re.compile(r"([\w_\-]+?<<.*?(?=(?:[\w_\-]+?<<)|$))", re.DOTALL)

#: These are all the option that contain flags of each step update accordingly
STEP_FLAGS_OPTIONS = [ENABLE_SECTION_OPTION, LOG_RETURNED_VALUE_OPTION, PASS_RETURNED_VALUE_OPTION]

STEP_RUNNER_STR_FORMAT = "%(step)s, %(step_type)s, Next Steps: %(all_next_steps)s, Extra Info %(extra_info)s"
NEXT_STEP_FORAMT = r"\r\n\t%(step)s If %(condition)s"

BLACKBOX_RUN_STR_FORMAT = r"%(module)s %(params)s"
BLACKBOX_RUN_STR_PARAMS_FORMAT = r"%s %s"
BLACKBOX_RET_CODE_NAME = "return_code"

PIPELINE_RUN_STR_FORMAT = r"%(module)s %(params)s"
PIPELINE_RUN_STR_PARAMS_FORMAT = r"%s %s"
PIPELIE_RT_VAL_NAME_TMP_CNF = "TempConfFiles"
PIPELIE_RT_VAL_NAME_CREATED_FILES = "CreatedDirsAndFiles"


# =====================classes=====================#
#: This named tuple is used for a step returned value, the first three params are flags for possible actions,
#  (change if options ar added) the last is the value itself and a name to address it.
ReturnedValue = namedtuple("ReturnedValue", "pass_on remove_files log value val_name")


class SessionConfigException(Exception):
    pass


class SessionStepReturnedValues(object):
    """
    This class is a container for a return value of a session step.
    """

    def __init__(self, step_name, step_type, returned_values):
        """
        The constructor .
        :param step_name: The name of the step returning the data.
        :type step_name: C{str}
        :param step_type: The type of the step retruning the data.
        :type step_type: Implementer of C{SessionStepRunner}
        :param returned_values: A list of the returned values and what to do with them.
        :type returned_values: C{ReturnedValue}
        """
        assert step_type in STEP_MODULE_TYPES_ENUM, INVALID_STEP_TYPE_ERR % step_type
        assert all([isinstance(rv, ReturnedValue) for rv in returned_values]), STEP_RETURNED_VALUE_ERR % (step_name,
                                                                                                          step_type)
        super(SessionStepReturnedValues, self).__init__()
        self.step_type = step_type
        self.step_name = step_name

        self.passed_on = {}
        for val in returned_values:
            if val.log:
                logging.info(RET_VAL_LOG_MSG % dict(name=step_name, type=step_type, ret_val=pformat(val)))
            if val.remove_files:
                remove_files(val.value)
            if val.pass_on:
                self.passed_on[val.val_name] = val.value


class SessionStepRunner(object):
    """
    This class is an interface for a step in a session.
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, config, step_name):
        """
        The constructor which should call parse_and_assert_step. * DO NOT OVERRIDE!!! *
        :param config: the loaded configuration.
        :param step_name: The name of the step (for accessing)
        """
        try:
            self.__parse_args(config, step_name)
        except InterpolationError , e:
            raise SessionConfigException(SESSION_CONF_FORMAT_ERR % (step_name, e.message))
        self._parse_and_assert_step()
        self._returned_values = None

    @abc.abstractmethod
    def run_step(self, print_only=False, **kwargs):
        """
        When this function is called the step should be executed.
        :param print_only: If set will only print the run command of the step instead of running it.
        :param kwargs: optional args.
        :return: The step success or failure.
        :rtype  C{bool}
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def revert(self, **kwargs):
        """
        This function is called in case the step failed. It reverts the known data(and files) to what was before.
        :param kwargs: optional args.
        :return: Success or failure
        :rtype: C{bool}
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def get_returned_values(self, **kwargs):
        """
        This function retrieves the step's returned values.
        :param kwargs: optional args.
        :return: The returned values.
        :rtype: C{SessionStepReturnedValues}
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def _parse_and_assert_step(self):
        """
        This function loads the needed params from the step and performs the necessary sanity checks.
        :return: None
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def configure_step(self):
        """
        An interactive function for configuring steps on run time
        :return:  None
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def step_to_str(self, returned_values):
        """
        This function is for logging the steps. Should be called after run_step.
        :param returned_values: The returned values from the step.
        :return: None
        """
        raise NotImplementedError()

    def next_step(self, success=True):
        """
        This is a getter fir the next step. It's calculated upon calling for runtime evaluation of the conditions.
        :param success: The success status of running the step.
        :type success: C{bool}
        :return: The name of the next step.
        :rtype: C{str}
        """
        if not success:
            logging.warn(ERROR_STEP_WRN % {'step': self.error_step})
            return self.error_step

        for n_step, condition in self.all_next_steps:
            if None is condition:  # This is a default fall-to step or the error step.
                return n_step
            if Session.ConfigFunctions.condition(cond=condition):
                return n_step

    def __parse_args(self, config, step):
        """
        This function converts the options for the steps and set values in the instance according to the options name.
         * DO NOT OVERRIDE!!! *
        :param config: a loaded C{Config}
        :param step: The name of the step to parse.
        :return: the parsed arguments.
        :rtype None
        """
        # Set flags
        for flag in STEP_FLAGS_OPTIONS:
            flag_val = config.get(step, flag) if config.has_option(step, flag) else None
            self.__setattr__(flag.lower(), flag_val)

        # Set names
        self.step = step

        self.error_step = config.get_error_step(step)
        self.step_type = config.get(step, SESSION_STEP_TYPE_OPTION)
        self.name = config.get(step, STEP_NAME_OPTION)

        self.process_name = STEP_PROCESS_NAME_FORMAT % dict(timestamp=date.today().isoformat(),
                                                            type=self.step_type, name=self.step,
                                                            desc=self.name.replace(" ", "_"))

        self.all_next_steps = config.get_all_next_steps(step)

        self.enabled = config.enabled(step)

        self.accept_passed_vals = config.getboolean(step, SESSION_ACCEPT_PASSED_VALUES_OPTION)

        # TODO: add target specific option using prefixes and searching - step_options = config.options(step)
        args = {}
        # parse and set options containing multiple values
        for multi_param_opt in STEP_MULTI_PARAMS_OPTIONS:
            if config.has_option(step, multi_param_opt):
                opt = dict()
                for param in STEP_MULTI_PARAMS_RE.findall(config.get(step, multi_param_opt)):
                    field, val = param.split(SESSION_PARAMS_SEP)
                    field = field.strip("\r\n")
                    val = val.strip("\r\n")
                    try:
                        opt[field] = eval(val)
                    except Exception ,e:
                        logging.exception(SESSION_MULTI_ARGS_ERR % (self.step, self.step_type, field, e.message))
                        raise AssertionError(SESSION_MULTI_ARGS_ERR % (self.step, self.step_type, field,
                                                                       str(e) + e.message))
                    if isinstance(type(opt[field]), dict):
                        for key in opt[field].keys():
                            if key.endswith(SESSION_OVERRIDE_ONLY_PARAM_FORMAT):
                                tmp = opt[field].pop(key)
                                if "" != tmp:  # only of value is overridden send it.
                                    opt[field][key.replace(SESSION_OVERRIDE_ONLY_PARAM_FORMAT, '')] = tmp


            else:
                opt = None
            self.__setattr__(multi_param_opt.lower(), opt)

        # parse and set option containing target names to use (e.g. function name)
        for target_name_opt in STEP_MODULE_NAME_PARAMS_OPTIONS:
            target_val = config.get(step, target_name_opt) if config.has_option(step, target_name_opt) else None
            self.__setattr__(target_name_opt.lower(), target_val)


class BlackBoxRunner(SessionStepRunner):
    """
    This class runs a black box step
    """

    def revert(self, **kwargs):
        return False  # Cannot possibly revert the run of an unknown external program...

    def configure_step(self):
        super(BlackBoxRunner, self).configure_step()

    def get_returned_values(self, **kwargs):
        return self._returned_values.passed_on

    def run_step(self, print_only=False, **kwargs):
        run_params = getattr(self,SESSION_RUN_PARAM_OPTIONS.lower(), None)
        module = getattr(self,SESSION_STEP_MODULE_OPTION.lower(), None)

        run_str_params = " ".join([BLACKBOX_RUN_STR_PARAMS_FORMAT % (key, val) for key, val in run_params.iteritems()])
        run_str_params = run_str_params.replace(SESSION_BLACKBOX_EMPTY_PARAM_OPTION, "")

        if print_only:
            run_str = BLACKBOX_RUN_STR_FORMAT % {'module': module, 'params': run_str_params}
            print STEP_PRINT_ONLY_FORMAT % dict(step=self.step, name=self.name, type=self.step_type,
                                                run_str=run_str)
            return True

        try:
            ret_code = subprocess.call(BLACKBOX_RUN_STR_FORMAT % {'module': module, 'params': run_str_params},
                                              shell=True)

            log = getattr(self, LOG_RETURNED_VALUE_OPTION.lower(), False)
            pass_on = getattr(self, PASS_RETURNED_VALUE_OPTION.lower(), False)
            delete_files = False  # TODO: it is possible to add option where the user can write temp files to delete.
            val_name = BLACKBOX_RET_CODE_NAME
            ret_val = ReturnedValue(pass_on, delete_files, log, ret_code, val_name)
            self._returned_values = SessionStepReturnedValues(self.step, self.step_type, [ret_val, ])

        except:
            logging.exception(SESSION_STEP_FAILURE_ERR % (self.step, self.step_type))
            return False

        return True

    def step_to_str(self, returned_values):
        n_steps = ""
        for next_step, con in self.all_next_steps:
            n_steps += NEXT_STEP_FORAMT % dict(step=next_step, condition=con)
        return STEP_RUNNER_STR_FORMAT % dict(step_name=self.step, step_type=STEP_MODULE_TYPE_BLACKBOX,
                                             next_steps=n_steps, extra_info="Running: %s" %
                                                                            getattr(self,
                                                                                    SESSION_STEP_MODULE_OPTION.lower(), None))

    def _parse_and_assert_step(self):
        assert None is getattr(self,SESSION_INIT_PARAM_OPTIONS.lower(), None), GOT_INIT_PARAMS_ERRS % \
                                                                         (STEP_MODULE_TYPE_BLACKBOX, self.step)
        assert None is not getattr(self,SESSION_STEP_MODULE_OPTION.lower(), None), NO_MODULE_ERR % self.step

        run_params = getattr(self,SESSION_RUN_PARAM_OPTIONS.lower(), None)
        if run_params is None:
            self.__setattr__(SESSION_RUN_PARAM_OPTIONS.lower(), "")


class WhiteBoxRunner(SessionStepRunner):
    """
    This class runs a white box step
    """

    def revert(self, **kwargs):
        return False  # TODO: Each white box should have a revert func

    def get_returned_values(self, **kwargs):
        # TODO: These should be an attribute in WhiteBox class
        return self._returned_values.passed_on

    def step_to_str(self, returned_values):
        # TODO: add additional info according to data of the whitebox instance
        n_steps = ""
        for next_step, con in self.all_next_steps:
            n_steps += NEXT_STEP_FORAMT % dict(step=next_step, condition=con)
        return STEP_RUNNER_STR_FORMAT % dict(step_name=self.step, step_type=STEP_MODULE_TYPE_WHITEBOX,
                                             next_steps=n_steps, extra_info="Module: %s" %
                                                                            getattr(self,
                                                                                    SESSION_STEP_MODULE_OPTION.lower(),
                                                                                    None))

    def run_step(self, print_only=False, **kwargs):
        success = True
        #comm_pipe = multiprocessing.Pipe()#  todo add...
        step_process_kwargs = self.kwargs
        if kwargs:
            step_process_kwargs.update(kwargs)

        if print_only:
            run_str = ".".join([self.module.__name__, self.target_func.func_name]) + r"(%s)" % str(step_process_kwargs )
            print STEP_PRINT_ONLY_FORMAT % dict(step=self.step, name=self.name, type=self.step_type, run_str=run_str)
            return True

        try:
            self.target_func(**step_process_kwargs)
            # pipeline_ps = multiprocessing.Process(target=self.target_func, name=self.process_name,
            #                                       args=(step_process_kwargs))
            # pipeline_ps.daemon = False #  TODO: use this to kill https://stackoverflow.com/questions/40866576/run-atexit-when-python-process-is-killed
            # pipeline_ps.start()
        except:
            success = False
            logging.exception(SESSION_STEP_FAILURE_ERR % (self.step, self.step_type))

        #comm_pipe[0].close()

        # ret_val = None
        # while True:
        #     try:
        #         if comm_pipe.poll():
        #             ret_val = comm_pipe.recv()
        #             break
        #     except EOFError:
        #         break
        #
        # pipeline_ps.join()

        ret_success = ReturnedValue(pass_on=True, remove_files=False, log=True, value=success,
                                    val_name=STEP_SUCCESS_RT_NAME)
        self._returned_values = SessionStepReturnedValues(self.step, self.step_type, [ret_success, ])
        return success

    def configure_step(self):
        pass

    def _parse_and_assert_step(self):
        assert None is not getattr(self, SESSION_STEP_TARGET_FUNC_OPTION.lower(), None),\
            NO_TARGET_WHITE_BOX_ERR % self.step
        assert None is not getattr(self, SESSION_STEP_MODULE_OPTION.lower(), None), NO_MODULE_ERR % self.step
        try:  # TODO: This should be a WhiteBox class instance only
            self.module = importlib.import_module(getattr(self, SESSION_STEP_MODULE_OPTION.lower(), None))
        except ImportError, e:
            raise ValueError(INVALID_MODULE_ERR % (self.step, e))

        func = getattr(self, SESSION_STEP_TARGET_FUNC_OPTION.lower(), None)
        try:
            self.target_func = getattr(self.module, func)
        except AttributeError:
            raise ValueError(INVALID_FUNC_ERR % (self.step, func, self.module))
        self.kwargs = getattr(self, SESSION_RUN_PARAM_OPTIONS.lower(), {})


class PipelineRunner(SessionStepRunner):
    """
    This class runs a pipeline step.
    """

    def revert(self, **kwargs):
        remove_files(self._returned_values.passed_on[PIPELIE_RT_VAL_NAME_CREATED_FILES])

    def get_returned_values(self, **kwargs):
        return self._returned_values.passed_on

    def get_run_info(self):
        """
        This function retrieves the main run params (for logging and such)
        :return: the run param string.
        :rtype C{str}
        """
        return "Conf- %s, RootDir- %s LogDir- %s, OutputDir- %s" % (self.pipeline_args.get('config_file', ERROR_VAL),
                                                               self.pipeline_args.get('root_dir', ERROR_VAL),
                                                               self.pipeline_args.get('log_path', ERROR_VAL),
                                                               self.pipeline_args.get('output_dir', ERROR_VAL))

    def step_to_str(self, returned_values):
        n_steps = ""
        for next_step, con in self.all_next_steps:
            n_steps += NEXT_STEP_FORAMT % dict(step=next_step, condition=con)
        extra_info = self.get_run_info()
        return STEP_RUNNER_STR_FORMAT % dict(step_name=self.step, step_type=STEP_MODULE_TYPE_PIPELINE,
                                             next_steps=n_steps, extra_info=extra_info)

    def run_step(self, print_only=False, **kwargs):
        comm_pipe = multiprocessing.Pipe()  # Open a pipe for communication
        self.pipeline_args['pipe'] = comm_pipe  # Set the pipe under the name as defined in PipelineManger.main
        if print_only:
            run_str = "Pipline" + r"(%s)" % self.get_run_info()
            print STEP_PRINT_ONLY_FORMAT % dict(step=self.step, name=self.name,type=self.step_type, run_str=run_str)
            return True
        kwargs = dict(self.pipeline_args)
        pipeline_ps = multiprocessing.Process(target=PipelineManger.main, name=self.process_name,
                                              kwargs=kwargs)
        pipeline_ps.daemon = False #  TODO: use this to kill https://stackoverflow.com/questions/40866576/run-atexit-when-python-process-is-killed
        try:
            success = False
            temp_confs = []
            created_dirs_and_files = []

            pipeline_ps.start()
            comm_pipe[0].close()

            while True:
                try:
                    if comm_pipe[1].poll():
                        temp_confs, created_dirs_and_files, success = comm_pipe[1].recv()
                        break
                except EOFError:
                    break

            pipeline_ps.join()

        except Exception as e:
            logging.exception(SESSION_STEP_FAILURE_ERR % (self.step, self.step_type))

        ret_cof_val = ReturnedValue(pass_on=False, remove_files=True, log=True, value=temp_confs,
                                    val_name=PIPELIE_RT_VAL_NAME_TMP_CNF)
        ret_success = ReturnedValue(pass_on=True, remove_files=False, log=True, value=success,
                                          val_name=STEP_SUCCESS_RT_NAME)
        ret_created_files = ReturnedValue(pass_on=True, remove_files=False, log=True, value=created_dirs_and_files,
                                          val_name=PIPELIE_RT_VAL_NAME_CREATED_FILES)
        self._returned_values = SessionStepReturnedValues(self.step, self.step_type, [ret_cof_val,
                                                                                      ret_created_files,
                                                                                      ret_success])
        if not success:
            logging.exception(SESSION_STEP_FAILURE_ERR % (self.step, self.step_type))
            return False

        return True

    def configure_step(self):
        pass

    def _parse_and_assert_step(self):
        assert None is getattr(self,SESSION_INIT_PARAM_OPTIONS.lower(), None), GOT_INIT_PARAMS_ERRS % \
                                                                         (STEP_MODULE_TYPE_PIPELINE, self.step)
        assert None is getattr(self,SESSION_STEP_TARGET_FUNC_OPTION.lower(), None), GOT_TARGET_FUNC_ERRS% \
                                                                         (STEP_MODULE_TYPE_PIPELINE, self.step)
        assert None is not getattr(self,SESSION_RUN_PARAM_OPTIONS.lower(), None), MISSING_RUN_PARAMS_ERRS% \
                                                                         (STEP_MODULE_TYPE_PIPELINE, self.step)
        run_params = getattr(self, SESSION_RUN_PARAM_OPTIONS.lower(), None)

        self.pipeline_args = run_params


# =====================functions===================#


#: This is a dictionary to map between the config and the relevant C{SessionStepRunner}
MODULE_TYPES_MAPPING = {
    STEP_MODULE_TYPE_WHITEBOX: WhiteBoxRunner,
    STEP_MODULE_TYPE_BLACKBOX: BlackBoxRunner,
    STEP_MODULE_TYPE_PIPELINE: PipelineRunner,
}


def read(config):
    """
    This function reads the config and parses it.
    :param config: The config to parse.
    :type config: C{Config}
    :return: A dictionary of the parsed steps.
    :rtype: C{dict} of C{str} to C{SessionStepRunner}
    """
    steps = {}
    try:
        all_steps = config.get_all_steps()
    except RuntimeError:
        raise SessionConfigException(SESSION_CONF_STRUCT_ERR % config.FileName)

    for step in all_steps:
        step_type = config.get(step, SESSION_STEP_TYPE_OPTION)
        assert step_type in MODULE_TYPES_MAPPING, INVALID_STEP_TYPE_ERR % step
        steps[step] = MODULE_TYPES_MAPPING[step_type](config, step)
        assert isinstance(steps[step], SessionStepRunner)

    return steps


def run_steps(first_step_name, log_file, kill_flag, steps, print_only):
    """
    This function is a single linear subprocesses following the pipeline
    :param first_step_name: The name of the first step.
    :param log_file: The logging file started with (to go back to after every step in case of changes)
    :param kill_flag: The path of the kill flag (if exist will not continue pass the current step.
    :param steps: The parsed steps.
    :param print_only: If set will only print the run command of each step instead of running it.
    :type steps: C{dict} of C{str} to C{SessionStepRunner}
    :return: None
    """

    next_step = first_step_name
    prev_step_rt_val = {}
    while next_step in steps:
        step = steps[next_step]
        assert isinstance(step, SessionStepRunner)
        if os.path.exists(kill_flag):
            return -1
        logging.debug(multiprocessing.current_process().name)
        logging.debug(START_STEP_MSG % {'step': next_step,
                                    'name': step.name})

        success = run_step(step, prev_step_rt_val, log_file, print_only)
        if success == STEP_RT_VAL_SKIPPED:
            logging.info(SKIPPED_STEP_MSG % step.step)

        next_step = step.next_step(success)
        if not print_only:
            prev_step_rt_val = step.get_returned_values()
        logging.info(NEXT_STEP_MSG % (step.step, next_step))


def run_step(step, prev_step_rt_val, log_file, print_only):
    """
    This function wraps the run of a single step in the session.
    :param step: The step runner instance containing the step
    :type step: C{SessionStepRunner}
    :param prev_step_rt_val: The returned values of the prev step
    :type prev_step_rt_val: C{dict}
    :param log_file: The logging file started with (to go back to after every step in case of changes)
    :type log_file: C{str}
    :param print_only: If set will only print the run command of each step instead of running it.
    :type print_only: C{bool}
    :return: success status.
    """

    if not step.enabled:
        return STEP_RT_VAL_SKIPPED

    logging.debug(RUNNING_STEP_MSG % {'step': step.step})

    kwargs = dict(prev_step_rt_val) if step.accept_passed_vals else {}

    step_log_file = os.path.join(os.path.dirname(log_file), SESSION_STEP_LOG_FILE_FORMAT % (date.today().isoformat(),
                                 step.name))
    init_logging_dict(step_log_file)
    success = step.run_step(print_only=print_only, **kwargs)
    init_logging_dict(log_file=log_file)
    if not success:
        logging.warn(SESSION_STEP_FAILURE_REVERT_WRN % (step.step, step.step_type))
        step.revert()
    else:
        logging.info(END_STEP_MSG % step.step)

    return success


def run_session(input_dir, output_dir, config_path, logging_dir, log_file, c_args={}, print_only=False,
                kill_flag=PipelineManger.KILL_FLAG,
                first_step=FIRST_STEP):
    config = Config.Config(config_path)
    config.load_config()
    update_dict = {"DEFAULT": {"input_dir": input_dir, "output_dir": output_dir,
                                                                 'timestamp': date.today().isoformat(),
                                                                 'log_dir': logging_dir}}
    update_dict["DEFAULT"].update(c_args)
    config_copy = config.create_view(added_sections=update_dict,
                                     conf_file_dir=logging_dir)

    steps = read(config_copy)

    run_steps(first_step, log_file, kill_flag, steps, print_only)

    logging.info(SESSION_END_MSG)

    os.remove(config_copy.file_name)

if __name__ == '__main__':

    desc = "A Session manager enabling running of multiple pipelines and processing steps easily"

    parser = argparse.ArgumentParser(prog='SessionManager', description=desc)

    parser.add_argument('-i', '--input_dir', metavar="input path", dest='runon', nargs='?', required=True,
                        help='The input dir to run on.')
    parser.add_argument('-o', '--output_dir', metavar="output path", dest='out_path', nargs='?', required=True,
                        help='Outputer dir for the output')
    parser.add_argument('-c', '--config_path', metavar="config path", dest='config_path', nargs='?', required=True,
                        help='The path to the session configuration file')
    parser.add_argument('-l', '--logging_dir', metavar="log dir path", dest='logging_dir', nargs='?', required=False,
                        help='Outputer dir for logging', default=".")
    parser.add_argument('-k', '--kill_flag', metavar="kill_flag", dest='kill_flag', nargs='?', required=False,
                        default="./killme", help="The path of the kill flag")
    parser.add_argument('-a', '--args', metavar="config extra args", dest='args', required=False, default={}, nargs='*',
                        help='named args (in the format of <var>=\"<val>\",<var>=\"<val>\") to set in the config.s')
    parser.add_argument('-p', '--print_only', dest="print_only", action='store_true',
                        required=False, help="If set will just print the run "
                                             "commands for each step rather than running them")

    options = parser.parse_args()
    c_args = convert_args_to_dict(options.args)

    log_file = os.path.join(options.logging_dir, SESSION_LOG_FILE_FORMAT % date.today().isoformat())


    init_logging_dict(log_file)

    logging.info(SESSION_START_MSG % dict(input_dir=options.runon, output_dir=options.out_path,
                                          logging_dir=options.logging_dir, config=options.config_path, c_args=c_args))

    run_session(input_dir=options.runon, output_dir=options.out_path, config_path=options.config_path,
                logging_dir=options.logging_dir, log_file=log_file, c_args=c_args,
                print_only=options.print_only, kill_flag=options.kill_flag)

