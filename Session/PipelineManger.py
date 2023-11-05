import sha

__author__ = 'Hillel'
"""
This module is based on Michal Barak's script - Pipeline.
"""
# =====================imports=====================#
# region Builtin Import
import logging.config
import multiprocessing
import argparse
import subprocess
from ConfigParser import NoOptionError
from datetime import datetime
from re import compile
import sys
import os
import operator

# endregion

# region Internal Imports
if __name__ == '__main__':
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from Config import *
from ConfigConsts import MAX_STEP_PROCESSES_OPTION, OVERALL_MAX_PS_OPTION, STEP_TYPE_OPTION, \
    STEP_NAME_OPTION, PROGRAM_NAME_OPTION, PROGRAM_PARAMS_OPTION, RECOVER_SECTION, FIRST_STEP
from Commons.general_functions import init_logging_dict, convert_args_to_dict, add_get_paths_function_to_argparser, \
    get_paths, get_path_depth, make_recursive_output_dir, get_sample_name, convert_params_to_bool_dict

from Commons.consts import ALL

# endregion

# =====================constants===================#
# A string format for the subprocesses names.
PIPELINE_SUB_PS_NAME = "Pipline_%s"

# ---defaults---
# The default value of the maximal number of processes running simultaneously.
MAX_PROCESSES_NUM = 10
# Default values for kill flag
KILL_FLAG = os.path.join(".", "end")
# log file's default path
LOG_DIR = os.path.join(".", "PipelineLogs", "pipeline %s.log")

# ---logging string formats---
START_STEP = "Process: %(ps_name)s; Started Step: %(step)s - %(name)s"
RUNNING_STEP = "Process: %(ps_name)s; Running Step: %(step)s"
ERROR_STEP_MSG = "Process: %(ps_name)s; Going To Error Step: %(step)s"
END_STEP = "Process: %(ps_name)s; Finished Step: %(step)s"
STARTED_RUNING = "Started Running"
FILE_FLOW_START = "Started Processing File %(file)s"
PS_END = "Process: %(ps_name)s; Exited With Code: %(ret_code)s"
RAN_ON_MSG = "Ran (Or Running) On %s Files So Far"

KW_ARGS_RE = compile("(?:^ *\w*=.*$)|(?:^$)")


# =====================functions===================#
def RunSteps(config, sema_dict, first_step_name=FIRST_STEP, kill_flag=KILL_FLAG, disabled_steps={}):
    # type: (Config, dict, str , str, dict) -> None

    """
    This function is a single linear subprocesses following the pipeline
    :param config: The configuration obejct.
    :type config: L{Config.Config}
    :param sema_dict: The semaphore to step dict.
    :param first_step_name: The name of the first step.
    :param kill_flag: The path to kill flag.
    :return: None
    """
    step_name = first_step_name
    while True:
        try:
            if not config.has_section(step_name):
                logging.info("Next Step %s Doesn't Exist - Assuming Exit Step!" % step_name)
                break

            if os.path.exists(kill_flag):
                return -1
            logging.debug(multiprocessing.current_process().name)
            sema_dict[step_name].acquire()
            logging.debug(START_STEP % {'ps_name': multiprocessing.current_process().name, 'step': step_name,
                                        'name': config.get(step_name, STEP_NAME_OPTION)})
            retcode = run_step(config, step_name, disabled_steps)
            sema_dict[step_name].release()

            if retcode == 0:
                config.set(step_name, RECOVER_SECTION, "True")
                logging.debug(END_STEP % {'ps_name': multiprocessing.current_process().name, 'step': step_name})
                step_name = config.get_next_step(step_name)
                logging.debug("going to " + step_name + "\n")
            else:
                step_name = config.get(step_name, ERROR_STEP_OPTION)
                logging.error(ERROR_STEP_MSG % {'ps_name': multiprocessing.current_process().name, 'step': step_name})
        except Exception, e:
            err_step = config.get(step_name, ERROR_STEP_OPTION)
            logging.exception("Pipeline Encountered An Exception On %s! Going To Error Step %s" % (step_name, err_step))
            logging.error(ERROR_STEP_MSG % {'ps_name': multiprocessing.current_process().name, 'step': err_step})
            step_name = err_step

    sema_dict[OVERALL_MAX_PS_OPTION].release()
    logging.debug(PS_END % {"ret_code": str(retcode), "ps_name": multiprocessing.current_process().name})


def run_for_a_time(cmd_line, wall_time=36000):
    delay = 1.0
    timeout = int(wall_time / delay)

    task = subprocess.Popen(cmd_line, shell=True)

    # while the process is still executing and we haven't timed-out yet
    ret_code = task.poll()
    while ret_code is None and timeout > 0:
        time.sleep(delay)
        timeout -= delay
        ret_code = task.poll()

    if ret_code is None:
        task.kill()
        return 15
    return ret_code


def run_step(config, step_name, disabled_steps={}):
    """
    This function runs a signle step in the pipeline.
    :param config:The configuration obejct.
    :type config: L{Config.Config}
    :param step_name:
    :return:
    """
    retcode = 0
    if not config.enabled(step_name) or step_name in disabled_steps:
        return retcode
    logging.debug(RUNNING_STEP % {'ps_name': multiprocessing.current_process().name, 'step': step_name})
    if config.get(step_name, STEP_TYPE_OPTION) == "cmd":
        program_name = config.get(step_name, PROGRAM_NAME_OPTION)
        program_paras = config.get(step_name, PROGRAM_PARAMS_OPTION)
        try:
            wall_time = config.get(step_name, WALL_TIME_OPTION)
        except NoOptionError:
            wall_time = 60 * 60 * 12  # == 12h
        logging.debug(" ".join([step_name + ":", program_name, program_paras]) + '\n')
        retcode = run_for_a_time(program_name + " " + program_paras, wall_time)
        logging.debug("Step %s returned: %s\n" % (step_name, str(retcode)))
    elif config.get(step_name, STEP_TYPE_OPTION) == "mail":
        '''text = config.get(step, "text")
        sender = config.get(step, "from")
        to = config.get(step, "to")
        logging.debug("sending mail\n")
        # SendEmail(config, text, sender, to)'''

    return retcode


def Pipeline_mt(config, sema_dict, log_dir, disabled_steps={}):
    logging.basicConfig(level=logging.DEBUG,
                        format='[%(asctime)s] %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%y-%m-%d %H:%M:%S',
                        filename=os.path.join(log_dir, multiprocessing.current_process().name + ".log"),
                        filemode='w')
    outdir = config.get("DEFAULT", "output_dir")
    if config.has_option("DEFAULT", "kill_flag"):
        kill_flag = config.get("DEFAULT", "kill_flag")
    else:
        kill_flag = KILL_FLAG

    os.chdir(outdir)
    RunSteps(config=config, sema_dict=sema_dict, kill_flag=kill_flag, disabled_steps=disabled_steps)
    del config


def run_on_dir(config_file, root_path, outdir=None, include_paths=None, include_paths_operator=operator.or_,
               exclude_paths=None, exclude_paths_operator=operator.or_, recursion_depth=100, follow_links=False,
               log_dir=None, kill_flag=KILL_FLAG, files_suffix="", truncs=(0, 2), n_suff_trunc=-1,
               trunc_after_join=False,
               config_args={}, disabled_steps={}, no_extension=False):

    temp_conf_files = []  # for cleanup
    created_dirs_and_files = []  # for revert
    exception = None  # for exception handling

    config = Config(config_file)
    config.load_config()

    config_copy = config.create_view(added_sections=dict(DEFAULT=config_args))
    temp_conf_files.append(config_copy.file_name)

    semaphore_dict = create_semaphores(config_copy, config_copy.get_all_steps())
    predicates_dict = {ALL: lambda path: path.endswith(files_suffix)}
    files = get_paths(root_path=root_path, must_include_paths=include_paths,
                      must_include_operator=include_paths_operator, exclude_paths=exclude_paths,
                      exclude_operator=exclude_paths_operator, follow_links=follow_links,
                      recursion_depth=recursion_depth, predicates_dict=predicates_dict)[ALL]
    file_counter = 0

    if not os.path.isdir(os.path.join(log_dir, "flags")):
        os.makedirs(os.path.join(log_dir, "flags"))
        created_dirs_and_files.append(os.path.join(log_dir, "flags"))
    try:
        success = True
        for full_file_name in files:

            file_name = os.path.basename(full_file_name)
            logging.info(RAN_ON_MSG % file_counter)

            if not file_name.endswith(".flg"):

                if os.path.exists(kill_flag):
                    return -1
                short_file_name = get_sample_name(file_name=file_name,
                                                  n_suff_to_keep=truncs[0],
                                                  chars_to_trunc=truncs[1],
                                                  suff_to_trunc=n_suff_trunc,
                                                  after_joining=trunc_after_join)
                # Check if the file is already processed in another run.
                flag_name = os.path.join(log_dir, "flags", short_file_name + ".flg")
                if not os.path.exists(flag_name):  # for more than one run parallel
                    try:
                        open(flag_name, "w")
                    except:
                        logging.warn("Could Not Create Flag For %s, This May Cause Down-Stream Problems!" % file_name)
                else:
                    logging.info("Flag For %s Exists, Skipping." % file_name)
                    continue  # The file is already processed in another run.

                logging.debug(FILE_FLOW_START % {'file': file_name})
                root = os.path.dirname(full_file_name)
                new_o_dir = make_recursive_output_dir(outdir=outdir,
                                                      root=root,
                                                      extension=short_file_name if not no_extension else "",
                                                      depth=get_path_depth(full_file_name) - get_path_depth(root_path))

                if not os.path.exists(new_o_dir):
                    os.makedirs(new_o_dir)
                    created_dirs_and_files.append(new_o_dir)

                updated_options = {"DEFAULT": {"input_dir": root,
                                               "file_name": short_file_name,
                                               "output_dir": new_o_dir,
                                               "input_file": full_file_name}
                                   }
                updated_options["DEFAULT"].update(config_args)
                config_copy = config.create_view(added_sections=updated_options)
                temp_conf_files.append(config_copy.file_name)

                semaphore_dict[OVERALL_MAX_PS_OPTION].acquire()

                tr = multiprocessing.Process(target=Pipeline_mt,
                                             args=(config_copy, semaphore_dict, log_dir, disabled_steps),
                                             name=PIPELINE_SUB_PS_NAME % short_file_name)

                tr.daemon = True
                tr.start()
                file_counter += 1
            active = multiprocessing.active_children()
            for tr in active:
                tr.join(timeout=1)

        active = multiprocessing.active_children()
        for tr in active:
            tr.join()

        logging.info("Delete The Following Files:\n%s" % "\t\n".join(temp_conf_files))
        '''ls -l --author /tmp|grep <user>|grep cnfcopy|awk '{print "/tmp/" $10}'|xargs rm'''
    except Exception as e:
        success = False
        exception = sys.exc_info()
    finally:
        return temp_conf_files, created_dirs_and_files, success, exception


def create_semaphores(config, steps):
    """
    This function creates semaphores for all the steps according to what was defined on each step's configuration.
    :param config: The L{Config} instace used for this run.
    :param steps: A list of all the steps in the config.
    :return: A dictionary containing all the semaphores mapped to the step name.
    :rtype: C{dict} of C{str} to C{multiprocessing.Semaphore}
    """
    sema_dict = dict()

    for step in steps:
        try:
            max_ps_num = config.get(step, MAX_STEP_PROCESSES_OPTION)  # get max num. of processes for the step
        except NoOptionError:
            max_ps_num = MAX_PROCESSES_NUM
        sema_dict[step] = multiprocessing.Semaphore(int(max_ps_num))  # init a semaphore with that number.
    try:
        max_ps_num = config.get(step, OVERALL_MAX_PS_OPTION)  # get max num. of processes for the step
    except NoOptionError:
        max_ps_num = MAX_PROCESSES_NUM
    sema_dict[OVERALL_MAX_PS_OPTION] = multiprocessing.Semaphore(int(max_ps_num))  # init a semaphore with that number.

    return sema_dict


# doesn't work yet
# def SendEmail(config, text, sender, to):
#     try:
#         if success:
#             success_str = "succsesfuly"
#         else:
#             success_str = "with errors"
#         msg = MIMEText("shalom")
#         msg['Subject'] = text
#         msg['From'] = sender
#         msg['To'] = to
#         # Send the message via our own SMTP server, but don't include the
#         # envelope header.
#         s = smtplib.SMTP("localhost")
#         s.sendmail(sender, to, msg.as_string())
#         s.quit()
#     except:
#         outfile.write("error sending mail\n")


def main(config_file, root_dir, output_dir=None, include_paths=None, include_paths_operator=operator.or_,
         exclude_paths=None, exclude_paths_operator=operator.or_, recursion_depth=100, follow_links=False,
         log_path=None, kill_flag=KILL_FLAG, files_suffix="", truncs=(0, 2), n_suff_trunc=-1, trunc_after_join=False,
         config_args={}, pipe=None, disabled_steps={}, no_extension=False):
    # type: (str, str, str, list, operator, list, operator, int, bool, str, str, str, list, int, bool, dict, multiprocessing.Pipe, dict, bool) -> None

    l_path = log_path if log_path else LOG_DIR % datetime.now().strftime("%y-%m-%d")
    to_del_files = []
    if not os.path.isdir(l_path):
        os.makedirs(l_path)

    timestamp = datetime.today().isoformat()

    init_logging_dict(os.path.join(l_path, multiprocessing.current_process().name + "." + timestamp + ".log"))

    logging.info(STARTED_RUNING)

    kill_flag = kill_flag if kill_flag else KILL_FLAG
    if config_file:
        if root_dir:
            to_del_files, created_dirs_and_files, success, exception = run_on_dir(root_path=root_dir,
                                                                                  outdir=output_dir,
                                                                                  include_paths=include_paths,
                                                                                  include_paths_operator=include_paths_operator,
                                                                                  exclude_paths=exclude_paths,
                                                                                  exclude_paths_operator=exclude_paths_operator,
                                                                                  recursion_depth=recursion_depth,
                                                                                  follow_links=follow_links,
                                                                                  config_file=config_file,
                                                                                  log_dir=log_path,
                                                                                  kill_flag=kill_flag,
                                                                                  files_suffix=files_suffix,
                                                                                  truncs=truncs,
                                                                                  n_suff_trunc=n_suff_trunc,
                                                                                  trunc_after_join=trunc_after_join,
                                                                                  config_args=config_args,
                                                                                  disabled_steps=disabled_steps,
                                                                                  no_extension=no_extension)

            logging.info("Finished Running Pipleline %s" % str(locals()))
    else:
        logging.warn("Failed Running Pipleline %s" % str(locals()))

    if exception:
        logging.error("Fatal Exection Occured While Running Pipeline!", exc_info=exception)
    if pipe:
        pipe[0].send((to_del_files, created_dirs_and_files, success))
        pipe[0].close()
        pipe[1].close()


# def convert_args_to_dict(arg):
#     assert all([KW_ARGS_RE.match(a) for a in arg.split(",")])
#     return eval("dict(" + arg + ")")

if __name__ == "__main__":
    desc = """Run on each file in a given directory a set of steps."""


    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(prog='Pipeline', description=desc, formatter_class=MyFormatter)

    add_get_paths_function_to_argparser(parser=parser)
    parser.add_argument('-f', '--files_suffix', metavar="files_suffix", dest='files_suffix', nargs='?', required=False,
                        default="", help="A suffix of the files to run on.")
    parser.add_argument('-l', '--log_path', metavar="log_path", dest='log_path', nargs='?', required=True,
                        help="The path where the logs (and flags) will be written.")
    parser.add_argument('-o', '--output_dir', metavar="output_dir", dest='output_dir', nargs='?', required=False,
                        default=None, help="The root directory for the output. "
                                           "(if not specified will be  taken from the config)")

    parser.add_argument('-c', '--config_file', metavar="config_file", dest='config_file', nargs='?', required=True,
                        help="The configuration file's path.")

    parser.add_argument('-k', '--kill_flag', metavar="kill_flag", dest='kill_flag', nargs='?', required=False,
                        default=KILL_FLAG, help="The path of the kill flag")

    parser.add_argument('-t', '--truncs', metavar="truncs", dest='truncs', nargs='?', required=False,
                        default="0,2", help="The truncation of the file name "
                                            "<suffixes to retain>,<chars to drop from filename>")

    parser.add_argument('-tn', metavar="n suffix to trunc", dest='n_suff_trunc', nargs='?', required=False,
                        default=-1, type=int, help="The index of the suffix to truncate to get the sample name.")

    parser.add_argument('--trunc_after_join', dest='trunc_after_join', required=False, action='store_true',
                        help='If set, will truncate file name *after* joining back the suffixes, this functionally nullifies -tn.')

    parser.add_argument('--no_extensions', dest='no_extension', action='store_true', required=False,
                        help="")
    parser.add_argument('-a', '--args', metavar="config extra args", dest='args', required=False,
                        default={}, nargs='*',
                        help='named args (in the format of <var>=\"<val>\",<var>=\"<val>\") to set in the config.s')
    parser.add_argument('--disbale', metavar="disable steps", dest='disabled_steps', required=False,
                        default=[], nargs='*', help='a list of steps to disable')

    options = parser.parse_args()

    args = convert_args_to_dict(options.args)
    truncs = [int(options.truncs.split(",")[0]), int(options.truncs.split(",")[1])]
    main(root_dir=options.root_dir,
         output_dir=options.output_dir,
         include_paths=options.include_prefixes,
         include_paths_operator=options.include_operator,
         exclude_paths=options.exclude_prefixes,
         exclude_paths_operator=options.exclude_operator,
         recursion_depth=options.recursion_depth,
         follow_links=options.follow_links,
         config_file=options.config_file,
         log_path=options.log_path,
         kill_flag=options.kill_flag,
         files_suffix=options.files_suffix,
         truncs=truncs,
         n_suff_trunc=options.n_suff_trunc,
         trunc_after_join=options.trunc_after_join,
         config_args=args,
         disabled_steps=convert_params_to_bool_dict(options.disabled_steps),
         no_extension=options.no_extension
         )

# =====================classes=====================#
