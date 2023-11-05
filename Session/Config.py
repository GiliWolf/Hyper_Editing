"""
config file for the running of the Pipeline procedure
Michal Barak
13-12-2011
Shalom Hillel Roth
17-10-2017
"""

import ConfigParser
import logging
import time
from tempfile import mkstemp

from ConfigConsts import *
from ConfigFunctions import condition


class Config(ConfigParser.SafeConfigParser):
    def __init__(self, file_name=r"Pipeline.conf"):
        ConfigParser.SafeConfigParser.__init__(self)
        self.__file_name = file_name

    def set_config_file_name(self, file_name):
        self.__file_name = file_name

    @property
    def file_name(self):
        return self.__file_name

    def load_config(self, update_time=True):
        read_ok = self.read(self.__file_name)
        assert read_ok, "Could Not Load Config File!"
        if update_time:
            timeadd = time.strftime("%d-%m-%Y-%H", time.gmtime())
            self.set_config_file_name(timeadd + ".cnf")

    def create_view(self, dropped_sections=[], added_sections={}, conf_file_dir=None):
        """
        This function creates a "view" copy of the config with the given changes, first dropped then added, overriding.
        :param dropped_sections: A list of section to drop from the copy.
        :type dropped_sections: C{list} of C{str}
        :param added_sections: A dictionary of the the wanted values to add in the format of section=>option=>value.
        :type added_sections: C{dict} of C{str} to C{dict} of C{str} to C{str}
        :param conf_file_dir: If given, where to save the temp conf.
        :type conf_file_dir: C{str}
        :return: The view copy created.
        :rtype: L{Config}
        """
        copy = self.deepcopy(conf_file_dir)
        for dropped_section in dropped_sections:
            copy.remove_section(dropped_section)

        for section in added_sections:
            for option in added_sections[section]:
                copy.set(section=section, option=option, value=added_sections[section][option])


        copy.WriteConfig()
        copy.load_config(False)
        return copy

    def deepcopy(self, dir=None):
        s_fname = self.__file_name
        handle, path = mkstemp(dir=dir, prefix=self.__file_name + "copy" + time.strftime("%d-%m-%Y-%H", time.gmtime()))
        self.set_config_file_name(path)
        self.WriteConfig()
        self.set_config_file_name(s_fname)

        copy = Config(file_name=path)
        copy.load_config(False)

        return copy

    def WriteConfig(self):
        with open(self.__file_name, 'w') as configfile:
            self.write(configfile)

    def getboolean_option_soft_case(self, step_name, option_name):
        """
        wraps the getboolean function and enables ignoring of case in option name
        :param str step_name: The name of the step (section)
        :param str option_name: The name of the option
        :return: The value stored
        :rtype: bool
        """
        try:
            as_is = self.getboolean(step_name, option_name)
        except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
            as_is = False

        try:
            upper = self.getboolean(step_name, option_name.upper())
        except(ConfigParser.NoSectionError, ConfigParser.NoOptionError):
            upper = False
        try:
            lower = self.getboolean(step_name, option_name.lower())
        except(ConfigParser.NoSectionError, ConfigParser.NoOptionError):
            lower = False

        return as_is or upper or lower

    def enabled(self, step_name):
        if not self.getboolean_option_soft_case(step_name, ENABLE_SECTION_OPTION):
            logging.info("Step %s Was Not Enabled, Skipping." % step_name)
            return False
        if self.has_option(step_name, RECOVER_SECTION) and self.getboolean(step_name, RECOVER_SECTION):
            return False
        if self.has_option(step_name, CONSTRAIN_COND_OPTION) and condition(self.get(step_name, CONSTRAIN_COND_OPTION)):
            logging.info("Step %s's Constraint (%s) Was Met, Skipping." % (
                step_name, self.get(step_name, CONSTRAIN_COND_OPTION)))
            return False
        return True

    def get_error_step(self, curr_step):
        """
        This function retrieves the error step for a step.
        :param curr_step: The step to retrieve for.
        :return: The name of the error step
        """
        try:
            error_step = self.get(curr_step, ERROR_STEP_OPTION)
        except ConfigParser.NoOptionError:
            # if error step wasn't specified "inject" the default value (even if the step
            # doesn't exist, this still have logging importance).
            error_step = ERROR_STEP

        return error_step

    def get_next_step(self, curr_step):
        """
        This function retrieves the next step for a given step.
        :param curr_step: The current step.
        :return: The next step according to the condition(s)
        """
        error_step = self.get_error_step(curr_step)

        # get the ordered pairs of next steps and conditions.
        next_s_pairs = [s.strip() for s in self.get(curr_step, NEXT_STEP_OPTION).split(NEXT_STEPS_SEPARATOR)]

        for step_condition_pair in next_s_pairs:
            split_res = [s.strip() for s in step_condition_pair.split(STEP_CONDITION_SEPARATOR)]
            if len(split_res) == 1:
                return split_res[0]
            step = split_res[0]
            con = self.get(curr_step, split_res[1])
            if condition(con):
                return step

        return error_step

    def get_all_next_steps(self, curr_step):
        """
        This function return all possible following steps for a given step (and their matching conditions).
        :param curr_step: The step to get the next steps for.
        :return: A list of matched pairs of a next step and its corresponding condition.
        """
        error_step = self.get_error_step(curr_step)

        # get the ordered pairs of next steps and conditions.
        next_s_pairs = [s.strip() for s in self.get(curr_step, NEXT_STEP_OPTION).split(NEXT_STEPS_SEPARATOR)]

        next_steps = []
        conditions = []
        for step_condition_pair in next_s_pairs:
            split_res = [s.strip() for s in step_condition_pair.split(STEP_CONDITION_SEPARATOR)]
            next_steps.append(split_res[0])
            if len(split_res) > 1:
                try:
                    conditions.append(self.get(curr_step, split_res[1]))
                except ConfigParser.InterpolationError, e:
                    logging.warning("%s\nWhile pre-loading the configuration an option was missing!, please check your "
                                    "config file (if it's overridden default value(s)the pipeline will run)" % e.message)
            else:
                conditions.append(STEP_SEMI_CONDITION)

            if error_step != curr_step:
                next_steps.append(error_step)
                conditions.append(STEP_SEMI_CONDITION)

        return zip(next_steps, conditions)

    def get_all_steps(self, root_step=FIRST_STEP):
        """
        This recursive function retrieves all the steps (unlike just sections) in a config.
        :param root_step: The step to start the recursion with.
        :return: A list of all the possible steps.
        """
        if self.has_section(root_step):
            steps = [root_step]
            next_steps = self.get_all_next_steps(root_step)
            for step, con in next_steps:
                if self.has_section(step):
                    steps.append(step)
                    steps.extend(self.get_all_steps(step))

            return set(steps)


