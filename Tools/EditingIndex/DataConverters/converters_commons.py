__author__ = 'Hillel'

# =====================imports=====================#
import os
import subprocess
# =====================functions=====================#


def pyi_safe_subprocess(popen_cmd):
    env = dict(os.environ)  # make a copy of the environment
    lp_key = 'LD_LIBRARY_PATH'  # for Linux and *BSD.
    lp_orig = env.get(lp_key + '_ORIG')
    if lp_orig is not None:
        env[lp_key] = lp_orig  # restore the original, unmodified value
    else:
        # This happens when LD_LIBRARY_PATH was not set.
        # Remove the env var as a last resort:
        env.pop(lp_key, None)

    return subprocess.check_output(popen_cmd, env=env, shell=True).strip().split("\n")


