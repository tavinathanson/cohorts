import os
import subprocess

def call_R(r_script_path, r_command_line_args=""):
    """
    Parameters
    _____
    R_script_path : str
    R_command_line_args : str
    """
    if os.name == "posix":
        command = "R -f {} {}".format(r_script_path, r_command_line_args)
        return subprocess.check_output(command, shell=True)
    else:
        proc_list = ['R', '-f', r_script_path] + r_command_line_args.split(" ")
        theproc = subprocess.Popen(proc_list)
        return theproc.communicate()

def write_r_template_to_file(r_template, path_to_write_r_script=None):
    """
    Writes the R template to file at the given path
    TO DO: any way to make this more efficient so as not to overwrite every
     time(check hashes?)

    Parameters
    _______________
    r_template : str
    path_to_write_r_script : str
    """
    if path_to_write_r_script is None:
        home_dir = os.path.expanduser('~')
        output_dir = "{}/tmp_R".format(home_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        final_path = "{}/r_template.R".format(output_dir)
    else:
        final_path = path_to_write_r_script

    with open(final_path, "w") as f:
        f.write(r_template)
    print "Wrote to {}".format(final_path)
    return final_path