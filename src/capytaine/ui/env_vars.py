import os


def look_for_boolean_var(var_name, default):
    """Return the value of the environment variable, else default."""
    if var_name in os.environ:
        env_var = os.environ[var_name].lower()
        if env_var in {'true', '1', 't'}:
            return True
        elif env_var in {'false', '0', 'f'}:
            return False
        else:
            raise ValueError("Invalid value '{}' for the environment variable {}.".format(os.environ[var_name], var_name))
    else:
        return default
