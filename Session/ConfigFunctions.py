from Commons.consts import LOGICAL_OPERATORS
from ConfigConsts import IN_CONDITION_SEPARATOR, CONDITIONS


def condition(cond):
    parts = cond.split(IN_CONDITION_SEPARATOR)
    i = 0
    res = None
    while i < len(parts):
        if parts[i] == IN_CONDITION_SEPARATOR:
            continue  # it's a typo
        if parts[i] in LOGICAL_OPERATORS:
            op = LOGICAL_OPERATORS[parts[i]]
            if parts[i + 1] in LOGICAL_OPERATORS:
                op2 = LOGICAL_OPERATORS[parts[i + 1]]
                con = CONDITIONS[parts[i + 2]]
                param = parts[i + 3]
                i += 4
                if res is None:
                    res = op(op2(con(param)))
                else:
                    res = op(res, op2(con(param)))
            else:
                con = CONDITIONS[parts[i + 1]]
                param = parts[i + 2]
                i += 3
                if res is None:
                    res = op(con(param))
                else:
                    res = op(res, con(param))
        else:
            con = CONDITIONS[parts[i]]

            while (parts[i + 1] == '') and i < len(parts) - 1:
                i += 1

            param = parts[i + 1]
            res = con(param)
            i += 2

    return res
