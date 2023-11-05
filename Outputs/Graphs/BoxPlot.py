__author__ = 'Hillel'
# =====================imports=====================#
import plotly.offline as py
import plotly.graph_objs as go
import numpy

py.init_notebook_mode()
# =====================constants===================#
#<Change Here According To Your CSV Outputer>
PATH = r"C:\Users\User\Documents\dropbox\Dropbox\Genomics\AutoImmun\Several\AluIndex.csv"
TITLE = 'Editing In Several AID'
GROUP_COL_NAME = "group"
SUB_GROUPS = True
EXCLUDE = ["united_group"]
SELECTED = []
SUB_GROUP_SEP = "_"
# =====================functions===================#

def get_arrays(path, group_col, exclude_cols=EXCLUDE, selected_cols=SELECTED):
    with open(path) as input_f:
        res = {}
        array = numpy.genfromtxt(input_f, names=True, delimiter=",", dtype=None)
        run_on_cols = selected_cols if selected_cols else array.dtype.names
        for vals_col in run_on_cols:
            if vals_col in exclude_cols or vals_col == group_col:
                continue
            groups = {}
            for group in set(array[group_col]):
                groups[group] = [c[0] for c in zip(array[vals_col].tolist(),array[group_col]) if c[1] == group]
            res[vals_col] = groups
    return res

def get_sub_groups(data):
    res = {}
    sub_groups = {}
    for sub_group in data[data.keys()[0]]:
        sub = sub_group.split(SUB_GROUP_SEP)
        res[sub_group] = SUB_GROUP_SEP.join(sub[:-1])
        sub_groups[sub_group] = sub[-1]

    return res, sub_groups


def generate_boxpilot_from_csv(path, sub_groups=False):
    res = {}
    my_data = get_arrays(PATH, GROUP_COL_NAME,exclude_cols=EXCLUDE,selected_cols=SELECTED)

    if sub_groups:
        sub_groups_d, sub_gs = get_sub_groups(my_data)
        groups = list(set(sub_groups_d.values()))
        sub_g_i = sorted(sub_gs.values())
        colors = ['hsl('+str(h)+',50%'+',50%)' for h in numpy.linspace(0, 255, len(groups))]
        c_dict = {}
        for i,c in enumerate(colors):
            c_dict[groups[i]] = c

    for val_col in my_data:
        boxes = []
        scatters = []
        title = TITLE + "(%s)" % val_col
        if sub_groups:
            runon = sorted(my_data[val_col],key=lambda e: [sub_g_i.index(sub_gs[e]), sub_groups_d[e]])
            for group in runon:
                boxes.append(go.Box(
                    name=group,
                    y=my_data[val_col][group],
                    boxpoints='all',
                    whiskerwidth=0.2,
                    marker=dict(size=4, color=c_dict[sub_groups_d[group]]),
                    line=dict(width=1),
                    legendgroup=sub_groups_d[group]

                ))
            for s_group in sub_g_i:
                x = []
                y  =[]
                for group in runon:
                    if sub_gs[group] == s_group:
                        x.append(group)
                        y.append(numpy.median(my_data[val_col][group]))
                scatters.append(go.Scatter(
                    name=s_group,
                    y=y,
                    x=x,
                    legendgroup=sub_groups_d[group],
                    showlegend=False,
                    hoverinfo="none",
                    line=dict(width=0.5)
                ))

            size_factor = len(boxes)
            graphs = boxes + scatters

        else:
            for group in sorted(my_data[val_col]):
                boxes.append(go.Box(
                    name=group,
                    y=my_data[val_col][group],
                    boxpoints='all',
                    whiskerwidth=0.2,
                    marker=dict(size=4,),
                    line=dict(width=1),
                ))
            size_factor = len(boxes)
            graphs = boxes



        layout = go.Layout(
            title=title,
            yaxis=dict(
                autorange=True,
                showgrid=True,
                zeroline=True,
                gridcolor='rgb(30, 30, 30)',
                gridwidth=1,
                zerolinecolor='rgb(30, 30, 30)',
                zerolinewidth=4,
            ),
            xaxis=dict(autorange=True,),

            margin=dict(
                l=40,
                r=30,
                b=0,
                t=100,
            ),
            paper_bgcolor='rgb(250, 250, 250)',
            plot_bgcolor='rgb(240, 255, 250)',
            showlegend=True,
            width = size_factor*40 if size_factor*40 > 1200 else 1200,
            height = size_factor*25 if size_factor*25 > 350 else 350
        )

        fig = go.Figure(data=graphs, layout=layout)
        py.iplot(fig)




generate_boxpilot_from_csv(PATH,SUB_GROUPS)