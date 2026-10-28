"""Figure 1 panels A and B (carried from figures/make_headline.py:render)."""


def render(pipeline, panels):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FuncFormatter, FixedLocator, NullLocator
    plt = pipeline.plt
    blue, orange, grey, black = (pipeline.BLUE, pipeline.ORANGE, pipeline.GREY, pipeline.BLACK)
    green = '#009E73'
    colors = {'retained': blue, 'second_pass_only': green, 'first_pass_only_in_library': orange, 'first_pass_only_outside_library': grey}
    fig = plt.figure(figsize=(12.2, 11.4))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.27, 1], width_ratios=[2.45, 1.15])
    ax = fig.add_subplot(grid[0, :])
    minor_labels = []
    for i, row in enumerate(panels['A']):
        counts = row['counts']
        for j, phase in enumerate(['first', 'final']):
            x = 3 * i + j
            start = 0
            categories = ['retained', 'first_pass_only_in_library', 'first_pass_only_outside_library'] if phase == 'first' else ['retained', 'second_pass_only']
            for category in categories:
                value = counts[category]
                ax.bar(x, value, bottom=start, width=0.77, color=colors[category], edgecolor='white', linewidth=0.6, hatch='//' if category == 'first_pass_only_outside_library' else None)
                if category == 'first_pass_only_in_library':
                    ax.text(x, start + value / 2, f'{value:,}', color='white', ha='center', va='center', fontsize=10)
                start += value
            ax.text(x, counts[phase] + 7500, f'{counts[phase]:,}', ha='center', va='bottom', fontsize=12, fontweight='bold', color=black)
            represented = counts[phase + '_represented_precursors']
            phase_label = 'First pass' if phase == 'first' else 'Final'
            minor_labels.append(f'{phase_label}\nRepresented\nprecursors:\n{represented:,}')
    ax.set(ylim=(0, 500000), xlim=(-0.7, 7.7), ylabel='Accepted precursor-cell observations')
    ax.set_title('A  More accepted observations, fewer represented precursors', loc='left', fontsize=15)
    ax.set_xticks([0.5, 3.5, 6.5], [row['label'].replace(' (', '\n(') for row in panels['A']])
    ax.set_xticks([0, 1, 3, 4, 6, 7], minor_labels, minor=True)
    ax.tick_params(axis='x', which='major', length=0, pad=75, labelsize=12)
    ax.tick_params(axis='x', which='minor', length=0, pad=5, labelsize=10)
    ax.set_yticks([0, 100000, 200000, 300000, 400000])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: '0' if not v else f'{v / 1000:g}k'))
    ax.legend(handles=[Patch(facecolor=colors[k], edgecolor='white', hatch='//' if 'outside' in k else None, label=v) for k, v in [('retained', 'Retained'), ('second_pass_only', 'Second-pass-only'), ('first_pass_only_in_library', 'First-pass-only: in library'), ('first_pass_only_outside_library', 'First-pass-only: outside library')]], frameon=False, ncol=2, loc='upper left', columnspacing=1.8, fontsize=11.5)
    forest = fig.add_subplot(grid[1, 0])
    values = fig.add_subplot(grid[1, 1], sharey=forest)
    values.axis('off')
    for i, row in enumerate(panels['B']):
        for j, entry in enumerate(row['contrasts']):
            y = i + (-0.17 if j == 0 else 0.17)
            color, marker = (green, 'o') if j == 0 else (orange, 's')
            pipeline.interval(forest, y, entry, 'point', 'low95', 'high95', color=color, marker=marker)
            text = f"{float(entry['point']):.3f} [{float(entry['low95']):.3f}, {float(entry['high95']):.3f}]"
            values.text(0, y, text, color=color, va='center', fontsize=11.5)
    labels = [row['label'].replace(': ', ':\n', 1).replace(' (pooled F3)', '\n(pooled F3)') for row in panels['B']]
    forest.set(xscale='log', xlim=(0.3, 30), ylim=(3.65, -0.9), yticks=list(range(4)), yticklabels=labels, xlabel='Wrong-organism fraction ratio (class / retained)')
    forest.set_title('B  Class-relative fractions (two-way 95% intervals)', loc='left', fontsize=15, pad=43)
    forest.axvline(1, color=grey, linestyle='--', linewidth=1)
    forest.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10, 20]))
    forest.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
    forest.xaxis.set_minor_locator(NullLocator())
    forest.tick_params(axis='y', labelsize=12)
    forest.legend(handles=[Line2D([], [], color=green, marker='o', linestyle='none', label='Second-pass-only / retained'), Line2D([], [], color=orange, marker='s', linestyle='none', label='First-pass-only / retained')], loc='lower left', bbox_to_anchor=(-0.01, 1.0), frameon=False, fontsize=11, ncol=1, borderpad=0)
    values.text(0, -0.72, 'Ratio [two-way 95% interval]', fontsize=11.5, fontweight='bold', va='center')
    return fig
