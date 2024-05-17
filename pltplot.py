def pltplot(x, y, title="", xlabel="", ylabel="", legend_label="", color="blue", show_legend=False, show_grid=False, grid_style="dotted"):
    """
    Creates a plot with extensive customization options including line color.

    Parameters:
    x (iterable): Data for the x-axis.
    y (iterable): Data for the y-axis.
    title (str): Title of the plot.
    xlabel (str): Label for the x-axis.
    ylabel (str): Label for the y-axis.
    legend_label (str): Label for the legend.
    color (str or tuple): Color of the plot line. Can be a string or RGB tuple.
    show_legend (bool): Whether to show the legend.
    show_grid (bool): Whether to show grid lines.
    grid_style (str): Style of the grid lines.
    """
    import matplotlib.pyplot as plt    # Create the plot

    plt.plot(x, y, label=legend_label, color=color)

'''
# Set title and labels
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Optional legend
    if show_legend:
        plt.legend()

    # Optional grid
    if show_grid:
        plt.grid(True, linestyle=grid_style)
    plt.show()
'''
