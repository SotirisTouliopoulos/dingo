# Pre- and post-sampling features to leverage flux sampling at both the strain and the community level

#### A contribution for the Google Summer of Code 2024 program


## Overall

#### A summary of the implemented methods, merged into the dingo library:

- preprocess for the reduction of metabolic models.
- inference of pairwise correlated reactions.
- visualization of a steady-states correlation matrix.
- construction of a weighted graph of the model's reactions with the correlation coefficients as weights.
- annotation of these weights to the metabolic model and extraction to an annotated SBML file.


## Preprocess

- Large metabolic models contain many reactions and metabolites. 
- Sampling in the flux space of such models requires an increased computational time due to the higher number of dimensions.
- Preprocessing of models can deal with this problem, as it removes reactions and thus decreases dimensional space.

#### The concept of reduction of metabolic models is illustrated in the following figure:

![Network Reduction Concept](/img/reduction.png)

I implemented a `PreProcess` class that identifies and removes 3 types of reactions:

- Blocked reactions: cannot carry a flux in any condition.
- Zero-flux reactions: cannot carry a flux while maintaining at least 90% of the maximum growth rate.
- Metabolically less-efficient reactions: require a reduction in growth rate if used.

Example of calling the `reduce` function of the `PreProcess` class to the E.coli core model:

    cobra_model = load_json_model("ext_data/e_coli_core.json")
    obj = PreProcess(cobra_model, tol=1e-6, open_exchanges=False)
    removed_reactions, final_dingo_model = obj.reduce(extend=False)

Explaining the parameters and the returned objects:

- `open_exchanges` is a parameter to the `find_blocked_reactions` function of the cobra library. It controls whether or not to open all exchange reactions to very high flux ranges.
- `removed_reactions` is a list that contains the names of the removed reactions. `final_dingo_model` is a reduced model with bounds of removed reactions set to 0.
- Users can decide if they want to remove an additional set of reactions, by setting the `extend` parameter to `True`.
  These reactions are the ones that do not affect the value of the objective function, when removed.

Reduction with the `PreProcess` class has been tested with the E. coli core and with other models used in dingo's publication article. 
A figure with results follows, which presents the number of remained reactions, after calling the `reduce` function with `extend` set both to `False` and `True`.
Reduction with `extend` set to `True` is based on sampling, so a slight change in the number of remained reactions may appear from time to time. 
It is recommended to use `extend` set to `False` in large models, due to the higher computational time.

![Reduction_Results_Plot](/img/reduction_results_plot.png)


## Correlated reactions

- Reactions in biochemical pathways, can be positive, negative or non correlated.
- Positive correlation between a set of reactions A and B means that if reaction A is active then reaction B is also active and vice versa.
- Negative correlation means that if reaction A is active, then reaction B is deactive.
- Zero correlation means that reaction A can have any status, regardless of the reaction's B status and vice versa.

I implemented a `correlated_reactions` function that calculates reactions steady states using dingo's `PolytopeSampler` class and creates a correlation matrix based on pearson correlation coefficient between pairwise set of reactions. This function also calculates a copula indicator for a specific set of reactions, to filter false positive correlations.

Example of calling the `correlated_reactions` function to the E.coli core model:

    dingo_model = MetabolicNetwork.from_json('ext_data/e_coli_core.json')

    sampler = PolytopeSampler(dingo_model)
    steady_states = sampler.generate_steady_states()
    corr_matrix = correlated_reactions(steady_states,
                                       pearson_cutoff = 0.99,
                                       indicator_cutoff=2,
                                       cells = 10,
                                       cop_coeff = 0.3,
                                       lower_triangle = False)

Explaining the parameters and the returned objects:

- `steady_states` are the steady states returned from dingo's `generate_steady_states` function.
- `pearson_cutoff` is a cutoff to replace all lower values with 0.
- `indicator_cutoff` is a cutoff to filter false positive correlations for correlations greater than the cutoff and replace them with 0.
- `cells` is a variable that defines the number of cells in the computed copulas
- `cop_coeff` is a variable that defines the width of copula's diagonal.
- `lower_triangle` is a `boolean` value that when `True` keeps only the lower triangular matrix. This can be useful for visualization purposes.
- `corr_matrix` is the calculated correlation matrix with dimensions equal to the number of reactions of the given model.

I also implemented a `plot_corr_matrix` function that visualizes the correlation matrix given as an input with a heatmap plot. In reduced models, there is an option to plot only the remained reactions names, if they are provided in a list.

Example of calling the `plot_corr_matrix` function:

    plot_corr_matrix(corr_matrix, reactions, format="svg")

Explaining the parameters:

- `corr_matrix` is a correlation matrix produced from the `correlated_reactions` function.
- `reactions` is a list with the reactions names that will appear as labels in the axes of the heatmap plot.
- `format` is a variable that defines the desired image saving format.

First we will examine the capabilities of the `correlated_reactions` function from the produced heatmap plots in the E. coli core model.

A heatmap from a symmetrical correlation matrix without pearson and indicator filtering:
![corr_matrix_no_cutoffs](/img/corr_matrix_no_cutoffs.png)

A heatmap from a symmetrical correlation matrix with `pearson_cutoff = 0.7` and without indicator filtering:
![corr_matrix_pearson](/img/corr_matrix_pearson.png)

A heatmap from a symmetrical correlation matrix with `pearson_cutoff = 0.7` and with `indicator_cutoff = 100`:
![corr_matrix_pearson_indicator](/img/corr_matrix_pearson_indicator.png)

A heatmap from a triangular correlation matrix with `pearson_cutoff = 0.7` and with `indicator_filtering = 100`:
![corr_matrix_pearson_indicator_triangle](/img/corr_matrix_pearson_indicator_triangle.png)

Now we will examine heatmap plots from reduced E. coli core models.

A heatmap from a symmetrical correlation matrix with `pearson_cutoff = 0.7` and with `indicator_filtering = 100`, from a reduced E. coli core model with `extend` set to `False`:
![corr_matrix_extend_false](/img/corr_matrix_extend_false.png)

A heatmap from a symmetrical correlation matrix with `pearson_cutoff = 0.7` and with `indicator_filtering = 100`, from a reduced E. coli core model with `extend` set to `True`:
![corr_matrix_extend_true.png](/img/corr_matrix_extend_true.png)

We can see that the additional reaction removed with `extend` set to `True` is `FRD7` which is the least correlated one across the matrix.
`FRD7` is a fumarate reductase, which in the E. coli core model appears in a loop with `SUCDi`, as seen in the following figure obtained from `ESCHER`:
![escher_frd7.png](/img/escher_frd7.png)


![dendrogram_no_cutoffs.png](/img/dendrogram_no_cutoffs.png)
![graph_both_no_cutoffs.png](/img/graph_both_no_cutoffs.png)

![graph_negative_pearson.png](/img/graph_negative_no_cutoffs.png)
![graph_positive_pearson.png](/img/graph_positive_no_cutoffs.png)

![dendrogram_pearson.png](/img/dendrogram_pearson.png)
![graph_both_pearson.png](/img/graph_both_pearson.png)

![subgraph_both_pearson.png](/img/subgraph_both_pearson.png)



