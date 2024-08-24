# Pre- and post-sampling features to leverage flux sampling at both the strain and the community level

#### A contribution for the Google Summer of Code 2024 program


## Overall

#### A summary of the implemented methods, merged into the dingo library:

- preprocess for the reduction of metabolic models.
- inference of pairwise correlated reactions.
- visualization of a steady-states correlation matrix.
- clustering of a steady-states correlation matrix.
- construction of a weighted graph of the model's reactions with the correlation coefficients as weights.


## Preprocess

- Large metabolic models contain many reactions and metabolites. 
- Sampling in the flux space of such models requires an increased computational time due to the higher number of dimensions.
- Preprocessing of models can deal with this problem, as it removes certain reactions and thus decreases dimensional space.

#### The concept of reduction in metabolic models is illustrated in the following figure created with BioRender <a href="#ref-3">[1]</a>:

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
- `removed_reactions` is a list that contains the names of the removed reactions. `final_dingo_model` is a reduced model with the bounds of removed reactions set to 0.
- Users can decide if they want to remove an additional set of reactions, by setting the `extend` parameter to `True`.
  These reactions are the ones that do not affect the value of the objective function, when removed.

Reduction with the `PreProcess` class has been tested with various models [2](#ref-2). Some of them are used in dingo's [3](#ref-3) publication article too. 
A figure with reduction results follows, which presents the number of remained reactions, after calling the `reduce` function with `extend` set both to `False` and `True`.
Reduction with `extend` set to `True` is based on sampling, so a slight change in the number of remained reactions may appear from time to time. 
It is recommended to use `extend` set to `False` in large models, due to the higher computational time.

![Reduction_Results_Plot](/img/reduction_results_plot.png)
![Reduction_Results_Plot2](/img/reduction_results_plot2.png)


## Correlated reactions

- Reactions in biochemical pathways, can be positive, negative or non correlated.
- Positive correlation between a set of reactions A and B means that if reaction A is active then reaction B is also active and vice versa.
- Negative correlation means that if reaction A is active, then reaction B is deactive and vice versa.
- Zero correlation means that reaction A can have any status, regardless of the reaction's B status and vice versa.

I implemented a `correlated_reactions` function that calculates reactions steady states using dingo's `PolytopeSampler` class and then creates a correlation matrix based on pearson correlation coefficient between pairwise set of reactions. This function also calculates a copula indicator for filtering correlations greater than the pearson cutoff.

Example of calling the `correlated_reactions` function to the E.coli core model:

    dingo_model = MetabolicNetwork.from_json('ext_data/e_coli_core.json')

    sampler = PolytopeSampler(dingo_model)
    steady_states = sampler.generate_steady_states()
    corr_matrix, indicator_dict = correlated_reactions(
                                  steady_states,
                                  pearson_cutoff = 0.99,
                                  indicator_cutoff=2,
                                  cells = 10,
                                  cop_coeff = 0.3,
                                  lower_triangle = False)

Explaining the parameters and the returned objects:

- `steady_states` are the reactions steady states returned from dingo's `generate_steady_states` function.
- `pearson_cutoff` is a cutoff to filter and replace all lower correlation values with 0.
- `indicator_cutoff` is a cutoff that corresponds to the copula's indicator. It filters correlations greater than the pearson cutoff.
- `cells` is a variable that defines the number of cells in the computed copulas.
- `cop_coeff` is a variable that defines the width of the copula's diagonal.
- `lower_triangle` is a `boolean` value that when `True` keeps only the lower triangular matrix. This can be useful for visualization purposes.
- `corr_matrix` is the calculated correlation matrix with dimensions equal to the number of reactions of the given model.
- `indicator_dict` is a dictionary containing filtered reactions combinations with the pearson's cutoff, alongside their copula's indicator value and a classification for the correlation.

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
`FRD7` is a fumarate reductase, which in the E. coli core model appears in a loop with `SUCDi`, as seen in the following figure obtained from `ESCHER` [4](#ref-4):
![escher_frd7.png](/img/escher_frd7.png)


## Clustering

- Clustering based on a dissimilarity matrix, reveals group of reactions with similar correlation values across the matrix.
- Reactions that belong to the same cluster might contribute to the same pathways.

I implemented a `cluster_corr_reactions` function that takes as an input a correlation matrix and calculates a dissimilarity matrix.
The dissimilarity matrix can be calculated in 2 ways, either substrating each value from 1, or substracting each absolute value from 1.

Example of calling the `cluster_corr_reactions` function:

    dissimilarity_matrix, labels, clusters = cluster_corr_reactions(
                                             corr_matrix,
                                             reactions=reactions,
                                             linkage="ward",
                                             t = 10.0,
                                             correction = True)

Explaining the parameters and the returned objects:

- `corr_matrix` is a correlation matrix produced from the `correlated_reactions` function.
- `reactions` is a list with the reactions names.
- `likage` defines the type of linkage. `single`, `average`, `complete`, `ward` are available linkage types.
- `t` defines a height threshold to cut the dendrogram and return the created clusters.
- `correction` is a boolean variable that if True, the dissimilarity matrix is calculated by substracting absolute values from 1.
- `dissimilarity_matrix` is the distance matrix created from the correlation matrix.
- `labels` are index labels that correspond to a specific cluster.
- `clusters` is a nested list, containing sublists with reactions that belong to the same cluster.

I also implemented a `plot_dendrogram` function, that plots a dendrogram given a dissimilarity matrix created from the `cluster_corr_reactions` function.

Example of calling the `plot_dendrogram` function:

    plot_dendrogram(dissimilarity_matrix, reactions , plot_labels=True, t=10.0, linkage="ward")

Explaining the parameters:

- `dissimilarity_matrix` is a distance matrix created from the `cluster_corr_reactions` function.
- `reactions` is a list with the reactions names.
- `plot_labels` is a variable that defines, whether the reaction names will appear in the x-axis.
- `t` defines a threshold that cuts the dendrogram at a specific height and colors the occuring clusters accordingly.
- `linkage` defines the type of linkage. Available linkage types are: `single`, `average`, `complete`, `ward`.

All the dendrograms and graphs in this page are created from correlation matrices of the E. coli core model.

This is a dendrogram created from the `plot_dendrogram` function from a correlation matrix without pearson filtering and with `correction = True`:
![dendrogram_no_cutoffs.png](/img/dendrogram_no_cutoffs.png)

This is a dendrogram from the same correlation matrix with `pearson_cutoff = 0.9999` and with `correction = True`:
![dendrogram_pearson.png](/img/dendrogram_pearson.png)

Far distinct clusters are observed. Graphs will reveal if these clusters interact with other clusters or reactions.


## Graphs

- Graph creation can reveal reactions networks.
- Subgraphs creation can reveal reactions subnetworks.

I implemented a `graph_corr_matrix` function that takes as an input a correlation matrix and creates network graphs from it.
Except from the initial graph, this function splits the graph into subgraphs.

Example of calling the `graph_corr_matrix` function:

    graphs, layouts = graph_corr_matrix(corr_matrix, reactions, correction=True, clusters=clusters)

Explaining the parameters and the returned objects:

- `corr_matrix` is a correlation matrix produced from the `correlated_reactions` function.
- `reactions` is a list with the reactions names that will appear as nodes in the graphs.
- `correction` is a boolean type variable that if set to True, it transforms the correlation matrix to the absolute correlation matrix.
- `clusters` is a nested list, containing sublists with reactions that belong to the same cluster, created from the `cluster_corr_reactions` function.
- `graphs` is a list containing created graph objects.
- `layouts` is a list containing graphs' layouts. Each layout has the same index with its corresponding graph from `graphs` list.

I also implemented a `plot_graph` function, that takes as an input a graph object with its corresponding layout and plots the resulting graph.

Example of calling the `plot_graph` function:

    plot_graph(graph, layout)

Explaining the parameters:

- `graph` is a graph object returned from the `graph_corr_matrix` function.
- `layout` is a layout that corresponds to the given graph object, also created from the `graph_corr_matrix` function.

Example of calling the `plot_graph` function recursively for every subgraph returned:

    for i in range(len(graphs)):
        graph = graphs[i]
        layout = layouts[i]
        plot_graph(graph, layout)

This is a graph created from the `plot_graph` function from a correlation matrix without pearson filtering and with `correction = True`:
![graph_no_cutoffs.png](/img/graph_no_cutoffs.png)

This is a graph created from a correlation matrix with `pearson_cutoff = 0.9999` and with `correction = True`:
![graph_pearson.png](/img/graph_pearson.png)

This is one of the subgraphs created from a correlation matrix with `pearson_cutoff = 0.9999` and with `correction = True`:
![subgraph_pearson.png](/img/subgraph_pearson.png)

This subgraph has 9 nodes that correspond to 9 reactions close to each other in the topology of the E. coli core model.
These reactions are: `PGI, G6PDH2R, PGL, GND, RPE, RPI, TKT1, TALA, TKT2`.
Their topology can be seen in the following figure obtained from `ESCHER`:
![escher_subgraph.png](/img/escher_subgraph.png)

We can see that `PGI` seems to contribute to a different pathway. However it shares a common metabolite with `G6PDH2r`.
If we apply a stricter pearson cutoff (equal to 0.99999), this reaction is removed from this subgraph and only the rest of the reactions remain.
This is an important observation. Looser cutoffs lead to wider sets of connected reactions that form larger metabolic pathways.


## Conclusion

- `dingo` is a python package for the analysis of metabolic networks. I expanded `dingo` by incorporating pre- and post- sampling features.
- Flux sampling is an unbiased method that can advance research in metabolic models both at the strain and the community level.
- From model reduction to pathways identification, the developed features incorporated to `dingo`, show the increased statistical value of flux sampling.


## References

- <a id="ref-1"></a> [1]: Created with BioRender.com.
- <a id="ref-2"></a> [2]: Zachary A. King, Justin Lu, Andreas Dräger, Philip Miller, Stephen Federowicz, Joshua A. Lerman, Ali Ebrahim, Bernhard O. Palsson, Nathan E. Lewis, BiGG Models: A platform for integrating, standardizing and sharing genome-scale models, Nucleic Acids Research, Volume 44, Issue D1, 4 January 2016, Pages D515–D522, https://doi.org/10.1093/nar/gkv1049.
- <a id="ref-3"></a> [3]: Chalkis A, Fisikopoulos V, Tsigaridas E, Zafeiropoulos H. dingo: a Python package for metabolic flux sampling. Bioinform Adv. 2024 Mar 22;4(1):vbae037. doi: 10.1093/bioadv/vbae037. PMID: 38586119; PMCID: PMC10997433.
- <a id="ref-4"></a> [4]: Zachary A. King, Andreas Dräger, Ali Ebrahim, Nikolaus Sonnenschein, Nathan E. Lewis, and Bernhard O. Palsson (2015) Escher: A web application for building, sharing, and embedding data-rich visualizations of biological pathways, PLOS Computational Biology 11(8): e1004321. doi:10.1371/journal.pcbi.1004321.
