# Pre- and post-sampling features to leverage flux sampling at both the strain and the community level

#### A contribution for the Google Summer of Code 2024 program

## Overall

#### A summary of the implemented methods, merged into the dingo library:

- preprocess for the reduction of metabolic models
- inference of pairwise correlated reactions
- visualization of a steady-states correlation matrix
- construction of a weighted graph of the model's reactions with the correlation coefficients as weights
- annotation of these weights to the metabolic model and extraction to an annotated SBML file

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

`cobra_model = load_json_model("ext_data/e_coli_core.json")`

`obj = PreProcess(cobra_model, tol=1e-6, open_exchanges=False)`

`removed_reactions, final_dingo_model = obj.reduce(extend=False)` 

Explaining the parameters and the returned objects:

- `open_exchanges` is a parameter to the `find_blocked_reactions` function of the cobra library. It controls whether or not to open all exchange reactions to very high flux ranges.
- `removed_reactions` is a list that contains the names of the removed reactions. `final_dingo_model` is a reduced model with bounds of removed reactions set to 0.
- Users can decide if they want to remove an additional set of reactions, by setting the `extend` parameter to `True`.
  These reactions are the ones that do not affect the value of the objective function, when removed.

