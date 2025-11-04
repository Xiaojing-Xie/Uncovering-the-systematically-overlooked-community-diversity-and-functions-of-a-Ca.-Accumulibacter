We implemented the dynamic EBPR framework in three small modules: a configuration/utility layer, a two-phase simulation core, and a batch runner.
The configuration layer holds all model-specific IDs and helper functions.
Here we define the metabolite and reaction IDs for O₂, inorganic phosphate, polyphosphate, the PolyP hydrolysis and synthesis reactions, the biomass reaction, and the exchange reactions for amino acids and carbon sources.
This layer also provides helpers to find exchange reactions from metabolite IDs, reset bounds, run pFBA, identify the biomass reaction, open amino acid uptake according to a simple rule, and compute the carbon number of a substrate from its formula.
Users only need to adapt this configuration to their own SBML models; the rest of the code does not depend on model-specific hard-coding.

The second layer is the two-phase simulation core.
It takes a model, a chosen carbon source exchange, the O₂ exchange, and a biomass reaction, together with a small set of timing and dosing parameters (anaerobic and aerobic durations, time step, initial finite dose of carbon source).
Internally, it adds simple “storage” reactions to allow selected internal metabolites to accumulate and be released.
For each time step of the anaerobic and aerobic windows, it solves a steady-state LP (with a phase-specific objective) and updates external and internal pools by Euler integration.
At the end of one cycle, it returns the time course (time, phase label, biomass, remaining substrate, etc.), total biomass increase, total substrate used, and biomass yields per mmol and per C-mmol of anaerobically consumed substrate.

The third layer is a batch driver that loops over all models and all candidate carbon sources defined in the configuration.
For each combination, it runs a “background” cycle with no extra carbon source (only amino acids and base medium).
It then runs a full cycle with the selected carbon source.
After that, it subtracts the background growth to obtain a net biomass increase.
It computes net yields (per mmol and per C-mmol), writes time-course CSV files for each run, and collects all summary statistics into long and wide tables for further plotting or analysis.

In our work, we used the following IDs for soluble carbon sources in the configuration layer: 
Glucose (cpd00027), Acetate (cpd00029), Pyruvate (cpd00020), Lactate (cpd00159), Succinate (cpd00036), Formate (cpd00047), Fructose (cpd00082), Glycerol (cpd00100), Citrate (cpd00137), Ethanol (cpd00363), Propionate (cpd00141), Malate (cpd00130), Methanol (cpd00116), Formaldehyde (cpd00055), Cellulose (cpd11746), Butyrate (cpd00211), Benzoate (cpd00153), Phenol (cpd00127), and Toluene (cpd01034).
For the 20 proteinogenic amino acids we used: Ala (cpd00035), Val (cpd00156), Leu (cpd00107), Ile (cpd00322), Pro (cpd00129), Cys (cpd00084), Met (cpd00060), Gly (cpd00033), Trp (cpd00065), Phe (cpd00066), Lys (cpd00039), Arg (cpd00051), His (cpd00119), Tyr (cpd00069), Thr (cpd00161), Glu (cpd00023), Gln (cpd00053), Asp (cpd00041), Asn (cpd00132), and Ser (cpd00054).
Users working with other models can replace these IDs in the configuration layer and keep the rest of the workflow unchanged.
