README: How ecological opportunity drives asymmetric phenotypic diversity in the gartersnakes, watersnakes, and allies (Natricidae:Thamnophiini)
https://doi.org/10.5061/dryad.ghx3ffc03

Description of the data and file structure
Files and variables
File: All_Thamno_Traits.csv
Description: Dataset of trait data for all Thamnophiini taxa.

Variables
Name_in_Tree: Scientific binomial included in the Thamnophiini phylogeny.
Lifestyle: Habitat categorization for each taxon.
Diet: Diet categorization for each taxon.
Clade: Clade categorization for each taxon.
Side: Geography categorization for each taxon.
Genus: Genus for each taxon.
File: Landmarks_Bilat_Curves.csv
Description: Datset of landmarks to be read into downstream geometric morphometrics scripts.

Variables
num: Number across all total landmarks.
lm: Number split across anatomical landmarks and curve semilandmarks.
Bone: Bone on which the landmark is placed.
type: If it is either anatomical landmark or cuver semilandmark.
Name: Name of the landmark as defined from json files exported from 3D Slicer.
Description: Abbreviation of position of the landmark.
Position: Whether landmark is lateral, medial, or on the right bone.
File: Specimens_Used.csv
Description: List of specimens CT-scanned for this study.

Variables
species: Taxon scientific binomial.
catnum: Museum catalog number of specimen.
sex: Sex of the specimen.
File: All_Thamnophiini_Ana_Curves_UCE_CatNum.tps
Description: Concatenated landmarks across all bones for all specimens.

File: Thamno_Morphology_Tree.tre
Description: Time-calibrated phylogeny of Thamnophiini used for all downstream comparative analyses in this study.

File: Braincase_PACA.nex
Description: Nexus file of the first four phylogentically aligned components for braincase shape to be used for BayesTraits and MuSSCRat analyses.

File: Diet.nex
Description: Nexus file of diet categorizations to be used for MuSSCRat analyses.

File: Geography.nex
Description: Nexus file of geography categorizations to be used for MuSSCRat analyses.

File: Habitat.nex
Description: Nexus file of habitat categorizations to be used for MuSSCRat analyses.

File: Maxilla_PACA.nex
Description: Nexus file of the first four phylogentically aligned components for maxilla shape to be used for BayesTraits and MuSSCRat analyses.

File: Palatopterygoid_PACA.nex
Description: Nexus file of the first four phylogentically aligned components for palatopterygoid arch shape to be used for BayesTraits and MuSSCRat analyses.

File: Skull_PACA.nex
Description: Nexus file of the first four phylogentically aligned components for skull shape to be used for BayesTraits and MuSSCRat analyses.

File: Snout_PACA.nex
Description: Nexus file of the first four phylogentically aligned components for snout shape to be used for BayesTraits and MuSSCRat analyses.

File: Thamnophiini_UCE_51_Taxa.nex
Description: Nexus file of time-calibrated Thamnophiini phylogeny to be used for BayesTraits and MuSSCRat analyses.

File: Mandible_PACA.nex
Description: Nexus file of the first four phylogentically aligned components for mandible shape to be used for BayesTraits and MuSSCRat analyses.

File: Suspensorium_PACA.nex
Description: Nexus file of the first four phylogentically aligned components for suspensorium shape to be used for BayesTraits and MuSSCRat analyses.

Code/software
Geometric Morphometrics Scripts

Sliding_Landmarks.R: R script to prepare landmark data for downstream geometric morphometric analyses.
PCA_Script.R: R script to conduct principal components analyses on shape data.
Convergence_Analyses.R: R script to conduct convergent evolution analyses on shape data.
Supplementary Functions

read.markups.json.R: Supplementary functions to read into Sliding_Landmarks.R to read JSON files exported from 3D Slicer. Taken from Rolfe et al. (2021).
MatchedLocalSuperimpositions.R: Supplementary functions to read into PCA_Script.R to perform local Procrustes superimposition. Taken from Rhoa et al. (2021).
MuSSCRat Scripts

**Geography_Braincase_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of geography on evolution of braincase shape.
**Geography_Mandible_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of geography on evolution of mandible shape.
**Geography_Maxilla_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of geography on evolution of maxilla shape.
**Geography_Palatopterygoid_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of geography on evolution of palatopterygoid arch shape.
**Geography_Skull_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of geography on evolution of skull shape.
**Geography_Snout_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of geography on evolution of snout shape.
**Geography_Suspensorium_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of geography on evolution of suspensorium shape.
**Habitat_Braincase_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of habitat on evolution of braincase shape.
**Habitat_Mandible_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of habitat on evolution of mandible shape.
**Habitat_Maxilla_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of habitat on evolution of maxilla shape.
**Habitat_Palatopterygoid_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of habitat on evolution of palatopterygoid arch shape.
**Habitat_Skull_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of habitat on evolution of skull shape.
**Habitat_Snout_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of habitat on evolution of snout shape.
**Habitat_Suspensorium_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of habitat on evolution of suspensorium shape.
**Diet_Braincase_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of diet on evolution of braincase shape.
**Diet_Mandible_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of diet on evolution of mandible shape.
**Diet_Maxilla_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of diet on evolution of maxilla shape.
**Diet_Palatopterygoid_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of diet on evolution of palatopterygoid shape.
**Diet_Skull_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of diet on evolution of skull shape.
**Diet_Snout_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of diet on evolution of snout shape.
**Diet_Suspensorium_relaxed_state_dependent_hypo_test.Rev: **MuSSCRat script to test state-dependence of diet on evolution of suspensorium shape.
