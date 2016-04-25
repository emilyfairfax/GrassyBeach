# GrassyBeach

#Introduction to the Model
<p><b>Navigating the Files</b><p>
<p>There are three main versions of the model: Beneficial_Relationship_Grassy_Beach.m, Neutral_Relationship_Grassy_Beach.m, and Harmful_Relationship_Grassy_Beach.m. There is an additional version that compiles the three relationships to run at once for easy comparison called Relationship_Comparison_Grassy_Beach.m. Also uploaded into this repository are movies and images showing the final results of the models.</p>

<p><b> Motivation </b> </p>
<p>The motivation for this model is to demonstrate how a plant's physiological response to burial by sand affects both it's growth pattern through time as well as the morphology of obstacle dunes which form downwind of the plant. It takes into account sand falling out of suspension behind the plant in the plant's "windshadow," as well as hillslope diffusion of the sand as it falls out of suspension and is added to the bed. </p>

<p><b>Relationship to Burial in Sand </b></p>
<p>In this particular version of the model, the three relationships with burial (beneficial, neutral, and harmful) are compared visually against one another. </p>
<p><i>Beneficial Relationship</i></p>
<p> In this particular version of the model, the plant has a beneficialrelationship with burial, such as that found in Sea Oats and American Dunegrass. This means that the growth rate of the plant is based on the following equation: G = c*gB(1-B), where G is the growth rate in vertical height/time, g is a scaling factor of a base growth rate, and B is the scaled biomass of the plant - in this case B is represented by the plant's height as a percentage of a max possible growth height. The "growth factor" term, c, accounts for how much of the plant is buried in the sand. There is an optimal burial percentage defined in the model - for burials less than or equal to this percentage, the amount of benefit the plant receives from burial increases with increasing burial. Above this percentage and increasing burial begins to decrease the benefit the plant receives from being buried. In this model the plant is never buried enough to cause negative growth.</p>
<p><i>Neutral Relationship</i></p>
<p>In this particular version of the model, the plant has a neutral relationship with burial. This means that the growth rate of the plant is based only on the following general equation: G = gB(1-B), where G is the growth rate in vertical height/time, g is a scaling factor of a base growth rate, and B is the scaled biomass of the plant - in this case B is represented by the plant's height as a percentage of a max possible growth height. There is no "growth factor" term to account for how much of the plant is buried in the sand since in this case it is unrelated to growth rate.</p>
<p><i>Harmful Relationship</i></p>
<p>In this particular version of the model, the plant has a harmful relationship with burial, such as that in any plant that has not evolved to thrive specifically in dunes. This means that the growth rate of the plant is based on the following equation: G = c*gB(1-B), where G is the growth rate in vertical height/time, g is a scaling factor of a base growth rate, and B is the scaled biomass of the plant - in this case B is represented by the plant's height as a percentage of a max possible growth height. The "growth factor" term, c, accounts for how much of the plant is buried in the sand. The more the plant is buried, the more it hurts the growth rate. There is a plant death burial percentage, at which point the plant growth rate goes to zero.  Once the growth rate stalls, the burial amount continues to increase due to diffusion of the accumulated dune. This starts irreversible the death of the plant, where the growth rate becomes negative until the whole plant is gone.</p>
<p><b>Model Output </b></p>
<p>This model will output a movie of the dune forming as it runs, and then upon completion produces figures showing the following values through time: total accumulated sand in the dune, plant height, percentage of plant buried in the sand, growth rate and aboveground exposed plant biomass. It also produces a plot of the growth factor term as it related to the buried biomass percentage. </p>
<p><b>Notes on Running the Model</b></p>
<p> Since this model includes diffusion, many of the parameters are fairly sensitive - in particular the time step and spacial steps, dx and dt. Changing these will likely result in an unstable model.</p>
<p><b>Enjoy</b></p>
<p>Feel free to email me with any questions at emily.fairfax@gmail.com</p>
