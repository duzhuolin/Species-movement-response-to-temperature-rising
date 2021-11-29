# Species_movement_response_to_temperature_rising
This model has two entities, species and landscape. Every species has respective thermal tolerance determined by initial temperature of the patch where it resides at the beginning. Initial temperature of each patch is caculated by its elevation and latitude depending on a selected point from the real world map. Temperature changes progress at varied anual rates. When the temperature of the patch rises to  the  tolerance  temperature  of  the  species  that  resides  on  this  patch,  the  species  will move to an adjacent patch possessing a satisfactory temperature. Step sizes are adjustable. The basic principle of this model is climate-driven changes in the geographical distribution of species, simulating the emergence of poleward and uphill shifts of species. The species are able to sense the change of temperature and make adaptive behavior, e.g., emigration. Likelihood of species moving to coolest adjecent patch is set to mount the stochasticity. The simulation will run for 300 years.
