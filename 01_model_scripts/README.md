### Models

##### Purpose

The purpose of the model(s) is to simulate a fishery overtime that is using FISHE to guide their management decisions. The model allows us to test the effectiveness of different management actions on the health of a fishery over time. 

##### Differences

There is one difference between the models. The *sim_fishery* model does not close if the fishery goes over the limit. Assuming a full closure is unrealistic, this model maintains some small level of fishing activity. If the fishery is over the limit in year *y*, the model will cut back the catch in year *y* by 95%. This cutback lasts a single year and the fishing pressure (*f*) in year *y+1* will be the same as year *y*. 

### Model Assumptions

The models are based on the following assumptions:

- Fisheries represented in the models have no current management, therefore initial biomass is due to some equilibrium that the fishery has reached over the years based on a fishing mortality rate (*f*).
- The current fishing mortality rate (*f*) is reached when catch is equivalent to the surplus production of the fishery at that biomass. 
- The fishing mortality rate (*f*) calculated at the beginning of the simulation is the true fishing mortality rate of the operating fishery.
- Fishing mortality (*f*) and *f/fmsy* are calculated at the beginning of each year and the change in *f* made by the management decision is applied to *f* in that same year.
- Climate change is affecting the productivity of fish stocks by influencing the growth term (*r*).
- Participants of the fishery adhere to the management decisions and the fishing mortality rate (*f*) that results from the management decision is the true fishing mortality rate of the fishery. 

### Performance Indicators/Reference Points

The model uses the following performance indicator: *f / fmsy* also called the *f ratio*

which is the ratio of the current fishing mortality rate (*f*) to the fishing mortality that would produce biomass of maximum sustainable yield (*fmsy*).

The reference points used in the model are a target ratio of 1 and a limit ratio of 2. A target ratio of 1 means the fishery is currently fishing at exactly the fishing pressure that gives the biomass of maximum sustainable yield, while a ratio of 2 means the fishery is fishing twice as hard as the fishing pressure that would give the biomass of maximum sustainable yield. The decision tree used in the model is as follows: 

| Condition                   | Action                                                       |
| --------------------------- | ------------------------------------------------------------ |
| *f / fmsy* > limit          | Close the fishery (sim_closure) or maintain *f* but reduce the catch (sim_fishery) |
| target < *f / fmsy* < limit | Reduce fishing pressure (*f*) by the harvest control rule    |
| target = *f / fmsy*         | Fishing pressure *f* remains the same                        |
| target > *f / fmsy*         | Increase fishing pressure *f* by 5%                          |

### Model Inputs

The two models have the same model inputs. Model inputs are also described in the `.R` files. 

-  `b` - initial biomass of the fishery
- `r` - intrinsic growth rate of species 
- `r_s` - pace and magnitude of climate change as a percent decline in productivity per year (refer to the [climate scenarios](./reference/climate_scenarios.md) document for examples)
- `r_p_s` - a proxy estimate for the pace and magnitude of climate change as a percent decline in productivity per year (values can also be chosen using the [climate scenarios](./reference/climate_scenarios.md) document as a reference) 
- `error` - the amount of sampling error as a decimal (i.e. 0.1 for 10% error)
- `p` - the shape parameter of the surplus production model; default is 0.2 for the Pella Tomlinson model 
- `k` - the estimated carrying capacity relative to the initial biomass; default is 10000
- `years` - the length of the simulation; default is 100 years
- `hcr` - the harvest control rule used if the fishery falls between the target and limit expressed as the decimal 1-hcr (i.e. if you will reduce fishing pressure by 10% the hcr is 0.9)

### Model Outputs

The two models have the same model outputs. Model outputs are also described in the `.R` files. 

- `b` - biomass in year (*y*)
- `c` - catch in year (*y*)
- `year` - year of the simulation
- `r` - growth rate in year (*y*)
- `r_p` - the perceived growth rate in year (*y*) based on a proxy estimate of climate change influences from the input `r_p_s`
- `f` - the fishing mortality rate in year (*y*) 
- `f_msy` - the fishing mortality rate that gives biomass for maximum sustainable yield; calculated using the current growth rate `r`
- `f_msy_p` - the perceived fishing mortality rate that gives biomass for maximum sustainable yield; calculated based on if and how managers are accuounting for the influence of climate change on the underlying growth of the species (`r_p`)
- `f_ratio` - the ratio of fishing mortality (`f`) to fishing mortality that would give the biomass of maximum sustainable yield; based on actual growth rates (`r`) in year (*y*)
- `f_ratio_p` - the ratio of fishing mortality (`f`) to fishing mortality that would give the biomass of maximum sustainable yield; based on the perceived growth rate (`r_p`) in year (*y*)
- `f_ratio_err` - an f / fmsy ratio used to make a management decision (as shown in the table above). This ratio is drawn from a lognormal distribution with a mean equal to the `f_ratio_p` and a coefficient of variation equal to the `error` input

