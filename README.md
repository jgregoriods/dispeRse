<h1>dispeRse</h1>
<h2>Simulation of demic diffusion with environmental constraints</h2>

**Jonas Gregorio de Souza**<br/>
jonas.gregorio@gmail.com<br/>
[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0001--7879--4531-green.svg)](https://orcid.org/0000-0001-6032-4443)<br/>

<h3>Motivation</h3>



<h3>Submodels</h3>

<h3>Growth</h3>

Growth is modelled as a density-dependent process by setting an upper boundary or saturation point, defined as the carrying capacity $K$. As the population approaches $K$, the growth rate $r$ decreases. It is assumed that population at the saturation point will remain in equilibrium as the birth and death rates cancel each other out. There is more than one proposal on how to model logistic growth; here, we use the Verhulst-Pearl equation for its simplicity:

$$
\frac{d{N}}{d{t}} = r N \left( 1 - \frac{N}{K} \right)
$$

Notice that, for economy, we write $r$, $N$ and $K$ for the growth rate, population and carrying capacity of a cell. We could also write $r_{i,j}$, $N_{i,j}$ etc. for each cell indexed by $(i,j)$. The solution to the equation above gives the new population after growth for $t$ years:

$$
N(t) = \frac{KN_0}{(K-N_0) e^{-rt} + N_0}
$$

Where $N_0$ is the initial population.

<h3>Emigration</h3>

Models of emigration as a density-dependent process usually assume that pressure to emigrate increases with population density due to diminishing returns. There are many competing models, most of which postulate a population threshold above which emigration happens. We will use one of the simplest, the asymptotic threshold model. Past a given threshold $\phi$, the probability of emigration grows asymptotically towards 100% as the population approaches carrying capacity.

By rewriting the equation of ..., we can model the number of migrants at each time step as:

$$
M =
    \begin{cases}
        0 & \text{if } \phi \leq \frac{N}{K} \\
        N - \phi K & \text{otherwise}
    \end{cases}
$$

Where $N$ is the population after growth.

<h3>Mobility</h3>

In normal terrain, the migrants are redistributed to the eight cells in the Moore neighborhood, provided they have positive environmental values and their population is below the emigration threshold. If a terrain layer is included with values for corridors and barriers, dispersal is affected in the following way: 1) population in cells marked as barriers will not disperse; 2) population in cells marked as corridors are allowed to disperse beyond the Moore neighborhood. In the latter case, dispersal can be as far as the distance defined by the parameter acceleration, and only if the more distant cells are also in a corridor. Redistribution of the population to the available cells occurs proportionally to the inverse of the square distance.

<p align="center"><img src="man/img/disp.png" width="300" /></p>
<h5><p align="center"><b>Figure 1. a.</b> Cells considered for migration in normal terrain. <b>b.</b> Cells considered for migration in a corridor (light grey cells) with acceleration=3. Values are the inverse of the square distance to each cell.</p></h5>

<h3>Effect of the environment</h3>

While the terrain layer affects mobility, the environment layer controls carrying capacity and growth rate. The environment layer represents the maximum population density that can be achieved in each cell as a fraction of the maximum absolute density in the study area, expressed in the interval [0,1]. Negative values can be used to represent oceans, glaciers or other areas that cannot be settled.

Each cell's carrying capacity is taken directly from the environment layer. The dependence of $r$ on the environment, however, is allowed to be nonlinear, with a shape controlled by the parameter $\gamma$. To be explicit, for each cell indexed by $(i,j)$, the local growth rate $r_{i,j}=K_{i,j}^\gamma$. Thus, when $\gamma=0$, the growth rate is the same for all cells; when $\gamma=1$, it varies linearly with the environment and so on.

As a useful procedure for simulating the effect of climate change, the environment and terrain layers can be updated at time steps defined by the user.

<h3>Example</h3>

We can test the model's performance on a well-known example, the Neolithic migration from the Near East to Europe.

```R
library(dispeRse)
library(raster)
library(rnaturalearth)
library(viridisLite)
```

```R
# for plotting
borders <- ne_download(scale=50, type="coastline", category="physical")
```

The package <code>dispeRse</code> includes the following datasets:
* euro_npp: a RasterStack with net primary production (NPP) in Europe, North Africa and Near East between 11 ka and 4 ka in 1 ka intervals. Scaled between 0 and 1. NPP has been calculated using the Miami formula on paleoclimate simulations from ... downloaded with the package <code>pastclim</code>.
* terr_npp: a Raster with elevation > 1750 m marked as barriers and major rivers and coastline marked as corridors. Elevation was reclassified from SRTM data provided by bioclim. For the rivers, we rasterized those classified as levels 1 and 2 in ... Coastline was taken from natural earth.
* ppnb: coordinates and earliest dates (median cal BP) for Near Eastern sites of the Late Pre-Pottery Neolithic B period (taken from Pinhasi).
* euro_dates: coordinates and earliest dates (median cal BP) for European Neolithic sites (taken from Pinhasi).

Since the environment will be updated at the time steps defined by the NPP stack, we need to create an equivalent stack for the terrain:

```R
terr <- raster::stack(replicate(8, euro_terr))
```

The following parameters will be used for the simulation: an emigration threshold of ca. 1/3 of $K$; an annual growth rate of 2.5%; and a generation time of 30 years. The default migration distance of 50 km per generation will be used, with an acceleration of 3 (150 km) for corridors (rivers and coastline). We will run the simulation for 5500 years, which is enough for the entire Europe to be settled. These parameters are similar, though not identical to, the ones used by Fort et al.

```R
sim <- simulate_dispersal(euro_npp, terr, ppnb, 5500, phi=0.33, r=0.025, t=30, updates=seq(10000,4000,-1000))
[1] "Preparing rasters..."
[1] "Running model..."
[1] "Done."
```

Let's visualize the simulated arrival times:

```R
plot(sim, col=viridis(10), xlim=c(-12, 46), ylim=c(28, 66))
contour(sim, nlevels=5, add=TRUE)
plot(borders, add=TRUE)
```

<p align="center"><img src="man/img/test.png" width="400" /></p>

We can evaluate the model on the European Neolithic dates by taking the mean absolute error:

```R
sim_dates <- extract(sim, euro_dates)
print(mean(abs(euro_dates$date - sim_dates), na.rm=TRUE))
[1] 528.0772
```

In conclusion, our model has a similar error to the one reported by Fort et al. This is a reasonable fit, but we can probably achieve a smaller error by fine-tuning the parameters.