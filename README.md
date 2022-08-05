<h1>dispeRse</h1>
<h2>Simulation of demic diffusion with environmental constraints</h2>

**Jonas Gregorio de Souza**<br/>
jonas.gregorio@gmail.com<br/>
[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0001--7879--4531-green.svg)](https://orcid.org/0000-0001-6032-4443)<br/>

<p>The growth rate is given by:</p>

$$
\frac{\partial{N}}{\partial{t}} = N r \left( 1 - \frac{N}{K} \right)
$$

<p>This is a test. Example formula:</p>

$$
N_t = \frac{KN_0}{(K-N_0) e^{-rt} + N_0}
$$

<p>The probability of fission:</p>

$$
m =
    \begin{cases}
        0 & \text{if } \phi \leq \frac{N}{K} \\
        N - \phi K & \text{otherwise}
    \end{cases}
$$

<img src="man/img/a.png" />