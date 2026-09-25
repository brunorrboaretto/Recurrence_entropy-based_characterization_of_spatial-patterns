# Recurrence_entropy-based_characterization_of_spatial-patterns
This repository contains the code and supporting material for the manuscript “Recurrence entropy-based characterization of spatial patterns,” which is currently under peer review.

## Paper overview
We introduce a recurrence-based framework for the characterization of two-dimensional spatial patterns by extending concepts from nonlinear time-series analysis to spatial data. The method operates directly on pixel arrays, generating spatial recurrence data that are decomposed into recurrence motifs. The statistical distribution of these motifs is then used to quantify the complexity and organization of spatial patterns through recurrence entropy.

By systematically varying the recurrence threshold, the method provides a compact representation of spatial organization in terms of the threshold and the corresponding recurrence entropy. This framework allows different types of spatial patterns and textures to be characterized and compared using a common recurrence-based approach.

## Repository Contents
The repository provides a self-contained example illustrating the main steps of the proposed framework:

- **Dataset generation:** code to generate fractional Brownian surfaces with different Hurst exponents;
- **Recurrence analysis:**  implementation of the recurrence microstate decomposition and recurrence entropy calculation;
- **Data output:** numerical results containing the Hurst exponent, random seed, system size, recurrence threshold, recurrence entropy, and mean recurrence rate;
- **Visualization:** generation of figures showing the spatial pattern together with the recurrence entropy and mean recurrence rate as functions of the recurrence threshold.
The example is designed to facilitate the reproduction of the main recurrence-based analysis and to provide a starting point for applying the method to other two-dimensional spatial patterns and textures.

## Python libraries:

- **NumPy**: Facilitates efficient handling and manipulation of large multi-dimensional arrays and provides a wide range of mathematical functions for numerical computations in Python;
- **matplotlib.pyplot** is a Python library commonly used for creating visualizations and plots, providing a high-level interface for generating a wide range of graphs and charts;

## Usage

Run the main script: 

<code> python example.py </code> 

The script generates fractional Brownian surfaces for a range of Hurst exponents and evaluates the recurrence entropy and mean recurrence rate over a set of recurrence thresholds. The numerical results are saved as a <code>.dat</code> file, while the corresponding figures are saved in the <code>results/</code> directory.

## Output

The numerical output contains the following columns:

<code> H    seed    N    threshold    entropy    recurrence_rate </code>

where:

- <code> H </code> is the Hurst exponent;
- <code> seed </code> is the random seed used to generate the surface;
- <code> N </code> is the linear size of the spatial pattern;
- <code> threshold </code> is the recurrence threshold;
- <code> entropy </code> is the normalized recurrence entropy;
- <code> recurrence_rate </code> is the mean recurrence rate.

### Citation

If you find this work helpful for your research, please consider reading and citing:

- Boaretto, B. R. R., Di Domenico, M., Prado, T. L., Lopes, S. R., Macau, E. E. N., & Masoller, C. (2026). "Recurrence entropy-based characterization of spatial patterns." *(under review).*  

A complete citation will be added once the manuscript is published.

--------------------------------------------------------------------------------------

Thank you for your interest in our research! </br>
We hope this repository will be a useful resource for researchers working on the characterization of spatial patterns, recurrence analysis, and image analysis in general. .</br>

Sincerely,</br>
Bruno R. R. Boaretto.

