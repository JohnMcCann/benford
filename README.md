# Benford's law of anomalous numbers

Benford's law can be used to analysis the numbers in a dataset. Caveats exist arond which datasets are likely to object Benford's law. Datasets that span multiple orders of magntiude, and whose values are the product of numerous distributions (random processes) will likely conform to Benford's law.

The project is unique in that a user can explore the distribution in any integer radix (base), with any starting digit location counting back from the leading digit, and for any length of consecutive digits. The analysis also comes complete with analysis if one can reject the null hypothesis.

See the `Benford_overview.ipynb` jupyter-notebook for a walk-through of the functionality of the project.

## visualizations

Distribution of leading digit of 50k draws from a log-normal distribution:

<p align="center">
<img src="https://github.com/JohnMcCann/benford/wiki/images/Benford_lognorm_50k.png" alt="Benford distribution of 50k log-normal draws" width="42.1%"/><img src="https://github.com/JohnMcCann/benford/wiki/images/Confidence_lognorm_50k.png" alt="dBenford confidence bounds for 50k log-normal distributionawing" width="49%"/>
</p>

Distribution of 2nd and 3rd digit in radix-5 of 1M draws from a (log-normal)x(log-normal) distribution:

<p align="center">
<img src="https://github.com/JohnMcCann/benford/wiki/images/Benford_radix5_digits2thru3.png" alt="Benford distribution of 1M 1Mx1M log-normal draws" width="42.1%"/><img src="https://github.com/JohnMcCann/benford/wiki/images/Confidence_radix5_digits2thru3.png" alt="Benford confidence bounds for 1Mx1M log-normal distribution" width="49%"/>
</p>
