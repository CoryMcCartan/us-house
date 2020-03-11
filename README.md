# U.S. House Election Prediction
This is a statistical model for the biannual elections of the United States
House of Representatives.  We use a two-stage Bayesian model to forecast voter
intent and election results, using past election data and generic congressional
ballot polling.  The model is implemented in [Stan](http://mc-stan.org) and R.

## Files

- `forecast.R` downloads new polling data and reruns the analysis.
- `model/` contains the model files, code, and data.
- `site/` contains the website that displays the analysis and model results.
- `docs/` contains a built copy of the website for serving.
