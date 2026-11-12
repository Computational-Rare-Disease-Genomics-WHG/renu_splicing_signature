# Per-sample subset refit weights

## Intent

Use the supplied insT and Stem III refit tables to score each sample with the
model fitted to its available junctions. Bundle the parameters so installed
`renu-score` needs neither the original Downloads files nor network access.

## Approach and scope

Bundle compact CSVs containing the intercept and gene weight columns only.
A non-missing weight identifies an available event, including coefficients
of zero. Preserve every supplied numeric value. Keeping the original CSVs
would also work but ships unused fitting and evaluation diagnostics. Encoding
subsets as numeric masks would save more space but make inspection harder.

Determine availability after the existing missing-count, zero-coverage,
minimum-depth, and matched-control-baseline checks. Select the exact subset
independently for every case and control. Use its coefficients and intercept
for weighted scores and probability, including for the complete set. Keep
uniform scores, missingness statuses, centering, and sample order unchanged.
No available events still produces unavailable scores rather than an
intercept-only prediction. Report the selected coefficients in event details.

The user approved using the new full-set rows too: this updates the complete
insT model, including its intercept to -8.33698995873104.

Bundled TOMLs explicitly opt into their bundled subset tables. Existing custom
TOMLs without that setting retain their supplied fixed weights and intercept.
No model fitting, classification thresholds, or new runtime dependencies.

## Verification

Check all 2,047 insT and 127 Stem III subsets and numeric values against the
source exports. Test independent numerical scores, samples with distinct
subsets, zero coefficients, each availability exclusion, absent events,
empty samples, details, custom-model compatibility, and installed wheel/source
distribution resources. Run the full existing suite and formatting checks.
