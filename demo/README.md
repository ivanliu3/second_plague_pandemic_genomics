## 3. Demo

The `demo/` directory contains example PCA data for ancient test individuals and present-day reference populations. The demo calculates PCA-based Mahalanobis distances using the first seven principal components.

### Instructions to run

From the top-level repository directory, run:

```bash
bash demo/run_demo.sh
```

Alternatively, from the `demo/` directory, run:

```bash
./run_demo.sh
```

The demo requires Python 3, NumPy, and SciPy.

### Expected output

The demo generates:

* `data.mahadist.csv`: Mahalanobis distances between ancient test individuals and reference populations
* `data.pvalue.csv`: p-values associated with the Mahalanobis-distance results

### Expected runtime

The demo completed in 3.85  seconds on a standard macOS desktop using Python 3.12. It is expected to complete in less than 1 minute on a normal desktop computer and requires no specialised hardware.
