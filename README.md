# LEVMpy

Python module and wrapper for the LEVM complex nonlinear least squares program.

[LEVM](http://jrossmacdonald.com/levmlevmw/) was written by J. Ross Macdonald.

Example of use:

```python
import levmpy as lv
expt = lv.Experiment("name-of-input-file", outsuffix="OUTFILE", path="path-name")
```

To fit the data use the command:

```python
expt.fit()
```

To write the output to a file use:

```python
expt.writetofile()
```


The following class variables of the Experiment class may be useful:

| Variable            | Description                                               |
|---------------------|-----------------------------------------------------------|
| `expt.freq`         | frequency data                                            |
| `expt.y`            | immittance data (can use `y1` and `y2` for complex data)  |
| `expt.r`            | uncertainty on y (can use `r1` and `r2` for complex data) |
| `expt.parameters`   | initial parameters                                        |
| `expt.celcap`       | empty cell capacitance                                    |
| `expt.x`            | fitted parameter estimates                                |
| `expt.rxsd`         | relative standard deviation on x                          |
| `expt.nfrei`        | number of free parameters                                 |
| `expt.jac`          | Jacobian matrix                                           |
| `expt.outputvals`   | model immittance data                                     |
| `expt.res`          | fit residual                                              |
| `expt.resmod`       | fit residual / model                                      |
| `expt.md`           | number of data points                                     |
| `expt.info`         | info flag from LMDER.f                                    |
| `expt.fnorm`        | L2 norm of residuals                                      |
| `expt.nfev`         | number of function evaluations                            |
| `expt.fqq`          | fit quality factor                                        |
| `expt.ndf`          | number of degrees of freedom                              |
| `expt.dattyp`       | data type (C, R or I)                                     |
| `expt.freqtyp`      | frequency type (F or None)                                |
