## Estimating test-retest reliability for a Stroop dataset

Gang Chen (twitter: @gangchen6)

### Preface ###

Properly estimating `test-retest reliability` has become a hot topic in recent years. Traditionally test-retest reliability is usually conceptualized and quantitatively estimated as `intraclass correlation` (ICC), whose definition can be found in, for example, this [wikipedia](https://en.wikipedia.org/wiki/Intraclass_correlation) page. The trational ICC formulation works fine when there is no or little measurement error. However, the adoption of ICC can be problematic when measurement accuracy becomes an issue in situations where the quantity under study is measured many times with substantial amount of variability. 

A new modeling framework is needed to handle measurement errors. In experiments where the effect is assessed through many trials, we need to construct a hierarchical/multilevel model that 1) chacterizes/maps the data structure as close to the data structure (or data generative mechanism) as possible, and 2) separates the cross-trial variability from the estimation process of test-retest reliability. It is this hierarchical modeling framework we would like to adopt and demonstrate in this blog. More detailed theorectical discussion from the modeling perspective can be found in the following papers.

* Haines et al., 2020. Learning from the Reliability Paradox: How Theoretically Informed Generative Models Can Advance the Social, Behavioral, and Brain Sciences (preprint). PsyArXiv. https://doi.org/10.31234/osf.io/xr7y3

* Rouder, J.N., Haaf, J.M., 2019. A psychometrics of individual differences in experimental tasks. Psychon Bull Rev 26, 452–467. https://doi.org/10.3758/s13423-018-1558-y

* Chen et al., 2021. Trial and error: A hierarchical modeling approach to test-retest reliability. NeuroImage 245, 118647. https://doi.org/10.1016/j.neuroimage.2021.118647

This blog intends to

1) lay out the hierarchical modeling structure;
2) demonstrate the implementation of the hierarhical model using the `R` package `brms` with a dataset from a Stroop experiment from [Hedge et al., 2018](https://doi.org/10.3758/s13428-017-0935-1).


### Hierarchical modeling framework ###

The modeling framework can be laid as below. Suppose that, in a test-retest experiment, the effect of interest (e.g., reaction time) $y_{crst}$ is measured at trial $t$ ($t=1,2,...,T$) during each of the two repetitions/sessions ($r=1,2$) for subject $s$ ($s=1,2,...,S$) under the condition $c$ ($c=1,2$). If one adopts the conventional ICC formulation, the data would have to be aggregated by collapsing trial dimension and obtain, for example, the average values $\overline{y}_{cs\cdot}$. However, test-retest reliability could be underestimated under some circumstances, and extent of underestimation depends on the relative magnitude of cross-trial variability compared to its cross-subject counterpart (Rouder et al., 2019; Chen et al., 2021). Here we build the following hierarchical model:

trial level: $y_{crst}~\sim ~\mathcal D (\mu_{crs},   \sigma^2);$\
subject level: $(\mu_{11s}, \mu_{21s}, \mu_{12s}, \mu_{22s})^T \sim ~\mathcal N(\boldsymbol 0_{4\times 1}, ~\boldsymbol S_{4\times 4}).$

Here the distribution $\mathcal D$ at the trial level can be any probability density that could properly capture the data generataive mechanism. The typical distributions for reaction time are Gaussian, Student's $t$, exponentially-modified Gaussian, (shifted) log-normal, etc. The variance-covariance matrix $\boldsymbol S_{4\times 4}$ captures the inter-relationships among the effects at the subject level. We know that, after scaling, the variance-covariance matrix $\boldsymbol S_{4\times 4}$ would show the correlation structure among the four components of $(\mu_{11s}, \mu_{21s}, \mu_{12s}, \mu_{22s})$. Later I will demonstrate how to extract the jewels in the crown from this matrix $\boldsymbol S_{4\times 4}$. 


I'd like to note the following few aspects-

The crucial aspect of the hierachical model above is that 
separation 




A few notes

$\sigma$

(I wish that the model could be expressed more elegantly using vector-matrix formulation, but the math notation support at gihub is quite limited at the moment.)


$$(\sigma_{11s}, \sigma_{21s}, \sigma_{12s}, \sigma_{22s})^T \sim ~\mathcal N(\boldsymbol 0_{4\times 1}, \boldsymbol S_{4\times 4})$$

$$(y_{1st}, ~y_{2st})^T ~\sim ~\mathcal D ((\mu_{c_1 r_1s},  \mu_{2s})^T, \Sigma_{2\times 2})$$
$$(\mu_{1s}, ~y_{2st})^T ~\sim ~\mathcal D ( $$

---

### This is a header

#### Some T-SQL Code

```tsql
SELECT This, [Is], A, Code, Block -- Using SSMS style syntax highlighting
    , REVERSE('abc')
FROM dbo.SomeTable s
    CROSS JOIN dbo.OtherTable o;
```

#### Some PowerShell Code

```powershell
Write-Host "This is a powershell Code block";

# There are many other languages you can use, but the style has to be loaded first

ForEach ($thing in $things) {
    Write-Output "It highlights it using the GitHub style"
}
```
