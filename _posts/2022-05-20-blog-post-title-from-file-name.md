## Estimating test-retest reliability for Stroop effect

Gang Chen (twitter: @gangchen6)

### Preface ###

Properly estimating `test-retest reliability` has become a hot topic in recent years. Traditionally test-retest reliability is usually conceptualized and quantitatively estimated as `intraclass correlation` (ICC), whose definition can be found in, for example, this [wikipedia](https://en.wikipedia.org/wiki/Intraclass_correlation) page. The trational ICC formulation works fine when there is no or little measurement error. However, the adoption of ICC can be problematic when measurement accuracy becomes an issue in situations where the quantity under study is measured many times with substantial amount of variability. 

A new modeling framework is needed to estimate test-retest reliabiltiy with measurement errors properly handled. In experiments where the effect is assessed through many trials, we need to construct a hierarchical/multilevel model that 

1) chacterizes/maps the data structure as close to the data structure (or data generative mechanism) as possible, and

2) separates the cross-trial variability from the estimation process of test-retest reliability. It is this hierarchical modeling framework we would like to adopt and demonstrate in this blog. More detailed theorectical discussion from the modeling perspective can be found in the following papers.

* Haines et al., 2020. Learning from the Reliability Paradox: How Theoretically Informed Generative Models Can Advance the Social, Behavioral, and Brain Sciences (preprint). PsyArXiv. https://doi.org/10.31234/osf.io/xr7y3

* Rouder, J.N., Haaf, J.M., 2019. A psychometrics of individual differences in experimental tasks. Psychon Bull Rev 26, 452–467. https://doi.org/10.3758/s13423-018-1558-y

* Chen et al., 2021. Trial and error: A hierarchical modeling approach to test-retest reliability. NeuroImage 245, 118647. https://doi.org/10.1016/j.neuroimage.2021.118647

This blog intends to

1) lay out the structure of a hierarchical modeling framework;

2) demonstrate the implementation of the hierarhical model using the `R` package `brms` with a dataset from a Stroop experiment ([Hedge et al., 2018](https://doi.org/10.3758/s13428-017-0935-1)).


### Hierarchical modeling framework ###

The modeling framework can be laid as below. Suppose that, in a test-retest experiment, the effect of interest (e.g., reaction time) $y_{crst}$ is measured at trial $t$ ($t=1,2,...,T$) during each of the two repetitions/sessions ($r=1,2$) for subject $s$ ($s=1,2,...,S$) under the condition $c$ ($c=1,2$). If one adopts the conventional ICC formulation, the data would have to be aggregated by collapsing trial dimension and obtain, for example, the average values $\overline{y}_{cs\cdot}$. However, test-retest reliability could be underestimated under some circumstances, and extent of underestimation depends on the relative magnitude of cross-trial variability compared to its cross-subject counterpart (Rouder et al., 2019; Chen et al., 2021). Here we build the following hierarchical model:

**trial** level: $y_{crst}~\sim ~\mathcal D (m_{cr}+\mu_{crs},   \sigma_{crs}^2);$\
**subject** level: $(\mu_{11s}, \mu_{21s}, \mu_{12s}, \mu_{22s})^T \sim ~\mathcal N(\boldsymbol 0_{4\times 1}, ~\boldsymbol S_{4\times 4}).$

Here the distribution $\mathcal D$ at the trial level can be any probability density that could properly capture the data generataive mechanism. The typical distributions for reaction time are Gaussian, Student's $t$, exponentially-modified Gaussian, (shifted) log-normal, etc. $m_{cr}$ is the population-level effect under condition $c$ during repetition $r$. The variance-covariance matrix $\boldsymbol S_{4\times 4}$ captures the inter-relationships among the subject-level effects $\mu_{crs}$. We know that, after scaling, the variance-covariance matrix $\boldsymbol S_{4\times 4}$ would show the correlation structure among the four components of $(\mu_{11s}, \mu_{21s}, \mu_{12s}, \mu_{22s})$. Later I will demonstrate how to extract the jewels in the crown from this matrix $\boldsymbol S_{4\times 4}$ and obtain test-retest reliability for various effects. (*I wish that the model could be expressed more elegantly using vector-matrix formulation, but the math notation support at gihub is quite limited at the moment.*)

The cross-trial variability $\sigma_{crs}$ can be further partitioned among the four combinations between the two factors of condition and session. Specifically, as shown in Haines (2020), the standard deviation $\sigma$ can be structured with three indices $c$, $r$, and $s$, and then assumed to be (mirroring the subject-level effects above):
    
$(\sigma_{11s}, \sigma_{21s}, \sigma_{12s}, \sigma_{22s})^T \sim ~\mathcal N((\gamma_{11}, \gamma_{21}, \gamma_{12}, \gamma_{22})^T, ~\boldsymbol G_{4\times 4}).$

Below we will adopt a hierarchical model with this fine-tuned structure for cross-trial variability.

Understanding the modeling formulation is important. Without jotting a model in mathematical formula, I would have trouble fully grasping a chunk of code (e.g., `brms` implementation). In fact, usually the model can be directly mapped to the numerical code. I'd like to note the following with regard to the hierarchical model for test-retest reliability -

* The crucial aspect of the hierachical model above is that the separation of cross-trial variability -- characterized by the variance $\sigma^2$ at the trial-level formulation -- from the test-retest relialiability (embedded in the variance-covariance matrix $\boldsymbol S_{4\times 4}$ at the subject level. It is the separation that allows the accurate estimation of test-retest reliability through the hierarchicaal model. It is also because of the lack of this seperation in the conventional ICC formulation that leads to the underestiation of test-retest reliability.

* The hierachical model formulated here is parameterized flatly with all the four combinations (indexed by $c$ and $r$). Personally I consider this parameterization is more intuitive and more generic. That is, the model considered here is differently from all the three papers cited above including my own (Chen et al., 2021).

    * Haines et al. (2020) adopted dummy coding for the two conditions with one condition coded as 1 while the other serves as the reference. Thus, each slope would correspond to the condition contrast (usually the effect of interest) and each intercept is associated with the reference condition per session. I might be wrong, but it seems that a strong assumption was made in Haines et al. (2020) that no correlation exists for the reference condition between the two repetitions/sessions. 
    
    * Rouder et al. (2019) use an indicator variable for the two conditions (0.5 for one condition and -0.5 for the other). Under this coding, each slope is the contrast between the two conditions per session, usually the effect of interest while each intercept is the average between the two conditions. One underlying assumption with the model in Rouder et al. (2019) was that no correlation exists between a slope (contrast) and an intercept (average). In addition, homogeneity of cross-trial variability was assumed.
    
    * Chen et al. (2021) utilized the same indicator and shared the same underlying assumptions as Rouder et al. (2019).
    
### Demo preparation ###

Here I'll use a dataset of Stroop task from Hedge et al. (2018) to demonstrate the hierarchical modeling approach to estimating test-retest reliability. The data were collected from 47 subjects who performed Stroop tasks with 2 conditions (congruent, incongruent), 2 sessions (about three weeks apart), 240 trial per condition per session. First, we download the Stroop data from this publicly accessible [site](https://osf.io/cwzds/) and put them in a directory called `stroop/`. Then, let's steal some `R` code (with slight modification) from Haines's nice [blog](http://haines-lab.com/post/2019-05-29-thinking-generatively-why-do-we-use-atheoretical-statistical-models-to-test-substantive-psychological-theories/thinking-generatively-why-do-we-use-atheoretical-statistical-models-to-test-substantive-psychological-theories/) and wrange the data a little bit:

```{r }
library(foreach); library(dplyr); library(tidyr)

data_path <- "stroop/"  # here I assume the download data are stored under directory 'stroop'
files_t1 <- list.files(data_path, pattern = "*1.csv")
files_t2 <- list.files(data_path, pattern = "*2.csv")

# Create a data frame in long format
long_stroop <- foreach(i=seq_along(files_t1), .combine = "rbind") %do% {
  # repetition/session 1
  tmp_t1 <- read.csv(file.path(data_path, files_t1[i]), header = F) %>% mutate(subj_num = i, time = 1)
  # repetition/session 2
  tmp_t2 <- read.csv(file.path(data_path, files_t2[i]), header = F) %>% mutate(subj_num = i, time = 2)
  # Condition (0: congruent, 1: neutral, 2: incongruent), Correct (1) or incorrect (0), RT (seconds)
  names(tmp_t1)[1:6] <- names(tmp_t2)[1:6] <- c("Block", "Trial", "Unused", "Condition", "Correct", "RT")
  rbind(tmp_t1, tmp_t2)
}

dat <- long_stroop[long_stroop$Condition!=1,]
```

Next, to apply the data to our hierarchical model, we flatten the two factors (condition and session) of $2\times 2$ structure into a dimension of four combinators with the following `R` code:

```{r}
dat[dat$Condition==2,'Condition'] <- 'inc' # incongruent
dat[dat$Condition==0,'Condition'] <- 'con' # congruent
dat$com <- paste0(dat$Condition, dat$time) # flatten the two factors
dat$sub <- paste0('s',dat$subj_num).       # subjects
write.table(dat[,c('sub', 'con', 'RT')], file = "stroop.txt", append = FALSE, quote = F, row.names = F)
```
Now, we have a data table called `stroop.txt` with the few lines like this:

```{r}
sub con RT
s1 con1 0.89078
s1 con1 0.72425
s1 con1 0.49442
s1 inc1 0.80221
s1 inc1 0.8327
...
```

With 47 subjects, 2 conditions, 2 sessions, 240 trial per condition per session, there are total 45120 rows in the data table. We purposefully flatten the two factor so that we have a factor with four levels:

```{r}
levels(dat$com)
```
These four levels of "con1", "con2", "inc1", and "inc2" correspond to congruent during session 1, congruent during session 2, incongruent during session 1, and incongruent during session 2. This coding will allow us to more intuitively/directly parameterize the correlations among all the four combinations. Our goal is to assess the test-retest reliability about the contrast between incongruent and congruent. Keep in mind that the conventional ICC only gave a lackluster reliability estimate of approximately 0.5 (Hedge et al., 2018; Haines et al., 2020; Chen et al., 2021).

OK, now we're ready for the next adventure.

### Estimating test-retest reliability using `brms` ###

Let's implement our hierarchical model with the newly obtained data `stroop.txt`. Run the following `R` code:

```{r}
dat <- read.table('stroop.txt', header=T)
library('brms')
options(mc.cores = parallel::detectCores())
m <- brm(bf(RT ~ 0+com+(0+com|sub), sigma ~ 0+com+(0+com|sub)), data=dat, 
         family=exgaussian, chains = 4, iter=1000)
save.image(file='stroop.RData')
```

You may be surprised to notice how simple the code is. The only line that codes our model using `brm` is quite straightforward (if you're familiar with the specification grammar used by the `R` package `lme4`) and it directly maps the data to our hierarchical model. Note that I fit the trial-level effects with an exponentially-modified Gaussian distribution for the probability density $\mathcal D$ in our hierarchicala model above. This implementation may take a few hours (within-chain parallelization would shorten the runtime), so leave your computer alone and come back later.

Let's check the results and make sure all the chains behaved properly. The following code
```{r}
load('stroop.RData')
summary(m)
```

reveals the summarized results:

```{r}
Group-Level Effects: 
~sub (Number of levels: 47) 
                                 Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS Tail_ESS
sd(comcon1)                         0.056     0.006    0.046    0.068 1.014      338      611
sd(comcon2)                         0.061     0.006    0.050    0.074 1.020      373      657
sd(cominc1)                         0.080     0.008    0.066    0.099 1.009      424      918
sd(cominc2)                         0.070     0.007    0.057    0.085 1.020      408      739
sd(sigma_comcon1)                   0.426     0.046    0.350    0.526 1.003      713      979
sd(sigma_comcon2)                   0.491     0.052    0.398    0.606 1.008      698     1059
sd(sigma_cominc1)                   0.415     0.046    0.336    0.515 1.003      719     1245
sd(sigma_cominc2)                   0.433     0.046    0.352    0.534 1.010      728     1163
cor(comcon1,comcon2)                0.738     0.071    0.584    0.854 1.017      390      628
cor(comcon1,cominc1)                0.927     0.027    0.865    0.966 1.005      656     1238
cor(comcon2,cominc1)                0.597     0.097    0.370    0.757 1.016      458      782
cor(comcon1,cominc2)                0.737     0.073    0.572    0.853 1.013      480      603
cor(comcon2,cominc2)                0.948     0.019    0.902    0.977 1.013      942     1341
cor(cominc1,cominc2)                0.695     0.082    0.504    0.826 1.012      531      829
cor(sigma_comcon1,sigma_comcon2)    0.815     0.061    0.676    0.915 1.005      745     1013
cor(sigma_comcon1,sigma_cominc1)    0.921     0.035    0.840    0.974 1.001      799     1156
cor(sigma_comcon2,sigma_cominc1)    0.727     0.082    0.534    0.856 1.002      949     1095
cor(sigma_comcon1,sigma_cominc2)    0.755     0.076    0.588    0.879 1.002      796      834
cor(sigma_comcon2,sigma_cominc2)    0.956     0.024    0.896    0.987 1.001     1168     1269
cor(sigma_cominc1,sigma_cominc2)    0.757     0.076    0.577    0.874 1.003      981     1235

Population-Level Effects: 
              Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS Tail_ESS
comcon1          0.642     0.008    0.626    0.657 1.009      176      416
comcon2          0.621     0.009    0.602    0.638 1.012      282      550
cominc1          0.705     0.011    0.682    0.727 1.007      182      457
cominc2          0.662     0.010    0.642    0.683 1.011      278      568
sigma_comcon1   -2.729     0.069   -2.866   -2.599 1.011      303      632
sigma_comcon2   -2.799     0.078   -2.953   -2.647 1.014      316      775
sigma_cominc1   -2.251     0.067   -2.383   -2.122 1.011      326      732
sigma_cominc2   -2.437     0.070   -2.575   -2.296 1.016      307      787

Family Specific Parameters: 
     Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS Tail_ESS
beta    0.187     0.001    0.184    0.189 1.000     3948     1820
```

In the above summary above, there is a lot of information to unpackage at various levels (population, condition, correlation, standard deviation, etc). But first, we should be happy that the four chains were well-behaved ($\hat R < 1.05$). In addition, we can use posterior predictive checks to verify the model quality:

```{r}
pp_check(m, ndraws = 100)
```

which shows that our model did a pretty good job - the simulated data (green cloud) based on the model fit well with the original RT data (black density curve):

<img alt="alt_text" width="360px" src="https://afni.nimh.nih.gov/sscc/staff/gangc/pub/ppc1.jpg" />

Is there any room for model improvement? Remeber that we used exponentially-modified Gaussian to fit the data (distribution $\mathcal D$ in the model) at the trial level. One may try other distributions such as shifted log-normal (as prefered in Haies et al. (2020)), Gaussian, inversge Gaussian Student's $t$, and (shifted) log-normal. Those alternative distributions could not compete with the exponentially-modified Gaussian as visually illustrated through posterior predictive checks as Fig. 5 in Chen et al. (2021). Model comparisons among these models can also be quantitively assessed through leave-one-out cross-validation using the function `loo` in `brms`. Keep in mind that even though the exponentially-modified Gaussian worked well for this particular dataset, a different disttribution might be appropriate for another dataset.

We should not forget our ultimate goal: estimating test-retest reliability! How to extract the information from the model output? Remember those four levels of "con1", "con2", "inc1", and "inc2" correspond to congruent during session 1, congruent during session 2, incongruent during session 1. Since in the current context, we are interested in the test-retest reliability about the contrast between incongruent and congruent. So we want to extract those model components of $(\mu_{11s}, \mu_{21s}, \mu_{12s}, \mu_{22s})$, and then obtain the correlation between $\mu_{21s}- \mu_{11s}$ and $\mu_{22s}- \mu_{21s})$. Here comes our finale:

```\{r}
ge <- ranef(m, summary = FALSE) # extract subject-Level effects
trr <- rep(0, 2000)
for(ii in 1:2000) trr[ii] <- cor(ge[["sub"]][ii, ,"cominc1"]-ge[["sub"]][ii, ,"comcon1"], 
                                 ge[["sub"]][ii, ,"cominc2"]-ge[["sub"]][ii, ,"comcon2"])
dens <- density(trr)
plot(density(trr), xlim=c(0.4,1), xlab='Test-Retest Reliability')
dens$x[which.max(dens$y)]  # show the peak of the density curve
```

The plot below shows the posterior distribution of test-retest reliability for cognitive inhibition effect (reaction time difference between incongruent and congruent tasks). Based on our hierarchical model, the mode (peak) for the test-retest reliability of the Stroop dataset is 0.82. This indicates that he underestimation by the conventional ICC(3,1) $\simeq 0.5$ is quite substantial. The "reliability crisis" in psychometrics and neuroimaging might be partly related to improper modeling. The reason for this large extent of underestimation is due to the substantial amount of cross-trial variablity compared to cross-subject variability. See more explanation in Chen et al. (2021) regarding the intriguing issue of cross-trial variablity as well as the crucial role of trial sample size relative to the subject sample size.

<img alt="alt_text" width="360px" src="https://afni.nimh.nih.gov/sscc/staff/gangc/pub/trr1.jpg" />

One nice aspect of our parameterization is the easy extraction for an effect of interest. Since the model is directly parameterized with the four combinations of "con1", "con2", "inc1", and "inc2", we could readily obtain other effects such as the average between the two conditions or each individual condition. Needless to say, the correlation structure among the four combinations is fully captured in the hierachical model.
