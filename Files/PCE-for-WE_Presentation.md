% Polynomial Chaos Expansions for Wind Energy: Tutorial
% **Juan P. Murcia:** PhD. student, jumu@dtu.dk \newline
  \ **Pierre-E. Réthoré:** Supervisor, Senior Scientist \newline

------------------

# Outline

\begin{enumerate}
\item (A very short) Introduction to Polynomial Chaos Expansions
\item A1. 2D Not correlated Normal
\item A2. 2D MvNormal
\item A3. 2D Conditionally Correlated
\item A4. 2D Copula Correlated
\item Conclusions
\end{enumerate}

------------------

# Uncertainty propagation problem

\begin{textblock*}{14cm}(-0.5cm,2cm)
\resizebox{\columnwidth}{!}{
\smartdiagramset{
    font=\Large,
    text width=4cm,
    module x sep=5,
    border color=white,
    set color list={white,red!50!gray,white},
    back arrow disabled=true,
    uniform arrow color=true,
    arrow color=black,
    arrow tip=latex
  }
\smartdiagram[flow diagram:horizontal]{
  Input variables \\ $\vx \in \mathbb{R}^M$,
  Model \\ $\model(\vx)$,
  Output variables \\ $\vy = \model(\vx) \in \mathbb{R}^L$}
}
\end{textblock*}

\begin{textblock*}{6cm}(1cm,6cm)
\includegraphics[height=0.4\textheight]{Figures/pdf_x.pdf}
\end{textblock*}

\begin{textblock*}{14cm}(-0.5cm,4.5cm)
\resizebox{\columnwidth}{!}{
\smartdiagramset{
    font=\Large,
    text width=4cm,
    module x sep=5,
    border color=white,
    set color list={white,red!50!gray,white},
    back arrow disabled=true,
    uniform arrow color=true,
    arrow color=black,
    arrow tip=latex
  }
\smartdiagram[flow diagram:horizontal]{
  \textbf{Random} inputs \\ $\pdf(\vx)$,
  Model \\ $\model(\vx)$,
  \textbf{Random} outputs \\ $\pdf(\vy)$}
}
\end{textblock*}

\begin{textblock*}{6cm}(10.5cm,6.5cm)
{\fontsize{3cm}{3.5cm}\selectfont ?}
\end{textblock*}


------------------

# (A very short) Introduction to PCE

**Monte-Carlo simulation**

- Obtain a response sample by evaluating the model in each input realization: $\vy_i = \model(\vx_i)$
- Input sample, $\vx_i$, can be generated using advanced sampling methods such as: Latin Hypercube sampling (LHS), Halton or Hammersley sequences.
- \color{dtured}{\textbf{Pros:}} \color{black} Very robust and easy to implement and parallelize
- \color{dtured}{\textbf{Cons:}} \color{black} Convergence is slow ($\propto N^{-1/2}$)

**Polynomial Chaos expansion**

- Build a polynomial surrogate of the model: $y(\vx) \approx \sum c_{l} \, \mpoly_{l} (\vx)$
- A polynomial basis, $\mpoly_{l} (\vx)$, is built with respect to $\pdf(\vx)$
- The model is evaluated, $\vy_i = \model(\vx_i)$, and \color{dtured}{\emph{projected/fitted}} \color{black} to the polynomial basis.
- The mean $\expect(\vy)$, variance $\variance(\vy)$ and Sobol's sensitivity index for each input $S_i$ are obtain from the polynomial coefficients $c_{l}$.
<!-- - A MC sample can be generated using the polynomial surrogate. -->
- \color{dtured}{\textbf{Pros:}} \color{black} Convergence is fast ($\propto N^{-m}$, $m>1$, $m$ is problem dependent)
- \color{dtured}{\textbf{Cons:}}  \color{black} How to define the order of the polynomials in each variable? How to avoid over fitting the model (Gibbs oscillations)?

------------------

# (A very short) Introduction to PCE

**Single uncertain variable** $x$
$$y(x) \approx \sum\limits_{l = 0}^{P} c_{l} \, \mpoly_{l} (x)$$

Define an inner product using $\pdf(x)$:
$$\langle  f, g \rangle  =\int f(x)\,g(x)\, \pdf(x)\,dx $$

The polynomial basis is constructed such that $\mpoly_{0}=1$ and:
$$\langle  \mpoly_{l},\mpoly_{k} \rangle = \begin{cases} 1 &\mbox{if } l=k \\ 0 &\mbox{if } l \neq k \end{cases}$$
$$\langle 1, \mpoly_{l} \rangle = 0 \quad \forall l>0 \quad \iff \quad \int \mpoly_{l}(x)\, \pdf(x)\,dx = 0 \quad \forall l>0 $$

------------------

# Methods to find the coefficients $c_l$

**Semi-Spectral projection (quadrature integration)**

- Use a quadrature rule  to approximate the integrals (nodes, $x_i$ and weights $\omega_i$). Gaussian quadrature is widely used.

\begin{center}
$c_l = \langle  y,  \mpoly_{l} \rangle  =\int y(x)\, \mpoly_{l}(x)\, \pdf(x)\,dx \approx \sum_{i=0}^N \omega_i \, y(x_i) \, \mpoly_{l}(x_i)$
\end{center}

- \color{dtured}{\textbf{Pros:}} \color{black} Very good for low number of dimensions
- \color{dtured}{\textbf{Cons:}} \color{black} Unstable for heavy tailed PDFs. Quadrature rules fail with most correlated variables

**Point collocation (polynomial fit)**

- Generate a small sample and fit the polynomial basis using Least squares or some other optimization method (e.g. LAR, LASSO).
- \color{dtured}{\textbf{Pros:}} \color{black} Very robust. Optimization algorithms are design to handle large number of dimensions (sparsity) and correlated inputs.
- \color{dtured}{\textbf{Cons:}} \color{black} Not as efficient as semi-spectral collocation.

------------------

# How to deal with correlated inputs?  \newline \color{dtured}{\normalsize{Rosenblatt Transformation [Rosenblatt 1952]}}

- Transforms the correlated input variables $(\vx)$ into a multi-dimensional uncorrelated uniform space $(\vw)$. Solve the propagation problem in the uncorrelated space: $\vy(\vx)=\vy(\mathbb{F}_Q^{-1}(\vw))$. Use Legendre polynomials for uniform variables.
- It consists in using the inverse of the CDF of each variable in sequence. Chaospy includes this transformation [Feinberg 2015]. Graph reproduced from Chaospy tutorials.

\begin{figure}[hc!]
  \includegraphics[width=0.8\textwidth]{Figures/Rosenblatt.pdf}
\end{figure}


------------------

# Outline

\begin{enumerate}
\item (A very short) Introduction to Polynomial Chaos Expansions
\item A1. 2D Not correlated Normal
\item A2. 2D MvNormal
\item A3. 2D Conditionally Correlated
\item A4. 2D Copula Correlated
\item Conclusions
\end{enumerate}

------------------

# Case A: A simple model \newline \color{dtured}{\normalsize{Model Description}}

\begin{textblock*}{14cm}(-0.5cm,3cm)
\resizebox{\columnwidth}{!}{
\smartdiagramset{
    font=\Large,
    text width=4cm,
    module x sep=5,
    border color=white,
    set color list={white,red!50!gray,white},
    back arrow disabled=true,
    uniform arrow color=true,
    arrow color=black,
    arrow tip=latex
  }
\smartdiagram[flow diagram:horizontal]{
  Input variables \\ $\vx =(I;a) \in \mathbb{R}^2$ ,
  Model \\ $\model(\vx; t) = I\,e^{-a \, t}$,
  Output variables \\ $\vy = \model(\vx;t) \in \mathbb{R}^L$}
}
\end{textblock*}

$$ \, $$
$$ \, $$
$$ \, $$
$$ \, $$

### Variables without uncertainty

- $t$: Location of the evaluation.

### Uncertain variables

- $I$: Initial condition.
- $a$: Rate of dissipation.

------------------

# Case A1: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Not correlated Normal}}

\begin{multicols}{2}

$ $

$ $

$y = I\,e^{-a \, t}$

$I \sim Normal(\mu_I=8., \sigma_I=2^{1/2})$

$a \sim Normal(\mu_a=0.8, \sigma_a=0.01^{1/2})$

$ $


\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/2_2xNormal_NoCorr_top.pdf}
\end{figure}

\begin{figure}[hc!]
  \includegraphics[height=0.45\textheight]{Figures/2_2xNormal_NoCorr_MC_sample.pdf}
\end{figure}
\end{multicols}

------------------

# Case A1: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Not correlated Normal}}

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/2_2xNormal_NoCorr_bot.pdf}\\
  \includegraphics[width=\textwidth]{Figures/2_2xNormal_NoCorr_Sens.pdf}
\end{figure}

------------------

# Case A1: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Not correlated Normal}}

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/2_2xNormal_NoCorr_bot.pdf}\\
  \includegraphics[width=0.4\textwidth]{Figures/2_Convergence_E_u.pdf}
  \includegraphics[width=0.4\textwidth]{Figures/2_Convergence_V_u.pdf}
\end{figure}

------------------

# Outline

\begin{enumerate}
\item (A very short) Introduction to Polynomial Chaos Expansions
\item A1. 2D Not correlated Normal
\item A2. 2D MvNormal
\item A3. 2D Conditionally Correlated
\item A4. 2D Copula Correlated
\item Conclusions
\end{enumerate}

------------------

# Case A2: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Correlated Normal}}

\begin{multicols}{2}

$ $

$y = I\,e^{-a \, t}$

$\begin{bmatrix}
    a \\
    I
\end{bmatrix} \sim Normal(\mu, \covar)$

$\mu=\begin{bmatrix}
    0.8 \\ 8.
\end{bmatrix} \,\, \covar=\begin{bmatrix}
    0.01 & 0.1 \\ 0.1 & 2.
\end{bmatrix}$

$ $

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/3_MultiNormal_top.pdf}
\end{figure}

\begin{figure}[hc!]
  \includegraphics[height=0.45\textheight]{Figures/3_MultiNormal_MC_sample_sp.pdf}
\end{figure}
\end{multicols}

------------------

# Case A2: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Correlated Normal}}

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/3_MultiNormal_bot.pdf}\\
  \includegraphics[width=\textwidth]{Figures/3_MultiNormal_Sens.pdf}
\end{figure}

------------------

# Case A2: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Correlated Normal}}

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/3_MultiNormal_bot.pdf}\\
  \includegraphics[width=0.4\textwidth]{Figures/3_Convergence_E_u.pdf}
  \includegraphics[width=0.4\textwidth]{Figures/3_Convergence_V_u.pdf}
\end{figure}

 ------------------

# Outline

\begin{enumerate}
\item (A very short) Introduction to Polynomial Chaos Expansions
\item A1. 2D Not correlated Normal
\item A2. 2D MvNormal
\item A3. 2D Conditionally Correlated
\item A4. 2D Copula Correlated
\item Conclusions
\end{enumerate}

------------------

# Case A3: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Conditionally Correlated}}

\begin{multicols}{2}

$ $

$y = I\,e^{-a \, t}$

$a \sim Uniform(0.6,0.8)$

$I \sim Weibull(k=2.,$

$\quad \quad  A=6(0.3+a)^4)$

$ $

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/4_Rosenblatt_top.pdf}
\end{figure}

\begin{figure}[hc!]
  \includegraphics[height=0.45\textheight]{Figures/4_Rosenblatt_MC_sample.pdf}
\end{figure}
\end{multicols}

------------------

# Case A3: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Conditionally Correlated}}

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/4_Rosenblatt_bot.pdf}\\
  \includegraphics[width=\textwidth]{Figures/4_Rosenblatt_Sens.pdf}
\end{figure}

------------------

# Case A3: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Conditionally Correlated}}

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/4_Rosenblatt_bot.pdf}\\
  \includegraphics[width=0.4\textwidth]{Figures/4_Convergence_E_u.pdf}
  \includegraphics[width=0.4\textwidth]{Figures/4_Convergence_V_u.pdf}
\end{figure}

-------------------

# Options: A better surrogate and MC

------------------

# Case A3: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Conditionally Correlated - A better surrogate and MC}}

\begin{multicols}{2}

$ $

$y = I\,e^{-a \, t}$

$a \sim Uniform(0.6,0.8)$

$I \sim Weibull(k=2.,$

$\quad \quad  A=6(0.3+a)^4)$

$ $

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/4_Rosenblatt_top.pdf}
\end{figure}

\begin{figure}[hc!]
  \includegraphics[height=0.45\textheight]{Figures/4b_Rosenblatt_MC_sample.pdf}
\end{figure}
\end{multicols}

------------------

# Case A3: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Conditionally Correlated - A better surrogate and MC}}

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/4b_Rosenblatt_bot.pdf}\\
  \includegraphics[width=\textwidth]{Figures/4b_Rosenblatt_Sens.pdf}
\end{figure}

------------------

# Case A3: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Conditionally Correlated - A better surrogate and MC}}

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/4b_Rosenblatt_bot.pdf}\\
  \includegraphics[width=0.4\textwidth]{Figures/4b_Convergence_E_u.pdf}
  \includegraphics[width=0.4\textwidth]{Figures/4b_Convergence_V_u.pdf}
\end{figure}

 ------------------

# Outline

\begin{enumerate}
\item (A very short) Introduction to Polynomial Chaos Expansions
\item A1. 2D Not correlated Normal
\item A2. 2D MvNormal
\item A3. 2D Conditionally Correlated
\item A4. 2D Copula Correlated
\item Conclusions
\end{enumerate}

------------------

# Case A4: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Copula Correlated}}

\begin{multicols}{2}

$ $

$y = I\,e^{-a \, t}$

$a \sim Uniform(0.6,0.8)$

$I \sim Weibull(k=2.,A=2)$

$Joe(CDF(a),CDF(I),theta=5.)$

$ $

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/5_Joe_copula_top.pdf}
\end{figure}

\begin{figure}[hc!]
  \includegraphics[height=0.45\textheight]{Figures/5_Joe_copula_MC_sample.pdf}
\end{figure}
\end{multicols}

------------------

# Case A4: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Copula Correlated}}

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/5_Joe_copula_bot.pdf}\\
  \includegraphics[width=\textwidth]{Figures/5_Joe_copula_Sens.pdf}
\end{figure}

------------------

# Case A4: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Copula Correlated}}

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/5_Joe_copula_bot.pdf}\\
  \includegraphics[width=0.4\textwidth]{Figures/5_Convergence_E_u.pdf}
  \includegraphics[width=0.4\textwidth]{Figures/5_Convergence_V_u.pdf}
\end{figure}

-------------------

# Options: A better surrogate and MC

------------------

# Case A4: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Copula Correlated - A better surrogate and MC}}

\begin{multicols}{2}

$ $

$y = I\,e^{-a \, t}$

$a \sim Uniform(0.6,0.8)$

$I \sim Weibull(k=2.,A=2)$

$Joe(CDF(a),CDF(I),theta=5.)$

$ $

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/5_Joe_copula_top.pdf}
\end{figure}

\begin{figure}[hc!]
  \includegraphics[height=0.45\textheight]{Figures/5b_Joe_copula_MC_sample.pdf}
\end{figure}
\end{multicols}

------------------

# Case A4: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Copula Correlated - A better surrogate and MC}}

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/5b_Joe_copula_bot.pdf}\\
  \includegraphics[width=\textwidth]{Figures/5b_Joe_copula_Sens.pdf}
\end{figure}

------------------

# Case A4: A simple model \newline \color{dtured}{\normalsize{Inputs: 2D Copula Correlated - A better surrogate and MC}}

\begin{figure}[hc!]
  \includegraphics[width=\textwidth]{Figures/5b_Joe_copula_bot.pdf}\\
  \includegraphics[width=0.4\textwidth]{Figures/5b_Joe_copula_Convergence_E_u.pdf}
  \includegraphics[width=0.4\textwidth]{Figures/5b_Joe_copula_Convergence_V_u.pdf}
\end{figure}

------------------

# Outline

\begin{enumerate}
\item (A very short) Introduction to Polynomial Chaos Expansions
\item A1. 2D Not correlated Normal
\item A2. 2D MvNormal
\item A3. 2D Conditionally Correlated
\item A4. 2D Copula Correlated
\item Conclusions
\end{enumerate}

------------------

# How to avoid over fitting and achieve sparsity? \newline \color{dtured}{\normalsize{Sparse linear model regression}}

- Least Absolute Shrinkage and Selection Operator problem (LASSO) is useful to avoid over-fitting and achieve sparsity in the PCE.
- LASSO is a least squares minimization problem with a $l_1$ penalization on the coefficients ($\vc$):
$$\min_\vc ||\vc \boldsymbol{\mpoly} - \text{\textbf{y}}_k||^2_2 + \alpha_k ||\vc||_1 = \min_\vc \sum_{i=0}^{N-1} \left[ \sum_{l=0}^{N_c-1} c_l \mpoly_l(\vw_i) - y_k(\vx_i) \right]^2 + \alpha_k \sum_{l=0}^{N_c-1} |c_l|$$

\begin{figure}[hc!]
  \includegraphics[width=0.3\textwidth]{Figures/overfitting.pdf}
  \includegraphics[width=0.3\textwidth]{Figures/overfitting_solved.pdf}
\end{figure}


------------------

# How to select the right sparsity? \newline \color{dtured}{\normalsize{k-fold cross validation}}

- It divides the dataset in k groups and uses k-1 groups ("folds") for training and the remaining for validation. Repeat this process until all the groups have been the validation set.
- k-fold cross validation is repeated for multiple values of the sparsity parameter.
- As a result it gives the optimal sparsity parameter ($\alpha_k$).

\begin{figure}[hc!]
  \includegraphics[width=0.3\textwidth]{Figures/overfitting_CV.pdf}
  \includegraphics[width=0.3\textwidth]{Figures/overfitting_solved_CV.pdf}
\end{figure}

------------------

# PCE of complex cases \newline \color{dtured}{\normalsize{Variable transformations steps}}

\begin{textblock*}{14cm}(-0.5cm,2.cm)
\resizebox{\columnwidth}{!}{
\smartdiagramset{
    font=\Large,
    text width=3cm,
    module x sep=4,
    border color=white,
    set color list={white,red!50!gray,white,red!50!gray,white,red!50!gray,white},
    back arrow disabled=true,
    uniform arrow color=true,
    arrow color=black,
    arrow tip=latex
  }
\smartdiagram[flow diagram:horizontal]{
  \textbf{Correlated} \\ inputs \\ $\pdf(\vx)$,
  Rosenblatt \\ transformation,
  Uniform \\ \textbf{independent} \\ $\pdf(\vw)$,
  PCE \\ $\vz=\hat{\model}(\vw)$,
  Overlimited \\ outputs \\ $\pdf(\vz)$,
  Logistic \\ transformation,
  \textbf{Correlated} \\ outputs \\ $\pdf(\vy)$}
}
\end{textblock*}

\vspace{2cm}

**Rosenblatt transformation**

- Used to decorrelate the variables
- Variables can be transformed to independent Uniform or Normal
- Inverse transformation used for MC sample. Use efficient sampling techniques in the unitary uniform uncorrelated space.

**PCE model surrogate**

- Polynomial chaos expansion working on the uncorrelated space.
- Trained using k-Fold validation to avoid over-fitting and prefer lower order polynomials (Least absolute shrinkage and selection operator - LASSO problem).

**Logistic transformation**

- Used to force fixed constrains in the outputs: i.e. to avoid overshoots.
- Can be used to smooth discontinuities and to impose only positive values.

------------------

# References {.allowframebreaks}
\fontsize{6}{6}\selectfont

[Rosenblatt 1952] Rosenblatt, M 1952 Remarks on a multivariate transformation. Annals of Mathematical Statistics Vol. 23, pp 470-472.]

[Feinberg 2015] Feinberg, J., & Langtangen, H. P. (2015). Chaospy: An open source tool for designing methods of uncertainty quantification. Journal of Computational Science, 11, 46-57.

------------------

# Questions?

\begin{textblock*}{14cm}(0cm,0cm)
\includegraphics[width=14cm]{dtu_bg_fiber.png}
\end{textblock*}

\begin{textblock*}{14cm}(2cm,0.5cm)
\color{white}{\large{\textbf{Questions?}}}
\hspace{5cm} \insertDTUWhiteLogo
\end{textblock*}

\begin{textblock*}{14cm}(7cm,5.5cm)
\begin{tikzpicture}[remember picture,overlay]
\node[fill=black, fill opacity=0.9,
      text=white, text opacity=1.0,
      rounded corners=5pt,
      font=\scriptsize,
      align=left] at (0, 0)
      { \begin{tabular}{l|l}
        J. P. Murcia & Technical University of Denmark (DTU) \\
        +45 2339 7790 & Building 101 \\
        jumu@dtu.dk & Ris{\o} Campus \\
        PhD Student & Frederiksborgvej 399\\
        DTU Wind Energy & 4000 Roskilde, Denmark \\
        \end{tabular}
      };
\end{tikzpicture}
\end{textblock*}
