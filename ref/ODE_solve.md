The dynamics of complex networks are often modeled by a system of coupled ordinary differential equations, such as \citet{barzel2013universality,cornelius2013realistic}, which use the ODE form as follows:

\begin{equation}\label{0ODE}
\frac{d x(t)}{d t}=\left(\begin{array}{c}
\frac{d x_1(t)}{d t} \\
\vdots \\
\frac{d x_p(t)}{d t}
\end{array}\right)=\left(\begin{array}{c}
f_1(x(t)) \\
\vdots \\
f_p(x(t))
\end{array}\right)=F(x(t)).
\end{equation}


Let \( x(t) = \big(x_1(t), \ldots, x_p(t)\big)^{\top} \in \mathbb{R}^p \) denote a system of \( p \) variables of interest, where \( F = \{f_1, \ldots, f_p\} \) represents the set of unknown functionals governing the regulatory relationships among \( x(t) \), and \( t \) indexes time over a standardized interval \( \mathscr{T} = [0, 1] \). In practice, the system \eqref{0ODE} is typically observed at discrete time points \( \{t_1, \ldots, t_n\} \), often with measurement error:

$$
y_i=x\left(t_i\right)+\epsilon_i, i=1, \ldots, n,
$$
where $y_i=\left(y_{i 1}, \ldots, y_{i p}\right)^{\top} \in \mathbb{R}^p$ denotes the observed data, $\epsilon_i=\left(\epsilon_{i 1}, \ldots, \epsilon_{i p}\right)^{\top} \in \mathbb{R}^p$ denotes the vector of random noise that are usually assumed to follow independent normal distribution with mean 0 and variance $\sigma_j^2, j=1, \ldots, p$. 



\subsection{Sparse Additive model}
In many biology systems, fully connected networks are rare due to potential vulnerabilities to external disturbances, and many current algorithms assume sparsity within the GRNs  \citep{allesina2012stability, mercatelli2020gene}.  Specifically, we say that $x_k$ regulates $x_j$ if $f_j$ is a functional of $x_k$. In other words, $x_k$ controls the change of $x_j$ through the functional $F_j$ on the derivative $d x_j / d t$. 
\citet{lu2011high} first proposed using the smoothly clipped absolute deviation (SCAD) penalty to select sparse edges in the above linear ODE form:

\begin{equation}\label{lODE}
X_j^{\prime}(t)=\sum_{i=1}^{p} \theta_{ji} X_i(t), \quad j=1,2, \ldots, p,
\end{equation}
where $\theta = \{\theta_{ji}\}$ quantifies the regulatory interactions and interrelations among genes within a network. In such networks, most of these $\theta_{ji}$ values are expected to be zero, reflecting the sparse connectivity characteristic of gene regulatory networks, and the question here really is how to perform model selection. 


Based on this, \citet{henderson2014network, wu2014sparse} extended their work to a more general additive nonparametric ordinary differential equation model, enabling the modeling of high-dimensional nonlinear gene regulatory networks:

\begin{equation}\label{nODE}
X_j^{\prime}(t)=\mu_j+\sum_{i=1}^{p} f_{ji}\left(X_i(t)\right), \quad j=1,2, \ldots, p.
\end{equation}

In this context, $\mu_i$ represents an intercept term, while $f_{j}$ defined in \eqref{0ODE} can be represented as the additive form of univarite functions $\sum_{i=1}^{p} f_{ji}\left(X_i(t)\right)$, which is used to quantify the nonlinear relationships.  The additive model formulation offers several key advantages. In the context of gene regulatory networks (GRNs), this approach enables: (1) identification of significant effects between individual genes, and (2) decomposition of dependent effect curves. These critical tasks necessitate the use of additive model structures. Moreover, based on the principle of sparsity commonly observed in gene regulatory networks and other biological systems, we typically assume that for each variable (gene) $X_j, j \in [p]$, the number of significant nonlinear effects $f_{ji}(\cdot)$ is small. Thus, for each $j \in [p]$ we consider the function space:

 \begin{align*}
     \mathcal{F}(\mathcal{H},s,p)\coloneqq &\Big\{f:\mathcal{X}^p\to \mathbb{R}\lvert\, f=\sum_{i=1}^{p}f_{i}(x_i),\\
     &\sum_{i=1}^{p} \mathbb{I}(f_{i} \neq 0) \le s,
 	{\rm\ and\ } f_i \in \mathcal{H}, \forall i\in [p]\Big\}.
 \end{align*}


To estimate the functions $f_{ji}$, we follow the idea of \citet{henderson2014network, wu2014sparse} and \citet{chen2017network} and approximate the unknown $f_{ji}$ with a truncated basis expansion. Consider a M-dimensional basis $\psi(x)=\left(\psi_1(x), \ldots, \psi_M(x)\right)^{\mathrm{T}}$, such that

\begin{equation}\label{f_ji}
f_{ji}\left(a_i\right)=\psi\left(a_i\right)^{\mathrm{T}} \theta_{ji}+\delta_{ji}\left(a_i\right), \quad \theta_{ji} \in \mathbb{R}^M,
\end{equation}
where $\delta_{ji}\left(a_i\right)$ denotes the residual. Although infinite-dimensional basis expansions have been employed in some literature, such as \citet{dai2022kernel}, and we posit that each functional component $f_{ji}$ resides within a reproducing kernel Hilbert space (RKHS), empirical evidence demonstrates that a finite-dimensional basis such as \cite{wu2014sparse, chen2017network} is highly advantageous for achieving enhanced sparsity and biological interpretability within our gene regulatory network framework: When modeling growth functions in ecological systems using polynomial approximations, Legendre polynomials offer a natural choice, which is the Gram-Schmidt orthogonalization of $[1,x,x^2,\cdots]$ on $x \in [0,1]$. However, to avoid overfitting while maintaining computational efficiency, we restrict the polynomial order to a finite degree $M$. This motivates the use of truncated Legendre polynomial expansions (up to order $M$) for the functional representation.

Combined with \eqref{nODE} and \eqref{f_ji}, the form of nonlinear additive ODE can be written as

\begin{equation}\label{nODE2}
X_j^{\prime}(t)=\theta_{j 0}+\sum_{i=1}^{p} \psi\left(X_i(t)\right)^{\mathrm{T}} \theta_{ji}+\sum_{i=1}^{p} \delta_{ji}\left(X_i(t)\right), \quad j=1, \ldots, p ,
\end{equation}
In practice, we can use common polynomial bases for $\psi(x)$. Here, we consider the first $ M $ basis functions. Under the RKHS assumption, in fact, when $ M $ is sufficiently large, the residual term becomes negligible. In a Biological network, we can assume that the highest order of interactions among variables is finite and relatively low. 
Thus, we can write the equation in multi-dimensional vector form as follows:


\begin{equation}\label{non-integral-based}
X_j^{\prime}(t)=\theta_{j 0}+\sum_{i=1}^{p} \psi\left(X_i(t)\right)^{\mathrm{T}} \theta_{ji}
\end{equation}

\begin{equation}\label{nODE3}
X^{\prime}(t) \equiv\left[\begin{array}{c}
\frac{d X_1(t)}{d t} \\
\vdots \\
\frac{d X_p(t)}{d t}
\end{array}\right] \approx \left[\begin{array}{c}
\theta_{10} \\
\vdots \\
 \theta_{p0}
\end{array}\right] + \left[\begin{array}{c}
  \sum_{i=1}^{p} \psi\left(X_i(t)\right)^{\mathrm{T}} \theta_{1i}  \\
    \vdots\\
    \sum_{i=1}^{p} \psi\left(X_i(t)\right)^{\mathrm{T}} \theta_{pi}
\end{array}\right], \quad t \in[0,1],
\end{equation}
where $\psi(X) = (\psi_{1}(X),\psi_2(X)\cdots,\psi_{M}(X))$ is M-dimensional known function, and the parameters denote $\{\theta_{jim}\}$, where $i \in [p], i\in [p], m \in [M]$. Moreover, departing from prior work, we assume the coefficients in the Karhunen-Loève (KL) expansion of these functions exhibit sparsity. which leads to an estimation problem that is "Sparse Additive model + Sparse coefficients". In the following, we will introduce a statistical problem related to \eqref{nODE3}. We will also discuss the topic of double sparsity and its application to non-linear ODEs and network structure.

Each regulatory function $f_{ji}(\cdot)$ includes an intercept term, which makes it difficult to identify individual intercepts. Any constant shift in one or more $f_{ji}(\cdot)$ functions can be absorbed into the intercept $\mu_j$ in equation~\eqref{nODE}, making the individual intercepts of $f_{ji}(\cdot)$ unidentifiable. In practice, only the sum of the intercepts (i.e., $\mu_j = \mu_{j0} + \sum_{i=1}^p \mu_{ji}$, where $\mu_{ji}$ is the intercept of $f_{ji}$) can be estimated. To resolve this issue and enhance biological interpretability, we impose the following constraint:
\begin{equation} \label{eq:constraint}
f_{ji}(0) = 0.
\end{equation}

This condition reflects the biologically meaningful assumption that a gene with zero expression level exerts no regulatory influence on other variables. Furthermore, this constraint removes the redundancy in intercept estimation and leads to a unique decomposition of $\mu_j$ and $f_{ji}(\cdot)$.


\subsection{Multi-task Learning and network structure}
We consider the multivariate regression problem that is simultaneously estimating $p$-dimensional vectors, and this regression form has an important related problem in statistics and machine learning, called multi-task learning \citet{tsy2011}. Specifically, we have $K$ regression tasks:

\begin{equation}\label{multi}
    \left\{
    \begin{aligned}
 y_{1} & =X_{11}\beta_{11}^*+X_{12}\beta_{12}^*+\cdots+X_{1L}\beta_{1L}^* +\xi_{1}, \\
& \vdots \\
 y_{K} & =X_{T1}\beta_{K1}^*+X_{K2}\beta_{K2}^*+\cdots+X_{KM}\beta_{KL}^* +\xi_{K}.
    \end{aligned}
    \right.
\end{equation}



Let $X_{kl} \in \mathbb{R}^n$ (with $k\in [K]$ and $l \in [L]$) denote the $l$th feature in the $k$th task, $\beta_1^*, \ldots, \beta_K^* \in \mathbb{R}^{L}$ be the parameters of interest, and $y_1, \ldots, y_K \in \mathbb{R}^{n}$ be the $n$-dimensional response vectors. The noise vectors $\xi_1, \ldots, \xi_K$ are assumed to be independent and identically distributed with zero mean. 

We assume that each coefficient vector $\beta^*_{.,1}, \ldots, \beta^*_{., L}$ can be regarded as representing $L$ workers, and that the total number of workers participating across all tasks does not exceed $s$, that is

\begin{equation}\label{group}
    \|\beta\|_{2,0} = \sum_{l=1}^L \mathbb{I}(\|\beta^*_{.,l}\|_2 \neq 0) \leq s.
\end{equation}

 Moreover, the number of task that a worker participate in is no more than $ s_0$, that is
 
\begin{equation}\label{element}
\|\beta_{.,l}\|_0 = \sum_{k=1}^{K}\mathbb{I}(\beta^*_{kl} \neq 0)  \leq s_0,~\forall l \in [L].
\end{equation}

Multitask learning has an equivalent form in linear regression. Let the response variable be $\mathbf{Y}^{\top} = [y_1^{\top}, \ldots, y_K^{\top}] \in \mathbb{R}^{nK}$, the block diagonal design matrix be $\mathbf{X} = \text{diag}(X_1, \ldots, X_K) \in \mathbb{R}^{nK \times KL}$, the parameter vector be $\mathbf{B}^{\top} = [\beta_1^{\top}, \ldots, \beta_K^{\top}] \in \mathbb{R}^{KL}$, and the random noise be $\Xi^{\top} = [\xi_1^{\top}, \ldots, \xi_K^{\top}] \in \mathbb{R}^{nK}$. Then, model \eqref{multi} transforms into a double sparse linear regression model.


$$
\mathbf{Y} = \mathbf{X}\mathbf{B}+ \Xi.
$$

It is worth noting that this framework has been extensively studied both theoretically (e.g., in multiclass classification \citet{abramovich2024classification}) and in practical applications, e.g., in Alzheimer’s disease diagnosis \citep{liu2018modeling, peng2019structured}.

Returning to the problem of gene network construction, the object we aim to estimate is a graph $G = (V, E)$, where $V$ and $E$ denote the set of vertices and edges, respectively. In equation \ref{nODE3}, we consider to estimate the parameters $\theta_{jim}$, $j \in [p], i \in [p], m \in [M]$. Let $ i $ and $ j $ denote the vertices of a graph, and let $ \theta_{ij} $ represent the edge between vertex $ i $ and vertex $ j $.
\begin{assumption}\label{assumption_sparse}
 We apply two sparsity assumptions \eqref{group} and \eqref{element}  to the estimation of sparse graphs：
    \begin{itemize}
    \item Let $\theta_{ji,\cdot}$ denote an $M$-dimensional vector representing the total effect of vertex $i$ on vertex $j$, which corresponds to the coefficients in a functional basis expansion from a nonparametric estimation perspective. We conclude that the edge $(j,i)$ does not exist when $\|\theta_{ji,\cdot}\|_2 = 0$. For graph estimation problems where sparsity of edges is assumed, this leads to:
   \begin{equation}\label{assm_sparse1}
      \sum_{i=1}^{p}\sum_{j=1}^{p}\mathbb{I}(\|\theta_{ji,\cdot}\|_2 \ne 0)<s.
   \end{equation}


    \item Let $\theta_{jim}$ represent the collection of all edges and their higher-order interactions. While we consider effects up to order $M$, we impose sparsity on the total number of significant effects, that is 
   \begin{equation}\label{assm_sparse2}
    \sum_{j = 1}^{p}\sum_{i = 1}^{p}\sum_{m = 1}^{M}\mathbb{I}(\theta_{jim} \ne 0) < s'.
     \end{equation}
    From a biological perspective, when higher-order effects are significant, their constituent linear effects are typically also significant, though the converse may not hold. Although our focus is on higher-order interactions, we penalize these effects when lower-order terms sufficiently capture the dynamical behavior, following the principle of Occam's razor.
\end{itemize}
\end{assumption}

Therefore, through Model \eqref{multi} and the sparsity assumptions \eqref{group} and \eqref{element}, we can simultaneously identify edges and vertices with significant effects in the network. We emphasize that based on the structural sparsity assumption of the parameter set $\{\theta_{jim}\}$ and its role in network reconstruction,  existing methods such as \citet{wu2014sparse,chen2017network}  have only considered group-wise sparsity for $\theta_{ji,\cdot}$, and such formulation naturally leads to the construction of a group lasso regression model \citet{Y2006}. We demonstrate that simultaneous consideration of both sparsity structures leads to significant performance improvement.  The conventional approach would employ the sparse group lasso framework \citep{friedman2010note, simon2013sparse, cai2019sparse}:
\begin{equation}
\min_{\theta} \left\{ \mathcal{L}(\theta) + \lambda_1 \sum_{g \in \mathcal{G}} \|\theta_g\|_2 + \lambda_2 \|\theta\|_1 \right\}.
\end{equation}

However, our proposed estimator for doubly-sparse structures builds upon \citet{zhang2023minimax}'s innovative iterative hard-thresholding algorithm. Theoretical and empirical results show that hard-thresholding outperforms lasso-type estimators in graph recovery problems. However, unlike \citet{zhang2023minimax} that directly considers linear regression problems, we need to transform the original nonlinear ODE problem through a series of steps. We will detail our algorithmic procedure in the following section.

