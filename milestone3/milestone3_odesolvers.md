+++
title = "Milestone 3"
hascode = true
rss = "Description"
rss_title = "Milestone 3"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 3 - ODE Solvers

\toc

---

## Time dependent source term

As a next step, we want to incorporate the annual change of the insolation and consider a simplified EBM, where we only average the coefficients in space but not in time. We apply the area averaging to our coefficients to get
$$
\label{eq:ODE2}
\overline{C} \deriv{T_A}{t} + A(CO_2) + B\, T_A = \overline{S_{sol}}(t),
$$

where we keep the time dependence of the solar source term.

@@colbox-blue
**Remark:** We get $\overline{S_{sol}}(t)$ by computing the area average of $S_{sol}(x,t) = (1-\alpha(x))S(x,t)$ at
every single time step separately.
@@

As mentioned in the previous sections, due to the complexity of the solar forcing term, this model is somewhat difficult to solve analytically. However, this time, we take the chance to consider numerical approximation methods, so-called ODE solvers, to numerically calculate an (approximate) ODE solution to \eqref{eq:ODE2}.

Typically, in ODE literature, a problem of the following form is considered
\begin{align}
y' (t) &= f(y,t),~~ t \ge 0,
\\
y(t=0) &= y_0,
\end{align}
where $y(t)$ is the unknown function we want to compute with initial conditions $y_0$ and the right-hand-side operator $f(y,t)$. We can reformulate our problem \eqref{eq:ODE2} into this typical ODE notation for the choice
\begin{align}
y(t) &= T_A(t),
\\
f(y,t) &= \frac{1}{\overline{C}} \left( \overline{S_{sol}}(t) -A - B \, T_A  \right).
\end{align}

@@colbox-blue
**Remark:** There is extensive knowledge on the theory on existence and uniqueness of solutions for ordinary differential equations in the literature.
The existence and uniqueness of solutions for ODEs depends on the continuity properties of the right-hand-side $f(y,t)$.
The most famous results are the Lemma of Peano (existence of solutions) and the Lemma of Picard-Lindelöf (uniqueness of solutions).
We refer the interested reader to the standard literature on ODEs.

Some examples:
> Arnold, V. I. (1992). Ordinary differential equations. Springer Science & Business Media.

> Butcher, J. C. (2016). Numerical methods for ordinary differential equations. John Wiley & Sons.
@@


## ODE Solvers

**Comment:** We emphasize that there are typically full courses dedicated to the topic
of construction and analysis of numerical methods for ODEs. The focus is often on the state-of-the-art method class: the
Runge--Kutta Method (implicit and explicit versions).
In this course we focus on simpler methods, namely the forward and backward Euler schemes, that are sufficient in terms of efficiency for the problems we are interested in. However, we encourage the interested reader to look for more sophisticated methods and apply them for our problem and see if there is an efficiency gain possible.


### The Euler method

There are many ways to derive and motivate the Euler time-integration method.
We choose the idea of an approximation to the time integral when giving the integral form of the ODE solution.

Consider the ODE
$$
y'(t) = f(y,t).
$$
Formally, the solution is given by
$$
y(t) - y(t_0) = \int_{t_0}^t y'(t) \d t = \int_{t_0}^t f(y,t) \d t.
$$

The problem is of course that the solution $y(t)$ is needed to compute the right-hand-side time integral.

We first consider a discretization of the
time axis into small time steps $\Delta t$

\fig{/assets/milestone3/TimeGrid.png}

We start with the known initial value
$$
y(t=t_0) = y_0.
$$
To obtain the solution at the next time step, $y(t_1)$, we look at the integral formula
$$
y(t_1) - y_0 = \int_{t_0}^{t_1} f(y,t) \d t,
$$
and make the simplest approximation to the integral via a rectangular quadrature rule:

\fig{/assets/milestone3/EulerIntegral.png}

We can use the left-sided rule,
$$
y(t_1) \approx y_1 = y_0 + (t_1 - t_0) f(y_0,t_0),
$$
or the right-sided rule,
$$
y(t_1) \approx y_1 = y_0 + (t_1 - t_0) f(y_1,t_1).
$$

Now that we have an approximation for $y_1$, we can proceed analogously for $y_2$, $y_3$, $\ldots$
\begin{align}
\text{Euler~(explicit):}& &y_{k+1} &= y_k + (t_{k+1} - t_k) f(y_k,t_k)
\\
\text{Euler~(implicit):}& &y_{k+1} &= y_k + (t_{k+1} - t_k) f(y_{k+1},t_{k+1}).
\end{align}

The first variant is called _explicit_
as the new value $y_{k+1}$ can be directly
computed with the information that is known at the current time level $t_k$. The second variant is
_implicit_ as the new unknown solution $y_{k+1}$ also appears on the right-hand side inside the function $f(.,t)$.
Hence, an algebraic equation needs to be solved for every time step. This equation is non-linear or linear, depending on the function $f(.,t)$.

These two methods are the simplest variants of an ODE solver and are formally first-order accurate approximations (can be shown via Taylor expansion). This
means that in general the solutions of the Euler methods are **not** exact solutions of the ODE, but have an (numerical) error.
First-order accuracy means that the error of the numerical approximation depends **linearly** on the time-step size $\Delta t = t_{k+1} - t_k$, i.e., the error behaves like $error\sim \mathcal{O}(\Delta t^1)$. Hence, by decreasing the time-step size, the error gets smaller. However, at the same time, the number of
time steps $\ntime$ gets larger, i.e., the amount of computations increases. 
As a consequence, we aim to choose the time-step size as large as possible, but small enough to obtain
the desired accuracy in our numerical simulation.

@@colbox-blue
**Remark:** It is easy to generate a second-order time-integration method based on the Euler methods. For instance, the implicit Crank--Nicolson scheme is
$$
y_{k+1} = y_k + \frac{t_{k+1} - t_k}{2} \left(f(y_{k},t_{k}) + f(y_{k+1},t_{k+1}) \right),
$$
which results when approximating the integral with a trapezoidal rule, i.e., by considering half of the explicit plus half of the implicit Euler update.
@@

## Stability of the explicit and implicit Euler methods

When deciding between the explicit or implicit Euler variants, one might get tempted to use directly  the much simpler explicit
Euler method: every solution can be directly
computed without solving an algebraic (possibly
non-linear) equation.
However, the simplicity of the explicit nature
comes at a hefty price: the choice of the size
of $\Delta t$ is not arbitrary anymore, but needs
to be small enough due to stability reasons.

As an example, let us consider the ODE
\begin{align}\label{eq:simpleODE}
y'(t) &= \lambda y(t),& &\lambda \in \mathbb{C},  \\
y(0) &= 1,
\end{align}
with the analytical solution $y(t) = e^{\lambda t}$.

Since $\lambda$ is a complex number, we can write it as $\lambda = a + b i$. Therefore, we get
\begin{align}
y(t) &= e^{at} e^{bti}  \\
 &= e^{at} \left( \cos(bt) + i \sin(bt) \right),
\end{align}
which shows that the real part corresponds to the amplitude of the solution (which grows in time for $a>0$ and decreases for $a<0$) and that the imaginary part is the oscillatory part.

\fig{/assets/milestone3/ODEsolution.png}

We can consider a perturbed initial condition,
$$
y_{\delta}(0) = 1+ \delta,
$$
where $\delta$ is small compared to $y_0=1$.
The exact solution is then
$$
y_{\delta}(t) = (1 + \delta) e^{\lambda t}.
$$

If we compare the two solutions,
$$
|y_{\delta}(t) - y(t)| = |\delta e^{\lambda t}|,
$$
we see that only for non-positive real parts, $a \le 0$,
the difference stays small in time even
when considering a small $\delta$ perturbation.
The problem is called stable for $\text{Re}(\lambda) = a \le 0$.
We want to have numerical methods that
are numerically stable when solving stable
problems!

If we apply first the explicit Euler
formula to solve the problem, we get
\begin{align}
y_{k+1} &= y_k + \Delta t \lambda y_k
\\
&= (1 + \Delta t \lambda) y_k.
\end{align}

If we look at the amplitude of the numerical solution we obtain
$$
|y_{k+1}| = |1 + \Delta t \lambda| |y_k|,
$$
which shows that $|1 + \Delta t \lambda|$ is the amplification factor of the explicit Euler method.
We can investigate for which values of $\Delta t$ this factor is smaller than $1$ to obtain:

\fig{/assets/milestone3/Stability_explicit_Euler.png}

> We need that the product $(\Delta t \lambda)$ is within the marked red circle to keep the numerical amplification factor smaller than $1$!

The red circle defines the numerical stability
region of the explicit Euler method and
indicates that for guaranteed numerical
stability, the size of  $\Delta t$ is limited.

@@colbox-blue
**Example:** Consider the ODE
\begin{align}
y'(t) &= -1000\, y(t) \\
y(0) &= 1.
\end{align}

The red circle stability area gives a
maximum time step size of $\Delta t_{\max} = \frac{1}{500}$,
which shows that the time-step size can get
very low for large values of $\lambda$. Such problems
are sometimes also referred to as _stiff problems_.
@@

The advantage of the explicit Euler method is the
simplicity and the computationally cheap algorithm.
The downside is that the size of the time
step can be very low due to numerical 
stability issues.

If we consider now the **implicit** Euler method for our simple ODE \eqref{eq:simpleODE}, we obtain
\begin{align}
y_{k+1} &= y_k + \lambda \Delta t y_{k+1} \\
y_{k+1} &= \frac{y_k}{1- \lambda \Delta t}.
\end{align}
In this simple case, we can directly solve the implicit equation!

The amplification factor is now $\left| \frac{1}{1- \lambda \Delta t} \right|$, which gives the stability region of
the implicit Euler method:

\fig{/assets/milestone3/Stability_implicit_Euler.png}

This shows that the implicit Euler scheme
is stable for all values of $\lambda \Delta t$, except for the
ones inside the blue circle. In particular, the
scheme is _unconditionally_ stable for **all** choices of $\Delta t$ in
case $\text{Re}(\lambda) = a \le 0$.

Hence, the big advantage of the implicit
variant is that $\Delta t$ can be chosen arbitrarily
regarding the stability and only needs
to be chosen for accuracy reasons.

@@colbox-blue
**Remark:** Depending on the stiffness
of the ODE problem and the computational
complexity of the implicit equation, either
the explicit or implicit method is more
efficient.
@@

## Temporal equilibrium simulation
In conclusion, we have now all the tools
available to solve our time dependent
area-averaged EBM with either explicit
or implicit Euler in time.

The only open question is how we can reach equilibrium in time for our simulation.
We need to realize that the source term is periodic in time (with a period of one year).
Hence, we are looking for a numerical solution that changes within the year, but then repeats itself for every following year afterwards. Once we have reached such an annual state, we consider that our simulation of the climate system has reached _equilibrium_ and no further change of the annual temperature will occur unless the problem definition (e.g. the parametrizations) is changed again. 

In practice, we have to compare the yearly tempreature solutions with each other, until the difference from one year to the next is smaller than a _given tolerance_.  Among the different choices of norms to compute the yearly solutions, the simplest option is to compute the Euclidean norm of the data vectors in time.