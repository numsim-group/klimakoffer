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

As a next step, we want to incorporate the annual change of the insolation and consider a simplified EBM, where we only average in space
but not in time. We apply the area avenging to our coefficients to get
$$
\overline{C} \partialderiv{\overline{T}}{t} + A(CO_2) + B \overline{T} = \overline{Q_{\alpha} S}(t),
$$

where we keep the time dependence of the
solar source term

@@colbox-blue
**Remark:** We get $\overline{Q_{\alpha} S}(t)$ by computing the area average of $(1-\alpha)S(x,t)$ at
every single time step separately.
@@

As it is a cumbersome to compute the 
analytical solution, we take the chance to
consider numerical approximation methods,
so-called ODE solves, to numerically calculate 
an (approximate) ODE solution.

Typically, in ODE literature, we consider a
problem of the form
\begin{align}
y' (t) &= f(y,t),~~ t \ge 0
\\
y(t=0) &= y_0,
\end{align}
which we get for the choice
\begin{align}
y(t) &= \overline{T}(t)
\\
f(y,t) &= \frac{1}{\overline{C}} \left( \overline{Q_{\alpha} S}(t) -A - B \, \overline{T}  \right).
\end{align}

@@colbox-blue
**Remark:** There is extensive knowledge on the theory on existence and uniqueness of solutions for ordinary differential equations in the literature.
The existence and uniqueness of solutions for ODEs depends on the continuity properties of the right-hand-side $f(y,t)$.
The most famous results are the Lemma of Peano (existence of solutions) and the Lemma of Picard-Lindelöf (uniqueness of solutions).
We refer the interested reader to the standard literature on ODEs.
@@

Some examples:
> Arnold, V. I. (1992). Ordinary differential equations. Springer Science & Business Media.

> Butcher, J. C. (2016). Numerical methods for ordinary differential equations. John Wiley & Sons.

## ODE Solvers

**Comment:** We emphasize that there are
typically full courses dedicated to the topic
of construction and analysis of numerical
methods for ODEs. Typically, the focus is
on the state-of-the-art method class: the
Runge--Kutta Method.
In this course we focus on simpler methods,
namely the forward and backward Euler schemes, that are sufficient in terms of efficiency
for the problems we are interested in. However,
we encourage the interested reader to look for more efficient/accurate methods and apply them for our problem.


### The Euler method

There are many ways to derive and motivate
the Euler time-integration method.
We choose the idea of an approximation to the
time integral.

Consider the ODE
$$
y'(t) = f(y,t).
$$
The solution is given by
$$
y(t - y(t_0)) = \int_{t_0}^t y'(t) \d t = \int_{t_0}^t f(y,t) \d t.
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

We can use the left-sided rule:
$$
y(t_1) \approx y_1 = y_0 + (t_1 - t_0) f(y_0,t_0),
$$
the right-sided rule:
$$
y(t_1) \approx y_1 = y_0 + (t_1 - t_0) f(y_1,t_1).
$$

Now that we have an approximation for $y_1$, we can proceed analogously for $y_2$, $y_3$, $\ldots$,
\begin{align}
\text{Euler~(explicit):}& &y_{j+1} &= y_j + (t_{j+1} - t_j) f(y_j,t_j)
\\
\text{Euler~(implicit):}& &y_{j+1} &= y_j + (t_{j+1} - t_j) f(y_{j+1},t_{j+1}).
\end{align}

The first variant is called _explicit_
as the new value $y_{j+1}$ can be directly
computed with the information that is known at $t_j$. The second variant is
_implicit_ as the new unknown solution
$y_{j+1}$ also appears on the right-hand side.
Hence, an algebraic equation needs to be
solved for every time step.

These two methods are the simplest variants
of an ODE solver and are formally first-order accurate approximations. This
means that, in general, the solutions of
the Euler methods are **not** exact
solutions of the ODE, but have an error.
First-order accuracy means that the
error of the approximation depends **linearly**
on the time-step size $\Delta t = t_{j+1} - t_j$. Hence, by decreasing
the time-step size, the error gets smaller.
However, at the same time, the number of
time steps $\ntime$ gets larger, i.e., the
amount of computations increases. 
As a consequence, we aim to choose the time-step size as large as possible, but small enough to obtain
the desired accuracy in our numerical
simulation.

@@colbox-blue
**Remark:**
@@