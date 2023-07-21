<!--
Add here global page variables to use throughout your website.
-->
+++
author = "Gregor Gassner and Andrés Rueda-Ramírez"
mintoclevel = 2
@def hasplotly = false
# Add here files or directories that should be ignored by Franklin, otherwise
# these files might be copied and, if markdown, processed by Franklin which
# you might not want. Indicate directories by ending the name with a `/`.
# Base files such as LICENSE.md and README.md are ignored by default.
ignore = ["node_modules/"]

# RSS (the website_{title, descr, url} must be defined to get RSS)
generate_rss = false
website_title = "An Introduction to Climate Modeling"
website_descr = "In this course, you will learn how to implement an energy balance model for climate simulations"
website_url   = "https://www.mi.uni-koeln.de/IntroToClimateModeling/"
+++

<!--
Add here global latex commands to use throughout your pages.
-->
\newcommand{\R}{\mathbb R}
\newcommand{\d}{\mathrm{d}}
\newcommand{\nlong}{\texttt{n\_longitude}}
\newcommand{\nlat}{\texttt{n\_latitude}}
\newcommand{\ndof}{\texttt{NDOF}}
\newcommand{\ntime}{\texttt{n\_timesteps}}
\newcommand{\lat}{\theta}
\newcommand{\Nabla}{\vec{\nabla}}
\newcommand{\colat}{\tilde \theta}
\newcommand{\long}{\varphi}
\newcommand{\scal}[1]{\langle #1 \rangle}
\newcommand{\partialderiv}[2]{ \frac{\partial {#1}}{\partial {#2} } }
\newcommand{\mat}[1]{\underline{\mathbf{#1}}}
\newcommand{\deriv}[2]{ \frac{d {#1}}{d {#2} } }
\newcommand{\col}[2]{~~~<span style="color:~~~#1~~~">~~~!#2~~~</span>~~~}
\newcommand{\diffop}{L}
\newcommand{\diffcoeff}{\widetilde D}
<!-- Name of repository for GitHub pages -->
@def prepath = "IntroToClimateModeling"
